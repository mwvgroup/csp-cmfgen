#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.equivalent_width`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from sndata.csp import dr1

from analysis import equivalent_width as equiv_width
from analysis import models
from analysis import utils


class FeatureIdentification(TestCase):
    """Test the identification of feature boundaries"""

    @classmethod
    def setUpClass(cls):
        # Create dummy spectrum
        wave_range = (7000, 8001, 100)
        cls.wavelength = np.arange(*wave_range)
        cls.flux = np.ones_like(cls.wavelength)  # Define continuum

        cls.peak_wavelengths = (7100, 7500)
        for peak in cls.peak_wavelengths:
            cls.flux[cls.wavelength == peak] = 10  # Add delta function peak

    def test_peak_coordinates(self):
        """Test the correct peak wavelength is found for a single flux spike"""

        expected_peak = self.peak_wavelengths[0]
        recovered_peak = equiv_width._calc_ew.get_peak_coordinates(
            self.wavelength,
            self.flux,
            expected_peak - 10,
            expected_peak + 10
        )

        self.assertEqual(expected_peak, recovered_peak)

    def test_unobserved_feature(self):
        """Test an error is raise if the feature is out of bounds"""

        max_wavelength = max(self.wavelength)
        with self.assertRaises(equiv_width.UnobservedFeature):
            equiv_width._calc_ew.get_peak_coordinates(
                self.wavelength,
                self.flux,
                max_wavelength + 10,
                max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature is found """

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        recovered_lower_peak = equiv_width._calc_ew.get_peak_coordinates(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        recovered_upper_peak = equiv_width._calc_ew.get_peak_coordinates(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'max'
        )

        self.assertEqual(
            lower_peak_wavelength, recovered_lower_peak, 'Incorrect min peak')

        self.assertEqual(
            upper_peak_wavelength, recovered_upper_peak, 'Incorrect max peak')

    def test_feature_bounds(self):
        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        feature_dict = {
            'lower_blue': lower_peak_wavelength - 10,
            'upper_blue': lower_peak_wavelength + 10,
            'lower_red': upper_peak_wavelength - 10,
            'upper_red': upper_peak_wavelength + 10
        }

        feat_start, feat_end = equiv_width.get_feature_bounds(
            self.wavelength, self.flux, feature_dict)

        self.assertEqual(
            lower_peak_wavelength, feat_start, 'Incorrect min peak')

        self.assertEqual(
            upper_peak_wavelength, feat_end, 'Incorrect max peak')


class EWCalculation(TestCase):
    """Test the measurement of pew values"""

    @classmethod
    def setUpClass(cls):
        """Create a top-hat absorption feature with a depth of 1"""

        # Meta parameters for generating a fake spectrum
        cls.feat_width = 15
        cls.cont_slope = 2
        cls.cont_intercept = 1

        # Define the continuum
        cls.feat_wave = np.arange(-10, cls.feat_width + 10, 1)
        cls.feat_flux = cls.cont_slope * cls.feat_wave + cls.cont_intercept

        # Add an absorption feature from 0 to cls.feat_width
        feature_indices = (
                (0 <= cls.feat_wave) & (cls.feat_wave <= cls.feat_width)
        )

        cls.feat_flux[feature_indices] = 2 * cls.feat_wave[feature_indices]

    def test_continuum_func(self):
        """Test fitting the continuum returns correct slope and intercept"""

        cont_func = equiv_width.fit_continuum_func(
            self.feat_wave, self.feat_flux, 0, self.feat_width)

        fit_cont_params = np.polyfit(self.feat_wave, cont_func(self.feat_wave))
        self.assertEqual(self.cont_slope, fit_cont_params[0],
                         'Wrong continuum slope.')

        self.assertEqual(self.cont_intercept, fit_cont_params[1],
                         'Wrong continuum y-intercept.')

    def test_calc_ew(self):
        """Test the measured and simulated pew agree"""

        pew, feat_start, feat_end = equiv_width.calc_pew(
            self.feat_wave, self.feat_flux, feat_start=0,
            feat_end=self.feat_width)

        self.assertEqual(self.feat_width, pew)
        self.assertEqual(0, feat_start)
        self.assertEqual(self.feat_width, feat_end)


class Tabulation(TestCase):
    """Tests for the compilation of pew results into a single table"""

    @classmethod
    def setUpClass(cls):
        # Load CMFGEN models
        models.register_sources()
        m102 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.02))
        m104 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.04))
        m14 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.4))
        m17 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.7))
        cls.model_list = [m102, m104, m14, m17]
        for model in cls.model_list:
            model.add_effect(sncosmo.F99Dust(), 'ext', 'rest')

        # Load arbitrary data
        dr1.download_module_data()
        demo_id = '2004ef'
        demo_data = dr1.get_data_for_id(demo_id)
        obs_dates, wavelengths, fluxes = utils.parse_spectra_table(demo_data)
        cls.obs_date = obs_dates[0]
        cls.wavelength = wavelengths[0]
        cls.flux = fluxes[0]

    def test_create_pew_table(self):
        """Test the pew summary table has the correct length and columns"""

        table = equiv_width._calc_ew.create_pew_summary_table(self.model_list)
        self.assertSequenceEqual(['model', 'version'], table.colnames)

        # Expect one row per model plus one for observed data
        expected_len = len(self.model_list) + 1
        self.assertEqual(expected_len, len(table), 'Table has incorrect len')

    def test_boundary_fixing(self):
        """Test feature bounds do/don't vary per model when they
        are/aren't fixed
        """

        # Pick an arbitrary feature
        test_feature = 'pW3'

        # Determine pew for free and fixed bounds
        free_ew = equiv_width.tabulate_pew_spectrum(
            self.obs_date, self.wavelength, self.flux, self.model_list, False)

        fixed_ew = equiv_width.tabulate_pew_spectrum(
            self.obs_date, self.wavelength, self.flux, self.model_list, True)

        # Check that free boundaries are not all the same
        num_free_starts = len(set(free_ew[test_feature + '_start']))
        num_free_ends = len(set(free_ew[test_feature + '_start']))
        self.assertGreater(num_free_starts, 1, 'Free lower bounds do not vary')
        self.assertGreater(num_free_ends, 1, 'Free upper bounds do not vary')

        # Check that fixed boundaries are all the same
        num_fixed_starts = len(set(fixed_ew[test_feature + '_start']))
        num_fixed_ends = len(set(fixed_ew[test_feature + '_start']))
        self.assertEqual(num_fixed_starts, 1, 'Fixed lower bounds vary')
        self.assertEqual(num_fixed_ends, 1, 'Fixed upper bounds vary')


# Todo: Compare interpolated pew results and against published pew data
class PEWInterpolation(TestCase):
    pass
