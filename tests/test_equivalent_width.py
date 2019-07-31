#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.equivalent_width`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from sndata.csp import dr1, dr3

from analysis import equivalent_width as equiv_width
from analysis import models
from analysis import utils

models.register_sources()
dr3.download_module_data()
dr3.register_filters()


def create_test_spectrum(feat_start, feat_end, slope, cont_intercept, feat_intercept):
    """Create a dummy spectrum with an absorption feature and linear continuum

    The y-intercept of the feature is fixed at 1. The feature starts at 0 and
    ends at ``feat_end``.

    Args:
        feat_start       (int): Feature starting wavelength
        feat_end         (int): Feature ending wavelength
        slope          (float): Slope of the continuum
        cont_intercept (float): y-intercept of the continuum flux
        feat_intercept (float): y-intercept of the feature flux

    Returns:
        An array representing wavelength values
        An array representing flux values
        The equivalent width of the feature
    """

    # Define the continuum
    wave = np.arange(feat_start - 10, feat_end + 10, 1)
    flux = slope * wave + cont_intercept

    # Add an absorption feature from 0 to feat_width
    feature_indices = (
            (0 <= wave) & (wave <= feat_end)
    )

    flux[feature_indices] = slope * wave[feature_indices] + feat_intercept
    ew = (feat_end - feat_start) * abs(cont_intercept - feat_intercept)

    return wave, flux, ew


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
        recovered_peak = equiv_width._calc_ew.get_peak_wavelength(
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
            equiv_width._calc_ew.get_peak_wavelength(
                self.wavelength,
                self.flux,
                max_wavelength + 10,
                max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature wavelengths are found"""

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        recovered_lower_peak = equiv_width._calc_ew.get_peak_wavelength(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        recovered_upper_peak = equiv_width._calc_ew.get_peak_wavelength(
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

    def test_continuum_func(self):
        """Test fitting the continuum returns correct slope and intercept"""

        # Create a fake test spectrum
        feat_start = 0
        feat_end = 15
        slope=2
        intercept = 3
        test_wave, test_flux, ew = create_test_spectrum(
            feat_start=feat_start,
            feat_end=feat_end,
            slope=slope,
            cont_intercept=intercept,
            feat_intercept=1) # The feature intercept doesn't matter here

        # Fit the continuum
        cont_func = equiv_width.fit_continuum_func(
            test_wave, test_flux, feat_start, feat_end)

        fit_cont_params = np.polyfit(test_wave, cont_func(test_wave), deg=1)

        self.assertAlmostEqual(
            slope, fit_cont_params[0], places=10,
            msg='Wrong continuum slope.')

        self.assertAlmostEqual(
            intercept, fit_cont_params[1], places=10,
            msg='Wrong continuum y-intercept.')

    def test_calc_ew(self):
        """Test the measured and simulated pew agree"""

        model_start = 0
        model_end = 15
        test_wave, test_flux, ew = create_test_spectrum(
            feat_start=model_start,
            feat_end=model_end,
            slope=0,
            cont_intercept=2,
            feat_intercept=1)

        pew, feat_start, feat_end = equiv_width.calc_pew(
            test_wave, test_flux, feat_start=model_start, feat_end=model_end)

        self.assertEqual(model_start, feat_start, 'Wrong feature start')
        self.assertEqual(model_end, feat_end, 'Wrong feature end')
        self.assertEqual(ew, pew, 'Wrong feature width')


class Tabulation(TestCase):
    """Tests for the compilation of pew results into a single table"""

    @classmethod
    def setUpClass(cls):
        # Load CMFGEN models
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
