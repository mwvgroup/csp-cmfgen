#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.spectra_chisq`` module."""

from unittest import TestCase

import numpy as np
import sncosmo

from analysis import spectra_chisq
from analysis.exceptions import UnobservedFeature


def create_test_spectra(flux_diff, wave_start, wave_end, error=1):
    """Define arrays of constant observed flux, modeled flux, and errors

    Args:
        flux_diff  (float): Constant offset in flux between model and
        wave_start (float): Starting wavelength ofr spectrum
        wave_end   (float): Ending wavelength ofr spectrum
        error      (float): A constant error value

    Returns:
        An array of wavelengths
        An array of flux values
        An array of values equal to ``error``
        An array of model fluxes
    """

    wave = np.arange(wave_start, wave_end)
    flux = np.ones_like(wave) * 20
    eflux = np.ones_like(wave) * error
    model_flux = flux - flux_diff

    return wave, flux, eflux, model_flux


class ChiSquaredCalculation(TestCase):
    """Test the calculation of chi-squared values"""

    @classmethod
    def setUpClass(cls):
        """Define arrays of constant observed flux, modeled flux, and errors

        Wavelength range goes from 1000 to 1100. Model and observed flux
        values differ by 20.
        """

        cls.flux_diff = 10  # Constant offset between model and observed flux

        start_wave, end_wave = 1000, 1101
        cls.wave, cls.flux, cls.eflux, cls.model_flux = \
            create_test_spectra(cls.flux_diff, start_wave, end_wave)

    def test_correct_return_value(self):
        """Test the correct chi-squared is returned"""

        start, end = min(self.wave), max(self.wave)
        chisq = spectra_chisq.calc_chisq(
            self.wave, self.flux, self.eflux, self.model_flux, start, end)

        expected_chisq = len(self.wave) * self.flux_diff ** 2
        self.assertEqual(expected_chisq, chisq)

    def test_bounds_out_of_range(self):
        """Test exception is raised for out of range start/end wavelengths"""

        start, end = min(self.wave) + 10, max(self.wave) + 10
        args = (self.wave, self.flux, self.eflux, self.model_flux, start, end)
        self.assertRaises(UnobservedFeature, spectra_chisq.calc_chisq, *args)


class BandBoundaries(TestCase):
    """Test the determination of bandpass boundaries"""

    @classmethod
    def setUpClass(cls):
        """Register a fake top-hat bandpass with sncosmo

        Wavelength range goes from 1000 to 1100. Top-hat goes from 1010 to
        1990. Band is registered as 'test_tophat'
        """

        cls.top_hat_height = 10
        cls.wave = np.arange(1000, 1101)
        cls.transmission = np.zeros_like(cls.wave)

        start_index = 10
        end_index = -10
        cls.start_wave = cls.wave[start_index]
        cls.end_wave = cls.wave[end_index]
        cls.transmission[start_index: end_index + 1] = cls.top_hat_height

        band = sncosmo.Bandpass(cls.wave, cls.transmission, name='test_tophat')
        sncosmo.register(band)

    def test_correct_wavelength_bounds(self):
        """Test the function identifies the correct wavelengths"""

        returned_start, returned_end = \
            spectra_chisq.band_limits('test_tophat', self.top_hat_height)

        self.assertEqual(self.start_wave, returned_start)
        self.assertEqual(self.end_wave, returned_end)

    def test_error_unregistered_band(self):
        """Test an exception is raised for an invalid band name"""

        band_name = 'this is a not a valid name'
        args = (band_name, 10)  # second value is a dummy value
        self.assertRaises(Exception, spectra_chisq.band_limits, *args)

    def test_transmission_too_high(self):
        """Test a ValueError is raised for a transmission limit higher that
        the peak of the transmission filter
        """

        args = ('test_tophat', self.top_hat_height + 1)
        self.assertRaises(ValueError, spectra_chisq.band_limits, *args)


class OutputTable(TestCase):
    """Test the generation of an empty table for storing chi-squared data"""

    @classmethod
    def setUpClass(cls):
        """Create an empty table for testing"""

        cls.extra_col_names = ['test_col1', 'test_col2']
        cls.test_table = spectra_chisq._chisq.create_empty_chisq_table(
            cls.extra_col_names)

    def test_col_names(self):
        """Test table has correct column names"""

        base_columns = ['obj_id', 'source', 'version', 'time']
        expected_cols = base_columns + self.extra_col_names
        self.assertSequenceEqual(expected_cols, self.test_table.colnames)

    def test_dtypes(self):
        """Test table has correct data types assigned to each column"""

        base_columns = [np.dtype('U100'), np.dtype('U100'), np.dtype('U100'), float]
        expected_dtype = base_columns + [float for _ in self.extra_col_names]
        returned_dtype = [self.test_table.dtype[c] for c in self.test_table.colnames]
        self.assertSequenceEqual(expected_dtype, returned_dtype)

    def test_is_empty(self):
        """Test the returned table is empty"""

        self.assertEqual(0, len(self.test_table))


class NewTableRow(TestCase):
    """Test the creation of new row entries for a table of chisq results"""

    @classmethod
    def setUpClass(cls):
        # Define two features and bands to determine chi-squared in
        cls.bands = ['sdssu', 'sdssg']
        cls.features = {'feat1': None, 'feat2': None}

        # Define model flux that is within the range of the first band/feature
        # but not the second
        flux_offset, start_wave, end_wave = 10, 3200, 4000
        cls.wave, cls.flux, cls.eflux, cls.model_flux = \
            create_test_spectra(flux_offset, start_wave, end_wave)

    def test_error_for_missing_args(self):
        """Test ValueError raised if ``bands`` or ``features`` not specified"""

        kwargs = dict(
            obj_id='gummy_id',
            model=sncosmo.Model('salt2'),
            time=1,
            wave=[1000],
            flux=[1],
            eflux=[1],
            t0=1
        )

        self.assertRaises(
            ValueError, spectra_chisq._chisq.create_new_table_row, **kwargs)

    def test_row_for_band(self):
        """Test correct row length and mask is returned when bands are given"""

        # Define dummy arguments
        obj_id = 'dummy_id'
        model = sncosmo.Model('salt2')
        time, t0 = 10, 1
        trans_limit = .05

        new_row, mask = spectra_chisq._chisq.create_new_table_row(
            obj_id=obj_id,
            model=model,
            time=time,
            wave=self.wave,
            flux=self.flux,
            eflux=self.eflux,
            t0=t0,
            bands=self.bands,
            trans_limit=trans_limit
        )

        num_base_columns = 4
        expected_len = num_base_columns + len(self.bands)
        expected_mask = [False for _ in range(num_base_columns)] + [False, True]

        self.assertEqual(expected_len, len(new_row), 'New row is wrong length')
        self.assertSequenceEqual(expected_mask, mask, 'New row has wrong mask')
