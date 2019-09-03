#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.lc_colors`` module."""

from unittest import TestCase

import numpy as np
from astropy.table import Table, vstack

from analysis import lc_colors


class Regression(TestCase):
    """Test the gaussian regression returns a "correct" fit to the data"""

    @classmethod
    def setUpClass(cls):
        """Create dummy fluxes that are a sin wave in two bands"""

        x1 = np.arange(0, 4 * np.pi, .1)
        x2 = x1 + np.pi / 4

        data_1 = Table({'time': x1, 'flux': np.sin(x1)})
        data_1['fluxerr'] = 0
        data_1['band'] = 'sdssu'

        data_2 = Table({'time': x2, 'flux': np.sin(x2)})
        data_2['fluxerr'] = 0
        data_2['band'] = 'sdssg'

        cls.data = vstack([data_1, data_2])

    def test_band_prediction(self):
        """Test correct flux returned for a single band"""

        band_name = 'sdssu'
        band_data = self.data[self.data['band'] == band_name]

        gp = lc_colors.fit_gaussian_process(self.data)
        gp_flux, gp_err = lc_colors.predict_band_flux(
            gp, band_name, band_data['time'])

        is_close = all(np.isclose(band_data['flux'], gp_flux))
        self.assertTrue(is_close, 'Prediction values differ significantly')

    def test_light_curve_prediction(self):
        """Test correct flux returned for a multiple bands"""

        bands = set(self.data['band'])
        times = [self.data[self.data['band'] == b]['time'] for b in bands]

        gp = lc_colors.fit_gaussian_process(self.data)
        gp_flux, gp_err = lc_colors.predict_light_curve(gp, bands, times)

        for i, band_name in enumerate(bands):
            band_data = self.data[self.data['band'] == band_name]
            self.assertTrue(
                all(np.isclose(band_data['flux'], gp_flux[i])),
                f'Prediction values differ significantly for {band_name}'
            )


class ChiSquared(TestCase):

    @classmethod
    def setUpClass(cls):
        """Create dummy fluxes that are constants in two bands"""

        x1 = np.arange(0, 11)
        x2 = np.arange(5, 16)
        cls.color_time_range = (5, 10)

        data_1 = Table({'time': x1})
        data_1['flux'] = 1
        data_1['fluxerr'] = 1
        data_1['band'] = 'sdssu'

        data_2 = Table({'time': x2})
        data_2['flux'] = 2
        data_2['fluxerr'] = 1
        data_2['band'] = 'sdssg'

        cls.data = vstack([data_1, data_2])

    def test_observed_time_range(self):
        """Test for the correct return for the given data"""

        color_times = lc_colors.get_color_times(self.data, 'sdssu', 'sdssg')
        self.assertSequenceEqual(self.color_time_range, color_times)

    def test_create_empty_output_table(self):
        """Test summary table has correct columns and meta data"""

        band_combos = [('sdssu', 'sdssg'), ('sdssg', 'sdssr')]
        expected_cols = ['obj_id', 'source', 'version']
        for band_tuple in band_combos:
            col_name = '_'.join(band_tuple)
            expected_cols += [col_name, col_name + '_err']

        table = lc_colors._tabulate.create_empty_output_table(band_combos, '_err')
        self.assertEqual(expected_cols, table.colnames, 'Wrong column names.')
