#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.lc_colors`` module."""

from unittest import TestCase

import numpy as np
from astropy.table import Table, vstack

from analysis import regression


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

        gp = regression.fit_gaussian_process(self.data)
        gp_flux, gp_err = regression.predict_band_flux(
            gp, band_name, band_data['time'])

        is_close = all(np.isclose(band_data['flux'], gp_flux))
        self.assertTrue(is_close, 'Prediction values differ significantly')

    def test_light_curve_prediction(self):
        """Test correct flux returned for a multiple bands"""

        bands = set(self.data['band'])
        times = [self.data[self.data['band'] == b]['time'] for b in bands]

        gp = regression.fit_gaussian_process(self.data)
        gp_flux, gp_err = regression.predict_light_curve(gp, bands, times)

        for i, band_name in enumerate(bands):
            band_data = self.data[self.data['band'] == band_name]
            self.assertTrue(
                all(np.isclose(band_data['flux'], gp_flux[i])),
                f'Prediction values differ significantly for {band_name}'
            )
