#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.utils`` module"""

from unittest import TestCase, skip

import numpy as np
from sndata.csp import dr3

from analysis import utils

dr3.register_filters(force=True)


@skip('This test known to fail due to bug in external dependency sndata')
class TestEffectiveWavelength(TestCase):
    """Test that ``get_effective_wavelength`` returns the correct wavelength"""

    def runTest(self):
        band_names = np.array(dr3.band_names)
        returned_lam = [utils.get_effective_wavelength(b) for b in band_names]

        is_failing = ~np.equal(dr3.lambda_effective, returned_lam)
        delta_wave = np.subtract(dr3.lambda_effective, returned_lam)
        delta_wave_dict = {band: delta for band, delta in
                           zip(band_names[is_failing], delta_wave[is_failing])}

        if delta_wave_dict:
            self.fail(
                f'Returned effective wavelengths differ by: {delta_wave_dict}')


class TestGetCSPt0(TestCase):
    """Test that ``get_csp_t0`` returns the correct value"""

    def test_bad_ids(self):
        """Test the correct error is raised for fake ids and masked values"""

        self.assertRaises(utils.NoCSPData, utils.get_csp_t0, 'fake_id')
        self.assertRaises(utils.NoCSPData, utils.get_csp_t0, '2004dt')

    def test_correct_return(self):
        """Test the correct value is returned for a known target"""

        returned_val = utils.get_csp_t0('2005kc')
        expected_val = 2453698.11
        self.assertEquals(expected_val, returned_val)


class TestGetCSPEBV(TestCase):
    """Test that ``get_csp_ebv`` returns the correct wavelength"""

    def test_bad_ids(self):
        """Test the correct error is raised for fake ids and masked values"""

        self.assertRaises(utils.NoCSPData, utils.get_csp_ebv, 'fake_id')
        self.assertRaises(utils.NoCSPData, utils.get_csp_ebv, '2009ds')

    def test_correct_return(self):
        """Test the correct value is returned for a known target"""

        returned_val = utils.get_csp_ebv('2005kc')
        expected_val = 0.132
        self.assertEquals(expected_val, returned_val)


class TestConvertJD(TestCase):
    """Test various date formats are converted to JD correctly"""

    def test_single_value(self):
        snoopy_date = 500
        mjd_date = snoopy_date + 53000
        jd_date = mjd_date + 2400000.5
        expected_jd = np.array(jd_date)

        self.assertEqual(
            expected_jd, utils.convert_to_jd(snoopy_date),
            'Incorrect date for snoopy format')

        self.assertEqual(
            expected_jd, utils.convert_to_jd(mjd_date),
            'Incorrect date for MJD format')

        self.assertEqual(
            expected_jd, utils.convert_to_jd(jd_date),
            'Incorrect date for JD format')


class TestCSPDataFilter(TestCase):
    """Test ``filter_has_csp_data`` correctly filters data tables"""

    def runTest(self):
        no_t0_id = '2004dt'
        no_ebv_id = '2009ds'
        has_data_id = '2005kc'

        no_t0_data = dr3.get_data_for_id(no_t0_id)
        no_ebv_data = dr3.get_data_for_id(no_ebv_id)
        has_data = dr3.get_data_for_id(has_data_id)

        self.assertEqual(False, utils.filter_has_csp_data(no_t0_data))
        self.assertEqual(False, utils.filter_has_csp_data(no_ebv_data))
        self.assertEqual(True, utils.filter_has_csp_data(has_data))
