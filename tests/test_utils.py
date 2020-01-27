#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.utils`` module"""

from unittest import TestCase, skip

import numpy as np
import sncosmo
from astropy.table import Column, Table
from sndata.csp import dr3

from analysis import utils

dr3.register_filters(force=True)


class ConvertJD(TestCase):
    """Tests for the ``convert_to_jd`` function"""

    @classmethod
    def setUpClass(cls):
        """Define test dates"""

        cls.snoopy_date = 500
        cls.mjd_date = cls.snoopy_date + 53000
        cls.jd_date = cls.mjd_date + 2400000.5
        cls.expected_jd = np.array(cls.jd_date)

    def test_snoopy_format(self):
        """Test conversion of the snoopy date format to JD"""

        self.assertEqual(
            self.expected_jd, utils.convert_to_jd(self.snoopy_date),
            'Incorrect date for snoopy format')

    def test_mjd_format(self):
        """Test conversion of the MJD format to JD"""

        self.assertEqual(
            self.expected_jd, utils.convert_to_jd(self.mjd_date),
            'Incorrect date for MJD format')

    def test_jd_format(self):
        """Test conversion of the JD format to JD"""

        self.assertEqual(
            self.expected_jd, utils.convert_to_jd(self.jd_date),
            'Incorrect date for JD format')


class GetCSPt0(TestCase):
    """Tests for the ``get_csp_t0`` function"""

    def test_bad_ids(self):
        """Test the correct error is raised for fake ids and masked values"""

        self.assertRaises(utils.NoCSPData, utils.get_csp_t0, 'fake_id')
        self.assertRaises(utils.NoCSPData, utils.get_csp_t0, '2004dt')

    def test_correct_return(self):
        """Test the correct value is returned for a known target"""

        returned_val = utils.get_csp_t0('2005kc')
        expected_val = 2453698.11
        self.assertEqual(expected_val, returned_val)


class GetCSPEBV(TestCase):
    """Tests for the ``get_csp_ebv`` function"""

    def test_bad_ids(self):
        """Test the correct error is raised for fake ids and masked values"""

        self.assertRaises(utils.NoCSPData, utils.get_csp_ebv, 'fake_id')
        self.assertRaises(utils.NoCSPData, utils.get_csp_ebv, '2009ds')

    def test_correct_return(self):
        """Test the correct value is returned for a known target"""

        returned_val = utils.get_csp_ebv('2005kc')
        expected_val = 0.132
        self.assertEqual(expected_val, returned_val)


class CSPDataFilter(TestCase):
    """Tests for the ``filter_has_csp_data`` function"""

    def test_filters_bad_tables(self):
        no_t0_id = '2004dt'
        no_ebv_id = '2009ds'
        has_data_id = '2005kc'

        no_t0_data = dr3.get_data_for_id(no_t0_id)
        no_ebv_data = dr3.get_data_for_id(no_ebv_id)
        has_data = dr3.get_data_for_id(has_data_id)

        self.assertEqual(False, utils.filter_has_csp_data(no_t0_data))
        self.assertEqual(False, utils.filter_has_csp_data(no_ebv_data))
        self.assertEqual(True, utils.filter_has_csp_data(has_data))


@skip('This test known to fail due to bug in external dependency sndata')
class EffectiveWavelength(TestCase):
    """Tests for the ``get_effective_wavelength`` function"""

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


class ParseSpectraTable(TestCase):
    """Tests for the ``parse_spectra_table`` function"""

    @classmethod
    def setUpClass(cls):
        """Create a table of simulated spectra

        The simulated table includes two spectra, each with a different
        observed date and a constant flux value.
        """

        wavelength_range = 7_000, 10_000
        cls.spec1_wave = np.arange(*wavelength_range).tolist()
        cls.spec1_flux = np.full_like(cls.spec1_wave, 100).tolist()
        cls.spec1_date = np.full_like(cls.spec1_wave, 2453100.).tolist()

        cls.spec2_wave = np.arange(*wavelength_range, 2).tolist()
        cls.spec2_flux = np.full_like(cls.spec2_wave, 200).tolist()
        cls.spec2_date = np.full_like(cls.spec2_wave, 2453200.).tolist()

        date = cls.spec1_date + cls.spec2_date
        wavelength = cls.spec1_wave + cls.spec2_wave
        flux = cls.spec1_flux + cls.spec2_flux

        test_table = Table(
            names=['date', 'wavelength', 'flux'],
            data=[date, wavelength, flux]
        )

        cls.returned_date, cls.returned_wave, cls.returned_flux = \
            utils.parse_spectra_table(test_table)

        print(cls.returned_date)

    def test_sorted_by_date(self):
        """Test the returned dates are sorted"""

        self.assertListEqual(
            sorted(self.returned_date), self.returned_date.tolist()
        )

    def test_correct_returned_dates(self):
        """Test the returned dates match the dates from the input table"""

        for i, spectrum_date in enumerate(self.returned_date):
            original_date = [self.spec1_date, self.spec2_date][i][0]
            self.assertEqual(original_date, spectrum_date)

    def test_correct_returned_wavelengths(self):
        """Test returned wavelengths match wavelengths from the input table"""

        for i, spectrum_wave in enumerate(self.returned_wave):
            original_wave = [self.spec1_wave, self.spec2_wave][i]
            self.assertListEqual(original_wave, list(spectrum_wave))

    def test_correct_returned_flux(self):
        """Test returned fluxes match fluxes from the input table"""

        for i, spectrum_flux in enumerate(self.returned_flux):
            original_flux = [self.spec1_flux, self.spec2_flux][i]
            self.assertListEqual(original_flux, list(spectrum_flux))


class CalcModelChisq(TestCase):
    """Tests for the ``calc_model_chisq`` function"""

    # We define the data, model, and result properties to ensure a consistent,
    # un-mutated set of arguments when testing

    @property
    def model(self):
        """The model to use for testing"""

        model = sncosmo.Model('salt2')
        model.update(self.data.meta)
        return model

    @property
    def data(self):
        """The data to use for testing"""

        # Type cast the 'band' column to U15 so we can use longer band names
        data = sncosmo.load_example_data()
        data['band'] = Column(data['band'], dtype='U15')
        return data

    @property
    def results(self):
        """The results of fitting ``self.data`` with ``self.model``"""

        vparams = ['t0', 'x0', 'x1', 'c']
        results, fitted_model = sncosmo.fit_lc(self.data, self.model, vparams)
        return results

    # noinspection PyTypeChecker
    def test_arg_mutation(self):
        """Test the model, data table, and results object are not mutated"""

        model = self.model
        data = self.data
        results = self.results

        utils.calc_model_chisq(data, results, model)
        self.assertTrue(all(self.data == data), 'Data table was mutated')
        self.assertListEqual(
            list(self.results.parameters),
            list(results.parameters),
            'Results were mutated')

        self.assertListEqual(
            list(self.model.parameters),
            list(model.parameters),
            'Model parameters were mutated')

    def test_out_of_range_phase(self):
        """Test values outside the model's phase range are dropped"""

        data = self.data
        expected_chisq, expected_dof = utils.calc_model_chisq(
            data, self.results, self.model)

        # Add values that are out of the model's phase range
        # Columns: 'flux', 'fluxerr', 'zp', and 'zpsys' set as dummy values
        high_phase = self.model.mintime() - 1
        low_phase = self.model.maxtime() + 1
        data.add_row([high_phase, 'sdssu', 0, 0, 25, 'ab'])
        data.add_row([low_phase, 'sdssu', 0, 0, 25, 'ab'])

        chisq, dof = utils.calc_model_chisq(data, self.results, self.model)
        self.assertEqual(expected_chisq, chisq, 'Chisq values are not equal')
        self.assertEqual(expected_dof, dof, 'Degrees of freedom are not equal')

    def test_out_of_range_wave(self):
        """Test values outside the model's wavelength range are dropped"""

        data = self.data
        expected_chisq, expected_dof = utils.calc_model_chisq(
            data, self.results, self.model)

        # Define wavelength / transmission for bands above and below the
        # model's wavelength range
        min_wave = self.model.minwave()
        max_wave = self.model.maxwave()
        high_out_of_range = np.arange(min_wave, max_wave + 2)
        low_out_of_range = np.arange(min_wave - 2, max_wave)
        transmission = np.ones_like(low_out_of_range)

        # Register the bands with sncosmo
        high_band = sncosmo.Bandpass(high_out_of_range, transmission,
                                     name='high')
        low_band = sncosmo.Bandpass(low_out_of_range, transmission, name='low')
        sncosmo.register(high_band, name='high')
        sncosmo.register(low_band, name='low')

        # Add values that are out of the model's wavelength range
        data.add_row([1, 'high', 0, 0, 25, 'ab'])
        data.add_row([1, 'low', 0, 0, 25, 'ab'])

        chisq, dof = utils.calc_model_chisq(data, self.results, self.model)
        self.assertEqual(expected_chisq, chisq, 'Chisq values are not equal')
        self.assertEqual(expected_dof, dof, 'Degrees of freedom are not equal')

    def test_unregistered_bands(self):
        """Test unregistered band names in the data table cause an error
        instead of being dropped.
        """

        # Add values with an unregistered filter
        data = self.data
        data['band'][0] = 'made up band'

        func = utils.calc_model_chisq
        args = (data, self.model, self.results)
        self.assertRaises(Exception, func, *args)

    def test_empty_table(self):
        """Test an error is raised the data table is empty"""

        col_names = self.data.colnames
        empty_table = Table(names=col_names)
        args = empty_table, self.results, self.model
        self.assertRaises(ValueError, utils.calc_model_chisq, *args)

    def test_correct_chisq_calculation(self):
        """Test the correct chisq and dof are returned for simulated data"""

        # Define necessary information to simulate a table of flux data
        model = self.model
        t0 = model.parameters[1]
        zp_system = 'ab'
        zero_point = 25
        band_name = 'sdssg'
        flux_offset = 100  # Difference between "observed" and model flux
        flux_err_coeff = .1  # flux error / flux

        # Create table of simulated flux
        phase = np.arange(-10, 30) + t0
        model_flux = model.bandflux(band_name, phase, zero_point, zp_system)
        flux = model_flux + flux_offset
        flux_err = flux * flux_err_coeff
        band = np.full(len(phase), band_name)
        zpsys = np.full(len(phase), zp_system)
        zp = np.full(len(phase), zero_point)

        data = Table(
            [phase, band, flux, flux_err, zp, zpsys],
            names=['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']
        )

        result = sncosmo.utils.Result(
            {'vparam_names': ['t0', 'x0', 'x1', 'c']})

        expected_chisq = np.sum(((flux - model_flux) / flux_err) ** 2)
        chisq, dof = utils.calc_model_chisq(data, result, model)
        self.assertEqual(expected_chisq, chisq)
