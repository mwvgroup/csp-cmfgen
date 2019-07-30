from unittest import TestCase

import numpy as np

from analysis import equivalent_width


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
        recovered_peak = equivalent_width._calc_ew.get_peak_coordinates(
            self.wavelength,
            self.flux,
            expected_peak - 10,
            expected_peak + 10
        )

        self.assertEqual(expected_peak, recovered_peak)

    def test_unobserved_feature(self):
        """Test an error is raise if the feature is out of bounds"""

        max_wavelength = max(self.wavelength)
        with self.assertRaises(equivalent_width.UnobservedFeature):
            equivalent_width._calc_ew.get_peak_coordinates(
                self.wavelength,
                self.flux,
                max_wavelength + 10,
                max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature is found """

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        recovered_lower_peak = equivalent_width._calc_ew.get_peak_coordinates(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        recovered_upper_peak = equivalent_width._calc_ew.get_peak_coordinates(
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

        feat_start, feat_end = equivalent_width.get_feature_bounds(
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
        cls.continuum_slope = 2
        cls.continuum_intercept = 1

        # Define the continuum
        cls.feat_wave = np.arange(-10, cls.feat_width + 10, 1)
        cls.feat_flux = cls.continuum_slope * cls.feat_wave + cls.continuum_intercept

        # Add an absorption feature from 0 to cls.feat_width
        feature_indices = (0 <= cls.feat_wave) & (cls.feat_wave <= cls.feat_width)
        cls.feat_flux[feature_indices] = 2 * cls.feat_wave[feature_indices]

    def test_continuum_func(self):
        """Test fitting the continuum returns correct slope and intercept"""

        cont_func = equivalent_width.fit_continuum_func(
            self.feat_wave, self.feat_flux, 0, self.feat_width)

        fit_cont_params = np.polyfit(self.feat_wave, cont_func(self.feat_wave))
        self.assertEqual(
            self.continuum_slope, fit_cont_params[0], 'Wrong continuum slope.')

        self.assertEqual(
            self.continuum_intercept, fit_cont_params[1],
            'Wrong continuum y-intecept.')

    def test_calc_ew(self):
        """Test the measured and simulated pew agree"""

        pew, feat_start, feat_end = equivalent_width.calc_pew(
            self.feat_wave, self.feat_flux, feat_start=0,
            feat_end=self.feat_width)

        self.assertEqual(self.feat_width, pew)
        self.assertEqual(0, feat_start)
        self.assertEqual(self.feat_width, feat_end)


class Tabulation(TestCase):

    def test_create_empty_pew_table(self):
        pass

    def test_shift_models_to_data(self):
        pass

    def test_boundary_fixing(self):
        pass
