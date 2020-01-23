#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.lc_colors`` module."""

from unittest import TestCase

import numpy as np
from astropy.table import Table

from analysis import lc_colors


class CreateEmptyOutputTable(TestCase):
    """Tests for the ``create_empty_output_table`` function"""

    @classmethod
    def setUpClass(cls):
        # Arguments passed to the tested function
        cls.suffixes = ['_err']
        cls.color_bands = [
            ('survey_band1', 'survey_band2'),
            ('survey_band3', 'survey_band4')
        ]

        # Return of the tested function
        cls.test_table = lc_colors.create_empty_output_table(
            cls.color_bands, suffixes=cls.suffixes)

    def test_correct_column_names(self):
        """Test the returned table has the correct column names"""

        expected_cols = ['obj_id', 'source', 'version']
        for band_tuple in self.color_bands:
            col_name = '_'.join(b.split('_')[-1] for b in band_tuple)
            expected_cols += [col_name, col_name + '_err']

        self.assertSequenceEqual(self.test_table.colnames, expected_cols)

    def test_table_is_empty(self):
        """Test the returned table is empty"""

        self.assertEqual(0, len(self.test_table))

    def test_table_is_masked(self):
        """Test the returned table is a masked table"""

        self.assertTrue(self.test_table.masked)


class GetObservedColorTimes(TestCase):
    """Tests for the ``get_observed_color_times`` function """

    @classmethod
    def setUpClass(cls):
        """Create a dummy table describing observation times for two bands"""

        band1_times = np.arange(10, 21).tolist()
        band2_times = np.arange(15, 26, 5).tolist()
        combined_times = band1_times + band2_times

        # Record start and end values for the overlap of band1 and band2
        cls.simulated_overlap = 15, 20

        band1_names = np.full_like(band1_times, 'band1', dtype='U50').tolist()
        band2_names = np.full_like(band2_times, 'band2', dtype='U50').tolist()
        combined_names = band1_names + band2_names

        cls.test_table = Table(
            names=['time', 'band'],
            data=[combined_times, combined_names]
        )

    def test_missing_band(self):
        """Test a ValueError is raised when passed an unobserved bandpass"""

        args = (self.test_table, 'band1', 'nonexistent_band')
        self.assertRaises(
            ValueError, lc_colors.get_observed_color_times, *args)

    def test_non_overlapping_bands(self):
        """Test a ValueError is raised when passed an non-overlapping bands"""

        # Create a table where observations in band1 and band2 do not overlap
        new_test_table = self.test_table.copy()
        new_test_table['time'][new_test_table['band'] == 'band1'] *= 10

        args = (new_test_table, 'band1', 'band2')
        self.assertRaises(
            ValueError, lc_colors.get_observed_color_times, *args)

    def test_recovers_correct_times(self):
        """Test the returned start/end times match the simulated data"""

        returned_overlap = lc_colors.get_observed_color_times(
            self.test_table, 'band1', 'band2')

        self.assertEqual(
            self.simulated_overlap, returned_overlap,
            'Incorrect start or end time returned for band overlap')
