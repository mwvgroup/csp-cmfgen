#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``classification`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from sncosmo.utils import Result

from analysis import band_fitting


class CreateEmptyTable(TestCase):
    """Tests for the ``create_empty_table`` function"""

    def test_is_empty(self):
        """Test the returned table is empty by default"""

        num_rows = len(band_fitting.create_empty_table([]))
        self.assertEqual(0, num_rows, 'Table is not empty')

    def test_correct_columns(self):
        """Test the returned table has the correct columns"""

        parameters = ['z', 't0', 'x0', 'x1', 'c']
        returned = band_fitting.create_empty_table(parameters).colnames
        expected = [
            'obj_id', 'band', 'source', 'pre_max', 'post_max', 'vparams',
            'z', 't0', 'x0', 'x1', 'c',
            'z_err', 't0_err', 'x0_err', 'x1_err', 'c_err',
            'chisq', 'ndof', 'b_max', 'delta_15', 'message', ]

        self.assertSequenceEqual(expected, returned)

    def test_is_masked(self):
        """Test the returned table is a masked Table by default"""

        is_masked = band_fitting.create_empty_table([]).masked
        self.assertTrue(is_masked, 'Table has no mask')

    def test_accepts_kwargs(self):
        """Test that the table constructor uses the kwargs"""

        num_columns = len(band_fitting.create_empty_table([]).colnames)
        dummy_row = np.ones(num_columns)
        table = band_fitting.create_empty_table([], rows=[dummy_row])
        self.assertEqual(1, len(table), 'Table did not add specified rows')


class FitResultsToDict(TestCase):
    """Tests for the ``fit_results_to_dict`` function"""

    # Don't limit output messages on test failures
    maxDiff = None

    def runTest(self):
        """Run fits for data with a known result and check the results are
        correctly formatted as a table row.
        """

        data = sncosmo.load_example_data()
        params = ['z', 't0', 'x0', 'x1', 'c']
        param_values = [data.meta[p] for p in params]

        # Define mock fit results for a model fit for the last four parameters
        # I.e. for a fit where z is specified
        vparams = params[1:]
        result = Result({
            'param_names': params,
            'vparam_names': vparams,
            'parameters': param_values,
            'errors': {p: .1 * v for p, v in data.meta.items()}
        })

        model = sncosmo.Model('salt2')
        model.update(data.meta)
        data.meta['obj_id'] = 'dummy_id'

        row = band_fitting.fit_results_to_dict(
            data, 'dummy_id', 'dummy_band', result, model)

        expected_row = {
            'obj_id': 'dummy_id',
            'band': 'dummy_band',
            'source': 'salt2',
            'pre_max': 15,
            'post_max': 25,
            'vparams': ','.join(vparams),
            'z': result.parameters[0],
            't0': result.parameters[1],
            'x0': result.parameters[2],
            'x1': result.parameters[3],
            'c': result.parameters[4],
            'z_err': result.errors['z'],
            't0_err': result.errors['t0'],
            'x0_err': result.errors['x0'],
            'x1_err': result.errors['x1'],
            'c_err': result.errors['c'],
            'chisq': 36.44,
            'ndof': len(data) - len(result.vparam_names),
            'b_max': -19.5,
            'delta_15': 0.953,
            'message': 'NONE'
        }

        self.assertEqual(expected_row, row)


class PlotLc(TestCase):
    """Tests for the ``_plot_lc`` function"""

    def runTest(self):
        """Test the returned figure is not empty"""

        data = sncosmo.load_example_data()
        model = sncosmo.Model(source='salt2')
        result, fitted_model = sncosmo.fit_lc(
            data, model,
            ['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z': (0.3, 0.7)})

        fig = band_fitting._plot_lc(data, result, fitted_model)
        self.assertTrue(fig.get_axes(), 'Returned figure has no axes')


class RunBandFits(TestCase):
    """Tests for the ``fit_single_target`` function"""

    @classmethod
    def setUpClass(cls):
        """Run a set of fits on example data from sncosmo"""

        cls.data = sncosmo.load_example_data()
        cls.data.meta['obj_id'] = 'dummy_id'

        model = sncosmo.Model('salt2')
        priors = {p: cls.data.meta[p] for p in model.param_names}

        cls.fit_results = band_fitting.fit_single_target(
            fit_func=sncosmo.fit_lc,
            data=cls.data,
            model=model,
            priors=priors,
            kwargs=dict(vparam_names=['t0', 'x0', 'x1', 'c'])
        )

    def test_all_bands_were_fit(self):
        """Test fits were performed for all bands and each individual band"""

        expected_fits = ['all'] + sorted(set(self.data['band']))
        actual_fits = sorted(set(self.fit_results['band']))
        self.assertListEqual(expected_fits, actual_fits)
