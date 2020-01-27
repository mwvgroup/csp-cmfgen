#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``band_fitting`` module performs fits of individual bands of observed light-curves

Function Documentation
----------------------
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
import sncosmo
from astropy.table import Table
from matplotlib import pyplot
from tqdm import tqdm

from . import utils

DUST = sncosmo.F99Dust()


def create_empty_table(parameters, **kwargs):
    """Create an empty table for storing fit results

    Columns:
        - obj_id
        - band
        - source
        - pre_max
        - post_max
        - num_params
        - *parameters
        - *parameters + _err
        - chisq
        - ndof
        - b_max
        - delta_15
        - message

    Args:
        parameters (iter): List of parameter names to add columns for
        Any arguments to pass ``astropy.Table``

    Returns:
        A masked astropy Table
    """

    # Specify column names
    names = ['obj_id', 'band', 'source', 'pre_max', 'post_max', 'vparams']
    names += list(parameters) + [param + '_err' for param in parameters]
    names += ['chisq', 'ndof', 'b_max', 'delta_15', 'message']

    # Specify column data types
    dtype = ['U20', 'U100', 'U100', int, int, 'U100']
    dtype += [float for _ in range(2 * len(parameters))]
    dtype += [float, float, float, float, 'U10000']

    # Unless otherwise specified, we default to returning a masked table
    kwargs = deepcopy(kwargs)
    kwargs.setdefault('masked', True)
    return Table(names=names, dtype=dtype, **kwargs)


def fit_results_to_dict(data, obj_id, band_name, results, fitted_model):
    """Format sncosmo fit results so they can be appended to an astropy table

    See the ``create_empty_table`` function for information on the assumed
    table format.

    Args:
        data         (Table): The data used in the fit
        obj_id         (str): The id of the object that was fit
        band_name      (str): The name of the band that was fit
        results     (Result): Fitting results returned by ``sncosmo``
        fitted_model (Model): A fitted ``sncosmo`` model

    Returns:
        Fit results as a dictionary
    """

    new_row = {
        'obj_id': obj_id,
        'band': band_name,
        'source': fitted_model.source.name,
        'vparams': ','.join(results.vparam_names)
    }

    # Determine number of points pre and post maximum
    t0 = results.parameters[results.param_names.index('t0')]
    new_row['pre_max'] = sum(data['time'] < t0)
    new_row['post_max'] = sum(data['time'] >= t0)

    # Add parameters and their errors
    params = {p: v for p, v in zip(results.param_names, results.parameters)}
    new_row.update(params)
    for param, error in results.errors.items():
        new_row[param + '_err'] = error

    # Calc chi-squared
    chisq, ndof = utils.calc_model_chisq(data, results, fitted_model)
    new_row['chisq'] = np.round(chisq, 2)
    new_row['ndof'] = ndof

    # Determine peak magnitude and decline rate
    b_max = fitted_model.source_peakabsmag('bessellb', 'ab')
    new_row['b_max'] = np.round(b_max, 2)

    peak_phase = fitted_model.source.peakphase('bessellb')
    b_0 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase)
    b_15 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase + 15)
    delta_15 = b_15 - b_0
    new_row['delta_15'] = np.round(delta_15, 3)

    # Add fitting exit status message. Not all fitting routines include
    # this attribute, so we assign a default value of 'NONE'.
    message = getattr(results, 'message', 'NONE')
    new_row['message'] = message
    return new_row


def _plot_lc(data, result, fitted_model):
    """Plot fit results

    Args:
        data         (Table): The data used in the fit
        result      (Result): The fit results
        fitted_model (Model): Model with params set to fitted values
    """

    fig = sncosmo.plot_lc(data, fitted_model, errors=result.errors)
    xs, d = utils.calc_model_chisq(data, result, fitted_model)
    print(f'chisq / ndof = {xs} / {d} = {xs / d}', flush=True)
    return fig


def fit_single_target(
        fit_func, data, model, priors=None, kwargs=None,
        out_table=None, show_plots=False):
    """Run fits to individual bands of an observed light-curves

    Args:
        data      (Table): Table of photometric data
        fit_func   (func): Function to use to run fits (eg. ``sncosmo.fit_lc``)
        model     (Model): The model to use when fitting
        priors     (dict): Priors to use when fitting
        out_table    (Table): Append results to an existing table
        kwargs     (dict): Kwargs to pass ``fit_func`` when fitting
        show_plots (bool): Plot and display each individual fit

    Returns:
       A table with results each model / dataset combination
    """

    obj_id = data.meta['obj_id']
    model.update(priors)
    kwargs = deepcopy(kwargs)
    if out_table is None:
        out_table = create_empty_table(model.param_names)

    # Fit data in all bands
    all_result, all_fit = fit_func(data, model, **kwargs)
    all_row = fit_results_to_dict(data, obj_id, 'all', all_result, all_fit)
    out_table.add_row(all_row)

    if show_plots:
        _plot_lc(data, all_result, all_fit)
        pyplot.show()

    # Fix t0 and redshift during individual band fits
    kwargs['vparam_names'] = set(kwargs['vparam_names']) - {'t0', 'z'}

    # The amplitude from a fit to all data works as a better initial guess
    kwargs['guess_amplitude'] = False

    # Fit data in individual bands
    data = data.group_by('band')
    for band_name, band_data in zip(data.groups.keys['band'], data.groups):
        band_result, band_fit = fit_func(band_data, all_fit, **kwargs)
        band_row = fit_results_to_dict(
            band_data, obj_id, band_name, band_result, band_fit)

        out_table.add_row(band_row)

        if show_plots:
            _plot_lc(band_data, band_result, band_fit)

    return out_table


def _tabulate_fits_for_model(
        data_iter, model, config, fit_func, out_table, out_path=None):
    """Tabulate fit results for a collection of data tables and a single model

    Results are appended to the table specified by ``out_table``. Any objects
    with Ids already present in this table are skipped.

    Args:
        data_iter  (iter): Iterable of photometric data for different SN
        model     (Model): The model to use when fitting
        fit_func   (func): Function to use to run fits
        out_table (Table): Table to append results to
        out_path    (str): Optionally cache results to file in real time
    """

    for data in data_iter:
        # Get fitting priors and kwargs
        obj_id = data.meta['obj_id']
        if obj_id in out_table['obj_id']:
            continue

        try:
            fit_single_target(
                fit_func,
                data,
                model,
                priors=config[obj_id]['priors'],
                kwargs=config[obj_id]['kwargs'],
                out_table=out_table
            )

        except KeyboardInterrupt:
            raise

        except Exception as e:
            raise
            e_str = str(e).replace("\n", "")
            e_name = type(e).__name__
            out_table.add_row({
                'obj_id': obj_id,
                'message': f'{e_name}: {e_str}'
            })

        if out_path:
            out_table.write(out_path)


def tabulate_band_fits(
        data_release, models, fit_func, config=None, out_path=None):
    """Tabulate fit results for a collection of data tables

    Results already written to out_path are skipped.

    Args:
        data_release (module): An sndata data release
        models         (list): A list of sncosmo models
        fit_func       (func): Function to use to run fits
        config         (dict): Specifies priors / kwargs for fitting each model
        out_path        (str): Optionally cache results to file in real time

    Returns:
       An astropy table with fit results
    """

    # Set default kwargs
    config = deepcopy(config) or dict()

    # Add meta_data to output table meta data
    if Path(out_path).exists():
        out_table = Table.read(out_path)

    else:
        params = set()
        for m in models:
            params = params.union(m.param_names)

        out_table = create_empty_table(params)
        out_table.meta['fit_func'] = fit_func.__name__

    for model in tqdm(models, desc='Models'):
        data_iter = data_release.iter_data(
            verbose={'desc': 'Targets', 'position': 1},
            filter_func=utils.filter_has_csp_data)

        _tabulate_fits_for_model(
            data_iter, model, config, fit_func, out_table, out_path)

    return out_table
