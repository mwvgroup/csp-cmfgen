#!/usr/bin/env python
# coding: utf-8

"""This module fits determines the chi-squared values for colors modeled by
an SNCosmo model.
"""

import warnings

import numpy as np
from astropy.table import Table
from scipy import integrate
from sndata.csp import dr1, dr3

from ._lc_regression import fit_gaussian_process, predict_color
from ..utils import make_pbar

with warnings.catch_warnings():
    warnings.simplefilter("ignore")


def get_csp_t0(obj_id):
    """Get the t0 value published by CSP DR3 for a given object

    Args:
        obj_id (str): The object Id value

    Return:
        The published MJD of maximum minus 53000
    """

    params = dr3.load_table(3)
    if obj_id not in params['SN']:
        raise ValueError(f'No published t0 for {obj_id}')

    # Subtract 53000 to convert to zero point in DR3 photometry data tables
    return params[params['SN'] == obj_id]['T(Bmax)'][0] - 53000


def get_csp_ebv(obj_id):
    """Get the E(B - V) value published by CSP DR1 for a given object

    Args:
        obj_id (str): The object Id value

    Return:
        The published E(B - V) value
    """

    extinction_table = dr1.load_table(1)
    if obj_id not in extinction_table['SN']:
        raise ValueError(f'No published E(B-V) for {obj_id}')

    data_for_target = extinction_table[extinction_table['SN'] == obj_id]
    return data_for_target['E(B-V)'][0]


def get_color_times(data, band1, band2):
    """Return the time range for which observations were performed
     in a pair of band passes

    Args:
        data  (Table): Data returned by sndata
        band1   (str): The name of a bandpass in data['band']
        band1   (str): The name of a bandpass in data['band']

    Returns:
        The start time
        The end time
    """

    band1_times = data[data['band'] == band1]['time']
    band2_times = data[data['band'] == band2]['time']
    if (len(band1_times) == 0) or (len(band2_times) == 0):
        raise ValueError(f'Bands do not overlap: {band1} : {band2}')

    start = max(min(band1_times), min(band2_times))
    end = min(max(band1_times), max(band2_times))
    return start, end


def calculate_color_residual(gp, model, band1, band2, tstart, tend):
    """Calculate the integrated residual per unit wavelength for color

    Args:
        gp        (GP): A fitted gaussian process
        model  (Model): An sncosmo model
        band1    (str): Name of the first band in the magnitude difference
        band2    (str): Name of the second band in the magnitude difference
        tstart (float): The start time of the integration
        tend   (float): The end time of the integration

    Returns:
        (The integrated residual) / (tend - tstart)
        (The error in integrated residual) / (tend - tstart)
    """

    def residual(time):
        data_color, data_color_error = predict_color(gp, band1, band2, [time])
        model_color = model.color(band1, band2, 'ab', time)
        return (data_color - model_color) / data_color_error

    return integrate.quad(residual, tstart, tend)


def create_empty_output_table(band_combos):
    """Create an empty astropy table for storing residuals

    Args:
        band_combos (list): List of tuples with band names to determine color for

    Returns:
        An astropy table
    """

    def make_col_name(b1, b2):
        return f"{b1.replace('csp_dr3_', '')}_{b2.replace('csp_dr3_', '')}"

    names, dtype = ['obj_id'], ['U100']
    for b1, b2 in band_combos:
        col_name = make_col_name(b1, b2)
        names += [col_name, col_name + '_err']
        dtype += [float, float]

    return Table(names=names, dtype=dtype)


def tabulate_residuals(data_release, model, band_combos, verbose=True):
    """Tabulate color residuals for a given model

    Args:
        data_release (module): An sndata data release
        model         (Model): An sncosmo model
        band_combos    (list): List of tuples with band names
        verbose        (bool): Whether to display a progress bar
    """

    # Construct iterator over data tables for each SN in the given data release
    total_iters = len(data_release.get_available_ids())
    data_iterable = make_pbar(
        data_release.iter_data(format_sncosmo=True),
        verbose=verbose,
        desc='Targets',
        total=total_iters)

    out_table = create_empty_output_table(band_combos)
    for data_table in data_iterable:
        new_row = [data_table.meta['obj_id']]
        gp = fit_gaussian_process(data_table)

        for band1, band2 in make_pbar(band_combos, verbose, desc='Colors', position=1):
            try:
                tstart, tend = get_color_times(data_table, band1, band2)
                resid, resid_err = calculate_color_residual(
                    gp, model, band1, band2, tstart, tend)

            except ValueError:
                resid, resid_err = np.nan, np.nan

            new_row += [resid, resid_err]

        out_table.add_row(new_row)

    return Table(out_table)
