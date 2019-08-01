#!/usr/bin/env python
# coding: utf-8

"""This module fits determines the chi-squared values for colors modeled by
an SNCosmo model.
"""

from copy import deepcopy

import numpy as np
from astropy.table import Table
from scipy import integrate

from ._lc_regression import fit_gaussian_process, predict_color
from ..utils import get_csp_ebv, make_pbar


def get_color_times(data, band1, band2):
    """Return the time range for which observations were performed
     in a pair of band passes

    Args:
        data  (Table): Data returned by sndata
        band1   (str): The name of a bandpass in data['band']
        band2   (str): The name of a bandpass in data['band']

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
        data_color = predict_color(gp, band1, band2, [time])
        model_color = model.color(band1, band2, 'ab', time)
        return (data_color - model_color) / (tend - tstart)

    return integrate.quad(residual, tstart, tend)


def create_empty_output_table(model, band_combos):
    """Create an empty astropy table for storing model residuals

    Args:
        model      (Model): sncosmo model the table will be used for
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

    out_table = Table(names=names, dtype=dtype, masked=True)
    out_table.meta['source'] = model.source.name
    out_table.meta['version'] = model.source.version
    out_table.meta['bands'] = band_combos
    return out_table


def tabulate_residuals(data_release, model, band_combos, verbose=True):
    """Tabulate color residuals for a given model

    Args:
        data_release (module): An sndata data release
        model         (Model): An sncosmo model
        band_combos    (list): List of tuples with band names
        verbose        (bool): Whether to display a progress bar
    """
    model = deepcopy(model)

    # Construct iterator over data tables for each SN in the given data release
    total_iters = len(data_release.get_available_ids())
    data_pbar = make_pbar(
        data_release.iter_data(format_sncosmo=True),
        verbose=verbose,
        desc='Targets',
        total=total_iters)

    out_table = create_empty_output_table(model, band_combos)
    for data_table in data_pbar:
        gp = fit_gaussian_process(data_table)
        obj_id = data_table.meta['obj_id']
        new_row, mask = [obj_id], [False]

        band_pbar = make_pbar(band_combos, verbose, desc='Colors', position=1)
        for band1, band2 in band_pbar:
            try:
                tstart, tend = get_color_times(data_table, band1, band2)
                model.set(extebv=get_csp_ebv(data_table.meta['obj_id']))
                # noinspection PyTupleAssignmentBalance
                resid, resid_err = calculate_color_residual(
                    gp, model, band1, band2, tstart, tend)

                new_row += [resid, resid_err]
                mask += [False, False]

            except ValueError:
                new_row += [np.nan, np.nan]
                mask += [True, True]

        out_table.add_row(new_row, mask=mask)

    return Table(out_table)
