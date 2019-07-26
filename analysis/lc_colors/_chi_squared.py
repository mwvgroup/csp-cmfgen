#!/usr/bin/env python
# coding: utf-8

"""This module fits determines the chi-squared values for colors modeled by
CMFGEN.
"""

from copy import deepcopy

import numpy as np
from sndata.csp import dr1, dr3

from ._gauss_regression import fit_gaussian_process, predict_light_curve


def get_csp_t0(obj_id):
    """Get the t0 value published by CSP DR3 for a given object

    Args:
        obj_id (str): The object Id value

    Return:
        The published time of maximum minus 53000
    """

    params = dr3.load_table(3).to_pandas()
    params.set_index('SN', inplace=True)
    t0 = params.loc[obj_id]['T(Bmax)']
    if np.isnan(t0):
        raise ValueError(f'No published t0 for {obj_id}')

    return t0 - 53000  # Convert to CSP standard


def get_csp_ebv(obj_id):
    """Get the E(B - V) value published by CSP DR1 for a given object

    Args:
        obj_id (str): The object Id value

    Return:
        The published E(B - V) value
    """

    extinction_table = dr1.load_table(1)
    data_for_target = extinction_table[extinction_table['SN'] == obj_id]
    return data_for_target['E(B-V)'][0]


def get_color_times(data, band_combos, delta=1):
    """Create an array of times for which observations
    are available in a pair of band passes

    Args:
        data  (Table): Data returned by sncosmo
        band_combos (list[tuple[str]]): List of band pass pairs
        delta (float): Resolution of the returned array
    """

    out_times = []
    for (band1, band2) in band_combos:
        band1_times = data[data['band'] == band1]['time']
        band2_times = data[data['band'] == band2]['time']

        # Skip colors not observed for the given target
        if not (bool(len(band1_times)) and bool(len(band2_times))):
            out_times.append([])
            continue

        start = max(min(band1_times), min(band2_times))
        end = min(max(band1_times), max(band2_times))
        out_times.append(np.arange(start, end, delta))

    return out_times


def calculate_color_chisq(obj_id, model):
    data = dr3.get_data_for_id(obj_id, format_sncosmo=True)
    times = np.arange(min(data['time']), max(data['time']))
    times -= get_csp_t0(obj_id)
    ext = get_csp_ebv(obj_id)

    gp = fit_gaussian_process(data)
    reg_flux, reg_unc = predict_light_curve(gp, set(data['band']), times)

    model = deepcopy(model)
    model.set(extebv=ext)

    # Todo Calculate chi-squared