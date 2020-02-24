#!/usr/bin/env python
# coding: utf-8

"""The ``lc_colors`` module fits tabulates the chi-squared values for colors modeled by
an SNCosmo model. It also tabulates modeled and observed delta color over 15
days.

Function Documentation
----------------------
"""

from copy import deepcopy

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from . import utils
from .regression import fit_gaussian_process
from .regression import predict_c_15, predict_color


def create_empty_output_table(color_bands, suffixes):
    """Create an empty astropy table for storing color dependant values

    Returns an empty table with columns for each color defined by color_bands.
    Additional columns can be added with the same names but with an added
    suffix by specifying the ``suffixes`` argument.

    Args:
        color_bands (list[Tuple]): List of tuples with two band names
        suffixes      (list[str]): Include additional columns with suffixes

    Returns:
        An empty, masked astropy table
    """

    names, dtype = ['obj_id', 'source', 'version'], ['U100', 'U100', 'U100']
    for b1, b2 in color_bands:
        col_name = f"{b1.split('_')[-1]}_{b2.split('_')[-1]}"
        names += [col_name] + [col_name + s for s in suffixes]
        dtype += [float for _ in range(len(suffixes) + 1)]

    return Table(names=names, dtype=dtype, masked=True)


def get_observed_color_times(data, band1, band2):
    """Return the time range for which observations overlap in two band passes

    Args:
        data  (Table): Astropy table with columns 'band' and 'time'
        band1   (str): The name of a bandpass in data['band']
        band2   (str): The name of a bandpass in data['band']

    Returns:
        Tuple with the start and end times
    """

    for band_name in (band1, band2):
        if band_name not in data['band']:
            raise ValueError(f'No observations for band {band_name}')

    band1_times = data[data['band'] == band1]['time']
    band2_times = data[data['band'] == band2]['time']

    # Check that bands have overlapping observations
    b1_start, b1_end = min(band1_times), max(band1_times)
    b2_start, b2_end = min(band2_times), max(band2_times)
    if (b1_start > b2_end) or (b2_start > b1_end):
        raise ValueError(f'Bands do not overlap: {band1} : {band2}')

    return max(b1_start, b2_start), min(b1_end, b2_end)


def calc_color_chisq(data_table, model, band1, band2, interval=1, prange=None):
    """Calculate the chi-squared for observed and modeled color

    Assumes time values in the data and model are relative to the same zero
    point.

    Args:
        data_table (Table): Data table from sndata (format_sncosmo=True)
        model      (Model): An sncosmo model
        band1        (str):
        band2        (str):
        interval     (int): Spacing between phases when summing chisq
        prange     (tuple): Optional start and end phase for color evolution

    Returns:
        - The chi-square value
        - The degrees of freedom
    """

    gp = fit_gaussian_process(data_table)

    try:
        obs_start, obs_end = get_observed_color_times(data_table, band1, band2)

    except ValueError:  # No observations for given bands
        return np.nan, np.nan

    # Determine integration bounds
    phase_start, phase_end = obs_start, obs_end
    if prange is not None:
        phase_start, phase_end = prange[0], prange[-1]
        if (phase_start < obs_start) or (obs_end < phase_end):
            return np.nan, np.nan

    phase_range = np.arange(phase_start, phase_end + interval, interval)
    data_color, data_err = predict_color(gp, phase_range, band1, band2)
    model_color = model.color(band1, band2, 'ab', phase_range)
    chisq = np.sum((data_color - model_color) / data_err)
    return chisq, len(phase_range)


def tabulate_chisq(
        data_release, models, colors, interval=1, prange=None, out_path=None):
    """Tabulate color chi-squared for multiple models

    Integrate the chi-squared of the color evolution over the phase range
    ``prange``. The phase range should be specified relative to B-band max.
    If ``prange`` is not specified, use the largest integration range allowed
    by the data on a color by color basis.

    Args:
        data_release (module): An sndata data release
        models         (list): A list of sncosmo models
        colors         (list): List of tuples with band names
        interval        (int): Spacing between phases when summing chisq
        prange        (tuple): Optional start and end phases
        out_path        (str): Optionally write to path with each iteration
    """

    out_table = create_empty_output_table(colors, suffixes=['_dof'])
    out_table.meta['prange'] = prange
    out_table.meta['bands'] = colors

    for model in tqdm(models, desc='Models'):
        table_iterator = data_release.iter_data(
            verbose={'desc': 'Targets', 'position': 1},  # Put pbar on 2nd line
            filter_func=utils.filter_has_csp_data)

        for data_table in table_iterator:
            # Set model values so time values are relative to t0
            model = deepcopy(model)
            obj_id = data_table.meta['obj_id']
            z = data_table.meta['z']
            t0 = utils.get_csp_t0(obj_id)
            ebv = utils.get_csp_ebv(obj_id)
            model.set(t0=t0, z=z, mwebv=ebv)

            # Shift observed times to be relative to t0
            data_table = deepcopy(data_table)
            data_table['time'] = data_table['time'] - t0

            # Populate output table
            new_row = [
                data_table.meta['obj_id'],
                model.source.name,
                model.source.version
            ]
            mask = [False, False, False]

            for band1, band2 in colors:
                chisq, dof = calc_color_chisq(
                    data_table, model, band1, band2, interval, prange)

                new_row += [chisq, dof]
                mask += np.isnan([chisq, dof]).tolist()

            out_table.add_row(new_row, mask=mask)

            if out_path:
                out_table.write(out_path, overwrite=True)

    return out_table


def tabulate_delta_15(
        data_release, models, band_combos, out_path=None, t0_band='csp_dr3_B'):
    """Tabulate observed and modeled delta color over 15 days

    Args:
        data_release (module): An sndata data release
        models         (list): A list of sncosmo models
        band_combos    (list): List of tuples with band names
        out_path        (str): Optionally write to path with each iteration
        t0_band         (str): Band to use when setting model t0 to peak

    Returns:
        An astropy table
    """

    out_table = create_empty_output_table(band_combos, suffixes=['_err'])
    data_iter = data_release.iter_data(
        verbose={'desc': 'Targets', 'position': 0},
        filter_func=utils.filter_has_csp_data)

    for data_table in data_iter:
        # Convert observation times to phase values
        obj_id = data_table.meta['obj_id']
        t0 = utils.get_csp_t0(obj_id)
        data_table['time'] = data_table['time'] - t0
        gp = fit_gaussian_process(data_table)

        new_row = [obj_id, 'CSP', 'DR3']
        mask = [False for _ in new_row]
        for band1, band2 in band_combos:
            try:
                start_phase, end_phase = get_observed_color_times(
                    data_table, band1, band2)

                if start_phase > 0 or end_phase < 15:
                    raise ValueError

            except ValueError:
                new_row += [np.NAN, np.NAN]
                mask += [True, True]

            else:
                new_row.extend(predict_c_15(gp, band1, band2))
                mask += [False, False]

        out_table.add_row(new_row, mask=mask)
        if out_path:
            out_table.write(out_path, overwrite=True)

    for model in models:
        new_row = ['model', model.source.name, model.source.version]
        mask = [False for _ in new_row]

        model = deepcopy(model)
        t0 = model.source.peakphase(t0_band)
        for band1, band2 in band_combos:
            c0 = model.color(band1, band2, 'ab', t0)
            c15 = model.color(band1, band2, 'ab', t0 + 15)
            new_row += [c15 - c0, np.NAN]
            mask += [False, True]

        out_table.add_row(new_row, mask=mask)
        if out_path:
            out_table.write(out_path, overwrite=True)

    return out_table
