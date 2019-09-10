#!/usr/bin/env python
# coding: utf-8

"""This module fits tabulates the chi-squared values for colors modeled by
an SNCosmo model. It also tabulates modeled and observed delta color over 15
days.
"""

from copy import deepcopy

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from ._lc_prediction import predict_c_15, predict_color
from ._lc_regression import fit_gaussian_process
from .. import utils


def get_color_times(data, band1, band2):
    """Return the time range for which observations overlap in two band passes

    Args:
        data  (Table): Astropy table with columns 'band' and 'time'
        band1   (str): The name of a bandpass in data['band']
        band2   (str): The name of a bandpass in data['band']

    Returns:
        The start time
        The end time
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


def create_empty_output_table(band_combos, suffix):
    """Create an empty astropy table for storing chi-squareds for a given model

    Args:
        band_combos (list[Tuple]): List of tuples with band names

    Returns:
        An astropy table
    """

    names, dtype = ['obj_id', 'source', 'version'], ['U100', 'U100', 'U100']
    for b1, b2 in band_combos:
        col_name = f"{b1.split('_')[-1]}_{b2.split('_')[-1]}"
        names += [col_name, col_name + suffix]
        dtype += [float, float]

    return Table(names=names, dtype=dtype, masked=True)


def create_new_chisq_row(data_table, model, band_combos, interval, prange):
    """Fit a light curve and return a list of chi-squared results

    Results are returned as a list formatted to be a new row for a table
    returned by ``create_empty_output_table``.

    Assumes time values in the data and model are relative to the same zero
    point.

    Args:
        data_table (Table): Data table from sndata (format_sncosmo=True)
        model      (Model): An sncosmo model
        band_combos (list): List of tuples with band names
        interval     (int): Spacing between phases when summing chisq
        prange     (tuple): Optional start and end phase for color evolution

    Returns:
        A list with object id, model info, chi-squares, and integration errors
        A mask for the first return
    """

    gp = fit_gaussian_process(data_table)
    obj_id = data_table.meta['obj_id']
    new_row = [obj_id, model.source.name, model.source.version]
    mask = [False for _ in new_row]

    for band1, band2 in band_combos:
        try:
            obs_start, obs_end = get_color_times(data_table, band1, band2)

        except ValueError:  # No observations for given bands
            new_row += [np.nan, np.nan]
            mask += [True, True]
            continue

        # Determine integration bounds
        phase_start, phase_end = obs_start, obs_end
        if prange is not None:
            phase_start, phase_end = prange[0], prange[-1]
            if (phase_start < obs_start) or (obs_end < phase_end):
                new_row += [np.nan, np.nan]
                mask += [True, True]
                continue

        phase_range = np.arange(phase_start, phase_end + interval, interval)
        data_color, data_err = predict_color(gp, phase_range, band1, band2)
        model_color = model.color(band1, band2, 'ab', phase_range)
        chisq = utils.chisq(data_color, data_err, model_color)

        new_row += [chisq, len(phase_range)]
        mask += [False, False]

    return new_row, mask


def tabulate_chisq(data_release, models, band_combos, interval=1,
                   prange=None, out_path=None):
    """Tabulate color chi-squared for multiple models

    Integrate the chi-squared of the color evolution over the phase range
    ``prange``. The phase range should be specified relative to B-band max.
    If ``prange`` is not specified, use the largest integration range allowed
    by the data on a color by color basis.

    Args:
        data_release (module): An sndata data release
        models         (list): A list of sncosmo models
        band_combos    (list): List of tuples with band names
        interval        (int): Spacing between phases when summing chisq
        prange        (tuple): Optional start and end phases
        out_path        (str): Optionally write to path with each iteration
    """

    out_table = create_empty_output_table(band_combos, suffix='_dof')
    out_table.meta['prange'] = prange
    out_table.meta['bands'] = band_combos

    for model in tqdm(models, desc='Models'):
        table_iterator = data_release.iter_data(
            format_sncosmo=True,
            verbose={'desc': 'Targets', 'position': 1},  # Put pbar on 2nd line
            filter_func=utils.filter_has_csp_data)

        for data_table in table_iterator:
            # Set model values so time values are relative to t0
            model = deepcopy(model)
            obj_id = data_table.meta['obj_id']
            z = data_table.meta['redshift']
            t0 = utils.get_csp_t0(obj_id)
            ebv = utils.get_csp_ebv(obj_id)
            model.set(t0=t0, z=z, extebv=ebv)

            # Shift observed times to be relative to t0
            data_table = deepcopy(data_table)
            data_table['time'] = data_table['time'] - t0

            # Populate output table
            new_row, mask = create_new_chisq_row(
                data_table, model, band_combos, interval, prange)

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

    out_table = create_empty_output_table(band_combos, suffix='_err')
    data_iter = data_release.iter_data(
        format_sncosmo=True, verbose={'desc': 'Targets'},
        filter_func=utils.filter_has_csp_data)

    for data_table in data_iter:
        # Convert observation times to phase values
        obj_id = data_table.meta['obj_id']
        t0 = utils.get_csp_t0(obj_id)
        data_table['time'] = data_table['time'] - t0
        gp = fit_gaussian_process(data_table)

        new_row = [obj_id, 'CSP', 'DR3']
        mask = [False for _ in new_row]
        for band1, band2 in tqdm(band_combos, desc='Colors', position=1):
            try:
                start_phase, end_phase = get_color_times(data_table, band1, band2)
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
