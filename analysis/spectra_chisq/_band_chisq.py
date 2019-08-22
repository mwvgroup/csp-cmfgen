#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``spectra_chisq`` module calculates the chi-squared values for modeled
spectra.

**This module is incomplete, and is intended in it's current form as a
template / outline for future work. It is missing color warping**
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from tqdm import tqdm

from .. import utils
from ..equivalent_width import get_feature_bounds
from ..exceptions import UnobservedFeature


def band_limits(band_name, trans_limit):
    """Return wavelength range where a band is above a given transmission

    Args:
        band_name     (str): Name of an sncosmo registered band
        trans_limit (float): The transmission limit

    Returns:
        The wavelengths, fluxes, and flux errors inside the bandpass
    """

    if band_name.lower() == 'all':
        return -np.inf, np.inf

    band = sncosmo.get_bandpass(band_name)
    transmission_limits = band.wave[band.trans >= trans_limit]
    return np.min(transmission_limits), np.max(transmission_limits)


def calc_chisq(wave, flux, flux_err, model_flux, start, end):
    """Calculate the chi-squared for a spectrum within a wavelength range

    Args:
        wave       (ndarray): An array of wavelengths
        flux       (ndarray): An array of flux values
        flux_err   (ndarray): An array of error values for ``flux``
        model_flux (ndarray): An array of model flux values
        start        (float): The starting wavelength for the band
        end          (float): The ending wavelength for the band

    Returns:
        A dictionary with the chi_squared value in each band
    """

    if start < np.min(wave) or np.max(wave) < end:
        raise UnobservedFeature

    indices = np.where((start < wave) & (wave < end))[0]
    chisq_arr = (flux[indices] - model_flux[indices]) / flux_err[indices]
    return np.sum(chisq_arr ** 2)


def create_empty_chisq_table(col_names):
    """Create an empty astropy table for storing chi-squared results

    Returned table has columns 'obj_id', 'source', 'version', 'time' and
    *col_names.

    Args:
        col_names (list): Additional columns to add

    Returns:
        An astropy table
    """

    names = ['obj_id', 'source', 'version', 'time']
    dtype = ['U100', 'U100', 'U100', float]
    names.extend(col_names)
    dtype.extend((float for _ in col_names))

    out_table = Table(names=names, dtype=dtype, masked=True)
    return out_table


def create_new_table_row(obj_id, model, time, wave, flux, eflux, t0,
                         features=None, bands=None, trans_limit=None):
    """Calculate chi-squareds for a spectrum in multiple bands or features

    Results for formatted as a new row for tables returned by
    ``create_empty_chisq_table``. If ``bands`` and ``trans_limit`` are given
    then the chi-squared is determined per band. If ``features`` is given, then
    the chi-squared is determined per feature.

    Args:
        obj_id        (str): Id for the object being considered
        model       (Model): An sncosmo model
        time        (float): Observation time of the spectra
        wave        (array): Wavelengths for the spectrum
        flux        (array): Fluxes for the spectrum
        eflux       (array): Error in ``flux``
        t0          (float): Peak time of target
        features     (dict): Dictionary of features (optional)
        bands        (list): List of sncosmo registered band names (optional)
        trans_limit (float): Transmission limit for defining band wave range

    Returns:
        A list representing a new table row
        A mask for the new row
    """

    new_row = [obj_id, model.source.name, model.source.version, time]
    mask = [False, False, False, False]
    model_flux = model.flux(time - t0, wave)

    if bands and trans_limit:
        for band in bands:
            wave_start, wave_end = band_limits(band, trans_limit)

            try:
                chisq = calc_chisq(
                    wave, flux, eflux, model_flux, wave_start, wave_end)
                new_row.append(chisq)
                mask.append(False)

            except UnobservedFeature:
                new_row.append(np.NAN)
                mask.append(True)

    elif features:
        for feature in features:
            try:
                wave_start, wave_end = get_feature_bounds(wave, flux, feature)
                chisq = calc_chisq(wave, flux, eflux, model_flux, wave_start, wave_end)
                new_row.append(chisq)
                mask.append(False)

            except UnobservedFeature:
                new_row.append(np.NAN)
                mask.append(True)

    else:
        raise ValueError(
            'Must specify either ``features`` or ``bands`` and ``trans_limit``')

    return new_row, mask


def tabulate_chisq(data_release, models, err_estimate=.03, features=None,
                   bands=None, trans_limit=.1, out_path=None):
    """Tabulate band specific chi-squared values for spectroscopic observations

    If ``bands`` and ``trans_limit`` are given then the chi-squared is
    determined per band  over the wavelength range where the
    transmission is above ``trans_limit``. Specifying ``bands = 'all'``
    will calculate the chi-squared for the entire spectrum.

    If ``features`` is given, then the chi-squared is determined per feature.
    See ``analysis.equivalent_width.features`` for an example.

    Error in the flux is assumed to be a fraction of the observed flux.
    Defaults assume a 3% error in observed spectra.

    Args:
        data_release (module): An sndata data release
        models         (list): List of sncosmo models
        err_estimate  (float): Error estimate for spectra as fraction of flux
        bands          (list): A list of band names
        trans_limit   (float): Transmission limit for defining band wave range
        features       (dict): Dictionary of feature names and meta data
        out_path        (str): Optionally write results to file

    Returns:
        An astropy table of chi-squared values
    """

    out_table = create_empty_chisq_table(bands)
    data_iter = data_release.iter_data(
        verbose={'desc': 'Targets'}, filter_func=utils.filter_has_csp_data)

    for data_table in data_iter:
        obj_id = data_table.meta['obj_id']
        z = data_table.meta['redshift']
        ebv = utils.get_csp_ebv(obj_id)
        t0 = utils.get_csp_t0(obj_id)

        obs_time, wave, flux = utils.parse_spectra_table(data_table)
        flux_err = err_estimate * flux

        for model in tqdm(models, desc='Models', position=1):
            model = deepcopy(model)
            model.set(extebv=ebv, z=z)

            for t, w, f, fe in zip(obs_time, wave, flux, flux_err):
                new_row, mask = create_new_table_row(
                    obj_id, model, t, w, f, fe, t0,
                    features, bands, trans_limit)

                out_table.add_row(new_row, mask=mask)
                if out_table:
                    out_table.write(out_path, overwrite=True)

    return out_table
