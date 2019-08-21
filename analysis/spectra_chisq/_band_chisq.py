#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``spectra_chisq`` module calculates the chi-squared values for modeled
spectra.

**This module is incomplete, and is intended in it's current form as a
template / outline for future work.**
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from tqdm import tqdm

from .. import utils
from ..exceptions import UnobservedFeature


def band_limits(band_name, trans_limit):
    """Return wavelength range where a band is above a given transmission

    Args:
        wave      (ndarray): An array of wavelengths
        flux      (ndarray): An array of flux values
        flux_err  (ndarray): An array of error values for ``flux``
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


def chisq(wave, flux, flux_err, model_flux, start, end):
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


def create_empty_output_table(band_names):
    """Create an empty astropy table for storing chi-squared results

    Args:
        band_names  (list): List band names

    Returns:
        An astropy table
    """

    names = ['obj_id', 'source', 'version', 'time']
    dtype = ['U100', 'U100', 'U100', float]
    names.extend(band_names)
    dtype.extend((float for _ in band_names))

    out_table = Table(names=names, dtype=dtype, masked=True)
    return out_table


def tabulate_chisq(data_release, models, bands, err_estimate=.03,
                   trans_limit=.1, out_path=None):
    """Tabulate band specific chi-squared values for spectroscopic observations

    The chi-squared is calculated over the wavelength range where the
    transmission is above ``transmission_limit``. Specifying ``bands = 'all'``
    will calculate the chisquared for the entire spectrum.

    Defaults to a 3% error in observed spectra and a 10% transmission limit.

    Args:
        data_release (module): An sndata data release
        models         (list): List of sncosmo models
        bands          (list): A list of band names
        err_estimate  (float): Error estimate for spectra as fraction of flux
        trans_limit   (float): Transmission limit for defining band wave range
        out_path        (str): Optionally write results to file

    Returns:
        An astropy table of chi-squared values
    """

    out_table = create_empty_output_table(bands)
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
                new_row = [obj_id, model.source.name, model.source.version, t]
                mask = [False, False, False, False]
                model_flux = model.flux(t - t0, w)

                for band in bands:
                    band_start, band_end = band_limits(band, trans_limit)

                    try:
                        chisq = chisq(w, f, fe, model_flux, band_start,
                                      band_end)
                        new_row.append(chisq)
                        mask.append(False)

                    except UnobservedFeature:
                        new_row.append(np.NAN)
                        mask.append(True)

                out_table.add_row(new_row, mask=mask)
                if out_table:
                    out_table.write(out_path, overwrite=True)

    return out_table
