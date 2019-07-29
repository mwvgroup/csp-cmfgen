#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Calculate the chi-squared values for modeled spectra"""

from copy import deepcopy

import numpy as np
from astropy.table import Table
from scipy.ndimage.filters import gaussian_filter1d
from sndata.csp import dr1  # Todo: Switch to dr2

from .utils import chisq_sum, get_spectra_for_id, make_pbar


# Todo: Write the color_warp_spectrum function
def color_warp_spectrum(wave, flux, flux_err, **kwargs):
    """Color warp a spectrum

    Args:
        wave       (ndarray): An array of wavelengths
        flux       (ndarray): An array of flux values
        flux_err   (ndarray): An array of error values for ``flux``

    Returns:
        An array of color warped flux values
        An array of error values for the color warped fluxes
    """

    return flux, flux_err


def filter_spectrum(flux, width):
    """Reduce the resolution of a spectrum using a Gaussian filter

    Uses a filter with the given width and a sigma of 1.

    Args:
        flux (ndarray): An array of flux values
        width    (int): Width of the

    Returns:
        An array of flux values
        An array of error values for the fluxes
    """

    if not isinstance(width, int):
        raise TypeError('Width must be an int')

    # See https://stackoverflow.com/questions/25216382/gaussian-filter-in-scipy
    truncate = (width - 1) / 2
    flux_gauss = gaussian_filter1d(flux, 1, truncate=truncate)
    return flux_gauss


def band_chisq_spectrum(wave, flux, flux_err, model_flux, band_ranges):
    """Calculate the chi-squared in different band-passes for a single spectrum

    Args:
        wave       (ndarray): An array of wavelengths
        flux       (ndarray): An array of flux values
        flux_err   (ndarray): An array of error values for ``flux``
        model_flux (ndarray): An array of model flux values
        band_ranges   (dict): A dictionary of band pass limits.
             {<band name (str)>: (<min wave (float)>, <max wave (float)>)}

    Returns:
        A dictionary with the chi_squared value in each band
    """

    out_data = dict()
    for band, (min_wave, max_wave) in band_ranges.items():
        indices = (min_wave < wave) & (wave < max_wave)
        chi = chisq_sum(flux[indices], model_flux[indices], flux_err[indices])
        out_data[band] = np.array(chi)

    return out_data


def band_chisq_spectra(wave, flux, flux_err, model_flux, band_ranges):
    """Calculate the chi-squared in different band-passes for multiple spectra

    Args:
        wave       (ndarray): A 2d array of wavelengths
        flux       (ndarray): A 2d array of flux values
        flux_err   (ndarray): A 2d array of error values for ``flux``
        model_flux (ndarray): A 2d array of model flux values
        band_ranges   (dict): A dictionary of band pass limits.
             {<band name (str)>: (<min wave (float)>, <max wave (float)>)}

    Returns:
        A dictionary with an array of chi_squared values in each band
    """

    out_data = {key: [] for key in band_ranges}
    for w, f, fe, mf in zip(wave, flux, flux_err, model_flux):
        for band, (min_wave, max_wave) in band_ranges.items():
            indices = (min_wave < w) & (w < max_wave)
            chisq = chisq_sum(f[indices], mf[indices], fe[indices])
            out_data[band].append(chisq)

    return out_data


def create_empty_output_table(model, band_names):
    """Create an empty astropy table for storing chi-squared results

    Args:
        model      (Model): sncosmo model the table will be used for
        band_names  (list): List band names

    Returns:
        An astropy table
    """

    names, dtype = ['obj_id'], ['U100']
    names.extend(band_names)
    dtype.extend((float for _ in band_names))

    out_table = Table(names=names, dtype=dtype, masked=True)
    out_table.meta['source'] = model.source.name
    out_table.meta['version'] = model.source.version
    return out_table


def tabulate_chi_squared(data_release, model, band_ranges, verbose=True):
    """Tabulate color residuals for a given model

    Args:
        data_release (module): An sndata data release
        model         (Model): An sncosmo model
        band_ranges    (dict): A dictionary of band pass limits.
        verbose        (bool): Whether to display a progress bar
    """

    model = deepcopy(model)
    band_ranges = deepcopy(band_ranges)

    # Construct iterator over object ids
    id_pbar = make_pbar(
        data_release.get_available_ids(), verbose=verbose, desc='Targets')

    out_table = create_empty_output_table(model, band_ranges.keys())
    for obj_id in id_pbar:
        obs_data = zip(*get_spectra_for_id(dr1, obj_id))
        for time, wave, flux, flux_err in obs_data:
            obs_flux, obs_err = color_warp_spectrum(wave, flux, flux_err)

            # Todo: Define ebv and width
            model.set(ebv=ebv)
            model_flux = filter_spectrum(model.flux(time, wave), width)

            chi_sq = band_chisq_spectrum(
                wave, obs_flux, obs_err, model_flux, band_ranges)

            band_ranges['obj_id'] = obj_id
            out_table.add_row(band_ranges)

    return out_table
