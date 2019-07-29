#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Calculate the chi-squared values for modeled spectra"""

import numpy as np

from .utils import chisq_sum


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
