#!/usr/bin/env python
# coding: utf-8

"""This module computes photometric data (e.g. flux, color) from a gaussian
regression.
"""

import numpy as np
from sndata import get_zp

from ..utils import get_effective_wavelength


def predict_band_flux(gp, band_name, times):
    """Return the flux modeled by a gaussian regression in a single band

    Args:
        gp         (GP): A fitted gaussian process
        band_name (str): Name of band pass to return flux for
        times    (list): Times to predict flux for

    Returns:
        An array of flux values for the given times
        The errors in each flux value
    """

    effective_wavelength = get_effective_wavelength(band_name)
    wavelengths = np.ones(len(times)) * effective_wavelength
    predict_x_vals = np.vstack([times, wavelengths]).T
    return gp(predict_x_vals, return_var=True)


def predict_light_curve(gp, bands, times):
    """Return the flux modeled by a gaussian regression in multiple bands

    Times can either be a one dimensional list of times to be used for all
    bands or a two dimensional list specifying different times per band.

    Args:
        gp           (GP): A fitted gaussian process
        bands (list[str]): Name of band passes to return flux for
        times      (list): Times to predict flux for

    Returns:
        A 2d array of flux values for each band
        A 2d array of errors for the predicted fluxes
    """

    if np.ndim(times[0]) == 0:
        lc = np.array([predict_band_flux(gp, band, times) for band in bands])

    elif np.ndim(times[0]) == 1:
        lc = np.array(
            [predict_band_flux(gp, b, t) for b, t in zip(bands, times)])

    else:
        raise ValueError('Times must be a one or two dimensional list')

    return lc[:, 0], lc[:, 1]


def predict_color(gp, time, band1, band2):
    """Return the color value modeled by a gaussian regression

    Returns the band1 - band2 color. Assumes fluxes returned by the
    Gaussian regression are measured relative to the same zero point.

    Args:
        gp     (GP): A fitted gaussian process
        time (list): A 2d array of times for each band combination
        band1 (str): Name of the first band in the magnitude difference
        band2 (str): Name of the second band in the magnitude difference

    Returns:
        The predicted color
        The error in the predicted color
    """

    zp1 = get_zp(band1)
    zp2 = get_zp(band2)

    band1_pred, band1_err = predict_band_flux(gp, band1, time)
    band2_pred, band2_err = predict_band_flux(gp, band2, time)
    color = -2.5 * (np.log10(band1_pred) - np.log10(band2_pred)) + zp1 - zp2

    band1_err_term = ((-2.5 * band1_err) / (np.log(10) * band1_pred)) ** 2
    band2_err_term = ((-2.5 * band2_err) / (np.log(10) * band2_pred)) ** 2
    error = np.sqrt(band1_err_term + band2_err_term)
    return color, error


def predict_c_15(gp, band1, band2, t0=0):
    """Return the change in color over 15 days

    Args:
        gp     (GP): A fitted gaussian process
        band1 (str): Name of the first band in the magnitude difference
        band2 (str): Name of the second band in the magnitude difference
        t0  (float): Time of maximum (Default: 0)

    Returns:
        The change in color over 15 days
        The error in the change in color
    """

    c15, err15 = predict_color(gp, [15 + t0], band1, band2)
    c0, err0 = predict_color(gp, [0 + t0], band1, band2)
    return c15 - c0, err15 + err0
