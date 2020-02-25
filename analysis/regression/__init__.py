#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``regression`` module fits light curves using a gaussian regression and
compares the fitted colors with modeled color values.

Usage Example
-------------

.. code-block:: python
   :linenos:

   from sndata.csp import DR3

   from analysis import regression

   # Download demo data
   dr3 = DR3()
   dr3.download_module_data()

   # Load demo data as an astropy table
   demo_id = dr3.get_available_ids()[0]
   demo_data = dr3.get_data_for_id(demo_id)
   print(demo_data)

   # Fit a Gaussian process to the data
   gp = lc_colors.fit_gaussian_process(demo_data)
   gp_flux, gp_unc = regression.predict_light_curve(gp, bands, time)

   # Pick the bands and time values to evaluate the Gaussian process at
   all_bands = list(set(demo_data['band']))
   time = np.arange(min(demo_data['time']) - 1, max(demo_data['time']) + 1)

   # Interpolate flux in a single band
   band_flux, band_err = regression.predict_band_flux(gp, all_bands[0], time)

   # Interpolate flux in multiple bands
   lc_flux, lc_err = regression.predict_light_curve(gp, all_bands, time)

   # Predict the color at various time values
   band1 = all_bands[0]
   band2 = all_bands[1]
   color, color_err = predict_color(gp, time, band1, band2)

Function Documentation
----------------------
"""

import numpy as np
from scipy.optimize import curve_fit
from sndata import get_zp

from ._lc_regression import fit_gaussian_process
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
        lc = np.array([predict_band_flux(gp, b, t) for b, t in zip(bands, times)])

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


def neg_exponential(x, a, b, c):
    """Negative exponential

    Args:
        x: The values to evaluate the exponential for
        a: The scale factor of the exponential
        b: The average of the exponential
        c: The offset of the exponential

    Returns:
        -a * np.exp(-b * x) + c
    """

    return -a * np.exp(-b * x) + c


def inverse_exp(x, a, b, c):
    """Inverse of a negative exponential

    Args:
        x: The values to evaluate the exponential for
        a: The scale factor of the exponential
        b: The average of the exponential
        c: The offset of the exponential

    Returns:
        np.log(a / (c - x)) / b
    """

    return np.log(a / (c - x)) / b


def scale_errors(obs_err, reg_err, maxfev=10000):
    """Exponentially scale regression errors to match observed errors

    Fits a negative exponential (``neg_exponential``) to the observed error
    and then uses the optimized parameters to apply the inverse of a negative
    exponential (``inverse_exp``) to the regression errors.

    Args:
        obs_err (ndarray): The observational errors
        reg_err (ndarray): The regression errors
        maxfev      (int): Max number of iterations when fitting

    Returns:
        - The scaled errors
        - The optimized parameters
        - The covariance of the optimized parameters
    """
    exp_opt, exp_cov = curve_fit(neg_exponential, obs_err, reg_err, maxfev=maxfev)
    scaled_errors = inverse_exp(reg_err, *exp_opt)
    return scaled_errors, exp_opt, exp_cov
