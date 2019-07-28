#!/usr/bin/env python
# coding: utf-8

"""This module fits light curves using gaussian regression as outlined in
Boone et al. 2019. Some of the source code is ported from the avocado package,
which is described by the same paper.
"""

from functools import partial
from warnings import warn

import george
import numpy as np
import sncosmo
from astropy.table import Table
from george import kernels
from scipy.optimize import minimize


def get_effective_wavelength(band_name):
    """Get the effective wavelength for a given band

    Band name must be registered with SNCosmo

    Args:
        band_name (str): The name of a registered bandpass

    Returns:
        The effective wavelength in Angstroms
    """

    return sncosmo.get_bandpass(band_name).wave_eff


# Ported from avocado
def get_kernel(scale, fix_scale, length_scale):
    """Return a Matern 3/2 Kernel

    Args:
        scale        (float): Initial estimate of the scale
        fix_scale     (bool): Fix the scale to the initial estimate
        length_scale (float): Initial kernel length scale in days

    Returns
        A Matern32Kernel object
    """

    kernel = (
            (0.5 * scale) ** 2 *
            kernels.Matern32Kernel([length_scale ** 2, 6000 ** 2], ndim=2)
    )
    kernel.freeze_parameter('k2:metric:log_M_1_1')

    if fix_scale:
        kernel.freeze_parameter('k1:log_constant')

    return kernel


# Ported from avocado
def fit_gaussian_process(data, fix_scale=True, length_scale=20.):
    """Fit photometric observations with a gaussian regression

    Args:
        data         (Table): Data table from sndata (format_sncosmo=True)
        fix_scale     (bool): Whether to fix the scale while fitting (Default: True)
        length_scale (float): The initial length scale to use for the fit.

    Returns:
        A Gaussian process conditioned on the object's light-curve.
    """

    if isinstance(data, Table):
        data = data.to_pandas()

    fluxes = data['flux']
    flux_errors = data['fluxerr']
    wavelengths = data['band'].map(get_effective_wavelength)
    times = data['time']

    # Use the highest signal-to-noise observation to estimate the scale.
    # Include an error floor so that for very high
    # signal-to-noise observations we pick the maximum flux value.
    signal_to_noises = (
            np.abs(fluxes) /
            np.sqrt(flux_errors ** 2 + (1e-2 * np.max(fluxes)) ** 2)
    )
    scale = np.abs(fluxes[signal_to_noises.idxmax()])

    # Construct gaussian process
    kernel = get_kernel(scale, fix_scale, length_scale)
    gp = george.GP(kernel)

    # Compute the covariance matrix
    x_data = np.vstack([times, wavelengths]).T
    gp.compute(x_data, flux_errors)

    # Define log likelihood calculations
    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(fluxes)

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(fluxes)

    bounds = [(0, np.log(1000 ** 2))]
    if not fix_scale:
        bounds = [(-30, 30)] + bounds

    fit_result = minimize(
        neg_ln_like,
        gp.get_parameter_vector(),
        jac=grad_neg_ln_like,
        bounds=bounds,
    )

    if not fit_result.success:
        warn("GP fit failed! Using guessed GP parameters. ")

    gp.set_parameter_vector(fit_result.x)
    return partial(gp.predict, fluxes)


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

    Args:
        gp      (GP): A fitted gaussian process
        bands (list): Name of band passes to return flux for
        times (list): Times to predict flux for

    Returns:
        A 2d array of flux values for each band
        A 2d array of errors for the predicted flux
    """

    lc = np.array([predict_band_flux(gp, band, times) for band in bands])
    return lc[:, 0, :], lc[:, 1, :]


def predict_color(gp, band1, band2, time):
    """Return the color value modeled by a gaussian regression

    Returns the band1 - band2 color.

    Args:
        gp     (GP): A fitted gaussian process
        band1 (str): Name of the first band in the magnitude difference
        band2 (str): Name of the second band in the magnitude difference
        time (list): A 2d array of times for each band combination

    Returns:
        The predicted color
        The variance in the predicted color
    """

    band1_pred, band1_var = predict_band_flux(gp, band1, time)
    band2_pred, band2_var = predict_band_flux(gp, band2, time)
    color = -2.5 * (np.log10(band1_pred) - np.log10(band2_pred))
    error = (
            (2.5 / np.log(10)) ** 2 *
            (
                    (band1_var / band1_pred) ** 2 +
                    (band2_var / band2_pred) ** 2
            )
    )
    return color, error
