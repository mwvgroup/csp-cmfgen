#!/usr/bin/env python
# coding: utf-8

"""This module fits light curves using gaussian regression as outlined in
Boone et al. 2019. The source code is ported from the avocado package, which
is described by the same paper.
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

    band = sncosmo.get_bandpass(band_name)
    return band.wave_eff


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
            kernels.Matern32Kernel(
                [length_scale ** 2, 6000 ** 2], ndim=2)
    )
    kernel.freeze_parameter('k2:metric:log_M_1_1')

    if fix_scale:
        kernel.freeze_parameter('k1:log_constant')

    return kernel


def fit_gaussian_process(data, fix_scale=True, length_scale=20.):
    """Fit a gaussian regressor to a data table

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


def predict_light_curve(gp, bands, times):
    """Predict the flux of a light curve

    Args:
        gp             (GP): A fitted gaussian process
        bands   (list[str]): Name of band passes to fit for
        times (list[float]): Times to predict flux for

    Returns:
        The flux in each band
        The errors in each band
    """
    predictions = []
    prediction_uncertainties = []

    for band in bands:
        effective_wavelength = get_effective_wavelength(band)
        wavelengths = np.ones(len(times)) * effective_wavelength

        pred_x_data = np.vstack([times, wavelengths]).T
        band_pred, band_pred_var = gp(pred_x_data, return_var=True)

        prediction_uncertainties.append(np.sqrt(band_pred_var))
        predictions.append(band_pred)

    predictions = np.array(predictions)
    prediction_uncertainties = np.array(prediction_uncertainties)
    return predictions, prediction_uncertainties
