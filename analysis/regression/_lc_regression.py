#!/usr/bin/env python
# coding: utf-8

"""This module ports logic from the ``avocado`` package for fitting
light-curves with a gaussian regression. This module is meant as a place holder
for if / when ``avocado`` is upgraded to handle data other than the PLAsTiCC
data set. See Boone et al. 2019 for more details.
"""

from functools import partial
from warnings import warn

import george
import numpy as np
from astropy.table import Table
from george import kernels
from scipy.optimize import minimize

from ..utils import get_effective_wavelength


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
        fix_scale     (bool): Whether to fix the scale while fitting
        length_scale (float): The initial length scale to use for the fit.

    Returns:
        A Gaussian process conditioned on the object's light-curve.
    """

    if isinstance(data, Table):
        data = data.to_pandas()

    # Format data. Not included in original Avocado code
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

    guess_parameters = gp.get_parameter_vector()
    bounds = [(0, np.log(1000 ** 2))]
    if not fix_scale:
        bounds = (
                [(guess_parameters[0] - 10, guess_parameters[0] + 10)]
                + bounds
        )

    fit_result = minimize(
        neg_ln_like,
        gp.get_parameter_vector(),
        jac=grad_neg_ln_like,
        bounds=bounds
    )

    if not fit_result.success:
        warn("GP fit failed! Using guessed GP parameters. ")

    gp.set_parameter_vector(fit_result.x)

    # Wrap the output function to return standard deviation instead of variance
    # Not included in original Avocado code
    def out_func(*args, **kwargs):
        p, pv = partial(gp.predict, fluxes)(*args, **kwargs)
        return p, np.sqrt(pv)

    return out_func
