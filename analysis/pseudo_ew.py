#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits light curves using a gaussian regression and compares
the fitted colors with colors predicted by CMFGEN."""

import numpy as _np
from astropy.table import Table as _Table

feature_table = _Table({
    'feature_name': ['pW1', 'pW2', 'pW3', 'pW4', 'pW5', 'pW6', 'pW7', 'pW8'],
    'feature_id': ['Ca ii H&K', 'Si ii λ4130', 'Mg ii, Fe ii', 'Fe ii, Si ii',
                   'S ii λ5449, λ5622', 'Si ii λ5972', 'Si ii λ6355',
                   'Ca ii IR triplet'],
    'lower_blue': [3500, 3900, 3900, 4500, 5150, 5550, 5800, 7500],
    'upper_blue': [3800, 4000, 4150, 4700, 5300, 5700, 6000, 8000],
    'lower_red': [3900, 4000, 4450, 5050, 5500, 5800, 6200, 8200],
    'upper_red': [4100, 4150, 4700, 5550, 5700, 6000, 6600, 8900]
})


class UnobservedFeature(Exception):
    pass


# noinspection PyTypeChecker, PyUnresolvedReferences
def get_peak_coordinates(wavelength, flux, lower_bound, upper_bound):
    """Return coordinates of the maximum flux within given wavelength bounds

    Args:
        wavelength  (array): An array of wavelength values
        flux        (array): An array of flux values
        lower_bound (float): Lower wavelength boundary
        upper_bound (float): Upper wavelength boundary

    Returns:
        The wavelength for the maximum flux value
        The maximum flux value
    """

    feature_indices = (lower_bound < wavelength) & (wavelength < upper_bound)
    if not any(feature_indices):
        raise UnobservedFeature('Feature not in spectral wavelength range.')

    feature_flux = flux[feature_indices]
    feature_wavelength = wavelength[feature_indices]

    peak_index = _np.argmax(feature_flux)
    return feature_wavelength[peak_index], feature_flux[peak_index]


def get_feature_coordinates(wavelength, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wavelength (array): An array of wavelength values
        flux       (array): An array of flux values
        feature      (row): A row from the CSP feature table

    Returns:
        A list with the blueward and redward wavelengths
        A list with the blueward and redward flux values
    """

    blue_wave, blue_flux = get_peak_coordinates(
        wavelength, flux, feature['lower_blue'], feature['upper_blue'])

    red_wave, red_flux = get_peak_coordinates(
        wavelength, flux, feature['lower_red'], feature['upper_red'])

    return [blue_wave, red_wave], [blue_flux, red_flux]


def get_continuum_func(blue_wave, red_wave, blue_flux, red_flux):
    """Fit the pseudo continuum for a given feature

    Args:
        blue_wave (float): Starting wavelength of the feature
        red_wave  (float): Flux value for ``blue_wave``
        blue_flux (float): Ending wavelength of the feature
        red_flux  (float): Flux value for ``red_flux``

    Return:
        A linear function fit to the flux peaks bounding a feature
    """

    # Fit points to a line
    m = (red_flux - blue_flux) / (red_wave - blue_wave)
    b = blue_flux - m * blue_wave
    return lambda wave: m * _np.array(wave) + b


# noinspection PyTypeChecker, PyUnresolvedReferences
def calc_pew(wavelength, flux, feat_start, feat_end, cont_func):
    """Calculate the pseudo equivalent width for a given feature

    Args:
        wavelength (array): An array of wavelength values
        flux       (array): An array of flux values
        feat_start (float): The starting wavelength of the feature
        feat_end   (float): The ending wavelength of the feature
        cont_func   (func): Vectorized function returning continuum
            flux for wavelength

    Returns:
       The pseudo equivilent width as a float
    """

    indices = (feat_start < wavelength) & (wavelength < feat_end)
    feature_wave = wavelength[indices]
    feature_flux = flux[indices]
    continuum_flux = cont_func(feature_wave)

    normalized_flux = feature_flux / continuum_flux
    return _np.trapz(normalized_flux, feature_wave)
