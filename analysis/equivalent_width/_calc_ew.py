#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides generic functions for measuring the pseudo equivalent
width of arbitrary features for a one or more spectra.
"""

import numpy as np
from astropy.table import Table


class UnobservedFeature(Exception):
    pass


# noinspection PyTypeChecker, PyUnresolvedReferences
def _get_peak_coordinates(wavelength, flux, lower_bound, upper_bound):
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

    # Make sure the given spectrum spans the given wavelength bounds
    if (min(wavelength) >= lower_bound) or (upper_bound >= max(wavelength)):
        raise UnobservedFeature('Feature not in spectral wavelength range.')

    # Select the portion of the spectrum within the given bounds
    feature_indices = (lower_bound <= wavelength) & (wavelength <= upper_bound)
    feature_flux = flux[feature_indices]
    feature_wavelength = wavelength[feature_indices]

    peak_index = np.argmax(feature_flux)
    return feature_wavelength[peak_index], feature_flux[peak_index]


def get_feature_coordinates(wavelength, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wavelength (array): An array of wavelength values
        flux       (array): An array of flux values
        feature      (row): A row from the ``feature`` table

    Returns:
        A list with the blueward and redward wavelengths
        A list with the blueward and redward flux values
    """

    blue_wave, blue_flux = _get_peak_coordinates(
        wavelength, flux, feature['lower_blue'], feature['upper_blue'])

    red_wave, red_flux = _get_peak_coordinates(
        wavelength, flux, feature['lower_red'], feature['upper_red'])

    return [blue_wave, red_wave], [blue_flux, red_flux]


def get_continuum_func(blue_wave, red_wave, blue_flux, red_flux):
    """Fit the pseudo continuum for a given feature

    Args:
        blue_wave (float): Starting wavelength of the feature
        red_wave  (float): Flux value of spectrum at ``blue_wave``
        blue_flux (float): Ending wavelength of the feature
        red_flux  (float): Flux value of spectrum at ``red_flux``

    Return:
        A linear function fit to the flux peaks bounding a feature
    """

    m = (red_flux - blue_flux) / (red_wave - blue_wave)
    b = blue_flux - m * blue_wave
    return lambda wave: m * np.array(wave) + b


# noinspection PyTypeChecker, PyUnresolvedReferences
def calc_pew(wavelength, flux, feat_start, feat_end, cont_func):
    """Generic function for calculating the pseudo equivalent width

    The ``cont_func`` argument should be a vectorized function that accepts
    an array of wavelengths in the same units as ``wavelength`` and
    returns flux in the same units as ``flux``.

    Args:
        wavelength (array): An array of wavelength values
        flux       (array): An array of flux values
        feat_start (float): The starting wavelength of the feature
        feat_end   (float): The ending wavelength of the feature
        cont_func   (func): Function returning continuum flux for wavelength

    Returns:
       The pseudo equivalent width as a float
    """

    # Select the portion of the spectrum within the given bounds
    indices = (feat_start < wavelength) & (wavelength < feat_end)
    feature_wave = wavelength[indices]
    feature_flux = flux[indices]

    # Normalize the spectrum and calculate the EW
    continuum_flux = cont_func(feature_wave)
    normalized_flux = feature_flux / continuum_flux
    return np.trapz(normalized_flux, feature_wave)


def _feature_table_pew(wavelength, flux, feature_table):
    """Calculate the pseudo equivalent width for a table of features

    See the ``equivalent_width.feature`` table for an example of the
    ``feature_table`` argument.

    Args:
        wavelength    (array): An array of wavelength values
        flux          (array): An array of flux values
        feature_table (Table): A table defining spectral features

    Returns:
       The pseudo equivalent width as a float
    """

    ew_values = []
    for feature in feature_table:
        try:
            feat_wave, feat_flux = get_feature_coordinates(wavelength, flux,
                                                           feature)

        except UnobservedFeature:
            ew_values.append(-99)

        else:
            cont_func = get_continuum_func(*feat_wave, *feat_flux)
            ew = calc_pew(wavelength, flux, feat_wave[0], feat_wave[1],
                          cont_func)
            ew_values.append(ew)

    return ew_values


def tabulate_pew(time, wavelength, flux, feature_table):
    """Tabulate the pseudo equivalent widths for multiple spectra / features

    See the ``equivalent_width.feature`` table for an example of the
    ``feature_table`` argument.

    Args:
        time           (list): A list of observed MJD dates for each spectrum
        wavelength    (array): A 2d list of wavelength values for each date
        flux          (array): A 2d list of flux values for each date
        feature_table (Table): A table defining spectral features

    Returns:
       A table of equivalent widths over time
    """

    # noinspection PyTypeChecker
    ew_values = np.array([_feature_table_pew(w, f, feature_table) for w, f in zip(wavelength, flux)])

    out_data = Table(
        names=feature_table['feature_name'],
        rows=ew_values,
        masked=True)

    out_data.mask = np.transpose(ew_values < 0)
    out_data['time'] = time
    return out_data[['time'] + out_data.colnames[:-1]]
