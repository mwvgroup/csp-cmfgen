#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides generic functions for measuring the pseudo equivalent
width of arbitrary features for a one or more spectra.
"""

from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table, vstack

from ..utils import get_csp_ebv, get_csp_t0, make_pbar, parse_spectra_table

with open(Path(__file__).parent / 'features.yml') as infile:
    FEATURES = yaml.load(infile, Loader=yaml.FullLoader)


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

    # Make sure the given spectrum spans the given wavelength bounds
    if (min(wavelength) >= lower_bound) or (upper_bound >= max(wavelength)):
        raise UnobservedFeature('Feature not in spectral wavelength range.')

    # Select the portion of the spectrum within the given bounds
    feature_indices = (lower_bound <= wavelength) & (wavelength <= upper_bound)
    feature_flux = flux[feature_indices]
    feature_wavelength = wavelength[feature_indices]

    peak_index = np.argmax(feature_flux)
    return feature_wavelength[peak_index]


def get_feature_bounds(wavelength, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wavelength (array): An array of wavelength values
        flux       (array): An array of flux values
        feature      (row): A dictionary defining feature parameters

    Returns:
        The starting wavelength of the feature
        The ending wavelength of the feature
    """

    feat_start = get_peak_coordinates(
        wavelength, flux, feature['lower_blue'], feature['upper_blue'])

    feat_end = get_peak_coordinates(
        wavelength, flux, feature['lower_red'], feature['upper_red'])

    return feat_start, feat_end


def get_continuum_func(wavelength, flux, feat_start, feat_end):
    """Fit the pseudo continuum for a given feature

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        feat_start   (float): The starting wavelength of the feature
        feat_end     (float): The ending wavelength of the feature

    Return:
        A linear function fit to the flux values bounding a feature
    """

    blue_flux = flux[wavelength == feat_start]
    red_flux = flux[wavelength == feat_end]
    m = (red_flux - blue_flux) / (feat_end - feat_start)
    b = blue_flux - m * feat_start
    return lambda wave: m * np.array(wave) + b


# noinspection PyTypeChecker, PyUnresolvedReferences
def calc_pew(wavelength, flux, feature=None, feat_start=None, feat_end=None):
    """Generic function for calculating the pseudo equivalent width

    If ``feat_start`` and ``feat_end`` are given, use them to determine the
    continuum flux. If they aren't specified by ``feature`` is, estimate
    the boundaries of the feature.

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        feature       (dict): A dictionary defining feature parameters
        feat_start   (float): The starting wavelength of the feature
        feat_end     (float): The ending wavelength of the feature

    Returns:
       The pseudo equivalent width as a float
    """

    if (feat_start is None) or (feat_end is None):
        if feature is None:
            raise ValueError(
                'Must specify either feature bounds or feature dict.')

        feat_start, feat_end = get_feature_bounds(wavelength, flux, feature)

    # Select the portion of the spectrum within the given bounds
    indices = (feat_start <= wavelength) & (wavelength <= feat_end)
    feature_wave = wavelength[indices]
    feature_flux = flux[indices]

    # Normalize the spectrum and calculate the EW
    cont_func = get_continuum_func(wavelength, flux, feat_start, feat_end)
    continuum_flux = cont_func(feature_wave)
    normalized_flux = feature_flux / continuum_flux
    pew = feat_end - feat_start - np.trapz(normalized_flux, feature_wave)
    return pew, feat_start, feat_end


def create_empty_pew_table(models):
    """Create an astropy Table to store pew results

    Args:
        models (list): A list of sncosmo models

    Returns:
        An astropy table
    """

    # First row represents observed data
    model_names = [f'{m.source.name}' for m in models]
    model_versions = [f'{m.source.name}' for m in models]
    out_table = Table(
        data=[['obs'] + model_names, [''] + model_versions],
        names=['model', 'version'],
        dtype=['U100', 'U100'])

    return out_table


def tabulate_pew_spectrum(time, wave, flux, models, fix_boundaries):
    """Tabulate the observed and modeled pew for multiple features

    Args:
        time          (float): The time of the observed spectrum
        wave  (ndarray): An array of wavelength values
        flux        (ndarray): An array of flux values
        models         (list): A list of sncosmo models
        fix_boundaries (bool): Fix feature boundaries to observed values

    Returns:
        An astropy table
    """

    out_table = create_empty_pew_table(models)
    for feat_name, feature in FEATURES.items():
        # Calculate pew for observed data
        try:
            pew, feat_start, feat_end = calc_pew(wave, flux, feature)
            pew_data = [[pew, feat_start, feat_end]]

        except UnobservedFeature:
            continue

        # Reset feature boundaries
        if not fix_boundaries:
            feat_start, feat_end = None, None

        # Calculate pew for models
        for model in models:
            model_pew_results = calc_pew(
                wave, model.flux(time, wave), feature, feat_start, feat_end)

            pew_data.append(model_pew_results)

        new_columns = np.transpose(pew_data)
        out_table[feat_name] = new_columns[0]
        out_table[feat_name + '_start'] = new_columns[1]
        out_table[feat_name + '_end'] = new_columns[2]

    return out_table


def shift_models_to_data(obj_id, models):
    """Set models to observed extinction and t0 values

    Args:
        obj_id  (str): A CSP object id
        models (list): A list of sncomso models
    """

    for model in models:
        model.set(extebv=get_csp_ebv(obj_id))

        # Todo: This time shift doesn't work
        model.set(t0=get_csp_t0(obj_id))


def tabulate_pew(data_release, models, fix_boundaries, verbose=True):
    """Tabulate the pseudo equivalent widths for multiple spectra / features

    Args:
        data_release (module): An sndata data release
        models         (list): A list of sncosmo models
        fix_boundaries (bool): Fix feature boundaries to observed values
        verbose        (bool): Whether to display a progress bar

    Returns:
       A table of equivalent widths over time
    """

    pew_data = []
    total_targets = len(data_release.get_available_ids())
    data_iter = make_pbar(
        iterable=data_release.iter_data(),
        verbose=verbose,
        desc='Targets',
        total=total_targets)

    for data_table in data_iter:
        try:
            shift_models_to_data(data_table.meta['obj_id'], models)

        except ValueError:
            continue

        data_arrays = parse_spectra_table(data_table)
        spectra_iter = make_pbar(
            iterable=zip(*data_arrays),
            verbose=verbose,
            desc='Spectra',
            position=1,
            total=len(data_arrays[0]))

        pew_table = vstack(
            [tabulate_pew_spectrum(*s, models, fix_boundaries) for s in
             spectra_iter])

        pew_table['obj_id'] = data_table.meta['obj_id']
        pew_data.append(pew_table)

    return vstack(pew_data)
