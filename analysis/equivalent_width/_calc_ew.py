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
from ..exceptions import NoCSPData, UnobservedFeature

with open(Path(__file__).parent / 'features.yml') as infile:
    FEATURES = yaml.load(infile, Loader=yaml.FullLoader)


# noinspection PyTypeChecker, PyUnresolvedReferences
def get_peak_wavelength(
        wavelength, flux, lower_bound, upper_bound, behavior='min'):
    """Return wavelength of the maximum flux within given wavelength bounds

    The behavior argument can be used to select the 'min' or 'max' wavelength
    when there are multiple wavelengths having the same peak flux value. The
    default behavior is 'min'.

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        lower_bound  (float): Lower wavelength boundary
        upper_bound  (float): Upper wavelength boundary
        behavior       (str): Return the 'min' or 'max' wavelength

    Returns:
        The wavelength for the maximum flux value
        The maximum flux value
    """

    # Make sure the given spectrum spans the given wavelength bounds
    if (min(wavelength) > lower_bound) or (upper_bound > max(wavelength)):
        raise UnobservedFeature('Feature not in spectral wavelength range.')

    # Select the portion of the spectrum within the given bounds
    feature_indices = (lower_bound <= wavelength) & (wavelength <= upper_bound)
    feature_flux = flux[feature_indices]
    feature_wavelength = wavelength[feature_indices]

    peak_indices = np.argwhere(feature_flux == np.max(feature_flux))
    behavior_func = {'min': np.min, 'max': np.max}[behavior]
    return behavior_func(feature_wavelength[peak_indices])


def get_feature_bounds(wavelength, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        feature        (row): A dictionary defining feature parameters

    Returns:
        The starting wavelength of the feature
        The ending wavelength of the feature
    """

    feat_start = get_peak_wavelength(
        wavelength, flux, feature['lower_blue'], feature['upper_blue'], 'min')

    feat_end = get_peak_wavelength(
        wavelength, flux, feature['lower_red'], feature['upper_red'], 'max')

    return feat_start, feat_end


def fit_continuum_func(wavelength, flux, feat_start, feat_end):
    """Fit the pseudo continuum for a given feature

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        feat_start   (float): The starting wavelength of the feature
        feat_end     (float): The ending wavelength of the feature

    Return:
        A linear function fit to the flux values bounding a feature
    """

    if not ((feat_start in wavelength) and (feat_end in wavelength)):
        raise ValueError('Feature bounds not in wavelength array.')

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
    cont_func = fit_continuum_func(wavelength, flux, feat_start, feat_end)
    continuum_flux = cont_func(feature_wave)
    normalized_flux = feature_flux / continuum_flux
    pew = np.trapz(1 - normalized_flux, feature_wave)
    return pew, feat_start, feat_end


def create_pew_summary_table(models):
    """Create an astropy Table to store pew results

    Args:
        models (list): A list of sncosmo models

    Returns:
        An astropy table
    """

    # First row represents observed data
    model_names = [f'{m.source.name}' for m in models]
    model_versions = [f'{m.source.version}' for m in models]
    out_table = Table(
        data=[['OBSERVED'] + model_names, [''] + model_versions],
        names=['model', 'version'],
        dtype=['U100', 'U100'])

    return out_table


def tabulate_pew_spectrum(time, wave, flux, models=(), fix_boundaries=True):
    """Tabulate the observed and modeled pew for multiple features

    Args:
        time          (float): The time of the observed spectrum
        wave        (ndarray): An array of wavelength values
        flux        (ndarray): An array of flux values
        models         (list): A list of sncosmo models
        fix_boundaries (bool): Fix feature boundaries to observed values

    Returns:
        An astropy table
    """

    out_table = create_pew_summary_table(models)
    out_table['time'] = time

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
            # Shift time to beginning of explosion
            t0 = model.source.peakphase('csp_dr3_B')
            model_pew_results = calc_pew(
                wave, model.flux(time - t0, wave), feature, feat_start, feat_end)

            pew_data.append(model_pew_results)

        new_columns = np.transpose(pew_data)
        out_table[feat_name] = new_columns[0]
        out_table[feat_name + '_start'] = new_columns[1]
        out_table[feat_name + '_end'] = new_columns[2]

    return out_table


def tabulate_pew_spectra(
        data_release, models=(), fix_boundaries=True, verbose=True):
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
    id_iter = make_pbar(
        iterable=data_release.get_available_ids(),
        verbose=verbose,
        desc='Targets',
        total=total_targets)

    for obj_id in id_iter:
        data_table = data_release.get_data_for_id(obj_id)
        time, wavelength, flux = parse_spectra_table(data_table)

        try:
            for model in models:
                model.set(extebv=get_csp_ebv(obj_id))

            # Shift observed time to B-band peak
            time -= get_csp_t0(obj_id)

        except NoCSPData:
            continue

        spectra_iter = make_pbar(
            iterable=zip(time, wavelength, flux),
            verbose=verbose,
            desc='Spectra',
            position=1,
            total=len(time))

        pew_table = vstack([tabulate_pew_spectrum(*s, models, fix_boundaries) for s in spectra_iter])
        pew_table['obj_id'] = data_table.meta['obj_id']
        pew_data.append(pew_table)

    return vstack(pew_data)


def tabulate_peak_model_pew(models):
    """Tabulate the pew for each feature at time of B band maximum

    Args:
        models (list): A list of sncosmo models

    Returns:
       A table of equivalent widths for each feature
    """

    out_table = create_pew_summary_table(models)
    out_table.remove_row(0)  # Remove row for observed data

    for feat_name, feature in FEATURES.items():
        pew_data = []
        for model in models:
            time = model.source.peakphase('csp_dr3_B')
            wave = model.source.interpolated_model()[1]
            wave = wave[(wave > 3000) & (wave < 10000)]
            model_pew_results = calc_pew(wave, model.flux(time, wave), feature)
            pew_data.append(model_pew_results)

        new_columns = np.transpose(pew_data)
        out_table[feat_name] = new_columns[0]
        out_table[feat_name + '_start'] = new_columns[1]
        out_table[feat_name + '_end'] = new_columns[2]

    return out_table
