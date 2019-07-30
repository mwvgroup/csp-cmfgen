#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Utility functions used across the analysis package."""

import sncosmo
from tqdm import tqdm


def make_pbar(iterable, verbose, **kwargs):
    """A wrapper for tqdm.tqdm

    Args:
        iterable (iterator): Data to iterate over
        verbose      (bool): Whether to display the progress bar
        Any other arguments for tqdm.tqdm

    Returns:
        An iterable
    """

    if verbose:
        return tqdm(iterable, **kwargs)

    else:
        return iterable


def get_effective_wavelength(band_name):
    """Get the effective wavelength for a given band

    Band name must be registered with SNCosmo

    Args:
        band_name (str): The name of a registered bandpass

    Returns:
        The effective wavelength in Angstroms
    """

    return sncosmo.get_bandpass(band_name).wave_eff


def parse_spectra_table(data):
    """Convert spectral data from an astropy table to arrays

    Args:
        data (Table): A table of spectral data from sncosmo

    Returns:
        A list of observed MJD dates for each spectra
        A 2d list of wavelength values for each date
        A 2d list of flux values for each date
    """

    obs_dates = list(set(data['date']))

    wavelength, flux = [], []
    for date in obs_dates:
        data_for_date = data[data['date'] == date]
        wavelength.append(data_for_date['wavelength'])
        flux.append(data_for_date['flux'])

    return obs_dates, wavelength, flux
