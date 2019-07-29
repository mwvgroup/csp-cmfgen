#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Utility functions used across the analysis package."""

import numpy as np
from tqdm import tqdm


def get_spectra_for_id(data_release, obj_id):
    """Get spectral data for a given data id

    Args:
        data_release (module): An sndata data access module
        obj_id          (str): The ID of the desired object

    Returns:
        A list of observed MJD dates for each spectra
        A 2d list of wavelength values for each date
        A 2d list of flux values for each date
    """

    data = data_release.get_data_for_id(obj_id)
    obs_dates = list(set(data['date']))

    wavelength, flux, flux_err = [], [], []
    for date in obs_dates:
        data_for_date = data[data['date'] == date]
        wavelength.append(data_for_date['wavelength'])
        flux.append(data_for_date['flux'])
        flux_err.append(.1 * data_for_date['flux'])  # Todo: Parse flux errors

    return obs_dates, wavelength, flux, flux_err


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


def chisq_sum(obs, exp, error):
    return np.sum((obs - exp) ** 2 / (error ** 2))
