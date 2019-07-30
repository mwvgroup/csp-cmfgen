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
