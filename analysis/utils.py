#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""General utility functions used across the analysis package."""

import numpy as np
import sncosmo
from sndata.csp import dr1, dr3
from tqdm import tqdm

from .exceptions import NoCSPData


def filter_has_csp_data(data_table):
    """A filter function for an SNData table iterator

    Returns whether the object ID associated with a data table has an
    available t0 and E(B - V) value.

    Args:
        data_table (Table): A table from sndata

    Returns:
        A boolean
    """

    obj_id = data_table.meta['obj_id']
    try:
        get_csp_t0(obj_id)
        get_csp_ebv(obj_id)

    except NoCSPData:
        return False

    else:
        return True


@np.vectorize
def convert_to_jd(time):
    """Convert MJD and Snoopy dates into JD

    Args:
        time (float): Time stamp in JD, MJD, or Snoopy format

    Returns:
        The time value in JD format
    """

    # Snoopy time format
    if time < 53000:
        return time + 53000 + 2400000.5

    # Snoopy time format
    elif 53000 < time < 2400000.5:
        return time + 2400000.5

    else:
        return time


def get_csp_t0(obj_id):
    """Get the t0 value published by CSP DR3 for a given object

    Args:
        obj_id (str): The object Id value

    Returns:
        The published MJD of maximum minus 53000
    """

    dr3.download_module_data()
    params = dr3.load_table(3)
    params = params[~params['T(Bmax)'].mask]
    if obj_id not in params['SN']:
        raise NoCSPData(f'No published t0 for {obj_id}')

    return params[params['SN'] == obj_id]['T(Bmax)'][0]


def get_csp_ebv(obj_id):
    """Get the E(B - V) value published by CSP DR1 for a given object

    Args:
        obj_id (str): The object Id value

    Returns:
        The published E(B - V) value
    """

    dr1.download_module_data()
    extinction_table = dr1.load_table(1)
    if obj_id not in extinction_table['SN']:
        raise NoCSPData(f'No published E(B-V) for {obj_id}')

    data_for_target = extinction_table[extinction_table['SN'] == obj_id]
    return data_for_target['E(B-V)'][0]


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

    Band name must be registered with SNCosmo.

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
    data.sort('wavelength')
    for date in obs_dates:
        data_for_date = data[data['date'] == date]
        wavelength.append(data_for_date['wavelength'])
        flux.append(data_for_date['flux'])

    obs_dates = np.array(obs_dates)
    return obs_dates, np.array(wavelength), np.array(flux)
