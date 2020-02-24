#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``utils`` module provides a variety of general utility functions used across the analysis package.

Usage Example
-------------

.. code-block:: python
   :linenos:

   from sndata.csp import dr3

   from analysis import utils

   # Convert dates from the Snoopy or MJD formats to JD format
   dates = [350, 53350]  # snoopy and MJD respectively
   jd_dates = utils.convert_to_jd(dates)
   print(jd_dates)

   # Get certain CSP data for a given target
   t0 = utils.get_csp_t0(('2005kc'))
   ebv = utils.get_csp_ebv(('2005kc'))

   # Iterate over photometric data while only including targets with a published t0 and E(B-V)
   for data in dr3.iter_data(filter_func=utils.filter_has_csp_data)

Function Documentation
----------------------
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.time import Time
from sndata.csp import DR1, DR3

from .exceptions import NoCSPData

dr1 = DR1()
dr3 = DR3()


@np.vectorize
def convert_to_jd(date):
    """Convert MJD and Snoopy date formats into Julian date (JD) format

    Args:
        date (float): Time stamp in JD, MJD, or SNPY format

    Returns:
        The time value in JD format
    """

    snoopy_offset = 53000  # Snoopy date format is MDJ minus 53000
    mjd_offset = 2400000.5  # MJD date format is JD minus 2400000.5
    date_format = 'mjd'

    if date < snoopy_offset:
        date += snoopy_offset

    elif date > mjd_offset:
        date_format = 'jd'

    t = Time(date, format=date_format)
    t.format = 'jd'
    return t.value


def get_csp_t0(obj_id):
    """Get the t0 value published by CSP DR3 for a given object

    Args:
        obj_id (str): The object Id value

    Returns:
        The published date of maximum in JD format
    """

    dr3.download_module_data()
    params = dr3.load_table(3)
    params = params[~params['T(Bmax)'].mask]  # Drop masked values
    if obj_id not in params['SN']:
        raise NoCSPData(f'No published t0 for {obj_id}')

    # Return date in Julian date format
    return convert_to_jd(params[params['SN'] == obj_id]['T(Bmax)'][0])


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

    return True


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
        - A sorted list of observed MJD dates for each spectra
        - A 2d list of wavelength values for each date
        - A 2d list of flux values for each date
    """

    obs_dates = sorted(set(data['time']))

    wavelength, flux = [], []
    data.sort('wavelength')
    for date in obs_dates:
        data_for_date = data[data['time'] == date]
        wavelength.append(np.array(data_for_date['wavelength']))
        flux.append(np.array(data_for_date['flux']))

    return convert_to_jd(obs_dates), np.array(wavelength), np.array(flux)


def calc_model_chisq(data, result, model):
    """Calculate the chi-squared for a given data table and model

    Chi-squared is calculated using parameter values from ``model``. Degrees
    of freedom are calculated using the number of varied parameters specified
    is the ``result`` object.

    Args:
        data    (Table): An sncosmo input table
        model   (Model): An sncosmo Model
        result (Result): sncosmo fitting result

    Returns:
        The un-normalized chi-squared
        The number of data points used in the calculation
    """

    data = deepcopy(data)

    # Drop any data that is not withing the model's range
    min_band_wave = np.array([sncosmo.get_bandpass(b).minwave() for b in data['band']])
    max_band_wave = np.array([sncosmo.get_bandpass(b).maxwave() for b in data['band']])
    data = data[
        (data['time'] >= model.mintime()) &
        (data['time'] <= model.maxtime()) &
        (min_band_wave >= model.minwave()) &
        (max_band_wave <= model.maxwave())
        ]

    if len(data) == 0:
        raise ValueError('No data within model range')

    return sncosmo.chisq(data, model), len(data) - len(result.vparam_names)
