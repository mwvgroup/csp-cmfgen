#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module summarizes the pseudo equivalent widths of various features for
a given target and compares them against values predicted by CMFGEN.
"""


def get_model_spectra(time, source):
    """Return the model spectra for a given model at multiple times

    Args:
        time     (list): A list of observed MJD dates for each returned spectra
        source (Source): An sncosmo source object

    Returns:
        A 2d list of wavelength values for each time value
        A 2d list of flux values for each time valueeach date
    """

    pass


def compare_target_and_models(time, wavelength, flux, feature_table, sources):
    """Tabulate modeled and measured pseudo equivalent widths

    Args:
        time           (list): A list of observed MJD dates for each spectrum
        wavelength    (array): A 2d list of wavelength values for each date
        flux          (array): A 2d list of flux values for each date
        feature_table (Table): A table defining spectral features
        sources        (list): A list of sncosmo source objects

    Returns:
        A table of modeled and measured equivalent widths for each feature
    """

    pass
