#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module summarizes the pseudo equivalent widths of various features for
a given target and compares them against values predicted by CMFGEN.
"""

from astropy.table import hstack

from ._calc_ew import tabulate_pew


def get_model_spectra(time, wavelength, source):
    """Return the model spectra for a given model at multiple times

    Args:
        time       (list): List of times for each returned spectra
        wavelength (list): A 2d list of wavelength values for each date
        source   (Source): An sncosmo source object

    Returns:
        A 2d list of flux values for each time value
    """

    return [source.flux(t, w) for t, w in zip(time, wavelength)]


# noinspection PyTypeChecker
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

    out_tables = [tabulate_pew(time, wavelength, flux, feature_table)]
    for source in sources:
        source_flux = get_model_spectra(time, wavelength, source)
        ew_table = tabulate_pew(time, wavelength, source_flux, feature_table)

        for col_name in ew_table.colnames:
            ew_table.rename_column(col_name, col_name + f'_C{source.version}')

        out_tables.append(ew_table)

    return hstack(out_tables)
