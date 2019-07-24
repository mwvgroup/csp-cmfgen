#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module measures the pseudo equivalent widths of observed spectra and
compares them with CMFGEN models.
"""

from astropy.table import Table as _Table

from ._calc_ew import UnobservedFeature, calc_pew, get_continuum_func, \
    get_feature_coordinates, tabulate_pew
from ._compare_models import compare_target_and_models

feature_table = _Table({
    'feature_name': ['pW1', 'pW2', 'pW3', 'pW4', 'pW5', 'pW6', 'pW7', 'pW8'],
    'feature_id': ['Ca ii H&K', 'Si ii λ4130', 'Mg ii, Fe ii', 'Fe ii, Si ii',
                   'S ii λ5449, λ5622', 'Si ii λ5972', 'Si ii λ6355',
                   'Ca ii IR triplet'],
    'lower_blue': [3500, 3900, 3900, 4500, 5150, 5550, 5800, 7500],
    'upper_blue': [3800, 4000, 4150, 4700, 5300, 5700, 6000, 8000],
    'lower_red': [3900, 4000, 4450, 5050, 5500, 5800, 6200, 8200],
    'upper_red': [4100, 4150, 4700, 5550, 5700, 6000, 6600, 8900]
})
