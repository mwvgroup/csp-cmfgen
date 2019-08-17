#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``equivalent_width`` module measures the pseudo equivalent widths of
observed and modeled spectra.
"""

from ._calc_ew import (
    FEATURES as features,
    UnobservedFeature,
    calc_pew,
    fit_continuum_func,
    get_feature_bounds,
    tabulate_peak_model_pew,
    tabulate_pew_spectra,
    tabulate_pew_spectrum)

__all__ = [
    'features',
    'UnobservedFeature',
    'calc_pew',
    'fit_continuum_func',
    'get_feature_bounds',
    'tabulate_peak_model_pew',
    'tabulate_pew_spectra',
    'tabulate_pew_spectrum'
]
