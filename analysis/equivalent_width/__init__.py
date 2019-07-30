#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module measures the pseudo equivalent widths of observed spectra and
compares them with CMFGEN models.
"""

from ._calc_ew import (
    FEATURES as features,
    UnobservedFeature,
    calc_pew,
    fit_continuum_func,
    get_feature_bounds,
    tabulate_pew,
    tabulate_pew_spectrum)
