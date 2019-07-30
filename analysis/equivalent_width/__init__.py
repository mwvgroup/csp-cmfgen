#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module measures the pseudo equivalent widths of observed spectra and
compares them with CMFGEN models.
"""

from ._calc_ew import (
    UnobservedFeature,
    calc_pew,
    compare_target_and_models,
    features,
    get_continuum_func,
    get_feature_coordinates,
    tabulate_pew)
