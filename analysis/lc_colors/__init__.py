#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``lc_colors`` module fits light curves using a gaussian regression and
compares the fitted colors with modeled color values."""

from ._lc_prediction import (
    predict_band_flux,
    predict_c_15,
    predict_color,
    predict_light_curve)
from ._lc_regression import fit_gaussian_process
from ._tabulate import (
    get_color_times,
    tabulate_chisq,
    tabulate_delta_15
)

__all__ = [
    'get_color_times',
    'tabulate_chisq',
    'tabulate_delta_15',
    'predict_band_flux',
    'predict_c_15',
    'predict_color',
    'predict_light_curve',
    'fit_gaussian_process'
]
