#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``lc_colors`` module fits light curves using a gaussian regression and
compares the fitted colors with modeled color values."""

from ._chi_squared import get_color_times, tabulate_residuals
from ._lc_regression import (
    fit_gaussian_process,
    predict_band_flux,
    predict_color,
    predict_light_curve)

__all__ = [
    'get_color_times',
    'tabulate_residuals',
    'fit_gaussian_process',
    'predict_band_flux',
    'predict_color',
    'predict_light_curve'
]