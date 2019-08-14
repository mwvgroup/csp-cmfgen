#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits light curves using a gaussian regression and compares
the fitted colors with colors predicted by CMFGEN."""

from ._chi_squared import get_color_times
from ._chi_squared import tabulate_residuals
from ._lc_regression import fit_gaussian_process
from ._lc_regression import predict_band_flux
from ._lc_regression import predict_color
from ._lc_regression import predict_light_curve
