#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits light curves using a gaussian regression and compares
the fitted colors with colors predicted by CMFGEN."""

from ._gauss_regression import fit_gaussian_process, predict_light_curve
