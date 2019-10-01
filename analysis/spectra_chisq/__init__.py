#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``spectra_chisq`` module calculates the chi-squared values for modeled
spectra.

**This module is incomplete. It does not consider color warping.**
"""

from ._chisq import calc_chisq, band_limits, tabulate_chisq
from .photometry_comparison import photometry_to_spectra_time, tabulate_synthetic_photometry

__all__ = [
    'band_limits',
    'calc_chisq',
    'tabulate_chisq',
    #'photometry_to_spectra_time',
    #'tabulate_synthetic_photometry'
]
