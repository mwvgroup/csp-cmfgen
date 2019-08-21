#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``spectra_chisq`` module calculates the chi-squared values for modeled
spectra.

**This module is incomplete. It does not consider color warping.**
"""

from ._band_chisq import band_chisq, band_limits, tabulate_chi_squared

__all__ = [
    'band_limits',
    'band_chisq',
    'tabulate_chi_squared'
]
