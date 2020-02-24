#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``spectra_chisq`` module calculates the chi-squared values for modeled
spectra.

**This module is incomplete. It does not consider color warping.**

Function Documentation
----------------------
"""

from ._chisq import band_limits, calc_chisq, tabulate_chisq

__all__ = [
    'band_limits',
    'calc_chisq',
    'tabulate_chisq'
]
