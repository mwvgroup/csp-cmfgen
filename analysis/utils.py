#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Utility functions used across the analysis package."""

import numpy as np
from tqdm import tqdm


def make_pbar(iterable, verbose, **kwargs):
    """A wrapper for tqdm.tqdm

    Args:
        iterable (iterator): Data to iterate over
        verbose      (bool): Whether to display the progress bar
        Any other arguments for tqdm.tqdm

    Returns:
        An iterable
    """

    if verbose:
        return tqdm(iterable, **kwargs)

    else:
        return iterable


def chisq_sum(obs, exp, error):
    return np.sum((obs - exp) ** 2 / (error ** 2))
