#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script converts ascii model directories to a single npz file. It will
convert all directories in the current working path.
"""

import os

import numpy as np
from astropy.table import Table
from tqdm import tqdm


def save_model_to_npy(in_dir, out_path):
    """Convert a directory of ascii models files to a single npz file

    Args:
        in_dir   (str): Directory path of input ascii model
        out_path (str): Where to write the output file
    """

    # Read summary table for the model
    tab = Table.read(os.path.join(in_dir, 'AGE_FLUX_TAB'), format='ascii')

    # Populate arrays with model values
    phase = np.array([float(x.rstrip('D0')) for x in tab['col1']])
    flux = []
    wavelength = []
    for row in tqdm(tab):
        flux_table = Table.read(os.path.join(in_dir, row[1] + '.fl'),
                                format='ascii')
        wavelength.append(np.array(flux_table['col1']))
        flux.append(np.array(flux_table['col2']))

    np.savez(out_path, phase=phase, wavelength=wavelength, flux=flux)

    # Down sample onto a uniform grid
    wavelength_grid = wavelength[0]
    flux_grid = []
    for w, f in zip(wavelength, flux):
        flux_grid.append(np.interp(wavelength_grid, w, f))

    # Save data
    grid_path = out_path.rstrip('.npz') + '_grid.npz'
    np.savez(grid_path, phase=phase, wavelength=wavelength_grid,
             flux=flux_grid)


if __name__ == '__main__':
    model_dir = './'
    for sub_dir in os.listdir(model_dir):
        save_model_to_npy(os.path.join(model_dir, sub_dir), sub_dir)
