#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script converts CMFGEN models from a directory of ascii files to a
single npz file. It will convert all directories (i.e., models) in the current
working directory.
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

ASCII_MODEL_DIR = Path(__file__).resolve().parent
NPZ_MODEL_DIR = ASCII_MODEL_DIR / 'NPZ_models'


def save_model_to_npz(in_dir, out_path):
    """Convert CMFGEN models from ASCII to npz format

    Args:
        in_dir   (Path): Directory of ASCII files for a given model
        out_path (Path): Output path of the npz file
    """

    out_path = out_path.with_suffix('.npz')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ascii_summary_table = Table.read(str(in_dir / 'AGE_FLUX_TAB'),
                                     format='ascii')
    phase = np.array(
        [float(x.rstrip('D0')) for x in ascii_summary_table['col1']])

    # Iterate over rows with file names for SEDs at different phase values
    flux = []
    wavelength = []
    for row in tqdm(ascii_summary_table, desc=f'Formatting {out_path.stem}'):
        file_name = row[1] + '.fl'
        flux_table = Table.read(in_dir / file_name, format='ascii')
        wavelength.append(np.array(flux_table['col1']))
        flux.append(np.array(flux_table['col2']))

    np.savez(out_path, phase=phase, wavelength=wavelength, flux=flux)

    # Down sample onto a uniform grid
    wavelength_grid = wavelength[0]
    flux_grid = []
    for w, f in zip(wavelength, flux):
        flux_grid.append(np.interp(wavelength_grid, w, f))

    grid_path = out_path.with_name(out_path.stem + '_grid.npz')
    np.savez(grid_path, phase=phase, wavelength=wavelength_grid,
             flux=flux_grid)


def format_models(out_dir):
    """Format all CMFGEN models included with this distribution for use by
    this package

    Args:
        out_dir (Path): Directory where output files are written
    """

    for sub_dir in ASCII_MODEL_DIR.glob('*/'):
        save_model_to_npz(sub_dir, out_dir / sub_dir.stem)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert CMFGEN models into NPZ files.')

    parser.add_argument(
        '-o', '--out_path',
        type=str,
        required=True,
        help='Directory where output files are written.')

    out_path = Path(parser.parse_args().out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    format_models(out_path)
