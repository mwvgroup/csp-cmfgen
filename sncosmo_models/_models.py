#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines custom SNCosmo.Source objects for four different CMFGEN
based models.
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
import sncosmo
from astropy.table import Table
from scipy.interpolate import RectBivariateSpline
from tqdm import tqdm

FILE_DIR = Path(__file__).resolve().parent
ASCII_MODEL_DIR = FILE_DIR / 'ascii_models'
NPZ_MODEL_DIR = FILE_DIR / 'NPZ_models'

VERSION_PATHS = {
    '1.02': NPZ_MODEL_DIR / 'DDC10_grid.npz',
    '1.04': NPZ_MODEL_DIR / 'DDC0_grid.npz',
    '1.4': NPZ_MODEL_DIR / 'DDC15_grid.npz',
    '1.7': NPZ_MODEL_DIR / 'SCH05p5_grid.npz'
}


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


def format_models(force=True):
    """Format all CMFGEN models included with this distribution for use by
    this package
    """

    if not NPZ_MODEL_DIR.exists() or force:
        tqdm.write('Formatting models for use by package...')
        for sub_dir in ASCII_MODEL_DIR.glob('*/'):
            save_model_to_npz(sub_dir, NPZ_MODEL_DIR / sub_dir.stem)


class GenericSource(sncosmo.Source):
    """SNCosmo source class for a CMFGEN model"""

    _param_names = ['x0']
    param_names_latex = ['x0']  # used in plotting display

    def __init__(self, path, version):
        self._path = path
        self._parameters = [1]
        self.name = 'CMFGEN'
        self.version = version

        model_data = np.load(path)
        self._phase = model_data['phase']
        self._wave = model_data['wavelength']
        self._model_flux = model_data['flux']

        self.spline = RectBivariateSpline(
            self._phase,
            self._wave,
            self._model_flux,
            kx=3, ky=3)

    def interpolated_model(self):
        """Return the phase, wavelength, and flux values used for fitting

        Returns the data with flux values down-sampled onto the same wavelength
        values for all phase values. The returned model is the same as the data
        used when fitting light curves.

        Returns:
            A 1D array of phase values
            A 1D array of wavelength values
            A 2D array of flux values
        """

        return deepcopy([self._phase, self._wave, self._model_flux])

    def original_model(self):
        """Return the phase, wavelength, and flux values of the model

        Returns the un-gridded model that this class is based on without any
        down-sampling. Flux values span different wavelength ranges
        for different phases.

        Returns:
            A 1D array of phase values
            A 2D array of wavelength values
            A 2D array of flux values
        """

        path = str(self._path).replace('_grid.npz', '.npz')
        data = np.load(path)
        return data['phase'], data['wavelength'], data['flux']

    def _flux(self, phase, wave):
        # Returns the flux corresponding to the given phase and wavelength

        amplitude = self._parameters
        return amplitude * self.spline(phase, wave)


def get_model(name=None, version=None):
    """Return an SNCosmo CMFGEN source for a given model

    Args:
        name          (None): A dummy variable for compatibility with sncosmo
        version (str, float): The version (mass) of the CMFGEN model

    Returns:
        A source object for the specified model
    """

    return GenericSource(VERSION_PATHS[str(version)], version)


def register_sources(force=False):
    """Register CMFGEN models with sncosmo

    Versions include: ``salt2_phase``, ``color_interpolation``

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    for version in VERSION_PATHS.keys():
        sncosmo.register_loader(
            data_class=sncosmo.Source,
            name='CMFGEN',
            func=get_model,
            version=str(version),
            force=force)

        sncosmo.register_loader(
            data_class=sncosmo.Source,
            name='CMFGEN',
            func=get_model,
            version=float(version),
            force=force)
