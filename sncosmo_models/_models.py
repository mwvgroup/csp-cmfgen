#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines custom SNCosmo.Source objects for four different CMFGEN
based models.
"""

from copy import deepcopy
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline

FILE_DIR = Path(__file__).resolve().parent
path_1_02 = FILE_DIR / 'DDC10_grid.npz'
path_1_04 = FILE_DIR / 'DDC0_grid.npz'
path_1_4 = FILE_DIR / 'DDC15_grid.npz'
path_1_7 = FILE_DIR / 'SCH05p5_grid.npz'
VERSIONS = (1.02, 1.04, 1.4, 1.7)
PATHS = (path_1_02, path_1_04, path_1_4, path_1_7)

for _path in PATHS:
    if not _path.exists():
        print(f'Unzipping model: {_path}')
        with ZipFile(str(_path) + '.zip') as zip_ref:
            zip_ref.extractall(FILE_DIR)


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

    def gridded_model(self):
        """Return the phase, wavelength, and flux values used for fitting

        Returns the data with flux values down-sampled onto the same wavelength
        values for all phase values. The returned model is the same as the data
        used when fitting light curves.

        Returns:
            A 1D array of phase values
            A 1D array of wavelength values
            A 2D array of flux values
        """

        data = (self._phase, self._wave, self._model_flux)
        return (deepcopy(array) for array in data)

    def raw_model(self):
        """Return the phase, wavelength, and flux values of the model

        Returns the raw model that this class is based on without any
        down-sampling. Flux values span different wavelength ranges
        for different phases.

        Returns:
            A 1D array of phase values
            A 2D array of wavelength values
            A 2D array of flux values
        """

        path = self._path.rstrip('_grid.npz') + '.npz'
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

    version_paths = {str(v): p for v, p in zip(VERSIONS, PATHS)}
    return GenericSource(version_paths[str(version)], version)


def register_sources(force=False):
    """Register CMFGEN models with sncosmo

    Versions include: ``salt2_phase``, ``color_interpolation``

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    for version in (1.04, 1.02, 1.4, 1.7):
        sncosmo.register_loader(
            data_class=sncosmo.Source,
            name='CMFGEN',
            func=get_model,
            version=version,
            force=force)
