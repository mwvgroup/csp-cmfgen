#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``models`` module defines custom Source objects for using various
supernova models with the ``sncosmo`` python package.

Available Models
----------------

+--------+---------+--------------------------------------------------+
| Name   | Version | Description                                      |
+========+=========+==================================================+
| CMFGEN |  1.02   | Explosion Model for a 1.02 solar mass progenitor |
+--------+---------+--------------------------------------------------+
| CMFGEN |  1.04   | Explosion Model for a 1.04 solar mass progenitor |
+--------+---------+--------------------------------------------------+
| CMFGEN |  1.4    | Explosion Model for a 1.4  solar mass progenitor |
+--------+---------+--------------------------------------------------+
| CMFGEN |  1.7    | Explosion Model for a 1.7  solar mass progenitor |
+--------+---------+--------------------------------------------------+

Usage Example
-------------

.. code-block:: python
   :linenos:

   import sncosmo
   from matplotlib import pyplot as plt

   from analysis import models

   # Make sncosmo aware of the sncosmo models
   models.register_sources()

   # Initialize a CMFGEN model where the version is the model mass
   source = sncosmo.get_source('CMFGEN', version='1.04')
   model = sncosmo.Model(source=source)

   # run the fit
   data = sncosmo.load_example_data()
   result, fitted_model = sncosmo.fit_lc(
       data, model,
       ['z', 't0', 'x0'],  # parameters of model to vary
       bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)

   # Plot results
   fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
   plt.show()

Function Documentation
----------------------
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline

NPZ_MODEL_DIR = Path(__file__).resolve().parent / 'NPZ_models'
VERSION_PATHS = {
    '1.02': NPZ_MODEL_DIR / 'DDC10_grid.npz',
    '1.04': NPZ_MODEL_DIR / 'DDC0_grid.npz',
    '1.4': NPZ_MODEL_DIR / 'DDC15_grid.npz',
    '1.7': NPZ_MODEL_DIR / 'SCH05p5_grid.npz'
}


# noinspection PyMissingConstructor
class GenericSource(sncosmo.Source):
    """SNCosmo source class for a CMFGEN model"""

    _param_names = ['x0']
    param_names_latex = ['x0']  # used in plotting display

    def __init__(self, path, version):
        """Load a CMFGEN spectral template from file

        Args:
            path: The path of the spectral template
            version: The version number to assign the loaded template
        """

        super().__init__()

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
        data = np.load(path, allow_pickle=True)
        return data['phase'], data['wavelength'], data['flux']

    def _flux(self, phase, wave):
        # Returns the flux corresponding to the given phase and wavelength

        amplitude = self._parameters[0]
        flux = amplitude * self.spline(phase, wave)
        flux[(phase < self.minphase()) | (self.maxphase() < phase)] = 0
        return flux


# noinspection PyUnusedLocal
def _get_source(name=None, version=None):
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

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    for version in VERSION_PATHS.keys():
        # String version key
        sncosmo.register_loader(
            data_class=sncosmo.Source,
            name='CMFGEN',
            func=_get_source,
            version=str(version),
            force=force)

        # Float version key
        sncosmo.register_loader(
            data_class=sncosmo.Source,
            name='CMFGEN',
            func=_get_source,
            version=float(version),
            force=force)


def _unzip_models():
    """Decompress any models that haven't already been decompressed"""

    from zipfile import ZipFile

    file_paths = list(VERSION_PATHS.values())
    file_paths += [f.with_name(f.name.replace('_grid', '')) for f in
                   file_paths]

    for path in file_paths:
        if not path.exists():
            print(f'Unzipping model: {path}')
            with ZipFile(str(path) + '.zip') as zip_ref:
                zip_ref.extractall(NPZ_MODEL_DIR)


_unzip_models()
SubChandra_1 = _get_source(version=1.04)
SubChandra_2 = _get_source(version=1.02)
Chandra = _get_source(version=1.4)
SuperChandra = _get_source(version=1.7)
