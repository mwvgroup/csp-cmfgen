#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines custom SNCosmo.Source objects for four different CMFGEN
based models. Source objects include:

SubChandra_1: A Sub-Chandrasekhar model for 1.04 solar mass progenitors
SubChandra_2: A Sub-Chandrasekhar model for 1.02 solar mass progenitors
Chandra: A Chandrasekhar model for 1.4 solar mass progenitors
SuperChandra: A Super-Chandrasekhar model for 1.7 solar mass progenitors
"""

import os
from copy import deepcopy
from glob import glob
from zipfile import ZipFile

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
SUB_1_PATH = os.path.join(FILE_DIR, 'DDC0_grid.npz')
SUB_2_PATH = os.path.join(FILE_DIR, 'DDC10_grid.npz')
CHAN_PATH = os.path.join(FILE_DIR, 'DDC15_grid.npz')
SUP_PATH = os.path.join(FILE_DIR, 'SCH05p5_grid.npz')

for path in (SUB_1_PATH, SUB_2_PATH, CHAN_PATH, SUP_PATH):
    if os.path.exists(path):
        continue

    print('Unzipping models...')
    for path in glob(os.path.join(FILE_DIR, '*.zip')):
        with ZipFile(path) as zip_ref:
            zip_ref.extractall(FILE_DIR)

    break


class GenericSource(sncosmo.Source):
    _param_names = ['x0']
    param_names_latex = ['x0']  # used in plotting display

    def __init__(self, path, name, version):
        self._path = path
        self._parameters = [1]
        self.name = name
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


SubChandra_1 = GenericSource(SUB_1_PATH, 'SubChandra', 1.04)
SubChandra_2 = GenericSource(SUB_2_PATH, 'SubChandra', 1.02)
Chandra = GenericSource(CHAN_PATH, 'Chandra', 1.4)
SuperChandra = GenericSource(SUP_PATH, 'SuperChandra', 1.7)
