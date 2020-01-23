#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Tests for the ``analysis.models`` module."""

from unittest import TestCase

import numpy as np
import sncosmo

from analysis import models


class GetModel(TestCase):
    """Test correct models are run by ``_get_model``"""

    def assertCorrectSource(self, source_name, version):
        """Assert the returned model has the given name and version"""

        returned_model = models._get_model(source_name, version=version)
        self.assertEqual(source_name, returned_model.name)
        self.assertEqual(version, returned_model.version)

    def runTest(self):
        """Test the correct source is returned for different CMFGEN versions"""

        source_name = 'CMFGEN'
        version = '1.02', '1.04', '1.4', '1.7'
        for v in version:
            self.assertCorrectSource(source_name, v)


class RegisterSources(TestCase):
    """Tests for the ``register_sources`` function"""

    @classmethod
    def setUpClass(cls):
        """Register models with ``sncosmo``"""

        # Force prevents errors if other classes have already registered models
        models.register_sources(force=True)
        cls.source_name = 'CMFGEN'
        cls.versions_to_test = '1.02', '1.04', '1.4', '1.7'

    def assertIsRegistered(self, source_name, version):
        """Assert the returned model has the given name and version"""

        returned_model = sncosmo.get_source(source_name, version=version)
        self.assertEqual(source_name, returned_model.name)
        self.assertEqual(version, returned_model.version)

    def test_model_versions_as_str(self):
        """Test correct source is registered for different CMFGEN versions"""

        for v in self.versions_to_test:
            self.assertIsRegistered(self.source_name, str(v))

    def test_model_versions_as_float(self):
        """Test correct source is registered for different CMFGEN versions"""

        for v in self.versions_to_test:
            self.assertIsRegistered(self.source_name, float(v))


class BaseSourceTestingClass(TestCase):
    """Tests for an arbitrary sncosmo Source object"""

    @classmethod
    def setUpClass(cls):
        cls.source_name = 'CMFGEN'
        cls.source_version = '1.4'
        cls.source = models._get_model(cls.source_name, cls.source_version)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the model's phase range is zero"""

        min_phase = self.source.minphase()
        print(min_phase)
        early_phase_flux = self.source.bandflux('sdssg', min_phase - 1)
        early_phase_is_zero = np.isclose(0, early_phase_flux, atol=1e-6)
        self.assertTrue(early_phase_is_zero, f'Non-zero flux {early_phase_flux}')

        max_phase = self.source.maxphase()
        late_phase_flux = self.source.bandflux('sdssg', max_phase + 1)
        lat_phase_is_zero = np.isclose(0, late_phase_flux, atol=1e-6)
        self.assertTrue(lat_phase_is_zero, f'Non-zero flux {late_phase_flux}')
