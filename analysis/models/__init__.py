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

>>> import sncosmo
>>> from matplotlib import pyplot as plt
>>>
>>> from analysis import models
>>>
>>> # Make sncosmo aware of the sncosmo models
>>> models.register_sources()
>>>
>>> # Initialize a CMFGEN model where the version is the model mass
>>> source = sncosmo.get_source('CMFGEN', version=1.04)
>>> model = sncosmo.Model(source=source)
>>>
>>> # run the fit
>>> data = sncosmo.load_example_data()
>>> result, fitted_model = sncosmo.fit_lc(
>>>     data, model,
>>>     ['z', 't0', 'x0'],  # parameters of model to vary
>>>     bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)
>>>
>>> # Plot results
>>> fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
>>> plt.show()
"""

from ._models import get_model as _get_model, register_sources

__all__ = ['register_sources']


def _unzip_models():
    """Decompress any models that haven't already been decompressed"""

    from zipfile import ZipFile
    from ._models import VERSION_PATHS, NPZ_MODEL_DIR

    file_paths = list(VERSION_PATHS.values())
    file_paths += [f.with_name(f.name.replace('_grid', '')) for f in
                   file_paths]

    for path in file_paths:
        if not path.exists():
            print(f'Unzipping model: {path}')
            with ZipFile(str(path) + '.zip') as zip_ref:
                zip_ref.extractall(NPZ_MODEL_DIR)


_unzip_models()
SubChandra_1 = _get_model(version=1.04)
SubChandra_2 = _get_model(version=1.02)
Chandra = _get_model(version=1.4)
SuperChandra = _get_model(version=1.7)
