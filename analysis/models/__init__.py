#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines custom SNCosmo.Source objects for four different CMFGEN
based models.

To use the model:
    import sncosmo
    from matplotlib import pyplot as plt

    import sncosmo_models

    # Check available versions
    print(sncosmo_models.versions)

    # Make sncosmo aware of the sncosmo models
    sncosmo_models.register_sources()

    # Initialize a CMFGEN model where the version is the model mass
    source = sncosmo.get_source('CMFGEN', version=1.04)
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
"""

from ._models import get_model as _get_model, register_sources


def _unzip_models():
    """Decompress any models that haven't already been decompressed"""

    from zipfile import ZipFile
    from ._models import VERSION_PATHS, NPZ_MODEL_DIR

    file_paths = list(VERSION_PATHS.values())
    file_paths += [f.with_name(f.name.replace('_grid', '')) for f in file_paths]

    for path in file_paths:
        if not path.exists():
            print(f'Unzipping model: {path}')
            with ZipFile(str(path) + '.zip') as zip_ref:
                zip_ref.extractall(NPZ_MODEL_DIR)


_unzip_models()
versions = (1.04, 1.02, 1.4, 1.7)
SubChandra_1 = _get_model(version=1.04)
SubChandra_2 = _get_model(version=1.02)
Chandra = _get_model(version=1.4)
SuperChandra = _get_model(version=1.7)
