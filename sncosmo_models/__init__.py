#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines custom SNCosmo.Source objects for four different CMFGEN
based models. Source objects include:

SubChandra_1: A Sub-Chandrasekhar model for 1.04 solar mass progenitors
SubChandra_2: A Sub-Chandrasekhar model for 1.02 solar mass progenitors
Chandra: A Chandrasekhar model for 1.4 solar mass progenitors
SuperChandra: A Super-Chandrasekhar model for 1.7 solar mass progenitors

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

from ._models import ASCII_MODEL_DIR as _ascii_dir
from ._models import format_models, get_model, register_sources

format_models(force=False)
SubChandra_1 = get_model(version=1.04)
SubChandra_2 = get_model(version=1.02)
Chandra = get_model(version=1.4)
SuperChandra = get_model(version=1.7)
