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
    import matplotlib

    # create a model
    model = sncosmo.Model(source=Chandra())

    # run the fit
    data = sncosmo.load_example_data()
    result, fitted_model = sncosmo.fit_lc(
        data, model,
        ['z', 't0', 'x0'],  # parameters of model to vary
        bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)

    # Plot results
    fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
"""

from ._models import Chandra, SubChandra_1, SubChandra_2, SuperChandra
