# Supernova Model Comparison in the NIR

A comparison of [CMFGEN](http://kookaburra.phyast.pitt.edu/hillier/web/CMFGEN.htm) supernova simulations against data from the Carnegie Supernova Project (CSP),
particularly in the NIR. 



## Data Access

Access to CSP light curves is provided by the `sndata` package. Full documentation is available [here](https://sn-data.readthedocs.io/en/latest/index.html).




## Model Access

A copy of the CMFGEN models in ascii format is provided in the *asccii_models/* directory. A version of these models ported for use with **SNCosmo** is available via the `sncosmo_models` module. It is important to note that these two models are not identical. The models represent the flux of a supernova at 1Kpc for a given time and wavelength. The wavelength grid in the asccii models changes with phase, where as in the version ported to SNCosmo the flux for a given phase has been interpolated onto a common wavelength grid. 



#### Using SNCosmo Example:

```Python

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
```
