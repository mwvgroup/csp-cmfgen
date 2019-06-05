# Supernova Model Comparison in the NIR

A comparison of [CMFGEN](http://kookaburra.phyast.pitt.edu/hillier/web/CMFGEN.htm) supernova simulations against data from the Carnegie Supernova Project (CSP),
particularly in the NIR. 



## Data Access

Access to CSP light curves is provided by the `SNData` package. Full documentation is available [here](https://sn-data.readthedocs.io/en/latest/index.html).




## Model Access

A copy of the CMFGEN models in ascii format is provided in the *asccii_models/* directory. A version of these models ported for use with **SNCosmo** is available via the `sncosmo_models` module. It is important to note that these two models are not identical. The models represent the flux of a supernova at 1Kpc for a given time and wavelength. The wavelength grid in the asccii models changes with phase, where as in the version ported to SNCosmo the flux for a given phase has been interpolated onto a common wavelength grid. 



#### To retrieve a model:

```Python
from sncosmo_models import Chandra

source = Chandra()

# The model interpolated onto a common wavelength grid
# This is the information used when fitting light-curves with SNCosmo
phase_grid, wavelength_grid, flux_grid = source.gridded_model()

# The model un-interpolated
# This information is provided for convenience, and is not used
# by SNCosmo in any way.
phase, wavelength, flux = source.raw_model()
```



#### To fit a light curve:

```Python
import sncosmo
from SNData.csp import dr3

# Make sure data is available
dr3.download_module_data()
dr3.register_filters()

# Get a data table
demo_table = next(dr3.iter_data(format_sncosmo=True))

# create a model
model = sncosmo.Model(source=source)
model.set(z=demo_table.meta['redshift']

# run the fit
result, fitted_model = sncosmo.fit_lc(
    data=demo_table,
    model=model,
    vparam_names=['t0', 'x0'])

# Plot results
fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
fig.show()
```



