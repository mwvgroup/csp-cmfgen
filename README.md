# Supernova Model Comparison in the NIR

A comparison of [CMFGEN](http://kookaburra.phyast.pitt.edu/hillier/web/CMFGEN.htm) supernova simulations against data from the Carnegie Supernova Project (CSP), particularly in the NIR. 



## Data Access

CSP data is provided by the `csp` module. Any data that is not available on your local machine is downloaded automatically.

```python
import csp

# A table of target data from the third CSP data release
print(csp.master_table)

# The names and effective wavelengths for the CSP bandpasses
bands = csp.band_names
lambda_eff = csp.lambda_effective

# Get an astropy table with photometric data for a specific target
demo_data = csp.get_data_for_id('2004dt')
print(demo_data)

# You can also get the same data formatted for use with SNCosmo
demo_data_formatted = csp.get_data_for_id('2004dt')
print(demo_data_formatted)

# Iterate over photometric data tables for all targets
for data in csp.iter_sncosmo_input():
    print(data)
    break
```



## Model Access

A copy of the CMFGEN models in ascii format is provided in the *asccii_models/* directory. A version of these models ported for use with **SNCosmo** is available via the `sncosmo_models` module. It is important to note that these two models are not identical. The models represent the flux of a supernova at 1Kpc for a given time and wavelength. The wavelength grid in the asccii models changes with phase, where as in the version ported to SNCosmo the flux for a given phase has been interpolated onto a common wavelength grid. 



To retrieve a model:

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



To fit a light curve:

```Python
import sncosmo
import matplotlib

# create a model
model = sncosmo.Model(source=source)

# run the fit
data = sncosmo.load_example_data()
result, fitted_model = sncosmo.fit_lc(
    data, model,
    ['z', 't0', 'x0'],  # parameters of model to vary
    bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)

# Plot results
fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
fig.show()
```



