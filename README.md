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
