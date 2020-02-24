#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""The ``regression`` module fits light curves using a gaussian regression and
compares the fitted colors with modeled color values.

Usage Example
-------------

.. code-block:: python
   :linenos:

   from sndata.csp import DR3

   from analysis import regression

   # Download demo data
   dr3 = DR3()
   dr3.download_module_data()

   # Load demo data as an astropy table
   demo_id = dr3.get_available_ids()[0]
   demo_data = dr3.get_data_for_id(demo_id)
   print(demo_data)

   # Fit a Gaussian process to the data
   gp = lc_colors.fit_gaussian_process(demo_data)
   gp_flux, gp_unc = regression.predict_light_curve(gp, bands, time)

   # Pick the bands and time values to evaluate the Gaussian process at
   all_bands = list(set(demo_data['band']))
   time = np.arange(min(demo_data['time']) - 1, max(demo_data['time']) + 1)

   # Interpolate flux in a single band
   band_flux, band_err = regression.predict_band_flux(gp, all_bands[0], time)

   # Interpolate flux in multiple bands
   lc_flux, lc_err = regression.predict_light_curve(gp, all_bands, time)

   # Predict the color at various time values
   band1 = all_bands[0]
   band2 = all_bands[1]
   color, color_err = predict_color(gp, time, band1, band2)

Function Documentation
----------------------
"""

from ._lc_prediction import (
    predict_band_flux,
    predict_c_15,
    predict_color,
    predict_light_curve
)
from ._lc_regression import fit_gaussian_process

__all__ = [
    'predict_band_flux',
    'predict_c_15',
    'predict_color',
    'predict_light_curve',
    'fit_gaussian_process'
]
