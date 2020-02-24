# Supernova Model Comparison in the NIR

A comparison of [CMFGEN](http://kookaburra.phyast.pitt.edu/hillier/web/CMFGEN.htm)
supernova simulations against data from the Carnegie Supernova Project (CSP),
particularly in the NIR.

## Project Goals

Due to the increased observation time required for spectroscopic observations,
the number spectroscopically observed supernovae that are publicly available
is significantly smaller than those with photometric observations. This
imbalance is expected to worsen over the next decade with large scale
photometric surveys like the Large Synoptic Survey Telescope (LSST) coming
online in 2021. This is a particular challenge when modeling the astrophysical
channels that result in a Type Ia supernova (SNe Ia). Spectroscopic
observations are crucial in constraining SNe Ia model due to their increased
wavelength resolution.

We attempt to predict how well we will be able to constrain spectroscopic SNe
Ia models using the predominantly photometric data that will be available in
coming years. As a baseline, we use observations taken by the Carnegie
Supernova Project primarily due to its employment of simultaneous
spectroscopic and photometric observing campaigns. This includes spectroscopic
data published in the first CSP data release (DR1) and photometric observations
published in the third (DR3).

Our analysis is broken down into four key steps: two spectroscopic and one
photometric. We start by comparing modeled spectra against flux values from
DR1 observations. This serves as a baseline for understanding how well the
models we consider represent our observed targets on a physical level. Second,
we supplement this baseline by comparing the predicted pseudo equivalent
widths of different features considered in the DR1 release paper. Thirdly,
we perform light curve fits of DR3 photometry and compare our ability to fit
in the optical vs the NIR. Finally, we consider how well the color evolution
of each model matches the photometric data.

## **Running the Analysis Pipeline**

The analysis pipeline for this project can be accessed either from within a
python environment or from the command line. The python interface recommended
when trying to extend or reproduce individual steps of the analysis. For usage
examples see the *notebooks/* directory. The command line interface is intended
for running complete sections of the analysis. For more instructions see
`python run_analysis.py -h`.

## Notes on the Models

A copy of the CMFGEN models in ascii format is provided in the
*ascii_models/* directory. A version of these models ported for use with
**SNCosmo** is available via the `analysis.models` module (see the
[docs](https://mwvgroup.github.io/nir-comparison/build/html/api_ref/models.html)
for more info). It is important to note that these two models are not
identical. The wavelength grid in the asccii models changes with phase, where
as the SNCosmo models have been interpolated onto a common wavelength grid.
Both models represent the flux of a supernova at 1 kpc for a given time and
wavelength.
