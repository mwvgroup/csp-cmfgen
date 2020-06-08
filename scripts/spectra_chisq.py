#!/usr/bin/env python
# coding: utf-8

"""Calculate chi-squared values for each CMFGEN model against the CSP DR1
dataset. Determine results for each bandpass and spectroscopic feature.
"""

import sys
from pathlib import Path

from sndata.csp import DR1, DR3
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from analysis import (
    equivalent_width,
    models,
    spectra_chisq
)

dr1 = DR1()
dr1.download_module_data()

dr3 = DR3()
dr3.download_module_data()
dr3.register_filters()

# Runtime settings
out_dir = Path(__file__).resolve().parent / 'results'
model_list = models.get_models_with_ext()
feature_list = equivalent_width.features
band_list = dr3.band_names

if __name__ == '__main__':
    out_dir.mkdir(exist_ok=True, parents=True)

    tqdm.write('Calculating chi-squared for bands')
    spectra_chisq.tabulate_chisq(
        data_release=dr1,
        bands=band_list,
        models=model_list,
        out_path=out_dir / 'band_chisq.ecsv'
    )

    tqdm.write('\nCalculating chi-squared for features')
    spectra_chisq.tabulate_chisq(
        data_release=dr1,
        features=feature_list,
        models=model_list,
        out_path=out_dir / 'feature_chisq.ecsv'
    )
