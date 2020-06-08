#!/usr/bin/env python
# coding: utf-8

"""Runs light-curve fits using CMFGEN models on CSP targets. Fits are
performed on a band-by-band basis.
"""

import sys
from pathlib import Path

import yaml
from sndata.csp import DR3
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from analysis import band_fitting, models

dr3 = DR3()
dr3.download_module_data()
dr3.register_filters()

# Runtime settings
base_dir = Path(__file__).resolve().parent
out_dir = base_dir / 'results'
config_path = base_dir.parent / 'config' / 'csp_config.yml'
out_path = out_dir / 'band_fits.ecsv'
model_list = models.get_models_with_ext()

if __name__ == '__main__':
    out_dir.mkdir(exist_ok=True, parents=True)

    with config_path.open() as infile:
        config = yaml.load(infile, Loader=yaml.FullLoader)

    band_fitting.tabulate_band_fits(
        data_release=dr3,
        models=model_list,
        config=config,
        out_path=out_path
    )

    tqdm.write('\n')
