#!/usr/bin/env python
# coding: utf-8

"""Runs light-curve fits using CMFGEN models on CSP targets. Fits are
performed on a band-by-band basis.
"""

from pathlib import Path

import sncosmo
import yaml
from sndata.csp import DR3
from tqdm import tqdm

from analysis import band_fitting
from analysis import models

# Configure data access
dr3 = DR3()
dr3.download_module_data()
dr3.register_filters()

# I/O settings and configuration data
base_dir = Path(__file__).resolve().parent
out_path = base_dir / 'results' / 'band_fits.ecsv'
out_path.parent.mkdir(parents=True, exist_ok=True)

config_path = base_dir.parent / 'config' / 'config.yml'
with config_path.open() as infile:
    config = yaml.load(infile, Loader=yaml.FullLoader)

# Instantiate models for fitting
sources = [models.SubChandra_1, models.SubChandra_2, models.Chandra, models.SuperChandra]
models = []
for source in sources:
    model = sncosmo.Model(source)
    model.add_effect(sncosmo.F99Dust(), 'mw', 'obs')
    models.append(model)

tqdm.write('Fitting band-passes')
band_fitting.tabulate_band_fits(
    data_release=dr3,
    models=models,
    config=config,
    out_path=out_path
)
