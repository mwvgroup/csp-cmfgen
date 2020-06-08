#!/usr/bin/env python
# coding: utf-8

"""Compare modeled and observed color values between CSP and CMFGEN."""

import sys
from pathlib import Path

from sndata.csp import DR1, DR3
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from analysis import lc_colors, models

models.register_sources()

dr1 = DR1()
dr1.download_module_data()

dr3 = DR3()
dr3.download_module_data()
dr3.register_filters()

# Runtime settings
out_dir = Path(__file__).resolve().parent / 'results' / 'colors'
model_list = models.get_models_with_ext()
phase_ranges = (None, (-10, 20), (-5, 20), (5, 50))

colors = [
    ('csp_dr3_u', 'csp_dr3_g'),
    ('csp_dr3_g', 'csp_dr3_r'),
    ('csp_dr3_r', 'csp_dr3_i'),
    ('csp_dr3_B', 'csp_dr3_i'),
    ('csp_dr3_B', 'csp_dr3_V'),
    ('csp_dr3_B', 'csp_dr3_V0'),
    ('csp_dr3_B', 'csp_dr3_V1'),

    ('csp_dr3_Y', 'csp_dr3_J'), ('csp_dr3_Y', 'csp_dr3_Jdw'),
    ('csp_dr3_Ydw', 'csp_dr3_J'), ('csp_dr3_Ydw', 'csp_dr3_Jdw'),

    ('csp_dr3_J', 'csp_dr3_H'), ('csp_dr3_J', 'csp_dr3_Hdw'),
    ('csp_dr3_Jdw', 'csp_dr3_H'), ('csp_dr3_Jdw', 'csp_dr3_Hdw'),
]

if __name__ == '__main__':
    out_dir.mkdir(exist_ok=True, parents=True)

    lc_colors.tabulate_delta_15(
        data_release=dr3,
        models=model_list,
        band_combos=colors,
        out_path=out_dir / f'delta_c_15.ecsv')

    for trange in phase_ranges:
        tqdm.write('Tabulating color chisq (phase range = {})'.format(trange))

        if trange is None:
            file_name = 'no_limit.ecsv'

        else:
            file_name = f'{trange[0]}_{trange[1]}.ecsv'.replace('-', 'n')

        lc_colors.tabulate_chisq(
            data_release=dr3,
            models=model_list,
            colors=colors,
            prange=trange,
            out_path=out_dir / file_name)
