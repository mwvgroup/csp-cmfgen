#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis package."""

import argparse
import warnings
from pathlib import Path

import sncosmo
from sndata.csp import dr1, dr3
from tqdm import tqdm

from analysis import equivalent_width, lc_colors, models, spectra_chisq

warnings.filterwarnings('ignore')
models.register_sources()
dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters()


def get_models(model_id_list):
    """Return a list of models corresponding to a list of model ids

    Model ids should be of the form ``<source name>,<source version>

    Args:
        model_id_list (list): A list of model ids

    Returns:
        A list of sncosmo Models
    """

    model_list = []
    for model_id in model_id_list:
        name, version = model_id.split(',')
        model = sncosmo.Model(sncosmo.get_source(name, version=version))
        model.add_effect(sncosmo.F99Dust(), 'ext', 'rest')
        model_list.append(model)

    return model_list


def run_lc_color(cli_args):
    """Run color analysis

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    unique_bands = [
        'csp_dr3_u',
        'csp_dr3_g',
        'csp_dr3_r',
        'csp_dr3_i',
        'csp_dr3_B',
        'csp_dr3_V',
        'csp_dr3_Y',
        'csp_dr3_J'
        'csp_dr3_H',
    ]

    color_combos = [(unique_bands[i], unique_bands[i + 1]) for i in range(len(unique_bands) - 1)]
    out_dir = Path(cli_args.out_dir) / 'color_evolution'
    tqdm.write('Tabulating color evolution')
    for model in get_models(cli_args.models):
        tqdm.write(f'Fitting {model.source.name} {model.source.version}')
        file_name = f'{model.source.name}_{model.source.version}.ecsv'
        t = lc_colors.tabulate_residuals(dr3, model, color_combos)
        t.write(out_dir / file_name, overwrite=True)
        tqdm.write('\n')


def run_ew(cli_args):
    """Run equivalent width analysis

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    models = get_models(cli_args.models)
    out_dir = Path(cli_args.out_dir) / 'equivalent_width'
    out_dir.mkdir(parents=True, exist_ok=True)

    tqdm.write('Tabulating peak model pew...')
    peak_pew = equivalent_width.tabulate_peak_model_pew(models)
    peak_pew.write(out_dir / f'peak_model_pew.ecsv', overwrite=True)

    for fix_boundaries in (True, False):
        start_msg = 'Tabulating equivalent widths (Fixed bounds = {})'
        tqdm.write(start_msg.format(fix_boundaries))

        ew_results = equivalent_width.tabulate_pew_spectra(
            data_release=dr1,
            models=get_models(cli_args.models),
            fix_boundaries=fix_boundaries)

        out_path = out_dir / f'fixed_{fix_boundaries}.ecsv'
        ew_results.write(out_path, overwrite=True)


def run_spec_chisq(cli_args):
    """Calculate chi-squared for spectra

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    out_dir = Path(cli_args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    spectra_chisq.tabulate_chi_squared(
        data_release=dr1,
        bands=dr3.band_names,
        models=get_models(cli_args.models),
        err_estimate=cli_args.err_estimate,
        trans_limit=cli_args.trans_limit,
        out_path=out_dir / 'spec_chisq.ecsv'
    )


# Parse command line input
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare various SN models against CSP data.')
    subparsers = parser.add_subparsers(help='')

    parser.add_argument(
        '-o', '--out_dir',
        type=str,
        default=['./'],
        help='Output directory')

    color_parser = subparsers.add_parser('lc_color',
        help='Compare color evolution with models.')

    color_parser.set_defaults(func=run_lc_color)
    color_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    ew_parser = subparsers.add_parser('equivalent_width',
        help='Calculate pseudo equivalent width values')

    ew_parser.set_defaults(func=run_ew)
    ew_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    spec_chisq_parser = subparsers.add_parser('spec_chisq',
        help='Calculate chi-squared for spectra')

    spec_chisq_parser.set_defaults(func=run_spec_chisq)
    spec_chisq_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    spec_chisq_parser.add_argument(
        '-e', '--err_estimate',
        type=float,
        default=None,
        help='Error estimate as a fraction of the flux')

    spec_chisq_parser.add_argument(
        '-f', '--features',
        type=str,
        nargs='+',
        default=None,
        help='Features to tabulate chi-squared for.')

    spec_chisq_parser.add_argument(
        '-b', '--bands',
        type=str,
        nargs='+',
        default=None,
        help='Bands to tabulate chi-squared for.')

    spec_chisq_parser.add_argument(
        '-t', '--trans_limit',
        type=float,
        default=None,
        help='Transmission cutoff applied to each band')

    cli_args = parser.parse_args()
    cli_args.func(cli_args)
