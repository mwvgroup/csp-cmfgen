import argparse
import warnings
from pathlib import Path

import sncosmo
from sndata.csp import dr1, dr3
from tqdm import tqdm

from analysis import equivalent_width
from analysis import lc_colors
from analysis import models
from analysis import spectra_chisq

warnings.filterwarnings('ignore')
models.register_sources()
dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters()

BAND_COMBOS = [
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


def run_color_15(cli_args):
    """Tabulate change in color over 15 days
    Args:
        cli_args (Namespace): Command line arguments
    """

    out_dir = Path(cli_args.out_dir)

    tqdm.write('Tabulating delta color 15')
    lc_colors.tabulate_delta_15(
        data_release=dr3,
        models=get_models(cli_args.models),
        band_combos=BAND_COMBOS,
        out_path=out_dir / f'delta_c_15.ecsv')

    tqdm.write('\n')


def run_color_chisq(cli_args):
    """Tabulate chi-squares for color evolution
    Args:
        cli_args (Namespace): Command line arguments
    """

    out_dir = Path(cli_args.out_dir) / 'color_chisq'

    out_dir.mkdir(exist_ok=True, parents=True)
    if (cli_args.start is None) or (cli_args.end is None):
        trange = None
        out_path = out_dir / 'no_limit.ecsv'

    else:
        trange = (cli_args.start, cli_args.end)
        out_path = out_dir / f'{trange[0]}_{trange[1]}.ecsv'.replace('-', 'n')

    tqdm.write('Tabulating color chisq (phase range = {})'.format(trange))
    lc_colors.tabulate_chisq(
        data_release=dr3,
        models=get_models(cli_args.models),
        band_combos=BAND_COMBOS,
        prange=trange,
        out_path=out_path)

    tqdm.write('\n')


def run_ew(cli_args):
    """Run equivalent width analysis
    Args:
        cli_args (Namespace): Command line arguments
    """

    model_list = get_models(cli_args.models)
    out_dir = Path(cli_args.out_dir) / 'equivalent_width'
    out_dir.mkdir(parents=True, exist_ok=True)

    tqdm.write('Tabulating peak model pew')
    peak_pew = equivalent_width.tabulate_peak_model_pew(model_list)
    peak_pew.write(out_dir / f'peak_model_pew.ecsv', overwrite=True)

    fix_boundaries = bool(cli_args.fix_boundaries)
    start_msg = 'Tabulating equivalent widths (Fixed bounds = {})'
    tqdm.write(start_msg.format(fix_boundaries))

    ew_results = equivalent_width.tabulate_pew_spectra(
        data_release=dr1,
        models=get_models(cli_args.models),
        fix_boundaries=fix_boundaries)

    out_path = out_dir / f'fixed_{fix_boundaries}.ecsv'
    ew_results.write(out_path, overwrite=True)
    tqdm.write('\n')


def run_spec_chisq(cli_args):
    """Calculate chi-squared for spectra

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    out_dir = Path(cli_args.out_dir) / 'spec_chisq'
    out_dir.mkdir(parents=True, exist_ok=True)

    file_name = 'bands.ecsv'
    if cli_args.features is not None:
        file_name = 'features.ecsv'
        cli_args.features = \
            {f: equivalent_width.features[f] for f in cli_args.features}

    tqdm.write('Calculating chi-squared for spectra')
    spectra_chisq.tabulate_chisq(
        data_release=dr1,
        bands=cli_args.bands,
        features=cli_args.features,
        models=get_models(cli_args.models),
        err_estimate=cli_args.err_ratio,
        trans_limit=cli_args.trans_limit,
        out_path=out_dir / file_name
    )
    tqdm.write('\n')


def run_synthetic_photometry(cli_args):
    """Tabulate synthetic photometry

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    out_path = Path(cli_args.out_dir) / 'synth_phot.ecsv'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tqdm.write('Tabulating synthetic photometry')
    spectra_chisq.tabulate_photometry(
        dr1, dr3, cli_args.err_ratio, out_path)

    tqdm.write('\n')


def create_parser():
    """Create a command line argument parser

    Returns:
        An ``argparse.ArgumentParser`` object
    """

    parser = argparse.ArgumentParser(
        description='Compare various SN models against CSP data.')
    subparsers = parser.add_subparsers(help='')

    parser.add_argument(
        '-o', '--out_dir',
        type=str,
        default=['./'],
        help='Output directory')

    # For tabulating chi-squared of color evolution
    color_chisq_parser = subparsers.add_parser(
        'color_chisq', help='Compare color evolution with models.')

    color_chisq_parser.set_defaults(func=run_color_chisq)
    color_chisq_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    color_chisq_parser.add_argument(
        '-i', '--interval',
        default=1,
        type=int,
        help='Spacing between phases when summing chisq')

    color_chisq_parser.add_argument(
        '-s', '--start',
        type=float,
        help='Start of phase range to use in chi-squared integration')

    color_chisq_parser.add_argument(
        '-e', '--end',
        type=float,
        help='End of phase range to use in chi-squared integration')

    # For tabulating change in color over 15 days
    color_15_parser = subparsers.add_parser(
        'color_15', help='Compare color evolution with models.')

    color_15_parser.set_defaults(func=run_color_15)
    color_15_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    color_15_parser.add_argument(
        '-b', '--t0_band',
        type=str,
        help='Band to use when setting model t0 to peak')

    # For tabulating pseudo equivalent width values
    ew_parser = subparsers.add_parser(
        'equivalent_width', help='Calculate pseudo equivalent width values')

    ew_parser.set_defaults(func=run_ew)
    ew_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    ew_parser.add_argument(
        '-b', '--fix_boundaries',
        required=True,
        type=int,
        help='Fix feature boundaries to observed values')

    spec_chisq_parser = subparsers.add_parser(
        'spec_chisq', help='Calculate chi-squared for spectra')

    spec_chisq_parser.set_defaults(func=run_spec_chisq)
    spec_chisq_parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        required=True,
        help='Models to use')

    spec_chisq_parser.add_argument(
        '-e', '--err_ratio',
        type=float,
        default=.03,
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
        default=.1,
        help='Transmission cutoff applied to each band')

    synth_phot_parser = subparsers.add_parser(
        'synth_phot', help='Tabulate synthetic photometry')

    synth_phot_parser.set_defaults(func=run_synthetic_photometry)

    synth_phot_parser.add_argument(
        '-e', '--err_ratio',
        type=float,
        default=.03,
        help='Error estimate as a fraction of the flux')

    return parser


if __name__ == '__main__':
    # Parse command line input
    parser = create_parser()
    cli_args = parser.parse_args()
    cli_args.func(cli_args)
