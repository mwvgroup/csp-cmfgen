"""Create table with synthetic magnitudes and regressed magnitudes.

See Figure 2 of Folatelli et al 2013: Spectroscopy of Type 1a Supernova by
the Carnegie Supernova Project.

Data release papers: https://sn-data.readthedocs.io/en/latest/module_docs/csp.html
"""

import numpy as np
import sncosmo
import sndata
from astropy.table import Table, join, vstack
from matplotlib import pyplot as plt

from .. import lc_colors, utils

major, minor, _ = sndata.__version__.split('.')
if int(major) == 0 and int(minor) < 7:
    raise RuntimeError('Update to sndata 0.7.0 or higher')


def calc_synthetic_phot(spec_data, bands, err_ratio=.03):
    """Get synthetic magnitude in dr3 band passes for a table of spectra

    Output table columns:
        - obj_id
        - date
        - band
        - syn_mag
        - err_msg
    
    Args:
        spec_data (Table): Table of spectroscopic data from sndata
        bands (list[str]): The name of the bands to tabulate photometry for
        err_ratio (float): The fraction of the flux to take as the error

    Returns:
        An astropy Table synthetic photometry results
    """

    spec_data = spec_data.copy()
    out_table = Table(
        names=['obj_id', 'date', 'band', 'syn_mag', 'err_msg'],
        dtype=['U100', float, 'U100', float, 'U1000'],
        masked=True
    )

    # Add appropriate units and assume flux error is 3% of flux
    spec_data['fluxerr'] = err_ratio * spec_data['flux']
    obj_id = spec_data.meta['obj_id']

    for band in bands:
        for date in set(spec_data['date']):
            nrepeat = 5
            spectrum = spec_data[spec_data['date'] == date]
            flux = np.atleast_2d(spectrum['flux'])

            try:
                # We create a dummy model that represents our spectrum and use
                # it to determine the band magnitude
                sncosmo_source = sncosmo.TimeSeriesSource(
                    phase=np.linspace(-50., 50., nrepeat),
                    wave=spectrum['wavelength'],
                    flux=np.repeat(flux, nrepeat, axis=0)
                )

                sncmag = sncosmo.Model(sncosmo_source).bandmag(band, 'AB', 0)

            except ValueError as e:
                row = [obj_id, date, -99, -99, str(e).replace('\n', ' ')]
                mask = [False, False, True, True, False]
                out_table.add_row(row, mask=mask)
                continue

            out_table.add_row([obj_id, date, band, sncmag, ''])

    return out_table


def _photometry_to_spectra_time(spec_data, phot_release):
    """Regress photometry (CSP DR3) to the time of spectroscopic observations

    See ``tabulate_synthetic_photo()`` for the expected data model of
    ``spec_data``.

    Args:
        spec_data (table): Table of synthetic photometry from

    Returns:
        A Table with regressed magnitudes
    """

    out_table = Table(
        names=['obj_id', 'date', 'band', 'photo_mag'],
        dtype=['U100', float, 'U100', float], masked=True
    )

    for obj_id in set(spec_data['obj_id']):
        # Get spectroscopic and photometric data for particular object ID
        id_spec_data = spec_data[spec_data['obj_id'] == obj_id]
        id_photo_data = phot_release.get_data_for_id(obj_id)

        # Fit photometric data with a gaussian process
        gauss_process = lc_colors.fit_gaussian_process(data=id_photo_data)

        for band in phot_release.band_names:

            # Get times of spectra obs
            spec_times = id_spec_data[id_spec_data['band'] == band]['date']

            if band in id_photo_data['band']:
                # Interpolate band data to the spectral times
                pred_flux, pred_flux_err = lc_colors.predict_band_flux(
                    gp=gauss_process, band_name=band, times=spec_times)

                # Convert flux to magnitude
                zp = sndata.get_zp(band_name=band)
                pred_mag = utils.flux_to_mag(flux=pred_flux, zp=zp)

                for i, mag in enumerate(pred_mag):
                    out_table.add_row([obj_id, spec_times[i], band, mag])

            else:
                out_table.add_row([obj_id, 0.0, band, np.nan])

    return out_table


def tabulate_synthetic_photometry(
        spec_release, phot_release, err_ratio=0.3, out_path=None):
    """Tabulates synthetic and observed photometry

    Observed photometry is interpolated to spectral times using a gausian
    regression.

    Args:
        spec_release (module): An spectroscopic data release from sndata
        phot_release (module): A photometric data release from sndata
        out_path        (str): Optionally write results to file

    Returns:
        A Table of synthetic and observed photometry
    """

    # Get synthetic photometry for all spectroscopic data
    spec_data = spec_release.iter_data(verbose=True)
    bands = phot_release.band_names
    tables = [calc_synthetic_phot(tbl, bands, err_ratio) for tbl in spec_data]
    synthetic_photometry = vstack(tables)

    # Regress all photometric data
    obs_photometry = _photometry_to_spectra_time(
        synthetic_photometry, phot_release)

    # Create output table
    combined = join(synthetic_photometry, obs_photometry, join_type='left')
    if out_path:
        combined.write(out_path, overwrite=True)

    return combined


# Todo: Move into a jupyter notebook
def plot_synthetic_results(observed, synthetic):
    """Plot synthetic vs observed magnitudes

    Args:
        observed  (Table): A table of observed photometric for a single object
        synthetic (Table): A table of synthetic photometry for a single object

    Returns:
        A matplotlib figure
        A matplotlib axis
    """

    obj_id = observed.meta['obj_id']
    bands = set(synthetic['band']).intersection(set(observed['band']))
    fig, axis = plt.subplots(1, 1, figsize=(8, 8))
    for i, band in enumerate(bands):
        obj_band_data = observed[observed['band'] == band]
        syn_band_data = synthetic[synthetic['band'] == band]

        axis.scatter(
            obj_band_data['time'],
            utils.flux_to_mag(obj_band_data['flux'], sndata.get_zp(band)),
            label='',
            c=f'C{i}',
            s=5
        )

        axis.plot(
            obj_band_data['time'],
            utils.flux_to_mag(obj_band_data['flux'], sndata.get_zp(band)),
            c=f'C{i}',
            linestyle='--',
            label=''
        )

        axis.scatter(
            syn_band_data['date'],
            syn_band_data['syn_mag'],
            c=f'C{i}',
            marker='s',
            label=band.split('_')[-1]
        )

    axis.set_title(obj_id)
    axis.invert_yaxis()
    axis.legend()
    return fig, axis
