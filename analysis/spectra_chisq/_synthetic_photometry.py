"""The ``_synthetic_photometry`` module creates a table with synthetic magnitudes and regressed magnitudes.

Synthetic magnitudes are calculated using CSP DR1 spectra (flux per wavelength; erg/s/cm^2/Angstrom).
Photometric magnitudes (DR3) are regressed to the time of the spectra.

Data release papers: https://sn-data.readthedocs.io/en/latest/module_docs/csp.html
"""

import numpy as np
import pandas as pd
import sncosmo
import sndata
from astropy import constants as const
from astropy import units as u
from astropy.table import join, QTable, Table, unique, vstack
from matplotlib import pyplot as plt
import sys

from .. import lc_colors, utils

major, minor, _ = sndata.__version__.split('.')
if int(major) == 0 and int(minor) < 7:
    raise RuntimeError('Update to sndata 0.7.0 or higher')


def calc_synthetic_ab_mag(wavelength, flux_per_wave, transmission, err_ratio=0.03):
    """Calculate synthetic AB magnitude for a single spectra.

    Synthetic AB magnitude calculated using Eq 6 of Casagrande & VandenBerg 2014
    'Synthetic Stellar Photometry I. General considerations and new transformations for broad-band systems.'
    arXiv:1407.6095v2

    Args:
        wavelength (1d array-like): Wavelengths at which flux is observed. CSP DR1 units are Angstrom. Units must be
                                    consistent with `flux_per_wave`.
        flux_per_wave (1d array-like): Flux per wavelength. Units: erg/s/cm^2/Angstrom if units for `wavelength` are
                                        Angstrom.
        transmission (1d array-like): Transmission function. Transmission at each `wavelength` for specific band. Unitless
        err_ratio (float): The fraction of the flux to take as the error

    Returns:
        The synthetic AB magnitude
    """

    # Constant defined in Casagrande & VandenBerg 2014 (below Eq 4)
    flux_freq0 = 3.631 * 10 ** -20 * u.erg / u.s / u.cm / u.cm / u.Hz  # erg/s/cm^2/Hz

    # Add and test units
    # TODO: remove units via .value?
    if not wavelength.unit:
        wavelength.unit = u.angstrom
    if not flux_per_wave.unit:
        flux_per_wave.unit = u.erg / u.s / u.cm / u.cm / u.angstrom

    # Integrate dimension-full quantities (no errors about incompatible units)
    # AB mag in Casagrande & VandenBerg 2014 (Eq 6)
    numerator = np.trapz(wavelength * flux_per_wave * transmission, wavelength)
    denominator = np.trapz(const.c.cgs * flux_freq0 * transmission / wavelength, wavelength)
    mag_ab = -2.5 * np.log10(numerator / denominator)

    # TODO propagate error in flux to magnitude error. Try with scipy?

    # TODO handle blank cells

    return mag_ab


def tabulate_synthetic_photometry(spec_data, bands, zp_cgs=-48.60, zp_jy=8.9, err_ratio=0.03):
    """Create table with synthetic flux and magnitude for a single obj_id.

    Output table columns:
        - obj_id
        - time
        - band
        - syn_mag
        - syn_flux (Jy), is flux per frequency
        - err_msg TODO

    Args:
        spec_data (table): Table of spectroscopic data from sndata
        bands (list[str]): The name of the bands to tabulate photometry for
        err_ratio (float): The fraction of the flux to take as the error
        zp_cgs (float): Zero point for AB magnitude with F_nu in cgs units.
                        Given below Eq 4 of Casagrande & VandenBerg 2014
        zp_jy (float): Zero point for AB magnitude with F_nu in units of Jy

    Returns:
        An astropy Table with synthetic photometry results
    """

    # FIXME uncomment?
    #spec_data = spec_data.copy()

    # Create quantity table to correctly handle squaring units
    # http://docs.astropy.org/en/latest/table/mixin_columns.html#quantity-and-qtable
    spec_data = QTable(spec_data)

    # Create output table
    out_table = Table(
        names=['obj_id', 'time', 'band', 'synth_mag', 'synth_flux_Jy'],
        dtype=['U100', float, 'U100', float, float]  # , masked=True # FIXME ?
    )

    # Get SN object ID
    obj_id = spec_data.meta['obj_id']

    # Add appropriate units if necessary
    if not spec_data['wavelength'].unit:
        spec_data["wavelength"].unit = u.angstrom
    # 'flux' gives flux per wavelength, flux_lambda
    if not spec_data['flux'].unit:
        spec_data['flux'].unit = u.erg / u.s / u.cm / u.cm / u.angstrom

    # Assume flux error
    spec_data['fluxerr'] = err_ratio * spec_data['flux']

    # Fill output table
    for band in bands:
        for time in set(spec_data['time']):
            # Get a single spectrum
            spectrum = spec_data[spec_data['time'] == time]
            # Get the transmission function for `band`
            transmission = sncosmo.get_bandpass(band)(spectrum['wavelength'])
            # Get synthetic AB magnitude
            syn_mag = calc_synthetic_ab_mag(wavelength=spectrum['wavelength'], flux_per_wave=spectrum['flux'],
                                            transmission=transmission)
            # Convert magnitude to flux
            syn_flux_per_freq = utils.mag_to_flux(mag=syn_mag, zp=zp_jy)  # units: Jy
            # Add data to output table
            out_table.add_row([obj_id, time, band, syn_mag, syn_flux_per_freq])

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
        names=['obj_id', 'time', 'band', 'photo_mag'],
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
            spec_times = id_spec_data[id_spec_data['band'] == band]['time']

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


def tabulate_photometry(spec_release, phot_release, err_ratio=0.03, out_path=None):
    """Tabulates synthetic and observed photometry

    Observed photometry is interpolated to spectral times using a gaussian
    regression.

    Args:
        spec_release (module): An spectroscopic data release from sndata
        phot_release (module): A photometric data release from sndata
        out_path        (str): Optionally write results to file
        err_ratio (float): The fraction of the flux to take as the error

    Returns:
        A Table of synthetic and observed photometry
    """

    # Get synthetic photometry for all spectroscopic data
    spec_data = spec_release.iter_data(verbose=True)  # FIXME is module sndata.csp.dr1
    bands = phot_release.band_names
    tables = [tabulate_synthetic_photometry(tbl, bands, err_ratio) for tbl in spec_data]
    # tables = [calc_synthetic_phot(tbl, bands, err_ratio) for tbl in spec_data]
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
            syn_band_data['time'],
            syn_band_data['synth_mag'],
            c=f'C{i}',
            marker='s',
            label=band.split('_')[-1]
        )

    axis.set_title(obj_id)
    axis.invert_yaxis()
    axis.legend()
    return fig, axis


def get_synth_photo_for_id(obj_id, spec_release, phot_release, time, err_ratio=0.03, temp_flux=None, temp_spec_tbl=None, out_path=None):
    """Get synthetic magnitudes for one object id using the flux from sncosmo model with extinction.
    This function is used with ``_model_spectra`` module to fit extinction parameters using synthetic photometry and
    observed photometry.

    Args:

    Returns:
    """

    ''' Do the following once: get photo_mag (well it will be in the final table anyway...)'''


    if temp_flux is None:
        # Get synthetic photometry for object id
        spec_data = spec_release.get_data_for_id(obj_id).copy()
        #spec_data = unique(spec_data, keys='time')

    # TODO only one band is used in temp_flux.... change this. Want one fit for all bands?
    if temp_flux is not None and temp_spec_tbl is not None:
        # Add row with fitted flux so this can pass through synth_phot
        # Must have same dimensions
        spec_data=temp_spec_tbl
        #print('------->temp spectra table: ', len(temp_spec_tbl))
        #print('------->len temp flux: ', len(temp_flux))
        #spec_data['flux'] = temp_flux##FIXME?

    bands = phot_release.band_names
    ##tables = [tabulate_synthetic_photometry(tbl, bands, err_ratio) for tbl in spec_data]
    # tables = [calc_synthetic_phot(tbl, bands, err_ratio) for tbl in spec_data]
    synthetic_photometry = tabulate_synthetic_photometry(spec_data, bands, err_ratio)

    # Regress photometric data
    obs_photometry = _photometry_to_spectra_time(
        synthetic_photometry, phot_release)
    # TODO
    # Mask observed photometry by `time`... looks like it's not already masked by spec_data?
    if time is not None:
        tables = [obs_photometry[obs_photometry['time'] == t] for t in time]
        obs_photometry = vstack(tables)
    #print(obs_photometry)

    # Create output table
    combined = join(synthetic_photometry, obs_photometry, join_type='left')

    return combined