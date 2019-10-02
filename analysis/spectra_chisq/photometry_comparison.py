"""Create table with synthetic magnitudes and regressed magnitudes.
See Figure 2 of Folatelli et al 2013: Spectroscopy of Type 1a Supernova by the Carnegie Supernova Project.
Data release papers: https://sn-data.readthedocs.io/en/latest/module_docs/csp.html
"""

from astropy import constants as const
from astropy import units as u
from astropy.table import Column, join, QTable, Table, unique, vstack
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import sncosmo
import sndata
from sndata.csp import dr1, dr3
from .. import lc_colors
from .. import spectra_chisq
from .. import utils


def tabulate_synthetic_photometry(spec_data):
    """Get synthetic flux and magnitude.
    
    Args:
        spec_data: Table of spectroscopic data

    Returns:
        out_table: Table with synthetic magnitudes

    """
    
    # Create quantity table to correctly handle squaring units
    # http://docs.astropy.org/en/latest/table/mixin_columns.html#quantity-and-qtable
    spec_data = QTable(spec_data)
    
    out_table = Table(
        names=['obj_id', 'date', 'band', 'syn_mag'],
        dtype=['U100', float, 'U100', float]
    )

    obj_id = spec_data.meta['obj_id']
    
    # Add appropriate units
    spec_data['wavelength'] *= u.angstrom
    # Flux per wavelength
    spec_data['flux'] *= u.erg / u.s / u.cm / u.cm / u.angstrom
    # Assume flux error is 3% of flux
    spec_data['fluxerr'] = 0.03 * spec_data['flux']

    for band in dr3.band_names:
        for date in set(spec_data['date']):
            spectrum = spec_data[spec_data['date'] == date]
            transmission = sncosmo.get_bandpass(band)(spectrum['wavelength'])
            print('transmission', transmission)
            # Convert flux per wavelength to flux per frequency: F_freq = jacobian * F_wave
            jacobian = (spectrum['wavelength'] ** 2) / const.c 
            f_freq = spectrum['flux'] * jacobian
            # Integrate to get synthetic flux
            band_flux = np.sum(f_freq * transmission)
            # band_flux = np.trapz(f_freq * transmission, spectrum['wavelength']/u.angstrom)
            if band_flux > 0:
                mag = -2.5 * np.log10(band_flux.to(u.jansky).value) + 8.9 # 8.9 is zp
                out_table.add_row([obj_id, date, band, mag])
                print(mag)
                
    return out_table


def mag_to_flux(mag, zp, offset=0):
    """Convert magnitude to flux (F_nu).

    Args:
        mag: Photometric AB magnitude
        zp: Zeropoint of band
        offset:
    Returns:
        flux: F_nu; flux per frequency. Units: Jansky
    """

    flux = 10**(-0.4 * (mag - zp - offset))

    return flux


def flux_to_mag(flux, zp, offset=0):
    """Converts flux (F_nu) to magnitude.

    Args:
        flux: F_nu; flux per frequency. Units: Jansky
        zp: Zeropoint of band
        offset:
    Returns:
        mag: AB magnitude
    """
    
    mag = -2.5 * np.log10(flux) + zp + offset
    
    return mag


def get_offset_for_band(band):
    """Get offset from AB magnitude.
    https://csp.obs.carnegiescience.edu/data/filters

    Args:
        band (str):

    Returns:
        offset (float): Offset from AB magnitude
    """

    '''
    tables = pd.read_html("https://csp.obs.carnegiescience.edu/data/filters")
    print(tables)
    #print(tables[0])
    #name_col = tables[0][0] 
    #offset_col = tables[0][4]
    '''

    #TODO create this dict elsewhere?

    offsets = {
                "u" : -0.058, "g" : -0.017, "r" : -0.006, "i" : 0.001, "B" : -0.125, "V0" : -0.016, "V1" : -0.019,
                "V" : -0.02, "Y": 0.629, "J" : 0.912, "Jrc2" : 0.901, "H" : 1.343, "Ydw" : 0.629, "Jdw" : 0.912, "Hdw" : 1.342
            }
    offset = offsets.get(band.rsplit('csp_dr3_')[1])
    
    return offset


def photometry_to_spectra_time(spec_data):
    """Regress photometry (CSP DR3) to the time of spectroscopic (CSP DR1) observations

    Args:
        spec_data (table): astropy table from `tabulate_synthetic_photo()`
    Returns:
        out_table: Table with regressed magnitudes
    """
    
    # Get list of each object ID
    obj_ids = unique(spec_data, keys='obj_id')['obj_id']

    #TODO return column instead of table?
    out_table = Table(
        names=['obj_id', 'date', 'band', 'photo_mag'],
        dtype=['U100', float, 'U100', float], masked=True
    )

    out_table['photo_mag'].fill_value = -99.0
    
    for oid in obj_ids:
        # Get spectroscopic and photometric data for particular object ID
        id_spec_data = spec_data[spec_data['obj_id'] == oid]
        id_photo_data = dr3.get_data_for_id(oid, format_sncosmo=True)

        for band in dr3.band_names:
            # Filter spec_table by band
            # Get times of spectra obs
            spec_times = id_spec_data[id_spec_data['band'] == band]['date']
            photo_data = id_photo_data[id_photo_data['band'] == band] #TODO make sure obj id is the same
            # Test if photo data table is empty; TODO test with bool({list of column data}) ?
            if len(photo_data) >= 1:
                # Train photometric data on gaussian process
                gauss_process = lc_colors._lc_regression.fit_gaussian_process(data=photo_data)
                # Predict photometric data at the spectra times for specific band
                # Units of Jy?
                pred_flux, pred_flux_err = lc_colors._lc_prediction.predict_band_flux(gp=gauss_process, band_name=band, times=spec_times)
                # Convert flux to magnitude
                pred_mag = flux_to_mag(flux=pred_flux, zp=sndata.get_zp(band_name=band), offset=get_offset_for_band(band=band))
                for i, mag in enumerate(pred_mag):
                    out_table.add_row([oid, spec_times[i], band, mag])
            if len(photo_data) < 1:
                # TODO different placeholder?
                out_table.add_row([oid, 0.0, band, np.nan])

    return out_table
     
       
def make_table():
    """Create photometry comparison table
    """

    # Get synthetic photometry for all spectroscopic data
    synthetic_photo_table = vstack([tabulate_synthetic_photometry(data) for data in dr1.iter_data(verbose=True)])
    # Regress all photometric data 
    photo_table = photometry_to_spectra_time(spec_data=synthetic_photo_table)
    # Create output table
    tbl = join(synthetic_photo_table, photo_table, join_type='left')
    # TODO display filename at prompt?
    # pathname = os.path.join('..', '..', 'csp', 'test_output_data')
    # fn = os.path.join(pathname, 'comparision.fsits')
    fn = '../../csp/test_output_data/photo_comp.fits'
    ans = input('Overwrite joint table {}? (y/n)'.format(fn))
    if ans == 'y':
        tbl.write(fn, overwrite=True)
    if ans == 'n':
        fn_end = input('New filename:')
        fn = '../../csp/test_output_data/{}.fits'.format(fn_end)
        tbl.write(fn)
    print(fn)
    return None
