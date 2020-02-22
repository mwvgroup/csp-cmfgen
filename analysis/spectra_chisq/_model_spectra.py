"""The ``_model_spectra`` module models spectra by varying the dust extinction from the Milky Way.

The model used is ``sncosmo.CCM89Dust``. The parameters A_v and amplitude are fit to the observed spectra. The best fit  parameters are determined by ``scipy.optimize.curve_fit``
The fit parameters are saved to a fits file.
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from functools import partial
from scipy.optimize import curve_fit

import sncosmo
import sndata
from sndata.csp import DR1, DR3


dr1, dr3 = DR1(), DR3()
dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters(force=True)


def jd_to_mjd(date):
    """Convert JD date to MJD
    
    Args:
        date (float): JD date
        
    Returns:
        (float) MJD date
    """
    
    return date + 2400000.0


def get_t0(obj_id, tbl=dr3.load_table(3), mjd=False):
    """Get peak brightness of SN in B-band
    
    Args:
        obj_id (str): specific SN object id
        tbl (Table): 
        jd (bool): True if the date in `tbl` is JD. False if it is MJD.
    
    Returns:
        t0 (float): MJD date of max brightness
    """
    
    # Get data for specific SN
    obj_dat = tbl[tbl['SN'] == obj_id]
    # Get t0 for specific SN
    t0 = obj_dat['T(Bmax)']
    if not mjd:
        # Convert t0 to MJD date
        t0_mjd = jd_to_mjd(date=t0)
    
    return t0_mjd


def create_model(model, obj_id, a_v, redshift, amplitude, photo_data, t0_tbl=dr3.load_table(3),
                        r_v=3.1, mag_sys='ab'):
    """Model Milky Way dust extinction. Set parameters: EBV (color excess), t0 (time of max brightness in B band),
    and redshift.
    Note: DR1 neglects extinction from host galaxy, so only MW dust is modeled currently

    
    Note: typical Av ~ 0.42 mag via Av = Rv * E(B-V) = 3.1 * 0.134mag ~ 0.42mag
    where Milky Way ebv = 0.134 (from Table 8 of DR3/Folatelli et al. 2013)
    and Milkway Way Rv = 3.1 from DR1 at the end of pg 9
    
    Args:
        model (): sncosmo dust model
        obj_id (str):
        a_v (float): visual (V filter) extinction in units of magnitude
    
    Returns:
        model (): Updated sncosmo model
    """
    
    ebv = a_v / r_v # unit: mags?
    # Get time of peak brightness
    t0 = get_t0(obj_id=obj_id)
    model.set(z=redshift, t0=t0, mwebv=ebv, amplitude=amplitude)
    
    return model


def get_obs_wave_and_flux(obj_id, spectra, date):
    """ """
    
    spectrum = spec_data[spec_data['time'] == date]
    wave = spectrum['wavelength']
    flux = spectrum['flux']
    
    return wave, flux


def get_model_flux(wave, amplitude, a_v, obj_id, date, model, photo_data, redshift):
    """Function to use for fitting amplitude and Av of model.
    Note: order of args determines which parameters are fit
    
    Args:
        wave (1d array-like of floats): wavelengths in units of Angstrom
        amplitude (float): fit parameter
        a_v (float): fit parameter
        obj_id:
        date: MJD?
        model (sncosmo dust model):
        
    Returns:
        mdl_flux (1d array-like of floats): Fluxes computed using the `model.` Flux is computed at each wavelength
        in `wave.` Units: ergs/s/cm^2/Angstrom
    """
    
    model = create_model(model=model, obj_id=obj_id, a_v=a_v, redshift=redshift, amplitude=amplitude, photo_data=photo_data)

    # https://sncosmo.readthedocs.io/en/v1.1.x/models.html
    # Get spectrum at a given time and wavelengths (OBSERVER frame flux, in ergs/s/cm^2/Angstrom)
    mdl_flux = model.flux(date, wave) 
    
    return mdl_flux


def log_fit(tbl, obj_id, date, amp, av):
    """Logs one row of fit parameters in astropy table.
    Headers: oid, band, date, amplitude, av (unit: mag)
    
    Args:
        tbl (astropy Table): table to write to
    
    Returns:
        Updated table
    """
    #TODO input param: filename
    
    tbl.add_row((obj_id, date, amp, av))
    
    return tbl


def plotter(x1, y1, x2, y2):
    """Make two panel plot. Upper plot: data and model flux vs wavelength. Lower plot: residuals (data - model)
    
    Args:
        x1 (1d array-like of floats): Wavelength for observed
        y1 (1d array-like of floats): Observed flux
        x2 (1d array-like of floats): Wavelength for model
        y2 (1d array-like of floats): Modeled flux
    """
    
    # Set figure size
    plt.rcParams['figure.figsize'] = [12, 10] # units: inches
    # Create two stacked figures
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace':0, 'width_ratios':[1]})
    
    # Upper plot: flux vs wavelength
    # Plot redshifted wavelength vs obs flux
    ax1.plot(x1, y1, linestyle=' ', marker='.', label='obs')
    # Plot wavelength vs model flux
    ax1.plot(x2, y2, linestyle=' ', marker='.', label='fit')
    # Labels
    ax1.set_ylabel('Flux [ergs/(s cm^2 Angstrom)]', fontsize=14)
    ax1.legend(fontsize=14)
    #plt.suptitle('id: {}, band: {}'.format(OID,BND), fontsize=16)
    # Lower plot: residuals (data - model) in lower plot
    ax2.plot(x2, y1 - y2, color='green', linestyle=' ', marker='.')
    ax2.axhline(0, color='black', linestyle=':')
    ax2.set_xlabel('wavelength [Angstrom]', fontsize=14)
    ax2.set_ylabel('Observed - model flux \n [ergs/(s cm^2 Angstrom)]', fontsize=14)

    return None


def plotter2(x1, y1, x2, y2):
    """Make two panel plot. Upper plot: data and model flux vs wavelength. Lower plot: residuals (data - model)
    
    Args:
        x1 (1d array-like of floats): Wavelength for observed
        y1 (1d array-like of floats): Observed flux
        x2 (1d array-like of floats): Wavelength for model
        y2 (1d array-like of floats): Modeled flux
    """
    
    plt.clf()
    # Set figure size
    plt.rcParams['figure.figsize'] = [12, 6] # units: inches
    
    # Plot: flux vs wavelength
    # Plot redshifted wavelength vs obs flux
    plt.plot(x1, y1, linestyle=' ', marker='.', label='obs')
    # Plot wavelength vs model flux
    plt.plot(x2, y2, linestyle=' ', marker='.', label='fit')
    # Labels
    plt.xlabel('wavelength [AA]')
    plt.ylabel('Flux [ergs/(s cm^2 Angstrom)]', fontsize=14)
    plt.legend(fontsize=14)
    #plt.suptitle('id: {}, band: {}'.format(OID,BND), fontsize=16)
    # Lower plot: residuals (data - model) in lower plot

    plt.show()

    return None


def get_fit_parameters(date, flux, model, obj_id, wave, photo_data, redshift, fit_u_band=True):
    """Get fit parameters

    fit_u_band (bool): If True, u-band is included in the fits #TODO add False case
    
    Returns:
        fit_params (list of floats): (amplitude, av)
    """
        
    # Create partial function to use `curve_fit` with fixed parameters: obj_id, date, model, redshift, photo_data
    # Need `par_model_flux` function args to be parameters that will be fit
    # `get_model_flux` returns OBSERVER frame flux (https://sncosmo.readthedocs.io/en/v1.1.x/models.html)
    partial_model_flux = partial(get_model_flux, obj_id=oid, date=date, model=model, redshift=redshift, photo_data=photo_data)
    # Fits the last two parameters of `par_model_flux` note order of args matters in `get_model_flux`!
    if fit_u_band: 
        fit_params, cov_matrix = curve_fit(partial_model_flux, wave, flux)
    
    return fit_params



# Create table to save fit parameters to
table = Table(names=('obj_id', 'time', 'amplitude', 'a_v'), dtype=('S20', 'S20', 'f8', 'f8'))

# Create model for Milky Way dust
mw_dust = sncosmo.CCM89Dust()
# SN model with Milky Way dust
sn_model = sncosmo.Model(source='hsiao', effects=[mw_dust], effect_names=['mw'], effect_frames=['obs'])

for oid in dr3.get_available_ids()[8:9]: #FIXME remove [:#] after done testing
    # Get spectroscopic data
    spec_data = dr1.get_data_for_id(oid)
    # Get uniques dates. Note: duplicate dates with different values in 'wavelength' column
    dates = set(spec_data['time'])
    # Get object specific parameters: photometric data and redshift
    photo_data = dr3.get_data_for_id(oid)
    redshift = photo_data.meta['z']
    
    for date in set(dates):
        # Observered wavelength and flux (in observer frame)
        wave, flux = get_obs_wave_and_flux(obj_id=oid, spectra=spec_data, date=date)
        # Fit the parameters amplitude and Av in SN model
        amp, av = get_fit_parameters(date=date, obj_id=oid, model=sn_model, wave=wave, flux=flux, redshift=redshift, photo_data=photo_data)
        mdl_flux = get_model_flux(wave, obj_id=oid, date=date, model=sn_model, 
                                  amplitude=amp, a_v=av, redshift=redshift, photo_data=photo_data)
        # TODO / FIXME many fits return initial parameter guesses
        print('SN: {}, redshift: {}, date: {}, amplitude: {}, Av: {}'.format(oid, redshift, date, amp, av))
        
        t = log_fit(tbl=table, obj_id=oid, date=date, amp=amp, av=av)

t.write('spec_fit_param_log.fits', format='fits', overwrite=True)
        
# `mdl_flux` is in observer frame
# if `wave` and `flux` are also in observer frame, why is plot offset?
plotter2(x1=wave, y1=flux, x2=wave, y2=mdl_flux)
