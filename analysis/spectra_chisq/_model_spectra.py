"""The ``_model_spectra`` module models spectra by varying the dust extinction from the Milky Way.

The model used is ``sncosmo.CCM89Dust``. The parameters A_v and amplitude are fit to the observed spectra. The best fit  parameters are determined by ``scipy.optimize.curve_fit``
The fit parameters are saved to a fits file.
"""

#FIXME unique variable names


import matplotlib.pyplot as plt
import sncosmo
from astropy.table import Table
from functools import partial
from scipy.optimize import curve_fit
from sndata.csp import DR1, DR3

dr1, dr3 = DR1(), DR3()
dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters(force=True)


def mjd_to_jd(date):
    """Convert MJD date to JD
    JD = MJD + 2400000.5
    
    Args:
        date (float): JD date
        
    Returns:
        (float) MJD date
    """
    
    return date + 2400000.0


def get_jd_t0(obj_id, tbl_num=3, mjd=True):
    """Get time of peak brightness of SN in B-band
    
    Args:
        obj_id (str): specific SN object id
        tbl_num (int): `3` refers to Table 3 in DR3 paper (The Carnegie Supernova Project I:
            Third Photometry Data Release of Low-Redshift Type Ia Supernovae and Other White Dwarf Explosions)
        mjd (bool): True if the date given by `tbl_num` is MJD. False if it is JD
    
    Returns:
        The JD date of max brightness as a float
    """

    tbl = dr3.load_table(tbl_num)
    
    # Get data for specific SN
    obj_dat = tbl[tbl['SN'] == obj_id]
    # Get t0 for specific SN
    t0 = obj_dat['T(Bmax)'] # MJD, see DR3
    if mjd:
        # Convert t0 to JD date
        return mjd_to_jd(date=t0)
    if not mjd:
        return t0


def create_model(model, obj_id, a_v, redshift, amplitude, r_v=3.1):
    """Model Milky Way dust extinction.
    Set the ``z``, ``t0``, ``mwebv``, and ``amplitude`` parameters of a Model for a given CSP target`

    Note: DR1 neglects extinction from host galaxy, so only MW dust is modeled currently
    Note: typical Av ~ 0.42 mag via Av = Rv * E(B-V) = 3.1 * 0.134mag ~ 0.42mag
    where Milky Way ebv = 0.134 (from Table 8 of DR3/Folatelli et al. 2013)
    and Milky Way Rv = 3.1 from DR1 at the end of pg 9
    
    Args:
        model (Model): sncosmo dust model
        obj_id (str):
        redshift (float)
        amplitude (float):
        a_v (float): visual (V filter) extinction in units of magnitude
        r_v (float):
    
    Returns:
        An updated sncosmo model
    """
    
    ebv = a_v / r_v # unit: mags
    # Get time of peak brightness in JD
    t0 = get_jd_t0(obj_id=obj_id)
    # From the example of `model.set()` on https://sncosmo.readthedocs.io/en/v1.6.x/models.html,
    # it looks like t0 is in JD
    model.set(z=redshift, t0=t0, mwebv=ebv, amplitude=amplitude)
    
    return model


def get_model_flux(wave, amplitude, a_v, obj_id, date, model, redshift):
    """Function to use for fitting amplitude and Av of model.
    Note: order of args determines which parameters are fit
    
    Args:
        wave (1d array-like of floats): wavelengths in units of Angstrom
        amplitude (float): fit parameter
        a_v (float): fit parameter
        photo_data (Table):
        redshift (float):
        obj_id:
        date: MJD?
        model (sncosmo dust model):
        
    Returns:
        The fluxes (as a 1d array-like) at each wavelength in `wave` computed using `model.` Units: ergs/s/cm^2/Angstrom
    """
    
    model = create_model(model=model, obj_id=obj_id, a_v=a_v, redshift=redshift, amplitude=amplitude)

    # https://sncosmo.readthedocs.io/en/v1.1.x/models.html
    # Get spectrum at a given time and wavelengths (OBSERVER frame flux, in ergs/s/cm^2/Angstrom)
    mdl_flux = model.flux(date, wave) 
    
    return mdl_flux


def plot_spectra_and_residuals(x1, y1, x2, y2):
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


def plot_spectra(x1, y1, x2, y2):
    """Plot the spectra from the data and the model (flux vs wavelength).
    
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


def get_fit_parameters(obj_id, date, flux, model, wave, redshift, fit_u_band=True):
    """Get fit parameters

    fit_u_band (bool): If True, u-band is included in the fits #TODO add False case
    
    Returns:
        - The best-fit amplitude as a float
        - The best-fit Av as a float
    """
        
    # Create partial function to use `curve_fit` with fixed parameters: obj_id, date, model, redshift, photo_data
    # Need `par_model_flux` function args to be parameters that will be fit
    # `get_model_flux` returns OBSERVER frame flux (https://sncosmo.readthedocs.io/en/v1.1.x/models.html)
    partial_model_flux = partial(get_model_flux, obj_id=obj_id, date=date, model=model, redshift=redshift)
    # Fits the last two parameters of `par_model_flux` note order of args matters in `get_model_flux`!
    if fit_u_band: 
        fit_params, cov_matrix = curve_fit(partial_model_flux, wave, flux)
    
    return fit_params


def run():
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
        _photo_data = dr3.get_data_for_id(oid)
        _redshift = _photo_data.meta['z']

        for _date in set(dates):
            # Observered wavelength and flux (in observer frame)
            spectrum = spec_data[spec_data['time'] == _date]
            _wave = spectrum['wavelength'] # Units: Angstrom
            _flux = spectrum['flux'] # Units: ergs/s/cm^2/Angstrom (flux per wavelength)
            # Fit the parameters amplitude and Av in SN model
            amp, av = get_fit_parameters(obj_id=oid, date=_date, model=sn_model, wave=_wave, flux=_flux, redshift=_redshift)
            mdl_flux = get_model_flux(wave=_wave, obj_id=oid, date=_date, model=sn_model,
                                      amplitude=amp, a_v=av, redshift=_redshift)
            # TODO / FIXME many fits return initial parameter guesses
            print('SN: {}, redshift: {}, date: {}, amplitude: {}, Av: {}'.format(oid, _redshift, _date, amp, av))

            table.add_row((oid, _date, amp, av))

    #table.write('spec_fit_param_log.fits', format='fits', overwrite=False)

    # `mdl_flux` is in observer frame
    # if `wave` and `flux` are also in observer frame, why is plot offset?
    plot_spectra(x1=_wave, y1=_flux, x2=_wave, y2=mdl_flux)


run()