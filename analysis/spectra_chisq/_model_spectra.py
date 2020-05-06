"""The ``_model_spectra`` module models spectra by varying the dust extinction from the Milky Way.

The model used is ``sncosmo.CCM89Dust``. The parameters A_v and amplitude are fit to the observed spectra. The best fit
parameters are determined by ``scipy.optimize.curve_fit``
The fit parameters are saved to a fits file.

Currently, run this module at the command line using: "python _model_spectra.py"
"""


import matplotlib.pyplot as plt
import numpy as np
import sys
import sncosmo
from astropy.table import Table
from functools import partial
from scipy.optimize import curve_fit
from sndata.csp import DR1, DR3

from astropy.table import join, vstack #unique

#from ._synthetic_photometry import get_synth_photo_for_id # FIXME uncomment if fit performed using photometry

dr1, dr3 = DR1(), DR3()
dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters(force=True)


def mjd_to_jd(time):
    """Convert MJD time to JD
    JD = MJD + 2400000.5
    
    Args:
        time (float): JD time
        
    Returns:
        (float) MJD time
    """
    
    return time + 2400000.0


def get_jd_t0(obj_id, tbl_num=3, mjd=True):
    """Get time of peak brightness of SN in B-band
    
    Args:
        obj_id (str): specific SN object id
        tbl_num (int): `3` refers to Table 3 in DR3 paper (The Carnegie Supernova Project I:
            Third Photometry Data Release of Low-Redshift Type Ia Supernovae and Other White Dwarf Explosions)
        mjd (bool): True if the time given by `tbl_num` is MJD. False if it is JD
    
    Returns:
        The JD time of max brightness as a float
    """

    tbl = dr3.load_table(tbl_num)
    
    # Get data for specific SN
    obj_dat = tbl[tbl['SN'] == obj_id]
    # Get t0 for specific SN
    t0 = obj_dat['T(Bmax)'] # MJD, see DR3
    if mjd:
        # Convert t0 to JD time
        return mjd_to_jd(time=t0)
    if not mjd:
        return t0


def create_model(model, obj_id, a_v, redshift, amplitude, r_v=3.1):
    """Model Milky Way dust extinction.
    Set the ``z``, ``t0``, ``mwebv``, and ``amplitude`` parameters of a Model for a given CSP target

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


def get_model_flux(wave, amplitude, a_v, obj_id, time, model, redshift):
    """Function to use for fitting amplitude and Av of model.
    Note: order of args determines which parameters are fit
    
    Args:
        wave (1d array-like of floats): wavelengths in units of Angstrom
        amplitude (float): fit parameter
        a_v (float): fit parameter
        photo_data (Table):
        redshift (float):
        obj_id:
        time (float or array): JD
        model (sncosmo dust model):
        
    Returns:
        The fluxes (as a 1d array-like) at each wavelength in `wave` computed using `model.` Units: ergs/s/cm^2/Angstrom
    """

    model = create_model(model=model, obj_id=obj_id, a_v=a_v, redshift=redshift, amplitude=amplitude)

    # https://sncosmo.readthedocs.io/en/v1.1.x/models.html
    # Get spectrum at a given observer-frame time and wavelengths (spectrum contains OBSERVER frame flux)
    # Note from docs: can supply a list or array of times and get a 2-d array back,
    # representing the spectrum at each time
    mdl_flux = model.flux(time, wave) # Units: ergs/s/cm^2/Angstrom
    
    return mdl_flux


def plot_spectra_and_residuals(x1, y1, x2, y2):
    """Make two panel plot. Upper plot: data and model flux vs wavelength. Lower plot: residuals (data - model)
    
    Args:
        x1 (1d array-like of floats): Observer frame wavelength
        y1 (1d array-like of floats): Observed flux
        x2 (1d array-like of floats): Observer frame wavelength
        y2 (1d array-like of floats): Modeled flux
    """
    
    # Set figure size
    plt.rcParams['figure.figsize'] = [12, 10] # units: inches
    # Create two stacked figures
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace':0, 'width_ratios':[1]})
    
    # Upper plot: flux vs wavelength
    # Plot wavelength vs obs flux
    ax1.plot(x1, y1, linestyle=' ', marker='.', label='obs')
    # Plot wavelength vs model flux
    ax1.plot(x2, y2, linestyle=' ', marker='.', label='fit')
    # Labels
    ax1.set_ylabel('flux [ergs/(s cm^2 Angstrom)]', fontsize=14)
    ax1.legend(fontsize=14)
    #plt.suptitle('id: {}, band: {}'.format(OID,BND), fontsize=14)
    # Lower plot: residuals (data - model) in lower plot
    ax2.plot(x2, y1 - y2, color='green', linestyle=' ', marker='.')
    ax2.axhline(0, color='black', linestyle=':')
    ax2.set_xlabel('obs frame wavelength [Angstrom]', fontsize=14)
    ax2.set_ylabel('Observed - model flux \n [ergs/(s cm^2 Angstrom)]', fontsize=14)

    return None


def plot_spectra(x1, y1, x2, y2):
    """Plot the spectra from the data and the model (flux vs wavelength).
    
    Args:
        x1 (1d array-like of floats): Observer frame wavelength
        y1 (1d array-like of floats): Observed flux
        x2 (1d array-like of floats): Observer frame wavelength
        y2 (1d array-like of floats): Modeled flux
    """
    
    plt.clf()
    # Set figure size
    plt.rcParams['figure.figsize'] = [12, 6] # units: inches
    
    # Plot: flux vs wavelength
    # Plot  wavelength vs obs flux
    plt.plot(x1, y1, linestyle=' ', marker='.', label='obs')
    # Plot wavelength vs model flux
    plt.plot(x2, y2, linestyle=' ', marker='.', label='fit')
    # Labels
    plt.xlabel('obs frame wavelength [Angstrom]')
    plt.ylabel('flux [ergs/(s cm^2 Angstrom)]', fontsize=14)
    plt.legend(fontsize=14)
    #plt.suptitle('id: {}, band: {}'.format(OID,BND), fontsize=14)
    # Lower plot: residuals (data - model) in lower plot

    plt.show()

    return None


def get_fit_parameters(obj_id, time, flux, model, wave, redshift, fit_u_band=True):
    """Get fit parameters (amplitude and Av) ...

    fit_u_band (bool): If True, u-band is included in the fits TODO add False case
    
    Returns:
        - The best-fit amplitude as a float
        - The best-fit Av as a float
    """
        
    # Create partial function to use `curve_fit()` with fixed parameters: obj_id, time, model, redshift, photo_data
    # Need `partial_model_flux()` function args to be parameters that will be fit
    # `get_model_flux()` returns OBSERVER frame flux (https://sncosmo.readthedocs.io/en/v1.1.x/models.html)
    partial_model_flux = partial(get_model_flux, obj_id=obj_id, time=time, model=model, redshift=redshift)
    # Fits the 2nd and 3rd parameters of `get_model_flux()` using the first parameter
    # Note: order of args matters in `get_model_flux()`
    if fit_u_band: 
        fit_params, cov_matrix = curve_fit(partial_model_flux, wave, flux)
    
    return fit_params


def get_model_synth_table(time, amplitude, a_v, obj_id, wave, model, redshift, spec_release, phot_release, temp_spec_tbl, spec_data):
    """Get synthetic magnitude from modeled spectra. `time` is the first arg if `time` will be used in `curve_fit()`...
    FIXME base curve_fit() on time?

    Args:

    Returns:

    """

    # FIXME this function is called many times

    # Get flux. Flux will be used to calculate synthetic magnitude
    # Get flux at each time. Wavelengths must be monotonically increasing

    # all_flux = [get_model_flux(wave, amplitude, a_v, obj_id, t, model, redshift) for t in np.unique(time)] # time can be an array

    # TODO: is it better to use lists?
    # Get the wavelengths that were imaged at each time
    all_flux = np.empty(len(np.unique(time)), dtype=object)
    for i, _time in enumerate(np.unique(time)):
        spectrum = spec_data[spec_data['time'] == _time]
        _wave = spectrum['wavelength']
        flux = get_model_flux(_wave, amplitude, a_v, obj_id, _time, model, redshift)
        # Set each element with an array (dtype=object)
        all_flux[i] = flux

    # Flatten array, it will be used as a column in a table
    mdl_flux = np.concatenate(all_flux, axis=None)

    ##all_flux = get_model_flux(wave, amplitude, a_v, obj_id, np.unique(time), model, redshift)
    ##all_flux = get_model_flux(np.unique(wave), amplitude, a_v, obj_id, np.unique(time), model, redshift)

    #mdl_flux = vstack(all_flux)
    # Flatten array that contains all fluxes TODO try np.flatten() ?
    ##mdl_flux = np.concatenate(all_flux, axis=None)

    print('----->len model flux', len(mdl_flux))
    print('----->len spectra table', len(spec_data))
    print('----->len time', len(time))

    # Get synthetic magnitude using flux
    combo_tbl, spec_tbl = get_synth_photo_for_id(obj_id, spec_release, phot_release, err_ratio=0.03, temp_flux=mdl_flux, temp_spec_tbl=temp_spec_tbl,
                                                 out_path=None, time=np.unique(time)) #TODO unique time?
    # Filter by obj_id FIXME?
    #combo_tbl = combo_tbl[combo_tbl['obj_id'] == obj_id]
    # FIXME mdl_flux should be the same length as output table...check this
    # Remove nans; cannot use nans in `curve_fit()`  which fits using photometric and synthetic magnitudes
    combo_tbl = combo_tbl[np.logical_not(np.isnan(combo_tbl['photo_mag']))]
    combo_tbl = combo_tbl[np.logical_not(np.isnan(combo_tbl['synth_mag']))]

    #synth_mag = combo_tbl['synth_mag']

    # TODO check that the length of this synthetic magnitude table is the same as the initial table

    #return synth_mag, combo_tbl['photo_mag'] #FIXME remove second output? or split into two functions?

    return combo_tbl


def get_model_synth_mag(time, amplitude, a_v, obj_id, wave, model, redshift, spec_release, phot_release, temp_spec_tbl, spec_data):
    """
    :param time:
    :param amplitude:
    :param a_v:
    :param obj_id:
    :param wave:
    :param model:
    :param redshift:
    :param spec_release:
    :param phot_release:
    :param temp_spec_tbl:
    :return:
    """

    tbl = get_model_synth_table(time, amplitude, a_v, obj_id, wave, model, redshift, spec_release, phot_release,
                          temp_spec_tbl,  spec_data=spec_data)

    return tbl['synth_mag']


def get_photo_fit_parameters(obj_id, time, model, wave, redshift, photo_mag, phot_release, spec_release, temp_spec_tbl, synth_data, spec_data,
                             fit_u_band=True):
    """Get fit parameters (amplitude and Av) ...
    Fit using the photometric magnitude and the synthetic magnitude.
    This function is distinct from `get_fit_parameters()` because `curve_fit()` is performed using time, ... etc
    and other different args, but the returns are the same (best fit to amplitude and Av)

    Args:

    Returns:
        - The best-fit amplitude as a float
        - The best-fit Av as a float
    """

    # `partial_synth_mag()` must be a function to perform fit using `curve_fit()`
    # can only return synth_mag ...
    partial_synth_mag = partial(get_model_synth_mag, obj_id=obj_id, wave=wave, model=model, redshift=redshift,
                                spec_release=spec_release, phot_release=phot_release, temp_spec_tbl=temp_spec_tbl, spec_data=spec_data)

    # Want to fit synthetic magnitude to observed magnitude at several times (time should be an array)
    # and want to perform fits for each band? (all bands imaged at each time) TODO/FIXME
    print('-----> len time used in curve_fit: ', len(time))
    print('-----> len photo_mag used in curve_fit (should be more than time if not filtering by band...: ', len(photo_mag))

    # TODO check that the length of this is the same as original `synth_data`
    tabl = get_model_synth_table(obj_id=obj_id, wave=wave, model=model, redshift=redshift,
                                spec_release=spec_release, phot_release=phot_release, temp_spec_tbl=temp_spec_tbl, spec_data=spec_data,
                            time=time, amplitude=10**-8, a_v=0.3)
    print('-----> length of synthetic and photo mag: ', len(tabl), len(photo_mag))
    photo_mag = tabl['photo_mag'] #FIXME remove; guesses for amp, av not releveant do not affect photo mag

    # Note: `curve_fit()` can handle duplicate numbers in `time` array
    fit_params, cov_matrix = curve_fit(partial_synth_mag, time, photo_mag, p0=[10**-8, 0.42])

    print('______________________________ fit parameters: ', fit_params)

    return fit_params


def fit_to_obs_spectra():
    """Extinct model to fit to the observed spectra. Tabulate results.
    Plots the last result in the loop currently. TODO move plotting elsewhere

    Output table columns:
        - obj_id
        - time
        - amplitude
        - a_v

    TODO run this via `run_analysis.py`
    """

    # Create table to save fit parameters to
    table = Table(names=('obj_id', 'time', 'amplitude', 'a_v'), dtype=('S20', 'S20', 'f8', 'f8'))

    # Create model for Milky Way dust
    mw_dust = sncosmo.CCM89Dust()
    # SN model with Milky Way dust
    sn_model = sncosmo.Model(source='hsiao', effects=[mw_dust], effect_names=['mw'], effect_frames=['obs'])

    for oid in dr3.get_available_ids()[8:9]: #FIXME remove [:#] after done testing
        # Get spectroscopic data for specific SN
        spec_data = dr1.get_data_for_id(oid)
        # Get observation times
        times = spec_data['time']
        # Get object specific parameters: photometric data and redshift
        _photo_data = dr3.get_data_for_id(oid)
        _redshift = _photo_data.meta['z']

        # Loop over unique observation times
        for _time in set(times):
            # Get spectrum taken at specific `_time`
            spectrum = spec_data[spec_data['time'] == _time]
            # Observer-frame (?) flux TODO
            _flux = spectrum['flux'] # Units: ergs/s/cm^2/Angstrom (flux per wavelength)
            # From the README.txt in 'CSP_spectra_Folatelli+2013.tgz' on https://csp.obs.carnegiescience.edu/data:
            # "WAVELENGTH IS EXPRESSED IN THE REST FRAME OF THE SUPERNOVA IN ALL CASES."
            rest_frame_wave = spectrum['wavelength'] # Units: Angstrom
            # Get the observer-frame wavelength, because sncosmo `.flux()` returns observer frame flux...
            obs_frame_wave = rest_frame_wave * (1 + _redshift)
            # Fit the parameters: amplitude and Av, in SN model
            _amp, _av = get_fit_parameters(obj_id=oid, time=_time, model=sn_model, wave=obs_frame_wave, flux=_flux,
                                           redshift=_redshift)

            # TODO / FIXME many fits return initial parameter guesses
            print('SN: {}, redshift: {}, time: {}, amplitude: {}, Av: {}'.format(oid, _redshift, _time, _amp, _av))

            table.add_row((oid, _time, _amp, _av))

    #table.write('extinction_fit_parameters.fits', format='fits', overwrite=False)

    # Get flux from the SN model using the best-fit parameters
    mdl_flux = get_model_flux(wave=obs_frame_wave, obj_id=oid, time=_time, model=sn_model,
                              amplitude=_amp, a_v=_av, redshift=_redshift)
    plot_spectra(x1=obs_frame_wave, y1=_flux, x2=obs_frame_wave, y2=mdl_flux)



def fit_band_to_obs_photo():
    """TODO perform fit for each photometric band?"""

    return 0


def fit_to_obs_photo(spec_release=dr1, phot_release=dr3):
    """Extinct spectra to fit synthetic magnitude to observed magnitude.

    Output table columns:
        - obj_id
        - amplitude
        - a_v

    Args:
        spec_release (module): An spectroscopic data release from sndata
        phot_release (module): A photometric data release from sndata

    Returns:

    """

    # Create table to save fit parameters to
    # TODO add column describing how fit was obtained?
    # TODO include time in data table?
    table = Table(names=('obj_id', 'amplitude', 'a_v'), dtype=('S20', 'f8', 'f8'))

    # Create model for Milky Way dust
    mw_dust = sncosmo.CCM89Dust()
    # SN model with Milky Way dust
    sn_model = sncosmo.Model(source='hsiao', effects=[mw_dust], effect_names=['mw'], effect_frames=['obs'])

    for oid in dr3.get_available_ids()[11:12]:  # FIXME remove [:#] after done testing

        # Get initial table of synthetic photometry and regressed photometry
        # so that the photometry and spectra are at the same times (the spectra times)
        # This table contains magnitudes for one SN for all observation times; need regressed times and regressed fluxes
        # Need spectroscopic data (wavelengths and fluxes) for synthetic photometry
        # This table may have nans
        _synth_data, _spec_data = get_synth_photo_for_id(obj_id=oid, spec_release=spec_release, phot_release=phot_release)
        # spec_data = dr1.get_data_for_id(oid)

        # Wavelength must be increasing when using sncosmo `flux()`, so sort by wavelength
        _spec_data.sort('wavelength')

        # Get redshift
        _redshift = dr3.get_data_for_id(oid).meta['z']

        # Get all data at same band (wavelength range TODO)
        spectrum = _spec_data # TODO filter wavelength range by band? From DR3: Y-band filter (900 to 1100 nm)
        # Get wavelengths
        rest_frame_wave = spectrum['wavelength']  # Units: Angstrom
        obs_frame_wave = rest_frame_wave * (1 + _redshift) # Units: Angstrom

        #_time = np.unique(synth_data['time']) # Each of the (9?) bands are imaged at each time ...
        _time = _synth_data['time']
        _photo_mag = _synth_data['photo_mag']

        # Maybe some nans from initial synthetic photometry table will not be nans with new extinction parameters;
        # include all times TODO ?

        # Get time, wave, etc from the same table so each array is the same length TODO

        # Fit the parameters amplitude and Av in SN model
        # NOTE: one spectrum in `spectrum` for each time
        amp, av = get_photo_fit_parameters(synth_data=_synth_data, obj_id=oid, time=_time, model=sn_model, wave=obs_frame_wave,
                                           redshift=_redshift, photo_mag=_photo_mag, temp_spec_tbl=spectrum, spec_data=_spec_data,
                                           spec_release=spec_release, phot_release=phot_release)
        # TODO / FIXME fits return initial parameter guesses
        # FIXME time is column name?
        # What should this table contain? just SN object id and the best-fit parameters?
        # the other table has 'time'; include photometric band?
        print('SN: {}, redshift: {}, amplitude: {}, Av: {}'.format(oid, _redshift, amp, av))

        table.add_row((oid, amp, av))

    print(table)


# TODO move to `run_analysis.py`
fit_to_obs_spectra()