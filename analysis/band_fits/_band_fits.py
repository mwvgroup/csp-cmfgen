from analysis import models, utils
from astropy.io import ascii
from astropy.table import Table, unique
import numpy as np
import os.path
from pathlib import Path
import sfdmap
import sncosmo
from sndata.csp import dr3


models.register_sources(force=True)
dr3.download_module_data()
dr3.register_filters(force=True)

dustmap = sfdmap.SFDMap("../../sfddata-master")
dust = sncosmo.F99Dust()


# Load models for different masses, with dust extinction.
m102 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.02), effects=[dust, dust], effect_names=['host', 'mw'],
                     effect_frames=['rest', 'obs'])
m104 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.04), effects=[dust, dust], effect_names=['host', 'mw'],
                     effect_frames=['rest', 'obs'])
m14 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.4), effects=[dust, dust], effect_names=['host', 'mw'],
                    effect_frames=['rest', 'obs'])
m17 = sncosmo.Model(sncosmo.get_source('CMFGEN', version=1.7), effects=[dust, dust], effect_names=['host', 'mw'],
                    effect_frames=['rest', 'obs'])

model_list = [m102, m104, m14, m17]

unique_bands = [
    'csp_dr3_u',
    'csp_dr3_g',
    'csp_dr3_r',
    'csp_dr3_i',
    'csp_dr3_B',
    'csp_dr3_V',
    'csp_dr3_V0',
    'csp_dr3_Y',
    'csp_dr3_Ydw',
    'csp_dr3_J',
    'csp_dr3_H'
]

BAND_FITS_DIR = Path(__file__).resolve().parent.parent.parent / 'results' / 'band_fits'


def get_band_fits(bands, fix_t0=False):
    """ Fits band data from csp dr3 to 4 CMFGEN models

    Args:
        bands (List): list of band names (strings) specifying what photometric data we're fitting with
        fix_t0 (bool): whether we're interested in fixing the T0 value to be +/- 3 days of published values

    Returns:
        Astropy Table of the fit results.  If previously calculated a cached copy is returned.
    """

    # sort the bands so that the generated files are saved with names that are not dependent on the band ordering
    bands.sort(key=unique_bands.index)
    file_location = str(BAND_FITS_DIR) + '/' + '_'.join([x.replace('csp_dr3_', '') for x in bands]) \
                    + ('_T0' if fix_t0 else '') + '.ecsv'
    if os.path.exists(file_location):
        return ascii.read(file_location)
    else:
        results = Table(names=('object ID', 'model', 'version', 'success', 'chisq', 'ndof', 'z', 't0', 't0_err',
                               'x0', 'x0_err', 'hostebv', 'hostebv_err', 'mwebv'),
                        dtype=('U10', 'U12', 'U12', bool, 'f4', 'i4', float, float, float, float, float, float,
                               float, float))
        results.meta = {'bands': bands.copy()}

        bands_set = set(bands)
        for index, data in enumerate(dr3.iter_data()):
            print(index)
            # if we're constraining T0, find out if there actually is published T0 data.
            # If there isn't, continue to next object in the survey
            bounds = None
            if fix_t0:
                try:
                    t0 = utils.get_csp_t0(data.meta['obj_id'])
                    bounds = {'t0': (t0 - 3, t0 + 3)}
                except utils.NoCSPData as no_data:
                    print(f"no T0 data for {data.meta['obj_id']}")
                    continue

            # make sure object data has readings for all requested bands
            all_observed_bands = set(unique(data, keys='band')['band'])
            if not bands_set.issubset(all_observed_bands): continue

            data_subset = data[np.in1d(data['band'], bands)]
            for model in model_list:
                try:
                    out_data = np.array([data_subset.meta['obj_id'], model.source.name, str(model.source.version)])
                    model.set(z=data_subset.meta['redshift'])
                    model.set(mwebv=dustmap.ebv(data_subset.meta['ra'], data_subset.meta['dec']))
                    result, fitted_model = sncosmo.fit_lc(
                        data_subset, model,
                        ['t0', 'x0', 'hostebv'],
                        bounds=bounds)
                    param_dict = dict(zip(result.param_names, result.parameters))
                    out_data = np.append(out_data,
                                         [1 if result['success'] else 0, result['chisq'], str(int(result['ndof'])),
                                          param_dict['z'], param_dict['t0'], result.errors['t0'],
                                          param_dict['x0'], result.errors['x0'], param_dict['hostebv'],
                                          result.errors['hostebv'],
                                          param_dict['mwebv']])

                    results.add_row(out_data)
                except Exception as err:
                    print("ERROR: ", data.meta['obj_id'], " - ", model.source.name, "_", model.source.version, " -> ",
                          err)
        ascii.write(results, file_location, format='ecsv', fast_writer=False)
        return results

