"""This script generates a config file that specifies supernova specific priors
and fitting arguments for CSP. The config files is saved in yaml format and
values typically follow the following template:

    object_id:
      kwargs:
        bounds:
          t0:
            - ? # Typically t0 - 10
            - ? # Typically t0 + 10
        phase_range:
          - -20
          - 50
    object_id:
      priors:
        z: ?
        t0: ?
"""

import sys
from pathlib import Path

import sfdmap
import yaml
from sndata.csp import dr3

sys.path.insert(0, '../')
from analysis.utils import filter_has_csp_data

dr3.download_module_data()
DUST_MAP = sfdmap.SFDMap('./schlegel98_dust_map/')

csp_measurements = dr3.load_table(3)
csp_measurements['T(Bmax)'] += 2400000.5  # Convert MJD to JD
t0_dict = dict(zip(csp_measurements['SN'], csp_measurements['T(Bmax)']))

if __name__ == '__main__':
    config_dict = dict()
    for table in dr3.iter_data(filter_func=filter_has_csp_data):
        obj_id = table.meta['obj_id']
        z = table.meta['z']
        ra = table.meta['ra']
        dec = table.meta['dec']

        # Typecast avoids issues with yaml export
        mwebv = float(DUST_MAP.ebv(ra, dec))
        t0 = float(t0_dict[obj_id])

        config_dict[obj_id] = {
            'priors': {
                't0': t0,
                'z': z,
                'mwebv': mwebv
            },
            'kwargs': {
                'bounds': {
                    't0': [t0 - 10, t0 + 10]
                },
                'phase_range': [-20, 50],
                'vparam_names': ['t0', 'x0']
            }
        }

    with open('config.yml', 'w') as out_file:
        yaml.dump(config_dict, out_file)
