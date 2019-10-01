from astropy.table import join, vstack

from analysis import lc_colors
from analysis import spectra_chisq
from sndata.csp import dr1, dr3


dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters(force=True)


def make_table():
    """Create photometry comparison table
    """

    # Get synthetic photometry for all spectroscopic data
    synthetic_photo_table = vstack([spectra_chisq.photometry_comparison.tabulate_synthetic_photometry(data) for data in dr1.iter_data(verbose=True)])
    # Regress all photometric data 
    photo_table = spectra_chisq.photometry_comparison.photometry_to_spectra_time(spec_data=synthetic_photo_table)
    # Create output table
    tbl = join(synthetic_photo_table, photo_table, join_type='left')
    # TODO display filename at prompt?
    ans = input('Overwrite joint table? (y/n)')
    if ans == 'y':
        tbl.write('../../../../comp.fits', overwrite=True)
    if ans == 'n':
        # TODO also input path?
        fn = input('New filename (ex: table_name) :')
        tbl.write('../../../../{}.fits'.format(fn))
    
    return None

make_table()