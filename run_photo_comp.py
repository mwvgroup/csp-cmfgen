from astropy.table import join, vstack

from analysis import lc_colors
from analysis import spectra_chisq
from sndata.csp import dr1, dr3


dr1.download_module_data()
dr3.download_module_data()
dr3.register_filters(force=True)


def run_make_table():

    spectra_chisq.photometry_comparison.make_table()


# TODO add plotting code?

run_make_table()