"""
Given spectrum from PFS / Veloce / CHIRON / FEROS...
* Measure Li EW at 6708A.
* Vis Halpha at 6562.8A
* Measure Teff, vsini, logg.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy import units as u, constants as const
from copy import deepcopy
from cdips_followup.spectools import (
    read_pfs, read_veloce, viz_1d_spectrum,
    get_Li_6708_EW, inspect_pfs, specmatch_viz_compare, specmatch_analyze,
    plot_orders
)

# spectrum_path = '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/Veloce/20200130_837.01_Bouma_final_combined.fits'

if __name__ == "__main__":

    do_orders = 0
    do_inspect = 0      # inspect to figure out require rest-frame shift
    do_li_ew = 0        # once rest-frame shift is known
    do_sm_viz = 0       # specmatch-emp check
    do_sm_analysis = 1  # for vsini, Teff, Rstar

    if do_orders:

        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outdir = '../results/spec_analysis/PFS/spec_viz_orders/'
        idstring = 'TIC268301217_20200203_3556'

        plot_orders(spectrum_path, wvsol_path=wvsol_path, outdir=outdir,
                    idstring=idstring)

    if do_sm_viz:
        # 4990 through 6410 in the SpecMatch library. b/c HIRES needs different
        # echelle settings >~6200A (https://www2.keck.hawaii.edu/inst/hires/)
        specmatch_viz_compare(wavlim=[5160,5210])

    if do_inspect:
        for xlim in ['assign', 'fullorder']:
            for night in ['3555', '3556']:
                inspect_pfs(night, 'Halpha', xlim)
                inspect_pfs(night, 'LiI', xlim)
                inspect_pfs(night, 'Mgb1', xlim)

    if do_li_ew:
        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outpath = '../results/spec_analysis/PFS/Li_EW/3556_Li_EW_shift2.30.png'

        get_Li_6708_EW(spectrum_path, wvsol_path=wvsol_path, xshift=2.30,
                       outpath=outpath)

    if do_sm_analysis:
        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outdir = '../results/spec_analysis/PFS/specmatch/'
        idstring = 'TIC268301217_20200203_3556'

        specmatch_analyze(spectrum_path, wvsol_path=wvsol_path, region='Mgb1',
                          outdir=outdir, idstring=idstring)
