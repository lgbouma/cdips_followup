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
    get_Li_6708_EW,
    given_vsys_get_li_target_wv, viz_compare
)

# spectrum_path = '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/Veloce/20200130_837.01_Bouma_final_combined.fits'

def inspect_pfs(nightstr, targetline, xlim=None):
    # NIST has the Li I resonance doublet listed with one transition at 6707.76
    # and the other at 6707.91 A.

    targetid = 'TIC268301217.01'
    datestr = '20200203_{}'.format(nightstr)
    xshift = 2.30
    # note: wavelength soln might change by night...
    spectrum_path = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.{}'.
        format(nightstr)
    )
    wvsol_path = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
    )

    flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)

    if xlim == 'fullorder':
        xlim = None
        vlines = None
        names = None
        if targetline == 'Halpha':
            target_wav = 6562.8
        elif targetline == 'LiI':
            target_wav = 6707.85

    elif xlim == 'assign':
        if targetline == 'Halpha':
            target_wav = 6562.8
            vlines = [target_wav]
            names = ['Halpha']
        elif targetline == 'LiI':
            target_wav = 6707.85
            # FeI at 6703.58 and 6705.1 (Berger+18 Fig3).
            # FeI at 6707.44, Li I at 6707.76, 6707.91
            # CaI lambda at ~6718.
            vlines = [6703.58, 6705.1, 6707.44, 6707.76, 6707.91, 6718]
            names = ['FeI', 'FeI', 'FeI', 'Li', '', 'CaI$\lambda$']
        delta_wav = 5
        xlim = [target_wav-delta_wav, target_wav+delta_wav]

    # retrieve the order corresponding to target wavelength
    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    order = int(np.argwhere((_preorder != wav_2d.shape[1]-1) & (_preorder != 0)))

    flx, wav = flx_2d[order, :], wav_2d[order, :]
    shiftstr = ''
    if isinstance(xshift, (float, int)):
        wav = deepcopy(wav) - xshift
        shiftstr = '_shift{:.2f}'.format(float(xshift))

    xstr = ''
    if isinstance(xlim, list) or isinstance(xlim, tuple):
        xstr = '_{}_xlim{:.1f}-{:.1f}'.format(
            targetline, xlim[0], xlim[1]
        )

    outname = '{}_{}_order{}{}{}.png'.format(
        datestr, targetid, str(order).zfill(2), xstr, shiftstr
    )

    outpath = '../results/spec_viz/PFS/{}'.format(outname)

    viz_1d_spectrum(flx, wav, outpath, xlim=xlim, vlines=vlines, names=names)


if __name__ == "__main__":

    do_viz = 0
    do_li_ew = 1

    if do_viz:
        for xlim in ['assign', 'fullorder']:
            inspect_pfs('3555', 'Halpha', xlim)
            inspect_pfs('3555', 'LiI', xlim)
            inspect_pfs('3556', 'Halpha', xlim)
            inspect_pfs('3556', 'LiI', xlim)

    if do_li_ew:
        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outpath = '../results/spec_viz/PFS/3556_Li_EW_shift2.30.png'

        get_Li_6708_EW(spectrum_path, wvsol_path=wvsol_path, xshift=2.30,
                       outpath=outpath)
