import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from astropy import units as u

from cdips_followup import __path__

from cdips_followup.spectools import (
    read_veloce, read_pfs,
    fit_continuum,
    measure_vsini,
    given_deltawvlen_get_vsys,
    VELOCE_ORDERS, VELOCE_VSINI_ORDERS_VS_NEXTGEN
)

from stringcheese.plotutils import savefig, format_ax

##########
# config #
##########
DEBUG = 1
TEST = 0

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = os.path.join(os.path.dirname(__path__[0]), 'results/spec_analysis')
TESTOUTDIR = os.path.join(OUTDIR, 'tests')
if not os.path.exists(TESTOUTDIR):
    os.mkdir(TESTOUTDIR)

#########
# tests #
#########

def test_veloce_continuum_fit():

    specname = '20200131_837.01_Bouma_final_combined.fits'
    spectrum_path = os.path.join(
        DATADIR, 'Veloce', specname
    )

    flx_2d, wav_2d = read_veloce(spectrum_path, start=400, end=-400)

    for o in VELOCE_ORDERS:

        flx, wav = flx_2d[o, :], wav_2d[o, :]
        wav, flx = wav[::-1], flx[::-1] # wavelength increasing order
        cont_flx, cont_norm_spec = fit_continuum(flx, wav, instrument='Veloce')

        #
        # plot the results
        #
        outpath = os.path.join(TESTOUTDIR, specname.replace(
            '.fits', '_order{}_continuum_fit.png'.format(str(o).zfill(2))))

        plt.close('all')
        f,axs = plt.subplots(nrows=2, ncols=1, figsize=(6,4))

        axs[0].plot(wav, flx, c='k', zorder=3, lw=1)
        axs[0].plot(wav, cont_flx, c='r', zorder=2, lw=1)

        sel = (
            (cont_flx > 0) &
            (cont_norm_spec.flux < 2) &
            (cont_norm_spec.flux > -1)
        )
        axs[1].plot(cont_norm_spec.wavelength[sel], cont_norm_spec.flux[sel],
                    c='k', lw=1)

        axs[0].set_ylabel('flux')
        axs[1].set_ylabel('contnorm flux')
        axs[1].set_xlabel('wavelength [angstrom]')

        for ax in axs:
            format_ax(ax)
        savefig(f, outpath, writepdf=False)


def test_measure_veloce_vsini():

    specname = '20200131_837.01_Bouma_final_combined.fits'
    targetname = '20200131_837.01'
    teff = 6300

    vsini, shift, gamma = measure_veloce_vsini(specname, targetname, teff,
                                               TESTOUTDIR)


def main():

    if DEBUG:
        given_deltawvlen_get_vsys(deltawvlen=0.167*u.AA, wvlen_0=6500*u.AA)

    if TEST:
        test_measure_veloce_vsini()
        test_veloce_continuum_fit()

if __name__ == "__main__":
    main()
