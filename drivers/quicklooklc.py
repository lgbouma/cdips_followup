"""
Given a TICID, make diagnostic plots to help do an initial analysis of any
light-curve of interest. (E.g., phase-folding, river plots, different wide
figure sizes).
"""

from glob import glob
import os, multiprocessing
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from astrobase import lcmath
from astrobase.periodbase import kbls
from astrobase.checkplot.png import _make_phased_magseries_plot
from astrobase.lcmath import phase_bin_magseries
from copy import deepcopy

from cdips.plotting import vetting_pdf as vp
from cdips.lcproc import mask_orbit_edges as moe

from numpy.polynomial.legendre import Legendre

from astrobase.services.tesslightcurves import (
    get_two_minute_spoc_lightcurves,
    get_hlsp_lightcurves,
    get_eleanor_lightcurves
)


def get_data(ticid, cdips=1, spoc=0, eleanor=0, cdipspre=0):

    outdir = '../results/quicklooklc/TIC{}'.format(ticid)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if cdipspre:
        # e.g., 5996151172781298304_llc.fits
        lcfiles = glob(os.path.join(outdir,'*_llc.fits'))

    if cdips:
        lcfiles = glob(os.path.join(outdir,'hlsp_cdips*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_hlsp_lightcurves(tic_id, hlsp_products=['CDIPS'],
                                           download_dir=outdir, verbose=True)

    if spoc:
        raise NotImplementedError

    if eleanor:
        lcfiles = glob(os.path.join(outdir,'hlsp_eleanor*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_eleanor_lightcurves(ticid, download_dir=outdir)

    if lcfiles is None:
        return None

    data = []
    for f in lcfiles:
        hdul = fits.open(f)
        data.append(hdul[1].data)

    return data


def explore_mag_lightcurves(data, ticid):

    for yval in ['TFA1','TFA2','TFA3','IRM1','IRM2','IRM3']:

        times, mags= [], []
        for ix, d in enumerate(data):

            savpath = (
                '../results/quicklooklc/TIC{}/mag_lightcurve_{}_{}.png'.
                format(ticid, yval, ix)
            )
            if os.path.exists(savpath):
                print('found {}, continue'.format(savpath))
                continue

            plt.close('all')
            f,ax = plt.subplots(figsize=(16,4))

            ax.scatter(d['TMID_BJD'], d[yval], c='k', s=5)
            times.append(d['TMID_BJD'])
            mags.append( d[yval] - np.nanmedian(d[yval]) )

            ax.set_xlabel('time [bjdtdb]')
            ax.set_ylabel(yval)
            ylim = ax.get_ylim()
            ax.set_ylim((max(ylim), min(ylim)))

            ax.set_title(ix)

            f.savefig(savpath, dpi=300, bbox_inches='tight')
            print('made {}'.format(savpath))

        times = np.hstack(np.array(times).flatten())
        mags = np.hstack(np.array(mags).flatten())

        stimes, smags, _ = lcmath.sigclip_magseries(
            times, mags, np.ones_like(mags), sigclip=[20,3], iterative=True
        )

        savpath = (
            '../results/quicklooklc/TIC{}/mag_lightcurve_{}_allsector.png'.
            format(ticid, yval)
        )
        if os.path.exists(savpath):
            print('found {}, continue'.format(savpath))
            continue

        plt.close('all')
        f,ax = plt.subplots(figsize=(16,4))

        ax.scatter(stimes, smags, c='k', s=5)

        # period = 11.69201165
        # epoch = 2458642.44550000
        # tra_times = epoch + np.arange(-100,100,1)*period

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # ax.set_ylim((min(ylim), max(ylim)))
        # ax.vlines(tra_times, min(ylim), max(ylim), color='orangered',
        #           linestyle='--', zorder=-2, lw=2, alpha=0.3)
        # ax.set_ylim((min(ylim), max(ylim)))
        # ax.set_xlim(xlim)

        ax.set_xlabel('time [bjdtdb]')
        ax.set_ylabel('relative '+yval)
        ylim = ax.get_ylim()
        ax.set_ylim((max(ylim), min(ylim)))

        ax.set_title(ix)

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))


def explore_eleanor_lightcurves(data, ticid, period=None, epoch=None):

    yval = 'CORR_FLUX'

    times, fluxs= [], []
    for ix, d in enumerate(data):

        outdir = '../results/quicklooklc/TIC{}'.format(ticid)
        savpath = os.path.join(outdir, 'eleanor_lightcurve_{}.png'.format(ix))

        plt.close('all')
        f,ax = plt.subplots(figsize=(16,4))

        sel = (d['QUALITY'] == 0) & (d[yval] > 0)

        ax.scatter(d['TIME'][sel], d[yval][sel], c='k', s=5)
        times.append(d['TIME'][sel])
        fluxs.append( d[yval][sel] / np.nanmedian(d[yval][sel]) )

        ax.set_xlabel('time [bjdtdb]')
        ax.set_ylabel(yval)
        ylim = ax.get_ylim()

        ax.set_title(ix)

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

    times = np.hstack(np.array(times).flatten())
    fluxs = np.hstack(np.array(fluxs).flatten())

    stimes, smags, _ = lcmath.sigclip_magseries(
        times, fluxs, np.ones_like(fluxs), sigclip=[8,3], iterative=True,
        magsarefluxes=True
    )

    savpath = os.path.join(
        outdir, 'eleanor_lightcurve_{}_allsector.png'.  format(yval)
    )

    plt.close('all')
    f,ax = plt.subplots(figsize=(16,4))

    ax.scatter(stimes, smags, c='k', s=5)

    if not epoch is None:
        tra_times = epoch + np.arange(-1000,1000,1)*period - 2457000

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        ax.set_ylim((min(ylim), max(ylim)))
        ax.vlines(tra_times, min(ylim), max(ylim), color='orangered',
                  linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
        ax.set_ylim((min(ylim), max(ylim)))
        ax.set_xlim(xlim)

    ax.set_xlabel('time [bjdtdb]')
    ax.set_ylabel('relative '+yval)

    ax.set_title(ix)

    f.savefig(savpath, dpi=400, bbox_inches='tight')
    print('made {}'.format(savpath))



def main():

    ticid = '279819212' # '268016868' # '166527623' # '34488204' # '165722603' # '268301217'
    # ticid = '220322660'
    # optional #
    period = None #1.12041120
    epoch = None #2458715.57911000

    cdips = 0
    spoc = 0
    eleanor = 0
    cdipspre = 1

    do_mag_lcs = 1
    do_eleanor_lcs = 0
    do_flux_lcs = 0
    do_detrending = 0
    do_pf = 0
    do_riverplot = 0

    data = get_data(ticid, cdips=cdips, spoc=spoc, cdipspre=cdipspre,
                    eleanor=eleanor)

    if do_eleanor_lcs:
        explore_eleanor_lightcurves(data, ticid, period=period, epoch=epoch)

    if do_mag_lcs:
        explore_mag_lightcurves(data, ticid)

    if do_flux_lcs:
        explore_flux_lightcurves(data, iscdips=cdips)

    if do_pf:
        do_phasefolds(data)

    if do_detrending:
        #for wl in [99]:
        for wl in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            detrend_lightcurve_wotan(data, window_length=wl, iscdips=cdips)
    if do_riverplot:
        make_riverplot(data)


if __name__ == "__main__":
    main()
