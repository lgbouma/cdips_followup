"""
Contents:
    get_tess_data
    explore_flux_lightcurves
    explore_mag_lightcurves
    explore_eleanor_lightcurves
    make_periodogram
"""
from glob import glob
import os, multiprocessing
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits

from astrobase import lcmath
from astrobase.periodbase import kbls
from astrobase.checkplot.png import _make_phased_magseries_plot
from astrobase.lcmath import sigclip_magseries, phase_magseries, phase_bin_magseries
from copy import deepcopy

from cdips.plotting import vetting_pdf as vp
from cdips.lcproc import mask_orbit_edges as moe
from cdips.lcproc import detrend as dtr
from cdips.utils.lcutils import _given_mag_get_flux

from numpy.polynomial.legendre import Legendre

from astrobase.services.tesslightcurves import (
    get_two_minute_spoc_lightcurves,
    get_hlsp_lightcurves,
    get_eleanor_lightcurves
)

from cdips_followup.paths import RESULTSDIR

def get_tess_data(ticid, outdir=None, cdips=0, spoc=0, eleanor=0, cdipspre=0):
    """
    High-level wrapper to download available TESS data for `ticid` to `outdir`.
    A few different formats are implemented.  If nothing is found, returns
    None.

    Args:
        ticid (str)
        outdir (str or None): path
    Returns:
        data (list): list of FITS data tables, per sector.
    """

    if outdir is None:
        outdir = os.path.join(RESULTSDIR, 'quicklooklc', f'TIC{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    if (not cdips) and (not spoc) and (not eleanor) and (not cdipspre):
        raise ValueError(
            'Expected at least one format of TESS LC to be requested.'
        )

    if cdipspre:
        # e.g., 5996151172781298304_llc.fits
        lcfiles = glob(os.path.join(outdir,'*_llc.fits'))

    if cdips:
        lcfiles = glob(os.path.join(outdir,'hlsp_cdips*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_hlsp_lightcurves(ticid, hlsp_products=['CDIPS'],
                                           download_dir=outdir, verbose=True)

    if spoc:
        lcfiles = glob(os.path.join(outdir,'mastDownload','TESS','*','tess*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_two_minute_spoc_lightcurves(ticid, download_dir=outdir)

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


def explore_mag_lightcurves(data, ticid, period=None, epoch=None):

    for yval in ['TFA1','TFA2','TFA3','IRM1','IRM2','IRM3','PCA1','PCA2','PCA3']:

        times, mags= [], []
        for ix, d in enumerate(data):

            savpath = (
                '../results/quicklooklc/TIC{}/mag_lightcurve_{}_{}.png'.
                format(ticid, yval, ix)
            )
            if os.path.exists(savpath):
                print('found {}, rewriting'.format(savpath))

            plt.close('all')
            f,ax = plt.subplots(figsize=(16,4))

            ax.scatter(d['TMID_BJD'], d[yval], c='k', s=5)
            times.append(d['TMID_BJD'])
            mags.append( d[yval] - np.nanmedian(d[yval]) )

            ax.set_xlabel('time [bjdtdb]')
            ax.set_ylabel(yval)
            ylim = ax.get_ylim()
            ax.set_ylim((max(ylim), min(ylim)))

            if not epoch is None:
                if np.min(d['TMID_BJD']) < 2450000 and epoch > 2450000:
                    epoch -= 2457000
                if np.min(d['TMID_BJD']) > 2450000 and epoch < 2450000:
                    epoch += 2457000

                tra_times = epoch + np.arange(-1000,1000,1)*period

                xlim = ax.get_xlim()
                ylim = ax.get_ylim()

                ax.vlines(tra_times, max(ylim), min(ylim), color='orangered',
                          linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
                ax.set_ylim((max(ylim), min(ylim)))
                ax.set_xlim(xlim)

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
            print('found {}, rewriting'.format(savpath))

        plt.close('all')
        f,ax = plt.subplots(figsize=(16,4))

        ax.scatter(stimes, smags, c='k', s=5)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        if not epoch is None:
            if np.min(d['TMID_BJD']) < 2450000 and epoch > 2450000:
                epoch -= 2457000
            if np.min(d['TMID_BJD']) > 2450000 and epoch < 2450000:
                epoch += 2457000

            tra_times = epoch + np.arange(-1000,1000,1)*period

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ax.vlines(tra_times, max(ylim), min(ylim), color='orangered',
                      linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
            ax.set_ylim((max(ylim), min(ylim)))
            ax.set_xlim(xlim)

        ax.set_xlabel('time [bjdtdb]')
        ax.set_ylabel('relative '+yval)
        ylim = ax.get_ylim()
        ax.set_ylim((max(ylim), min(ylim)))

        ax.set_title(ix)

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))


def _get_ylim(y_obs):
    iqr = (np.nanpercentile(y_obs, 75) -
           np.nanpercentile(y_obs, 25))
    ylower = np.nanmedian(y_obs) - 2.5*iqr
    yupper = np.nanmedian(y_obs) + 2.5*iqr
    return (ylower, yupper)


def explore_flux_lightcurves(data, ticid, outdir=None, period=None, epoch=None,
                             isspoc=True, detrend=False, window_length=None,
                             do_phasefold=0, badtimewindows=None, get_lc=False,
                             require_quality_zero=1):
    """
    Given a list of SPOC 2 minute data FITS tables, stitch them across sectors
    and make diagnostic plots.

    Args:

        data (list): from `get_tess_data`.

        ticid (str): TIC ID.

        outdir (str): diagnostic plots are written here. If None, goes to
        cdips_followup results directory.

    Optional kwargs:

        period, epoch (float): optional

        detrend (bool): optional, will apply standard wotan rotation-flattening.

        badtimewindows (list): to manually mask out, [(1656, 1658), (1662,
            1663)], for instance.

        get_lc (bool): if True, returns time and flux arrays.

        require_quality_zero (bool): if True, sets QUALITY==0, throwing out
        lots of data.
    """

    if not isspoc:
        raise NotImplementedError

    yval = 'PDCSAP_FLUX'

    if outdir is None:
        outdir = os.path.join(RESULTSDIR, 'quicklooklc', f'TIC{ticid}')

    times, fluxs= [], []
    for ix, d in enumerate(data):

        savpath = os.path.join(outdir, f'spoc_lightcurve_{ix}.png')
        if detrend:
            savpath = os.path.join(outdir, f'spoc_lightcurve_detrended_{ix}.png')

        plt.close('all')
        f,ax = plt.subplots(figsize=(16*2,4*1.5))

        if require_quality_zero:
            sel = (d['QUALITY'] == 0) & (d[yval] > 0)
            print(42*'.')
            print('WRN!: omitting all non-zero quality flags. throws out good data!')
            print(42*'.')
        else:
            sel = (d[yval] > 0)
        if badtimewindows is not None:
            for w in badtimewindows:
                sel &= ~(
                    (d['TIME'] > w[0])
                    &
                    (d['TIME'] < w[1])
                )

        x_obs = d['TIME'][sel]
        y_obs = d[yval][sel] / np.nanmedian(d[yval][sel])

        if detrend:
            ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)
            y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs)

        if detrend:
            ax.plot(x_obs, y_trend, c='r', lw=0.5, zorder=3)
        else:
            ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)

        times.append( x_obs )
        fluxs.append( y_obs )

        ax.set_xlabel('time [bjdtdb]')
        ax.set_ylabel(yval)
        ylim = ax.get_ylim()

        ax.set_title(ix)

        if detrend:
            _ylim = _get_ylim(y_trend)
        else:
            _ylim = _get_ylim(y_obs)

        ax.set_ylim(_ylim)

        if not epoch is None:
            tra_times = epoch + np.arange(-1000,1000,1)*period

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ax.set_ylim((min(ylim), max(ylim)))
            ax.vlines(tra_times, min(ylim), max(ylim), color='orangered',
                      linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
            ax.set_ylim((min(ylim), max(ylim)))
            ax.set_xlim(xlim)

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

    times = np.hstack(np.array(times).flatten())
    fluxs = np.hstack(np.array(fluxs).flatten())

    stimes, smags, _ = lcmath.sigclip_magseries(
        times, fluxs, np.ones_like(fluxs), sigclip=[8,3], iterative=True,
        magsarefluxes=True
    )

    savpath = os.path.join(
        outdir, 'spoc_lightcurve_{}_allsector.png'.  format(yval)
    )
    if detrend:
        savpath = os.path.join(
            outdir, 'spoc_lightcurve_detrended_{}_allsector.png'.  format(yval)
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

    if do_phasefold:

        assert isinstance(period, float) and isinstance(epoch, float)

        # use times and fluxs, instead of the sigma clipped thing.

        plt.close('all')
        fig, ax = plt.subplots(figsize=(4,3))

        #
        # ax: primary transit
        #
        tdur_guess = 7
        phasebin = 6e-3
        minbinelems = 2
        tdur_by_period = (tdur_guess/24)/period
        plotxlim = (-0.5, 0.5)
        plotylim = (0.989, 1.004)

        _make_phased_magseries_plot(ax, 0, times, fluxs,
                                    np.ones_like(fluxs)/1e4, period, epoch,
                                    True, True, phasebin, minbinelems,
                                    plotxlim, '', xliminsetmode=False,
                                    magsarefluxes=True, phasems=0.8,
                                    phasebinms=4.0, verbose=True)
        ax.set_ylim(plotylim)

        ax.vlines(1/6, min(plotylim), max(plotylim), color='orangered',
                  linestyle='--', zorder=-2, lw=1, alpha=0.8)
        ax.set_ylim(plotylim)

        dstr = 'detrended' if detrend else ''
        savpath = os.path.join(
            outdir, f'spoc_lightcurve_{dstr}_{yval}_allsector_phasefold.png'
        )

        fig.savefig(savpath, dpi=400, bbox_inches='tight')
        print(f'made {savpath}')

        csvpath = savpath.replace('png','csv')
        pd.DataFrame({
            'time': times, 'flux': fluxs
        }).to_csv(csvpath, index=False)
        print(f'made {csvpath}')

    if get_lc:
        return times, fluxs


def make_periodogram(data, ticid, pipeline, period_min=0.1, period_max=20,
                     manual_peak=None, samples_per_peak=50):

    if pipeline == 'spoc':
        time = data['TIME']
        flux = data['PDCSAP_FLUX']
        err = data['PDCSAP_FLUX_ERR']
        qual = data['QUALITY']
        sel = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) # (qual == 0)
        time, flux, err = time[sel], flux[sel], err[sel]
        med_flux = np.nanmedian(flux)
        flux /= med_flux
        err /= med_flux
    elif pipeline == 'cdips':
        time = data['time']
        flux = data['flux']
        err = data['err']

    ##########################################

    from astropy.timeseries import LombScargle
    from scipy.signal import find_peaks

    plt.close('all')

    fig, ax = plt.subplots()

    nterms_0=6
    ls = LombScargle(time, flux, err, nterms=nterms_0)
    freq, power = ls.autopower(minimum_frequency=1/period_max,
                               maximum_frequency=1/period_min,
                               samples_per_peak=samples_per_peak)

    period = 1/freq
    ax.plot(period, power)

    ls_freq_0 = freq[np.argmax(power)]
    ls_period_0 = 1/ls_freq_0
    ax.axvline(ls_period_0, alpha=0.4, lw=1, color='C0')
    ax.axvline(2*ls_period_0, alpha=0.4, lw=1, color='C0')
    ax.axvline(0.5*ls_period_0, alpha=0.4, lw=1, color='C0')

    if manual_peak:
        ax.axvline(manual_peak, alpha=0.4, lw=1, color='C1')
        ax.axvline(2*manual_peak, alpha=0.4, lw=1, color='C1')
        ax.axvline(0.5*manual_peak, alpha=0.4, lw=1, color='C1')

    peaks, props = find_peaks(power, height=1e-1)
    print(period[peaks])

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlabel('period [d]')
    ax.set_ylabel('power')

    outdir = '../results/quicklooklc/TIC{}'.format(ticid)
    savpath = os.path.join(outdir, 'ls_periodogram.png')
    fig.savefig(savpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(savpath))

    ##########################################

    flux_fit_0 = ls.model(time, ls_freq_0)
    flux_r1 = flux - flux_fit_0

    nterms_1=6
    ls = LombScargle(time, flux_r1, err, nterms=nterms_1)
    freq, power = ls.autopower(minimum_frequency=1/period_max,
                               maximum_frequency=1/period_min)

    period = 1/freq

    ls_freq_1 = freq[np.argmax(power)]
    ls_period_1 = 1/ls_freq_1

    flux_fit_1 = ls.model(time, ls_freq_1)
    flux_r2 = flux_r1 - flux_fit_1

    ##########

    nterms_2=6
    ls = LombScargle(time, flux_r2, err, nterms=nterms_2)
    freq, power = ls.autopower(minimum_frequency=1/period_max,
                               maximum_frequency=1/period_min)

    period = 1/freq
    ls_freq_2 = freq[np.argmax(power)]
    ls_period_2 = 1/ls_freq_2

    flux_fit_2 = ls.model(time, ls_freq_2)
    flux_r3 = flux_r2 - flux_fit_2



    ##########################################

    plt.close('all')
    fig, axs = plt.subplots(nrows=4, figsize=(14,10), sharex=True)

    axs[0].scatter(time, flux, c='k', s=1, zorder=1)
    axs[0].plot(time, flux_fit_0, zorder=2, lw=1)
    axs[0].set_title(f'model: {nterms_0} fourier terms at {ls_period_0:.6f} days',
                     fontsize='x-small')

    axs[1].scatter(time, flux_r1, c='k', s=1, zorder=1)
    axs[1].plot(time, flux_fit_1, zorder=2, lw=1)
    axs[1].set_title(f'model: {nterms_1} fourier terms at {ls_period_1:.6f} days',
                     fontsize='x-small')

    axs[2].scatter(time, flux_r2, c='k', s=1, zorder=1)
    axs[2].plot(time, flux_fit_2, zorder=2, lw=1)
    axs[2].set_title(f'model: {nterms_2} fourier terms at {ls_period_2:.6f} days',
                     fontsize='x-small')

    axs[3].scatter(time, flux_r3, c='k', s=1, zorder=1)
    axs[3].plot(time, flux_fit_2-flux_fit_2, zorder=2, lw=1)
    axs[3].set_title('', fontsize='x-small')

    axs[0].set_ylabel('raw')
    axs[1].set_ylabel('resid 1')
    axs[2].set_ylabel('resid 2')
    axs[3].set_ylabel('resid 3')
    axs[3].set_xlabel('btjd')

    fig.tight_layout()

    savpath = os.path.join(outdir, 'ls_models_resids.png')
    fig.savefig(savpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(savpath))




    ##########################################
    show_tls = 0

    if show_tls:
        from transitleastsquares import transitleastsquares
        model = transitleastsquares(time, flux)
        results = model.power(period_min=period_min, period_max=period_max,
                              M_star_min=0.1, M_star_max=5, R_star=0.5, M_star=0.5)

        print(42*'-')
        print(
            't0: {}'.format(results.T0)+
            '\nperiod: {}'.format(results.period)+
            '\nperiod_unc: {}'.format(results.period_uncertainty)
        )
        print(42*'-')

        ##########################################

        plt.close('all')
        fig = plt.figure()
        ax = plt.gca()
        ax.axvline(results.period, alpha=0.4, lw=3)
        plt.xlim(np.min(results.periods), np.max(results.periods))
        for n in range(2, 10):
            ax.axvline(n*results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax.axvline(results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        plt.ylabel(r'SDE')
        plt.xlabel('Period (days)')
        plt.plot(results.periods, results.power, color='black', lw=0.5)
        plt.xlim(0, max(results.periods))

        savpath = os.path.join(outdir, 'tls_periodogram.png')
        fig.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

        ##########################################

        plt.close('all')
        fig = plt.figure()
        plt.scatter(results.folded_phase, results.folded_y, color='gray', s=2,
                    alpha=0.8, zorder=4, linewidths=0)

        pd = phase_bin_magseries(results.folded_phase, results.folded_y,
                                 binsize=0.01)
        plt.scatter(pd['binnedphases'], pd['binnedmags'], color='black', s=8,
                    alpha=1, zorder=5, linewidths=0)

        plt.title(f'period {results.period:.8f}d')
        plt.xlabel('Phase')
        plt.ylabel('Relative flux')
        plt.xticks([0, 0.25, 0.5, 0.75, 1])
        plt.grid(which='major', axis='both', linestyle='--', zorder=-3,
                 alpha=0.5, color='gray')

        pct_80 = np.percentile(results.model_folded_model, 80)
        pct_20 = np.percentile(results.model_folded_model, 20)
        center = np.nanmedian(results.model_folded_model)
        delta_y = (10/6)*np.abs(pct_80 - pct_20)
        plt.ylim(( center-0.7*delta_y, center+0.7*delta_y ))

        savpath = os.path.join(outdir, 'phasefold.png')
        fig.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))




