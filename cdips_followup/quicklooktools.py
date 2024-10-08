"""
Contents:
    get_tess_data
    get_kepler_data
    explore_flux_lightcurves
    explore_mag_lightcurves
    explore_eleanor_lightcurves
    make_periodogram
    make_riverplot
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
import cdips.utils.lcutils as lcu

from wotan import slide_clip

from os.path import join

from numpy.polynomial.legendre import Legendre

from astrobase.services.tesslightcurves import (
    get_two_minute_spoc_lightcurves,
    get_hlsp_lightcurves,
    get_eleanor_lightcurves
)

from cdips_followup.paths import RESULTSDIR

LCKEYDICT = {
    'spoc': {
        'flux': 'PDCSAP_FLUX', 'time': 'TIME', 'quality': 'QUALITY',
        'time_offset':2457000, 'prov': 'spoc', 'inst': 'tess'
    },
    'kepler': {
        'flux': 'PDCSAP_FLUX', 'time': 'TIME', 'quality': 'SAP_QUALITY',
        'time_offset':2454833, 'prov': 'spoc', 'inst': 'kepler'
    },
    #
    # QLP: per Huang+2020, "SAP_FLUX" is the optimal aperture
    # time is in BTJD
    #
    'qlp': {
        'flux': 'SAP_FLUX', 'time': 'TIME', 'quality': 'QUALITY',
        'time_offset':2457000, 'prov': 'mit', 'inst': 'tess'
    },
    'cdips': {
        'flux': 'PCA1', 'time': 'TMID_BJD', 'quality': 'IRQ1',
        'time_offset':0, 'prov':'cdips', 'inst':'tess'
    }
}



def get_tess_data(ticid, outdir=None, cdips=0, spoc=0, eleanor=0, cdipspre=0,
                  qlp=0, unpopular=0):
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
        outdir = join(RESULTSDIR, 'quicklooklc', f'TIC{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    if (
        (not cdips)
        and (not spoc)
        and (not eleanor)
        and (not cdipspre)
        and (not qlp)
        and (not unpopular)
    ):
        print(
            '..........'
            'WRN! Expected at least one format of TESS LC to be requested. '
            'Returning None'
            '..........'
        )
        return None

    if cdipspre:
        # e.g., 5996151172781298304_llc.fits
        lcfiles = glob(join(outdir,'*_llc.fits'))

    if cdips:
        lcfiles = glob(join(outdir,'hlsp_cdips*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_hlsp_lightcurves(ticid, hlsp_products=['CDIPS'],
                                           download_dir=outdir, verbose=True)

    if spoc:
        lcfiles = glob(join(outdir,'*','TESS',f'*tess*{ticid}*','tess*lc.fits'))
        if len(lcfiles) == 0:
            lcfiles = get_two_minute_spoc_lightcurves(ticid, download_dir=outdir)

    if eleanor:
        lcfiles = glob(join(outdir,'hlsp_eleanor*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_eleanor_lightcurves(ticid, download_dir=outdir)

    if qlp:
        # manually downloaded from MAST portal.
        lcfiles = glob(join(outdir, 'QLP', 'hlsp_qlp*fits'))
        if len(lcfiles) == 0:
            lcfiles = get_hlsp_lightcurves(ticid, hlsp_products=['QLP'],
                                           download_dir=outdir, verbose=True)

    if unpopular:
        from astrobase.services.tesslightcurves import get_unpopular_lightcurve
        from cdips_followup.paths import FFICACHEDIR
        lcfiles = glob(join(outdir,f'*tess*{ticid}*lc.csv'))
        if len(lcfiles) == 0:
            lcfiles = get_unpopular_lightcurve(
                ticid, ffi_dir=FFICACHEDIR, overwrite=False, lc_dir=outdir
            )

    if lcfiles is None or lcfiles == []:
        return None, None

    data = []
    hdrs = []
    for f in lcfiles:
        if f.endswith(".fits"):
            hdul = fits.open(f)
            data.append(hdul[1].data)
            hdrs.append(hdul[0].header)
        elif f.endswith(".csv"):
            df = pd.read_csv(f)
            data.append(df)
            hdrs.append([])

    return data, hdrs


def get_kepler_data(ticid, outdir=None):
    """
    Reads in Kepler data for anything downloaded off MAST (manually).

    Args:
        ticid (str)
        outdir (str or None): path
    Returns:
        data (list): list of FITS data tables, per sector.
    """

    if outdir is None:
        outdir = join(RESULTSDIR, 'quicklooklc', f'TIC{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    lcfiles = glob(
        join(outdir,'MAST*','Kepler','kplr*lc*','kplr*llc.fits')
    )
    if len(lcfiles) == 0:
        raise NotImplementedError

    data = []
    for f in lcfiles:
        hdul = fits.open(f)
        data.append(hdul[1].data)

    return data


def explore_eleanor_lightcurves(data, hdrs, ticid, period=None, epoch=None,
                                require_quality_zero=0, detrend=0,
                                do_phasefold=0, do30mincadence=1):

    ykey = 'CORR_FLUX'

    times, fluxs= [], []
    for ix, d in enumerate(data):

        hdr = hdrs[ix]
        sectorstr = f"s{str(hdr['SECTOR']).zfill(4)}"

        plt.close('all')
        if not detrend:
            f,ax = plt.subplots(figsize=(16,4))
            axs = [ax]
        else:
            f,axs = plt.subplots(figsize=(16,8), nrows=2, sharex=True)

        sel = (d['QUALITY'] == 0) & (d[ykey] > 0)
        if require_quality_zero:
            if 'QUALITY' in d.names:
                qualkey = 'QUALITY'
            else:
                qualkey = 'SAP_QUALITY'
            sel = (d[qualkey] == 0) & (d[ykey] > 0)
            print(42*'.')
            print('WRN!: omitting all non-zero quality flags. throws out good data!')
            print(42*'.')
        else:
            sel = (d[ykey] > 0)

        x_obs = d['TIME'][sel]
        y_obs = d[ykey][sel]

        cadence_min = np.nanmedian(np.diff(d['TIME'][sel]))*24*60

        # Years >=3
        binnedflag = ''
        if np.round(cadence_min) < 25:
            from astrobase.lcmath import time_bin_magseries
            bd = time_bin_magseries(x_obs, y_obs, binsize=1800, minbinelems=3)
            x_obs = bd['binnedtimes']
            y_obs = bd['binnedmags']
            binnedflag = '_30minbinned'

        outdir = f'../results/quicklooklc/TIC{ticid}'
        savpath = join(
            outdir, f'eleanor_lightcurve_{sectorstr}{binnedflag}.png'
        )
        if detrend:
            savpath = join(
                outdir, f'eleanor_lightcurve_detrended_{sectorstr}{binnedflag}.png'
            )

        if detrend:
            axs[0].scatter(x_obs, y_obs/np.nanmedian(y_obs), c='k', s=4,
                           zorder=5)

            # # default pspline detrending
            if detrend=='pspline':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs)

            # in some cases, might prefer the biweight
            elif detrend == 'biweight':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs,
                                                  method='biweight', cval=5,
                                                  window_length=0.5,
                                                  break_tolerance=0.5)

            elif detrend == 'minimal':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs,
                                                  method='biweight', cval=5,
                                                  window_length=5,
                                                  break_tolerance=0.5)

            elif detrend == 'best':
                # NOTE: window_length plays an important role in what you get!!
                dtr_method, break_tolerance, window_length = 'best', 0.5, 0.40
                dtr_dict = {'method':dtr_method,
                            'break_tolerance':break_tolerance,
                            'window_length':window_length}

                starid = f"TIC{ticid}_{sectorstr}_w{window_length}_eleanor"
                from cdips.lcproc.find_planets import run_periodograms_and_detrend
                r, dtr_time, dtr_flux, dtr_stages_dict = run_periodograms_and_detrend(
                    starid, x_obs, y_obs, dtr_dict, return_extras=True,
                    magisflux=True, transit_depth_min=100e-6, orbitgap=0.5
                )
                trend_flux = dtr_stages_dict['trend_flux']
                trend_time = dtr_stages_dict['trend_time']

            else:
                raise NotImplementedError

        if detrend == 'best':
            axs[0].plot(trend_time, trend_flux, c='r', lw=0.5, zorder=3)
            axs[1].scatter(dtr_time, dtr_flux, c='k', s=4, zorder=4)
        elif detrend:
            axs[0].plot(x_obs, y_trend, c='r', lw=0.5, zorder=3)
        else:
            ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)

        if not epoch is None:
            tra_times = epoch + np.arange(-1000,1000,1)*period - 2457000

            for ax in axs:
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()

                ax.set_ylim((min(ylim), max(ylim)))
                ax.vlines(tra_times, min(ylim), max(ylim), color='orangered',
                          linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
                ax.set_ylim((min(ylim), max(ylim)))
                ax.set_xlim(xlim)

        for ix, ax in enumerate(axs):
            ax.set_ylabel(ykey)
            if ix == 0:
                ax.set_title(sectorstr)
            if len(axs) > 1:
                if ix != len(axs)-1:
                    continue
            ax.set_xlabel('time [bjdtdb]')

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

        times.append(x_obs)
        fluxs.append( y_obs / np.nanmedian(y_obs) )

    times = np.hstack(np.array(times).flatten())
    fluxs = np.hstack(np.array(fluxs).flatten())

    stimes, smags, _ = lcmath.sigclip_magseries(
        times, fluxs, np.ones_like(fluxs), sigclip=[5,2.5], iterative=True,
        magsarefluxes=True
    )

    savpath = join(
        outdir, f'eleanor_lightcurve_{ykey}_allsector.png'
    )
    if detrend:
        savpath = join(
            outdir, f'eleanor_lightcurve_detrended_{ykey}_allsector.png'
        )

    # do the sigma clipped
    x_obs, y_obs = stimes, smags

    plt.close('all')
    f,ax = plt.subplots(figsize=(16,4))

    ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)
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
    ax.set_ylabel('relative '+ykey)

    ax.set_title(ix)

    f.savefig(savpath, dpi=400, bbox_inches='tight')
    print('made {}'.format(savpath))

    if do_phasefold:

        assert isinstance(period, float) and isinstance(epoch, float)

        #
        # ax: primary transit
        #
        phasebin = 1e-2
        minbinelems = 2
        plotxlims = [(-0.5, 0.5), (-0.05,0.05)]
        xlimstrs = ['xwide','xnarrow']
        plotylim = None # (0.994, 1.005)
        do_vlines = False

        for plotxlim, xstr in zip(plotxlims, xlimstrs):

            plt.close('all')
            fig, ax = plt.subplots(figsize=(4,3))

            _make_phased_magseries_plot(ax, 0, x_obs, y_obs,
                                        np.ones_like(fluxs)/1e4, period, epoch,
                                        True, True, phasebin, minbinelems,
                                        plotxlim, '', xliminsetmode=False,
                                        magsarefluxes=True, phasems=0.8,
                                        phasebinms=4.0, verbose=True)
            if plotylim is not None:
                ax.set_ylim(plotylim)

            if do_vlines:
                ax.vlines(1/6, min(plotylim), max(plotylim), color='orangered',
                          linestyle='--', zorder=-2, lw=1, alpha=0.8)
                ax.set_ylim(plotylim)

            dstr = 'detrended' if detrend else ''
            savpath = join(
                outdir, f'eleanor_lightcurve_{dstr}_{ykey}_{xstr}_allsector_phasefold.png'
            )

            fig.savefig(savpath, dpi=400, bbox_inches='tight')
            print(f'made {savpath}')

        csvpath = savpath.replace('png','csv')
        # sigma clipped and detrended
        pd.DataFrame({
            'time': x_obs, 'flux': y_obs
        }).to_csv(csvpath, index=False)
        print(f'made {csvpath}')




def explore_mag_lightcurves(data, hdrs, ticid, period=None, epoch=None,
                            do_phasefold=False):

    for ykey in ['IRM1','IRM2','IRM3','PCA1','PCA2','PCA3','TFA1','TFA2','TFA3','PCA1']:

        times, mags= [], []
        for ix, d in enumerate(data):

            hdr = hdrs[ix]

            sectorstr = f"s{str(hdr['SECTOR']).zfill(4)}"

            savpath = (
                f'../results/quicklooklc/TIC{ticid}/'
                f'mag_lightcurve_{ykey}_{sectorstr}.png'
            )
            if os.path.exists(savpath):
                print('found {}, rewriting'.format(savpath))

            plt.close('all')
            f,ax = plt.subplots(figsize=(16,4))

            ax.scatter(d['TMID_BJD'], d[ykey], c='k', s=5)
            times.append(d['TMID_BJD'])
            mags.append( d[ykey] - np.nanmedian(d[ykey]) )

            ax.set_xlabel('time [bjdtdb]')
            ax.set_ylabel(ykey)
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
            format(ticid, ykey)
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
        ax.set_ylabel('relative '+ykey)
        ylim = ax.get_ylim()
        ax.set_ylim((max(ylim), min(ylim)))

        ax.set_title(ix)

        f.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

    x_obs = deepcopy(stimes)
    y_obs = deepcopy(smags)

    if do_phasefold:

        assert isinstance(period, float) and isinstance(epoch, (float,int))

        #
        # ax: primary transit
        #
        phasebin = 1e-2
        minbinelems = 2
        plotxlims = [(-0.5, 0.5), (-0.05,0.05)]
        xlimstrs = ['xwide','xnarrow']
        plotylim = None # (0.994, 1.005)
        do_vlines = False

        for plotxlim, xstr in zip(plotxlims, xlimstrs):

            plt.close('all')
            fig, ax = plt.subplots(figsize=(4,3))

            _make_phased_magseries_plot(ax, 0, x_obs, y_obs,
                                        np.ones_like(y_obs)/1e4, period, epoch,
                                        True, True, phasebin, minbinelems,
                                        plotxlim, '', xliminsetmode=False,
                                        magsarefluxes=True, phasems=0.8,
                                        phasebinms=4.0, verbose=True)
            if plotylim is not None:
                ax.set_ylim(plotylim)

            if do_vlines:
                ax.vlines(1/6, min(plotylim), max(plotylim), color='orangered',
                          linestyle='--', zorder=-2, lw=1, alpha=0.8)
                ax.set_ylim(plotylim)

            detrend = False
            dstr = 'detrended' if detrend else ''
            outdir = os.path.dirname(savpath)
            savpath = join(
                outdir, f'cdips_lightcurve_{dstr}_{ykey}_{xstr}_allsector_phasefold.png'
            )

            fig.savefig(savpath, dpi=400, bbox_inches='tight')
            print(f'made {savpath}')

        csvpath = savpath.replace('png','csv')
        # sigma clipped and detrended
        pd.DataFrame({
            'time': x_obs, 'flux': y_obs
        }).to_csv(csvpath, index=False)
        print(f'made {csvpath}')


def _get_ylim(y_obs):
    assert len(y_obs)>0
    iqr = (np.nanpercentile(y_obs, 75) -
           np.nanpercentile(y_obs, 25))
    ylower = np.nanmedian(y_obs) - 2.5*iqr
    yupper = np.nanmedian(y_obs) + 2.5*iqr
    return (ylower, yupper)


def explore_flux_lightcurves(
    data, hdrs, ticid, outdir=None, period=None, epoch=None, pipeline=None,
    detrend=False, window_length=None, do_phasefold=0, badtimewindows=None,
    get_lc=False, require_quality_zero=1, forceylim=None, normstitch=True,
    slideclipdict={'window_length':1, 'high':3, 'low':20},
    mask_orbit_edges=False, bintime=None
):
    """
    Given a list of SPOC 2 minute data FITS tables, stitch them across sectors
    and make diagnostic plots.

    Args:

        data (list): from `get_tess_data`, contents [hdulistA[1].data,
            hdulistB[1].data], etc..

        ticid (str): TIC ID.

        pipeline (str): one of ['cdips', 'spoc', 'eleanor', 'cdipspre',
        'kepler', 'qlp'].  This is used to access the flux, provenance, etc.

        outdir (str): diagnostic plots are written here. If None, goes to
        cdips_followup results directory.

    Optional kwargs:

        period, epoch (float): optional

        detrend (bool, or string): 'biweight' or 'pspline' accepted. Default
        parameters assumed for each.

        badtimewindows (list): to manually mask out, [(1656, 1658), (1662,
            1663)], for instance.

        get_lc (bool): if True, returns time and flux arrays.

        require_quality_zero (bool): if True, sets QUALITY==0, throwing out
        lots of data.

        normstitch (bool): normalize flux across sectors s.t. the relative
        amplitude remains fixed.

        slideclipdict (dict): e.g., {'window_length':1, 'high':3, 'low':8} for
        1 day sliding window, exclude +3MAD from median above and -8MAD from
        median below.
    """

    assert isinstance(data, list), 'Expected list of FITStables.'

    if pipeline not in ['spoc', 'kepler', 'qlp', 'cdips']:
        raise NotImplementedError

    if isinstance(epoch, float):
        if epoch < 2450000:
            raise ValueError(f'Expected epoch in BJDTDB. Got epoch={epoch:.6f}.')

    ykey = LCKEYDICT[pipeline]['flux']
    xkey = LCKEYDICT[pipeline]['time']
    qualkey = LCKEYDICT[pipeline]['quality']
    prov = LCKEYDICT[pipeline]['prov']
    inst = LCKEYDICT[pipeline]['inst']

    if outdir is None:
        outdir = join(RESULTSDIR, 'quicklooklc', f'TIC{ticid}')

    times, fluxs= [], []
    for ix, d in enumerate(data):

        hdr = hdrs[ix]
        sectorstr = f"s{str(hdr['SECTOR']).zfill(4)}"

        savpath = join(
            outdir, f'TIC{ticid}_{prov}_{inst}_lightcurve_{sectorstr}.png'
        )
        if detrend:
            savpath = join(
                outdir,
                f'TIC{ticid}_{prov}_{inst}_lightcurve_{detrend}_{sectorstr}.png'
            )

        plt.close('all')
        f,ax = plt.subplots(figsize=(16*2,4*1.5))

        if require_quality_zero:
            okkey = 0 if pipeline in 'spoc,kepler,qlp'.split(',') else 'G'
            sel = (d[qualkey] == okkey) & (d[ykey] > 0)
            print(42*'.')
            print('WRN!: omitting all non-zero quality flags. throws out good data!')
            print(42*'.')
        else:
            sel = (d[ykey] > 0)
        if badtimewindows is not None:
            for w in badtimewindows:
                sel &= ~(
                    (d[xkey] > w[0])
                    &
                    (d[xkey] < w[1])
                )

        # correct time column to BJD_TDB
        x_offset = LCKEYDICT[pipeline]['time_offset']
        x_obs = d[xkey][sel] + x_offset

        # get the median-normalized flux
        y_obs = d[ykey][sel]
        if pipeline == 'cdips':
            y_obs, _ = _given_mag_get_flux(y_obs, y_obs*1e-3)
        y_obs /= np.nanmedian(y_obs)

        if mask_orbit_edges:
            x_obs, y_obs, _ = moe.mask_orbit_start_and_end(
                x_obs, y_obs, raise_expectation_error=False, orbitgap=0.7,
                orbitpadding=12/(24),
                return_inds=True
            )

        # slide clip -- remove outliers with windowed stdevn removal
        if isinstance(slideclipdict, dict):
            y_obs = slide_clip(x_obs, y_obs, slideclipdict['window_length'],
                               low=slideclipdict['low'],
                               high=slideclipdict['high'], method='mad',
                               center='median')


        if detrend:
            ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)

            # # default pspline detrending
            if detrend=='pspline':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs)
                x_trend = deepcopy(x_obs)

            # in some cases, might prefer the biweight
            elif detrend == 'biweight':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs,
                                                  method='biweight', cval=5,
                                                  window_length=0.5,
                                                  break_tolerance=0.5)
                x_trend = deepcopy(x_obs)

            elif detrend == 'minimal':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs,
                                                  method='biweight', cval=2,
                                                  window_length=3.5,
                                                  break_tolerance=0.5)
                x_trend = deepcopy(x_obs)

            elif detrend == 'median':
                y_obs, y_trend = dtr.detrend_flux(x_obs, y_obs,
                                                  method='median',
                                                  window_length=0.6,
                                                  break_tolerance=0.5,
                                                  edge_cutoff=0.)
                x_trend = deepcopy(x_obs)

            elif detrend == 'best':
                from cdips.lcproc.find_planets import run_periodograms_and_detrend
                dtr_dict = {'method':'best', 'break_tolerance':0.5, 'window_length':0.5}
                lsp_options = {'period_min':0.1, 'period_max':20}

                # r = [source_id, ls_period, ls_fap, ls_amplitude, tls_period, tls_sde,
                #      tls_t0, tls_depth, tls_duration, tls_distinct_transit_count,
                #      tls_odd_even, dtr_method]

                r, search_time, search_flux, dtr_stages_dict = run_periodograms_and_detrend(
                    ticid, x_obs, y_obs, dtr_dict,
                    period_min=lsp_options['period_min'],
                    period_max=lsp_options['period_max'], dtr_method='best',
                    return_extras=True,
                    magisflux=True
                )
                y_trend, x_trend, dtr_method = (
                    dtr_stages_dict['trend_flux'],
                    dtr_stages_dict['trend_time'],
                    dtr_stages_dict['dtr_method_used']
                )
                x_obs, y_obs = deepcopy(search_time), deepcopy(search_flux)
                print(f'TIC{ticid} TLS results')
                print(f'dtr_method_used: {dtr_method}')
                print(r)

            else:
                raise NotImplementedError

        if detrend:
            ax.plot(x_trend, y_trend, c='r', lw=0.5, zorder=3)
        else:

            exptime_sec = np.nanmedian(np.diff(x_obs))*24*60*60

            if exptime_sec >= 600:
                ax.scatter(x_obs, y_obs, c='k', s=4, zorder=2)

            elif exptime_sec <= 600:
                ax.scatter(
                    x_obs, y_obs, c='darkgray', s=1, zorder=1, alpha=0.4
                )
                from astrobase.lcmath import time_bin_magseries

                if isinstance(bintime, (int,float)):
                    pass
                    s = 12
                else:
                    bintime = 1800.
                    s = 4

                bin_dict = time_bin_magseries(
                    x_obs, y_obs, binsize=bintime, minbinelems=2
                )
                ax.scatter(
                    bin_dict['binnedtimes'], bin_dict['binnedmags'], c='k',
                    s=s, zorder=2
                )

        times.append( x_obs )
        fluxs.append( y_obs )

        ax.set_xlabel('time [bjdtdb]')
        ax.set_ylabel(ykey)
        ylim = ax.get_ylim()

        ax.set_title(ix)

        if detrend:
            _ylim = _get_ylim(y_trend)
        else:
            _ylim = _get_ylim(y_obs)

        ax.set_ylim(_ylim)
        if isinstance(forceylim, list) or isinstance(forceylim, tuple):
            ax.set_ylim(forceylim)

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

    if normstitch:
        times, fluxs, _ = lcu.stitch_light_curves(
            times, fluxs, fluxs, magsarefluxes=True, normstitch=True
        )
    else:
        times = np.hstack(np.array(times).flatten())
        fluxs = np.hstack(np.array(fluxs).flatten())

    # NOTE: this call is deprecated
    stimes, smags, _ = lcmath.sigclip_magseries(
        times, fluxs, np.ones_like(fluxs), sigclip=[20,50], iterative=True,
        magsarefluxes=True
    )

    savpath = join(
        outdir, f'TIC{ticid}_{prov}_{inst}_lightcurve_{str(ykey).zfill(2)}_allsector.png'
    )
    if detrend:
        savpath = join(
            outdir, f'TIC{ticid}_{prov}_{inst}_lightcurve_{detrend}_{str(ykey).zfill(2)}_allsector.png'
        )

    plt.close('all')
    f,ax = plt.subplots(figsize=(16,4))

    exptime_sec = np.nanmedian(np.diff(stimes))*24*60*60

    if exptime_sec >= 600:
        ax.scatter(stimes, smags, c='k', s=1, zorder=2)

    elif exptime_sec <= 600:
        ax.scatter(
            stimes, smags, c='darkgray', s=1, zorder=1, alpha=0.4,
            edgecolors='darkgray', marker='.'
        )
        from astrobase.lcmath import time_bin_magseries
        if isinstance(bintime, (int,float)):
            pass
            s = 5
        else:
            bintime = 3600.
            s = 1

        bin_dict = time_bin_magseries(
            stimes, smags, binsize=bintime, minbinelems=2
        )
        ax.scatter(
            bin_dict['binnedtimes'], bin_dict['binnedmags'], c='k',
            s=s, zorder=2
        )

    if not epoch is None:
        tra_times = epoch + np.arange(-1000,1000,1)*period

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        ax.set_ylim((min(ylim), max(ylim)))
        ax.vlines(tra_times, min(ylim), max(ylim), color='orangered',
                  linestyle='--', zorder=-2, lw=0.5, alpha=0.3)
        ax.set_ylim((min(ylim), max(ylim)))
        ax.set_xlim(xlim)

    ax.set_xlabel('time [bjdtdb]')
    ax.set_ylabel('relative '+ykey)

    ax.set_title(ix)

    f.savefig(savpath, dpi=400, bbox_inches='tight')
    print('made {}'.format(savpath))

    csvpath = savpath.replace('.png','_sigclipped.csv')
    pd.DataFrame({
        'time': stimes, 'flux': smags,
    }).to_csv(csvpath, index=False)
    print(f'made {csvpath}')


    if do_phasefold:

        assert (
            isinstance(period, (float,int)) and isinstance(epoch, (float,int))
        )

        #
        # ax: primary transit
        #
        if inst == 'kepler':
            phasebin = 1e-3
        elif inst == 'tess':
            phasebin = 5e-3
        minbinelems = 2
        plotxlims = [(-0.5, 0.5), (-0.05,0.05)]
        xlimstrs = ['xwide','xnarrow']
        plotylim = [0.9, 1.08]#None #[0.9,1.1]
        do_vlines = False

        for plotxlim, xstr in zip(plotxlims, xlimstrs):

            plt.close('all')
            fig, ax = plt.subplots(figsize=(4,3))

            # use times and fluxs, instead of the sigma clipped thing.
            _make_phased_magseries_plot(ax, 0, times, fluxs,
                                        np.ones_like(fluxs)/1e4, period, epoch,
                                        True, True, phasebin, minbinelems,
                                        plotxlim, '', xliminsetmode=False,
                                        magsarefluxes=True, phasems=0.8,
                                        phasebinms=4.0, verbose=True)
            if isinstance(plotylim, (list, tuple)):
                ax.set_ylim(plotylim)
            else:
                plotylim = _get_ylim(fluxs)
                ax.set_ylim(plotylim)

            if do_vlines:
                ax.vlines(1/6, min(plotylim), max(plotylim), color='orangered',
                          linestyle='--', zorder=-2, lw=1, alpha=0.8)
                ax.set_ylim(plotylim)

            dstr = detrend if detrend else ''
            savpath = join(
                outdir, f'TIC{ticid}_{prov}_{inst}_lightcurve_{dstr}_{ykey}_{xstr}_allsector_phasefold.png'
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


def make_periodogram(data, ticid, pipeline, id_str=None, outdir=None, period_min=0.1,
                     period_max=20, manual_peak=None, samples_per_peak=50,
                     nterms_0=6):

    if pipeline in ['spoc', 'kepler']:
        time = data[LCKEYDICT[pipeline]['time']]
        flux = data[LCKEYDICT[pipeline]['flux']]
        err = flux*1e-4
        qual = data[LCKEYDICT[pipeline]['quality']]
        sel = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err) # (qual == 0)
        time, flux, err = time[sel], flux[sel], err[sel]
        med_flux = np.nanmedian(flux)
        flux /= med_flux
        err /= med_flux
    elif pipeline == 'cdips':
        time = data['time']
        flux = data['flux']
        err = data['err']
    elif pipeline == 'unpopular':
        time = data['time']
        flux = 1+data['dtr_flux']
        err = 1e-4*flux
    else:
        raise NotImplementedError

    ##########################################

    from astropy.timeseries import LombScargle
    from scipy.signal import find_peaks

    plt.close('all')

    fig, ax = plt.subplots()

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

    tstr = f'LS period = {ls_period_0:.6f} d'
    ax.set_title(tstr)
    ax.set_yscale('linear')
    ax.set_xscale('log')

    ax.set_xlabel('period [d]')
    ax.set_ylabel('power')

    if outdir is None:
        outdir = '../results/quicklooklc/TIC{}'.format(ticid)
    s = '' if id_str is None else f'_{id_str}'
    savpath = join(outdir, f'ls_periodogram{s}.png')
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

    s = '' if id_str is None else f'_{id_str}'
    savpath = join(outdir, f'ls_models_resids{s}.png')
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

        savpath = join(outdir, 'tls_periodogram.png')
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

        savpath = join(outdir, 'phasefold.png')
        fig.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))


def make_riverplot(data):

    msg = (
        'Implementation resides in complexrotators project'
    )
    raise NotImplementedError(msg)
