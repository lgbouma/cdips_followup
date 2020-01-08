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

def get_data(cdips=1, spoc=0):

    if cdips:
        # data: list of fits records
        lcfiles = glob('../results/PTFO_8-8695/*fits')
    if spoc:
        lcfiles = glob('../results/PTFO_8-8695/MAST_2020-01-06T1938/TESS/tess2018349182459-s0006-0000000264461976-0126-s/tess2018349182459-s0006-0000000264461976-0126-s_lc.fits')

    data = []
    for f in lcfiles:
        hdul = fits.open(f)
        data.append(hdul[1].data)

    return data


def explore_mag_lightcurves(data):

    for yval in ['TFA1','TFA2','TFA3','IRM1','IRM2','IRM3']:

        for ix, d in enumerate(data):

            title = 'sector_7' if ix == 0 else 'sector_9'
            savpath = (
                '../results/PTFO_8-8695/mag_lightcurve_{}_{}.png'.
                format(yval, title)
            )
            if os.path.exists(savpath):
                print('found {}, continue'.format(savpath))
                continue

            plt.close('all')
            f,ax = plt.subplots(figsize=(16,4))

            ax.scatter(d['TMID_BJD'], d[yval], c='k', s=5)

            ax.set_xlabel('time [bjdtdb]')
            ax.set_ylabel(yval)
            ylim = ax.get_ylim()
            ax.set_ylim((max(ylim), min(ylim)))

            ax.set_title(title)

            f.savefig(savpath, dpi=300, bbox_inches='tight')
            print('made {}'.format(savpath))


def explore_flux_lightcurves(data, iscdips=1):

    if iscdips:
        for yval in ['TFA1','TFA2','TFA3','IRM1','IRM2','IRM3']:

            for ix, d in enumerate(data):

                title = 'sector_6' if ix == 0 else None
                savpath = (
                    '../results/PTFO_8-8695/flux_lightcurve_{}_{}.png'.
                    format(yval, title)
                )
                if os.path.exists(savpath):
                    print('found {}, continue'.format(savpath))
                    continue

                plt.close('all')
                f,ax = plt.subplots(figsize=(16,4))

                flux = vp._given_mag_get_flux(d[yval])

                ax.scatter(d['TMID_BJD'], flux, c='k', s=5)

                ax.set_xlabel('time [bjdtdb]')
                ax.set_ylabel(yval)

                ax.set_title(title)

                f.savefig(savpath, dpi=300, bbox_inches='tight')
                print('made {}'.format(savpath))

    else:

        for ix, d in enumerate(data):

            title = 'sector_6_spoc' if ix == 0 else None
            savpath = (
                '../results/PTFO_8-8695/flux_lightcurve_{}.png'.
                format(title)
            )

            plt.close('all')
            f,ax = plt.subplots(figsize=(32,4))

            flux = d['PDCSAP_FLUX']/np.nanmedian(d['PDCSAP_FLUX'])

            ax.scatter(d['TIME'], flux, c='k', s=5)

            ax.set_xlabel('time [bjdtdb]')
            ax.set_ylabel('rel flux PDCSAP')

            ax.set_title(title)

            f.savefig(savpath, dpi=300, bbox_inches='tight')
            print('made {}'.format(savpath))


def get_light_detrended_data(data):
    """
    clip +/12 hours on every orbit
    fit out a quadratic to each orbit

    return a dict with keys:
        orbit1, orbit2, orbit3, orbit4, sector7, sector9, allsectors

        and each key leads to another dictionary with time and all available
        apertures+detrended stages.
    """

    data_dict = {}
    dtrtypes = ['IRM','PCA','TFA']
    apnums = [1,2,3]
    sectornums = [7,9]

    orbitgap = 1
    orbitpadding = 0.5

    for sectornum, sector_data in zip(sectornums, data):

        time = sector_data['TMID_BJD']

        flux_sector_dict = {}
        for dtrtype in dtrtypes:
            for apnum in apnums:
                k = dtrtype+str(apnum)
                this_flux = (
                    vp._given_mag_get_flux(sector_data[k])
                )

                # now mask orbit start and end
                time_copy = deepcopy(time)

                trim_time, trim_this_flux = moe.mask_orbit_start_and_end(
                    time_copy, this_flux, orbitgap=orbitgap,
                    expected_norbits=2, orbitpadding=orbitpadding
                )

                # now fit out quadratics from each orbit and rejoin them
                norbits, trim_time_groups = lcmath.find_lc_timegroups(
                    trim_time, mingap=orbitgap)
                if norbits != 2:
                    raise AssertionError('expected 2 orbits')

                save_flux = []
                for time_group in trim_time_groups:
                    trim_time_orbit = trim_time[time_group]
                    trim_this_flux_orbit = trim_this_flux[time_group]

                    order = 2 # fit out a quadtric! #FIXME CHECK
                    p = Legendre.fit(trim_time_orbit, trim_this_flux_orbit, order)
                    trim_this_flux_orbit_fit = p(trim_time_orbit)

                    save_flux_orbit = (
                        trim_this_flux_orbit / trim_this_flux_orbit_fit
                    )

                    save_flux.append(save_flux_orbit)

                flux_sector_dict[k] = np.concatenate(save_flux)

        # update time to be same length as all the trimmed fluxes
        time = trim_time

        sectorkey = 'sector{}'.format(sectornum)
        data_dict[sectorkey] = {}
        data_dict[sectorkey]['time'] = time
        data_dict[sectorkey]['fluxes'] = flux_sector_dict

    sectorkeys = ['sector{}'.format(s) for s in sectornums]

    # create the "allsectors" dataset
    data_dict['allsectors'] = {}
    data_dict['allsectors']['time'] = np.concatenate([
        data_dict[sectorkey]['time'] for sectorkey in sectorkeys
    ])

    data_dict['allsectors']['fluxes'] = {}
    for dtrtype in dtrtypes:
        for apnum in apnums:
            k = dtrtype+str(apnum)
            data_dict['allsectors']['fluxes'][k] = np.concatenate([
                data_dict[sectorkey]['fluxes'][k] for sectorkey in sectorkeys
            ])

    return data_dict


def do_phasefolds(data):

    data_dict = get_light_detrended_data(data)

    dtrtypes = ['IRM','PCA','TFA']
    apnums = [1,2,3]

    startp, endp = 1.96, 2.1
    maxtransitduration = 0.4
    nworkers = multiprocessing.cpu_count() - 2

    for plotxlim in [[-0.7, 0.7], [-0.25, 0.25], [-0.5, -0.2]]:

        for chunktype in list(data_dict.keys()):

            xmin = min(plotxlim)
            xmax = max(plotxlim)

            savpath = (
                '../results/PTFO_8-8695/{}_phasefold_xmin{}_xmax{}.png'.
                format(chunktype, xmin, xmax)
            )
            if os.path.exists(savpath):
                print('found {}, continue'.format(savpath))
                continue

            plt.close('all')
            fig,axs = plt.subplots(3,3,figsize=(12,12))

            for i, dtrtype in enumerate(dtrtypes):
                for j, apnum in enumerate(apnums):
                    flux_k = dtrtype+str(apnum)

                    time = data_dict[chunktype]['time']
                    flux = data_dict[chunktype]['fluxes'][flux_k]
                    err_flux = 1e-3 * flux

                    # get epoch, period
                    blsdict = kbls.bls_parallel_pfind(
                        time, flux, err_flux, magsarefluxes=True,
                        startp=startp, endp=endp,
                        maxtransitduration=maxtransitduration,
                        nworkers=nworkers, sigclip=None)
                    fitd = kbls.bls_stats_singleperiod(
                        time, flux, err_flux, blsdict['bestperiod'],
                        maxtransitduration=maxtransitduration,
                        magsarefluxes=True, sigclip=None, perioddeltapercent=5)

                    # nb. can access: fitd['phases'], fitd['phasedmags'],
                    # fitd['blsmodel'], fitd['period'], fitd['epoch'],

                    # get phased points and fluxes.
                    periodind = 0
                    phasebin = 0.01
                    minbinelems = 7

                    _make_phased_magseries_plot(
                        axs[i,j], periodind,
                        time, flux, err_flux,
                        2.029202, 2458562.47131,
                        #blsdict['bestperiod'], 2458562.47131,
                        #fitd['period'], fitd['epoch'],
                        True, True, phasebin, minbinelems, plotxlim, "bls",
                        lspmethodind=0, xliminsetmode=False, twolspmode=False,
                        magsarefluxes=True, verbose=True, phasems=2.0,
                        phasebinms=4.0, xticksize=None, yticksize=None,
                        titlefontsize=6,
                        makegrid=True, lowerleftstr=flux_k,
                        lowerleftfontsize=10
                    )

            fig.text(0.55,1, chunktype, ha='center')
            #fig.text(0.55,0, 'Phase', ha='center')

            for ax in axs.flatten():
                ax.get_yaxis().set_tick_params(which='both', direction='in')
                ax.get_xaxis().set_tick_params(which='both', direction='in')

            fig.tight_layout(h_pad=0.2, w_pad=0.2)

            fig.savefig(savpath, dpi=300, bbox_inches='tight')
            print('made {}'.format(savpath))


def make_riverplot(data):

    data_dict = get_light_detrended_data(data)

    period = 2.029202
    sampling_rate = 30/(60*24)

    #epoch = 2458562.47131
    time = data_dict['sector7']['time']
    offset = -period/2
    epoch = np.min(time) + offset

    # t_offset = 2458543

    f, axs = plt.subplots(nrows=2, ncols=1, figsize=(7,8))

    for _ix, ax in enumerate(axs):

        k = 'sector7' if _ix == 0 else 'sector9'

        time = data_dict[k]['time']
        flux = data_dict[k]['fluxes']['IRM1']

        transit_number = np.floor((time - epoch)/period)

        N_transits = int(np.max(transit_number) - np.min(transit_number))
        N_points_per_period = int(np.floor(period/sampling_rate))

        flux_arr = np.ones((N_transits, N_points_per_period))

        t_starts = [
            ix*period + epoch for ix in
            range(int(np.min(transit_number)), int(np.max(transit_number)))
        ]

        transit_index = [
            ix for ix in
            range(int(np.min(transit_number)), int(np.max(transit_number)))
        ]

        for ix, t_start in enumerate(t_starts):

            ind = (time > t_start) & (time < t_start + period)

            flux_to_insert = flux[ind]

            flux_arr[ix, :len(flux_to_insert)] = flux_to_insert

        im = ax.pcolor(np.array(range(N_points_per_period))/N_points_per_period,
                       transit_index,
                       flux_arr,
                       cmap='Blues_r')
                       #cmap='Greys_r') #looks alright too

    axs[0].set_xticks([])
    axs[1].set_xlabel('phase')
    axs[0].set_ylabel('transit number')
    axs[1].set_ylabel('transit number')

    f.tight_layout()

    outpath = '../results/PTFO_8-8695/riverplot.png'
    f.savefig(outpath, bbox_inches='tight', dpi=400)
    print('saved {}'.format(outpath))


def detrend_lightcurve_wotan(data, window_length=0.25, iscdips=True):

    #NOTE: junky. doesnt work. 
    # try two simple harmonic oscillators with periods separated by a factor of
    # two...

    if iscdips:
        sector_data = data[0]
        time = sector_data['TMID_BJD']
        flux = vp._given_mag_get_flux(sector_data['IRM1'])
    else:
        d = data[0]
        time = d['TIME']
        flux = d['PDCSAP_FLUX']/np.nanmedian(d['PDCSAP_FLUX'])

    from wotan import flatten

    _, trend_lc2 = flatten(
        time,                 # Array of time values
        flux,                 # Array of flux values
        method='biweight',
        robust=True,          # Iteratively clip 2-sigma outliers until convergence
        window_length=window_length,    # The length of the filter window in units of ``time``
        break_tolerance=0.5,  # Split into segments at breaks longer than that
        return_trend=True,    # Return trend and flattened light curve
    )
    if window_length == 99:
        trend_lc2 = np.ones_like(flux)*np.nanmedian(flux)

    # flatten_lc2, trend_lc2 = flatten(
    #         time,                  # Array of time values
    #         flux,                  # Array of flux values
    #         method='gp',
    #         kernel='periodic',     # GP kernel choice
    #         kernel_period=0.498818,  # GP kernel period
    #         kernel_size=50,        # GP kernel length
    #         break_tolerance=0.5,   # Split into segments at breaks longer than
    #         return_trend=True,     # Return trend and flattened light curve
    # )


    # plot!
    plt.close('all')
    if iscdips:
        f,axs = plt.subplots(nrows=2, ncols=1, figsize=(16,7))
    else:
        f,axs = plt.subplots(nrows=2, ncols=1, figsize=(32,7))

    axs[0].scatter(time, flux, c='k', s=5, zorder=2)
    axs[0].plot(time, trend_lc2, c='orange', lw=1, zorder=3)

    dtr_flux = flux-trend_lc2
    sel = ~np.isnan(dtr_flux)
    axs[1].scatter(time[sel], dtr_flux[sel], c='k', s=5)

    axs[1].set_xlabel('time [bjdtdb]')
    axs[0].set_ylabel('IRM1 flux')
    axs[1].set_ylabel('residual')

    axs[0].set_title('biweight detrend, {}d'.format(window_length))

    if not iscdips:
        axs[1].set_ylim([-0.1,0.05])
        if window_length == 99:
            axs[1].set_ylim([-0.1,0.1])

    isspoc = '' if iscdips else '_spoc2min'
    savpath = '../results/PTFO_8-8695/detrend_lc{}_{:.2f}d.png'.format(
        isspoc, window_length
    )
    f.savefig(savpath, dpi=300, bbox_inches='tight')
    print('made {}'.format(savpath))

    if not iscdips:

        from transitleastsquares import transitleastsquares
        model = transitleastsquares(
            time[sel],
            1 + dtr_flux[sel] + np.nanmedian(dtr_flux[sel])
        )
        results = model.power(period_min=0.3, period_max=0.6, M_star_min=0.1,
                              M_star_max=5, R_star=0.5, M_star=0.5)

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

        savpath = '../results/PTFO_8-8695/detrend_lc{}_{:.2f}d_periodogram.png'.format(
            isspoc, window_length
        )
        fig.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))

        plt.close('all')
        fig = plt.figure()
        plt.plot(results.model_folded_phase, results.model_folded_model,
                 color='orange', zorder=3)
        plt.scatter(results.folded_phase, results.folded_y, color='black', s=2,
                    alpha=0.9, zorder=4, linewidths=0)
        #plt.xlim(0.35, 0.65)
        plt.title('{}: {}d'.format(results.T0, results.period))
        plt.xlabel('Phase')
        plt.ylabel('Relative flux');

        savpath = '../results/PTFO_8-8695/detrend_lc{}_{:.2f}d_phasefold.png'.format(
            isspoc, window_length
        )
        fig.savefig(savpath, dpi=300, bbox_inches='tight')
        print('made {}'.format(savpath))




def main():

    cdips = 0
    spoc = 1

    do_mag_lcs = 0
    do_flux_lcs = 0
    do_detrending = 1
    do_pf = 0
    do_riverplot = 0

    data = get_data(cdips=cdips, spoc=spoc)

    if do_mag_lcs:
        explore_mag_lightcurves(data)
    if do_flux_lcs:
        explore_flux_lightcurves(data, iscdips=cdips)
    if do_detrending:
        #for wl in [99]:
        for wl in [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]:
            detrend_lightcurve_wotan(data, window_length=wl, iscdips=cdips)
    if do_pf:
        do_phasefolds(data)
    if do_riverplot:
        make_riverplot(data)

if __name__ == "__main__":
    main()
