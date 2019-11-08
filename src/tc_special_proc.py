from glob import glob
import os, multiprocessing
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from astrobase import lcmath
from astrobase.periodbase import kbls
from astrobase.checkplot.png import _make_phased_magseries_plot
from copy import deepcopy

from cdips.plotting import vetting_pdf as vp
from cdips.lcproc import mask_orbit_edges as moe

from numpy.polynomial.legendre import Legendre

def get_data():
    # data: list of fits records
    lcfiles = glob('../data/20191102_trojan_cand/*fits')

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
                '../results/tc_special_proc/mag_lightcurve_{}_{}.png'.
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


def explore_flux_lightcurves(data):

    for yval in ['TFA1','TFA2','TFA3','IRM1','IRM2','IRM3']:

        for ix, d in enumerate(data):

            title = 'sector_7' if ix == 0 else 'sector_9'
            savpath = (
                '../results/tc_special_proc/flux_lightcurve_{}_{}.png'.
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
                '../results/tc_special_proc/{}_phasefold_xmin{}_xmax{}.png'.
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


def main():

    data = get_data()
    explore_mag_lightcurves(data)
    explore_flux_lightcurves(data)
    do_phasefolds(data)

if __name__ == "__main__":
    main()
