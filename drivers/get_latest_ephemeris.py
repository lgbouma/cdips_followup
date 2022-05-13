"""
get_latest_ephemeris.py: given a source identifier, pull a) CDIPS light curves
and b) unpopular light curves for the object, and fit a default transit model.

TESS observes continually.  If you are slow enough (and I am), this means
you may sometimes want to get the latest ephemeris information for a given
star, without needing to say, wait for the next batch of CDIPS light curves.

(Run 2022/05/12, when S49 was max uplaod at
https://outerspace.stsci.edu/display/TESS/TESS+Holdings+Available+by+MAST+Service#TESSHoldingsAvailablebyMASTService-NorthernHemisphere:)

Contents:
    load_cdips_lc
    load_cpm_lc
    _plot_comparison
"""
#############
## LOGGING ##
#############

import logging
from astrobase import log_sub, log_fmt, log_date_fmt

DEBUG = False
if DEBUG:
    level = logging.DEBUG
else:
    level = logging.INFO
LOGGER = logging.getLogger(__name__)
logging.basicConfig(
    level=level,
    style=log_sub,
    format=log_fmt,
    datefmt=log_date_fmt,
)

LOGDEBUG = LOGGER.debug
LOGINFO = LOGGER.info
LOGWARNING = LOGGER.warning
LOGERROR = LOGGER.error
LOGEXCEPTION = LOGGER.exception

#############
## IMPORTS ##
#############

import os
from os.path import join
from glob import glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from astropy.io import fits

from astrobase.services.tesslightcurves import (
    get_hlsp_lightcurves, get_unpopular_lightcurve
)
from astrobase.services.identifiers import gaiadr2_to_tic

from cdips_followup.paths import PHOTDIR, RESULTSDIR, LOCALDIR

from cdips.utils import lcutils as lcu
from cdips.lcproc.find_planets import run_periodograms_and_detrend

def main():

    #
    # Objects for which to retrieve the latest and greatest light curves.
    #
    source_ids = [
        # PP0
        "2136224721950935168", # two new EM sectors!  and they look good.
        #"2072431332214824576", # got EM data: TOI-3561... not obviously good in either pipeline's LC though. wtf.
        #"4519168132218628992", # no EM (yet, S53,54)
        #"1935439776865267072", # no EM
        #"253885554213280896", # no EM
        #"521454526326028928", # no EM (yet, S52)
        #"542223510700499968", # no EM (yet, S52)
        #"1968083349388709120", # no EM (yet, S55)
        #"1970742930580834176", # no EM (yet, S55)
        #"1984154395449640832", # no EM
        #"2204810126994432384", # no EM
        #"4599686193337056256", # no EM (yet, S52)
        # PP1
        # "526763312065766016", # no EM (yet, S52)
        # "568619413331898240",  # no EM (yet, S52+53)
        # "1203566251429799680", # no EM (yet, S51)
        # "2056621793787000704", # no EM (yet, S55)
        # "2190239295819005440" # no EM (yet, S55)
    ]

    for source_id in source_ids:
        get_latest_ephemeris(source_id)


def load_cdips_lc(lcpath):
    # get time & flux. (follows procedure from do_initial_period_finding).
    modelid = 'simpletransit'
    starid = (
        os.path.basename(lcpath).
        replace('.fits','').
        replace('hlsp_cdips_tess_ffi_','')
    )

    fit_savdir = os.path.dirname(lcpath)
    lc_csvpath0 = join(fit_savdir, f'{starid}_{modelid}_rawlc.csv')
    lc_csvpath1 = join(fit_savdir, f'{starid}_{modelid}_detrendedlc.csv')

    if os.path.exists(lc_csvpath0) and os.path.exists(lc_csvpath1):
        return starid, pd.read_csv(lc_csvpath0), pd.read_csv(lc_csvpath1)

    hdul = fits.open(lcpath)
    hdr = hdul[0].header

    APNAME = 'PCA1'
    source_id, time, mag, xcc, ycc, ra, dec, _, _ = (
        lcu.get_lc_data(lcpath, mag_aperture=APNAME)
    )
    err_mag = hdul[1].data['IRE1']

    dtr_method, break_tolerance, window_length = 'best', 0.5, 0.5
    dtr_dict = {'method':dtr_method,
                'break_tolerance':break_tolerance,
                'window_length':window_length}

    r, dtr_time, dtr_flux, dtr_stages_dict = run_periodograms_and_detrend(
        source_id, time, mag, dtr_dict, return_extras=True
    )

    # save initial (orbit-edge-masked) time and flux, w/out detrending.
    outdf = pd.DataFrame({
        'time': dtr_stages_dict['time'],
        'flux_'+APNAME: dtr_stages_dict['flux'],
    })
    outdf.to_csv(lc_csvpath0, index=False)
    LOGINFO(f'Wrote {lc_csvpath0}')

    dtrdf = pd.DataFrame({
        'search_time': dtr_stages_dict['search_time'],
        'search_flux_'+APNAME: dtr_stages_dict['search_flux'],
    })
    dtrdf.to_csv(lc_csvpath1, index=False)
    LOGINFO(f'Wrote {lc_csvpath1}')

    return starid, outdf, dtrdf


def load_cpm_lc(source_id, lcpath):

    modelid = 'simpletransit'
    starid = (
        os.path.basename(lcpath).
        replace('.csv','')
    )

    fit_savdir = os.path.dirname(lcpath)
    lc_csvpath0 = join(fit_savdir, f'{starid}_{modelid}_rawlc.csv')
    lc_csvpath1 = join(fit_savdir, f'{starid}_{modelid}_detrendedlc.csv')

    if os.path.exists(lc_csvpath0) and os.path.exists(lc_csvpath1):
        return starid, pd.read_csv(lc_csvpath0), pd.read_csv(lc_csvpath1)

    df = pd.read_csv(lcpath)
    time = np.array(df.time)
    flux = np.array(df.dtr_flux)
    flux += 1
    flux = flux / np.nanmedian(flux)

    dtr_method, break_tolerance, window_length = 'best', 0.5, 0.5
    dtr_dict = {'method':dtr_method,
                'break_tolerance':break_tolerance,
                'window_length':window_length}

    r, dtr_time, dtr_flux, dtr_stages_dict = run_periodograms_and_detrend(
        source_id, time, flux, dtr_dict, return_extras=True, magisflux=True
    )

    # save initial (orbit-edge-masked) time and flux, w/out detrending.
    outdf = pd.DataFrame({
        'time': dtr_stages_dict['time'],
        'flux_cpm': dtr_stages_dict['flux'],
    })
    outdf.to_csv(lc_csvpath0, index=False)
    LOGINFO(f'Wrote {lc_csvpath0}')

    # with detrending
    dtrdf = pd.DataFrame({
        'search_time': dtr_stages_dict['search_time'],
        'search_flux_cpm': dtr_stages_dict['search_flux'],
    })
    dtrdf.to_csv(lc_csvpath1, index=False)
    LOGINFO(f'Wrote {lc_csvpath1}')

    return starid, outdf, dtrdf


def _plot_comparison(source_id, cdips_d, cpm_d, figpath, detrended=True):

    cpm_sectors = sorted([k.split("_")[1] for k in list(cpm_d.keys()) if not
                          k.endswith("_dtr")])
    cdips_sectors = sorted([k.split("-")[1] for k in list(cdips_d.keys()) if
                            not k.endswith("_dtr")])
    n_sectors = len(cpm_sectors)

    plt.close("all")
    fig, axs = plt.subplots(
        nrows=n_sectors, figsize=(12,4*n_sectors)
    )

    for ix in range(n_sectors):

        sector = cpm_sectors[ix]
        if detrended:
            cpm_key = [k for k in cpm_d.keys() if sector in k and
                       k.endswith("_dtr")]
            cdips_key = [k for k in cdips_d.keys() if sector in k and
                         k.endswith("_dtr")]
        else:
            cpm_key = [k for k in cpm_d.keys() if sector in k and not
                       k.endswith("_dtr")]
            cdips_key = [k for k in cdips_d.keys() if sector in k and not
                         k.endswith("_dtr")]
        assert len(cpm_key) == 1

        df_cpm = cpm_d[cpm_key[0]]

        if len(cdips_key) == 1:
            df_cdips = cdips_d[cdips_key[0]]
            have_cdips = True
            xkey, ykey = 'time', 'flux_PCA1'
            if detrended:
                xkey, ykey = 'search_time', 'search_flux_PCA1'
            ydiff = 1.4*(
                np.nanpercentile(df_cdips[ykey], 95)
                -
                np.nanpercentile(df_cdips[ykey], 5)
            )
        else:
            ydiff = 0
            have_cdips = False

        xkey, ykey = 'time', 'flux_cpm'
        if detrended:
            xkey, ykey = 'search_time', 'search_flux_cpm'
        axs[ix].scatter(df_cpm[xkey], df_cpm[ykey], c="C0", s=3, label="CPM",
                        zorder=2, rasterized=True, linewidths=0)

        if have_cdips:
            x0 = 2457000
            xkey, ykey = 'time', 'flux_PCA1'
            if detrended:
                xkey, ykey = 'search_time', 'search_flux_PCA1'
            axs[ix].scatter(df_cdips[xkey]-x0, df_cdips[ykey], c="k", s=3,
                            label="CDIPS(PCA1)", zorder=2, rasterized=True,
                            linewidths=0)

        txt = sector
        axs[ix].text(0.97,0.97,txt, transform=axs[ix].transAxes,
                     ha='right',va='top', color='k')


    axs[0].legend()

    fig.text(-0.01, 0.5, 'Relative flux', va='center', rotation=90)
    fig.text(0.5, -0.01, "Time - 2457000 [Days]", ha='center', va='center')
    fig.tight_layout()

    fig.savefig(figpath, bbox_inches="tight", dpi=300)


def get_latest_ephemeris(source_id):
    """
    Given a source identifier, pull a) CDIPS light curves and b) unpopular
    light curves for the object, and fit a default transit model using CDIPS
    when available, else unpopular.
    """

    tic_id = gaiadr2_to_tic(str(source_id))

    photdir = join(LOCALDIR, f'TIC{tic_id}')
    ephemdir = join(RESULTSDIR, f'latest_ephemerides')
    for d in [ephemdir, photdir]:
        if not os.path.exists(d):
            os.mkdir(d)

    lcfiles0 = get_hlsp_lightcurves(tic_id, hlsp_products=['CDIPS'],
                                    download_dir=photdir, verbose=True)

    lcfiles1 = get_unpopular_lightcurve(tic_id, download_dir=photdir)

    ephemcsv = join(ephemdir, f'{source_id}_ephem.csv')

    if len(lcfiles0) == len(lcfiles1):
        LOGINFO(70*'/')
        txt = f"No EM data for GDR2 {source_id} / TIC {tic_id}..."
        LOGINFO(f"{txt:/^70}")
        LOGINFO(70*'/')
        colorder = ['source_id', 'ticid', 'period',
                    'period_unc', 'epoch', 'epoch_unc', 'depth', 'depth_unc',
                    'duration', 'duration_unc', 'ephemeris_origin']
        colvals = [source_id, tic_id,
                   np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
                   'get_latest_ephemeris_no_extended_mission_data']
        df = pd.DataFrame({k:v for k,v in zip(colorder, colvals)}, index=[0])
        df.to_csv(ephemcsv, index=False)
        return -1

    # We have new extended mission data.  Load in the light curves.
    cdips_d = {}
    for lcpath in lcfiles0:
        starid, df, dtr_df = load_cdips_lc(lcpath)
        cdips_d[starid] = df
        cdips_d[starid+"_dtr"] = dtr_df

    cpm_d = {}
    for lcpath in lcfiles1:
        starid, df, dtr_df = load_cpm_lc(source_id, lcpath)
        cpm_d[starid] = df
        cpm_d[starid+"_dtr"] = dtr_df

    figpath = join(ephemdir, f'gdrtwo_{source_id}_cdips_cpm_comparison.png')
    if not os.path.exists(figpath):
        _plot_comparison(source_id, cdips_d, cpm_d, figpath, detrended=False)
    else:
        LOGINFO(f"Found {figpath}")

    figpath = join(ephemdir, f'gdrtwo_{source_id}_cdips_cpm_detrended_comparison.png')
    if not os.path.exists(figpath):
        _plot_comparison(source_id, cdips_d, cpm_d, figpath, detrended=True)
    else:
        LOGINFO(f"Found {figpath}")

    #
    # TODO here would be to fit the data, if you thought it was worthwhile.
    #

    import IPython; IPython.embed()
    # Fit everything. TODO: copy-paste from L953 of fit_models_to_gold
    #FIXME FIXME FIXME



if __name__ == "__main__":
    main()
