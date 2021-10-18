"""
Want RV vs flux
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os
from glob import glob
from cdips_followup.paths import DATADIR, RESULTSDIR
from astropy.io import fits
from aesthetic.plot import savefig, format_ax, set_style
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import curve_fit

def get_rvs():
    rvpath = os.path.join(
        DATADIR, 'HIP_67522',
        'HIP67522_minerva_chiron_merged_radvel_20211014.csv'
    )
    return pd.read_csv(rvpath)

def get_minerva_mean_rvs(df):

    df = df[df.tel.str.contains('MA_fiber')]
    utimes = np.unique(df.time)

    print('Beginning MINERVA averaging...')

    wmnvels, werrs = [], []
    notallfour = 0
    for t in utimes:

        try:
            assert len(df[df.time == t]) == 4
        except:
            print(f"WRN! {t} has {len(df[df.time == t])} time(s), not 4")
            notallfour += 1

        sdf = df[df.time == t]

        vel = np.array(sdf.mnvel)
        errvel = np.array(sdf.errvel)

        wmnvels.append(
            np.sum(vel/(errvel**2))/(np.sum(1/errvel**2))
        )
        werrs.append(
            np.sqrt(
                1/(np.sum(1/errvel**2))
            )
        )

    print(f"{notallfour-1} of {len(utimes)} MA points had 3 tels running")

    mdf = pd.DataFrame({
        'time': utimes,
        'mnvel': wmnvels,
        'errvel': werrs,
        'tel': 'MA_binned'
    })

    return mdf


def get_phot():
    lcpath = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/HIP_67522/'
        'TESS/MAST_2021-10-14T0045/TESS/tess2021118034608-s0038-0000000166527623-0209-s/'
        'tess2021118034608-s0038-0000000166527623-0209-s_lc.fits'
    )
    hl = fits.open(lcpath)
    return (
        hl[1].data['TIME'],
        hl[1].data['PDCSAP_FLUX'],
        hl[1].data['PDCSAP_FLUX_ERR'],
        hl[1].data['QUALITY']
    )


def plot_rv_flux_vs_time():

    rvdf = get_rvs()
    mrvdf = get_minerva_mean_rvs(rvdf)
    mrvdf = pd.concat([mrvdf, rvdf[rvdf.tel.str.contains("CHIRON")]])
    outdir = os.path.join(RESULTSDIR, 'HIP_67522')
    outpath = os.path.join(
        outdir, 'HIP67522_binnedminerva_chiron_merged_radvel_20211014.csv'
    )
    outmrvdf = mrvdf
    outmrvdf['mnvel'] *= 1e3
    sel0 = mrvdf.tel.str.contains("CHIRON")
    mrvdf.loc[sel0, "mnvel"] = mrvdf[sel0].mnvel - np.nanmean(mrvdf[sel0].mnvel)
    sel0 = mrvdf.tel.str.contains("MA_binned")
    mrvdf.loc[sel0, "mnvel"] = mrvdf[sel0].mnvel - np.nanmean(mrvdf[sel0].mnvel)
    outmrvdf['errvel'] *= 1e3
    outmrvdf.to_csv(outpath, index=False)
    print(f"made {outpath}")

    time,flux,flux_err,qual = get_phot()
    sel = (qual == 0) & np.isfinite(flux)
    stime,sflux,sflux_err = time[sel], flux[sel], flux_err[sel]
    mu = np.median(sflux)
    sflux = sflux / mu - 1
    sflux_err /= mu
    stime += 2457000 #BTJD to BJD


    #
    # make plot
    #
    plt.close('all')
    set_style()

    fig, axs = plt.subplots(nrows=3, figsize=(4,6), sharex=True)

    #
    # RVs
    #
    ax = axs[0]

    tels = np.unique(rvdf.tel)
    for ix, tel in enumerate(tels):
        sel = (rvdf.tel == tel)
        x = rvdf[sel].time
        y = rvdf[sel].mnvel
        y_err = rvdf[sel].errvel
        ax.errorbar(
            x, y, y_err, marker='o', elinewidth=0.5, capsize=4, lw=0, mew=0.5,
            color=f'C{ix}', markersize=3, zorder=5, label=tel
        )
    ax.legend(loc='best', fontsize='xx-small')
    ax.set_ylabel('RV [km/s]')

    ax = axs[1]

    tels = np.unique(mrvdf.tel)
    for ix, tel in enumerate(tels):
        sel = (mrvdf.tel == tel)
        x = mrvdf[sel].time
        y = mrvdf[sel].mnvel - np.nanmean(mrvdf[sel].mnvel)
        y_err = mrvdf[sel].errvel
        ax.errorbar(
            x, y, y_err, marker='o', elinewidth=0.5, capsize=4, lw=0, mew=0.5,
            color=f'C{ix}', markersize=3, zorder=-ix+5, label=tel
        )
    ax.legend(loc='best', fontsize='xx-small')
    ax.set_ylabel('RV-$\mu$ [km/s]')


    #
    # phot
    #
    ax = axs[2]

    ax.scatter(
        stime, sflux, c='k', s=0.5, rasterized=True, linewidths=0,
        zorder=1
    )

    # # NOTE not great
    # fn = interp1d(
    #     stime, sflux, kind='quadratic', bounds_error=True, fill_value=np.nan
    # )
    fn = UnivariateSpline(stime, sflux)
    fn.set_smoothing_factor(1e-2) # tune this number

    itime = np.linspace(stime.min(), stime.max(), 2000)
    ax.plot(
        itime, fn(itime), c='C0', zorder=5, lw=0.5, label='spline'
    )
    ax.legend(loc='best', fontsize='xx-small')


    ax.set_ylabel('Flux')
    ax.set_xlabel('BJDTDB')

    ax.set_xlim([2459330, 2459372])

    # set naming options
    outdir = os.path.join(RESULTSDIR, 'HIP_67522')
    outpath = os.path.join(outdir, 'rv_flux_vs_time.png')
    savefig(fig, outpath, dpi=400)


def plot_all_rvs():

    rvdf = get_rvs()
    mrvdf = get_minerva_mean_rvs(rvdf)
    mrvdf = pd.concat([mrvdf, rvdf[rvdf.tel.str.contains("CHIRON")]])

    time,flux,flux_err,qual = get_phot()
    sel = (qual == 0) & np.isfinite(flux)
    stime,sflux,sflux_err = time[sel], flux[sel], flux_err[sel]
    mu = np.median(sflux)
    sflux = sflux / mu - 1
    sflux_err /= mu
    stime += 2457000 #BTJD to BJD


    #
    # make plot
    #
    plt.close('all')
    set_style()

    fig, axs = plt.subplots(nrows=2, figsize=(4,6), sharex=True)

    #
    # RVs
    #
    ax = axs[0]

    tels = np.unique(rvdf.tel)
    for ix, tel in enumerate(tels):
        sel = (rvdf.tel == tel)
        x = rvdf[sel].time
        y = rvdf[sel].mnvel
        y_err = rvdf[sel].errvel
        ax.errorbar(
            x, y, y_err, marker='o', elinewidth=0.5, capsize=4, lw=0, mew=0.5,
            color=f'C{ix}', markersize=3, zorder=5, label=tel
        )
    ax.legend(loc='best', fontsize='xx-small')
    ax.set_ylabel('RV [km/s]')

    ax = axs[1]

    tels = np.unique(mrvdf.tel)
    for ix, tel in enumerate(tels):
        sel = (mrvdf.tel == tel)
        x = mrvdf[sel].time
        y = mrvdf[sel].mnvel - np.nanmean(mrvdf[sel].mnvel)
        y_err = mrvdf[sel].errvel
        ax.errorbar(
            x, y, y_err, marker='o', elinewidth=0.5, capsize=4, lw=0, mew=0.5,
            color=f'C{ix}', markersize=3, zorder=-ix+5, label=tel
        )
    ax.legend(loc='best', fontsize='xx-small')
    ax.set_ylabel('RV-$\mu$ [km/s]')

    ax.set_xlabel('BJDTDB')

    # set naming options
    outdir = os.path.join(RESULTSDIR, 'HIP_67522')
    outpath = os.path.join(outdir, 'all_rvs.png')
    savefig(fig, outpath, dpi=400)


def _linear_model(xdata, m, b):
    return m*xdata + b

def plot_rv_vs_flux():

    #
    # get data
    #
    rvdf = get_rvs()
    mrvdf = get_minerva_mean_rvs(rvdf)
    mrvdf = pd.concat([mrvdf, rvdf[rvdf.tel.str.contains("CHIRON")]])

    time,flux,flux_err,qual = get_phot()
    sel = (qual == 0) & np.isfinite(flux)
    stime,sflux,sflux_err = time[sel], flux[sel], flux_err[sel]
    mu = np.median(sflux)
    sflux = sflux / mu - 1
    sflux_err /= mu
    stime += 2457000 #BTJD to BJD

    #
    # fit spline used to predict flux at RV points
    #

    fn = UnivariateSpline(stime, sflux)
    fn.set_smoothing_factor(1e-2) # tuned in rv_flux_vs_time

    #
    # make plot
    #
    plt.close('all')
    set_style()

    fig, ax = plt.subplots(figsize=(4,4))

    tels = np.unique(mrvdf.tel)
    prvs,pflxs = [], []
    for ix, tel in enumerate(tels):

        sel = (
            (mrvdf.tel == tel)
            &
            (mrvdf.time < stime.max())
            &
            (mrvdf.time > stime.min())
        )
        t = mrvdf[sel].time

        rv = mrvdf[sel].mnvel - np.nanmean(mrvdf[sel].mnvel)
        rv_err = mrvdf[sel].errvel

        pred_flux = fn(t)

        prvs.append(rv)
        pflxs.append(pred_flux)

        ax.errorbar(
            pred_flux, rv, rv_err, marker='o', elinewidth=0.5, capsize=4, lw=0,
            mew=0.5, color=f'C{ix}', markersize=3, zorder=-ix+5, label=tel
        )

    prvs = np.hstack(prvs)
    pflxs = np.hstack(pflxs)

    slope_guess = -5e-2
    intercept_guess = 0
    p_opt, p_cov = curve_fit(
        _linear_model, pflxs, prvs,
        p0=(slope_guess, intercept_guess)#, sigma=prv_err
    )
    lsfit_slope = p_opt[0]
    lsfit_slope_err = p_cov[0,0]**0.5
    lsfit_int = p_opt[1]
    lsfit_int_err = p_cov[1,1]**0.5

    iflux = np.linspace(pflxs.min(), pflxs.max(), 1000)
    ax.plot(iflux, _linear_model(iflux, lsfit_slope, lsfit_int), '-',
            label='Best fit', color='k', alpha=0.5, lw=1)

    txt = (
        f'Slope: {lsfit_slope:.4f} +/- {lsfit_slope_err:.4f}\n'
        f'implies {abs(lsfit_slope/lsfit_slope_err):.2f}σ different from zero.\n'
        f'Intercept: {lsfit_int:.5f} +/- {lsfit_int_err:.5f}'
    )

    ax.set_ylabel('RV-$\mu$ [km/s]')
    ax.set_xlabel('Flux [spline interp]')
    ax.set_title('Does flux predict RV during overlap?\n'+txt, fontsize='x-small')

    ax.legend(loc='best', fontsize='xx-small')

    # set naming options
    outdir = os.path.join(RESULTSDIR, 'HIP_67522')
    outpath = os.path.join(outdir, 'rv_vs_flux.png')
    savefig(fig, outpath, dpi=400)



def main():
    plot_all_rvs()
    plot_rv_flux_vs_time()
    plot_rv_vs_flux()

if __name__ == "__main__":
    main()
