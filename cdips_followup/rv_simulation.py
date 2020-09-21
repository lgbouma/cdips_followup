import os, pickle
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import exoplanet as xo
import theano.tensor as tt
import pymc3 as pm
from aesthetic.plot import set_style, savefig

def get_simulated_rvs(sigma_obs, N_obs, start_time, end_time, xoparams, xo_df,
                      makeplot=True):

    #
    # a GP model is not defined without data. so _time, _rv, and _rv_err are
    # fake data, got by multiplying the photometry by the expected spot-induced
    # RV amplitude.
    #
    eps = 1e-5
    _time = np.ascontiguousarray(xo_df.time, dtype=np.float64)
    _rv = np.ascontiguousarray(xo_df.rv, dtype=np.float64)
    _rv_err = np.ones_like(_rv)*30

    t_train = np.linspace(_time.min(), _time.max(), int(1e4))
    t_true = np.linspace(start_time.jd, end_time.jd, int(1e4))
    t_long = np.linspace(_time.min(), end_time.jd, int(1e4))
    np.random.seed(42)
    t_obs = np.sort(np.random.uniform(start_time.jd, end_time.jd, N_obs))
    assert np.all(np.diff(t_obs) > 0)

    with pm.Model() as model:

        mean = pm.Normal('mean', mu=np.mean(_rv), sd=np.std(_rv))
        upper_amp = 13.8 # 14 works
        log_amp = pm.Bound(pm.Normal, lower=10, upper=upper_amp)(
            'log_amp', mu=xoparams['log_amp'], sd=10
        )
        P_rot = pm.Normal('P_rot', mu=xoparams['P_rot'], sd=0.01)
        log_Q0 = pm.Bound(pm.Normal, upper=7, lower=0.7)(
            'log_Q0', mu=xoparams['log_Q0'], sd=2
        )
        log_deltaQ = pm.Normal('log_deltaQ', mu=xoparams['log_deltaQ'], sd=0.1)
        mix = pm.Uniform('mix', lower=1e-4, upper=1)

        kernel = xo.gp.terms.RotationTerm(
            log_amp=log_amp, period=P_rot, log_Q0=log_Q0, log_deltaQ=log_deltaQ, mix=mix
        )

        gp = xo.gp.GP(kernel, _time, _rv_err**2, mean=mean)
        gp.marginal('gp', observed=_rv)

        pm.Deterministic('pred_train', gp.predict(t_train))
        pm.Deterministic('pred_long', gp.predict(t_long))
        pm.Deterministic('pred_obs', gp.predict(t_obs))
        pm.Deterministic('pred_true', gp.predict(t_true))

        map_soln = xo.optimize(start=model.test_point)

    # gp model doesn't really work because i don't understand it. period etc
    # are fine; scale is wonky.
    rv_train = map_soln['pred_train']
    rv_long = map_soln['pred_long']
    rv_obs_spot = map_soln['pred_obs']
    rv_true_spot = map_soln['pred_true']
    print(map_soln)

    orbit = xo.orbits.KeplerianOrbit(
        period=xoparams['period'], b=xoparams['b'], t0=xoparams['t0'],
        ecc=xoparams['ecc'], omega=xoparams['omega'],
        m_star=xoparams['m_star'], r_star=xoparams['r_star']
    )

    np.random.seed(42)

    rv_true_kep = orbit.get_radial_velocity(t_true, K=xoparams['K_orb']).eval()

    rv_obs_kep = orbit.get_radial_velocity(t_obs, K=xoparams['K_orb']).eval()
    rv_obs_noise = np.random.normal(loc=0, scale=sigma_obs, size=len(t_obs))

    rv_true = rv_true_spot + rv_true_kep
    rv_obs = rv_obs_spot + rv_obs_kep + rv_obs_noise

    rv_err_obs = np.ones_like(t_obs)*sigma_obs

    if makeplot:

        ####################
        # plot #1: training data

        from cdips_followup.paths import RESULTSDIR

        set_style()
        plt.close('all')

        f, ax = plt.subplots(figsize=(8,3))
        ax.scatter(_time, _rv, s=0.5, c='k', label='train data', zorder=3)
        ax.plot(t_train, rv_train, label='gp fit', lw=0.3, zorder=1)
        ax.legend()
        ax.set_xlabel('time')
        ax.set_ylabel('RV [m/s]')
        figpath = os.path.join(RESULTSDIR, 'HIP_67522',
                               f'training_data_short_Nobs{N_obs}.png')
        savefig(f, figpath, writepdf=0)

        ####################
        # plot #2: training + predict

        set_style()
        plt.close('all')

        f, ax = plt.subplots(figsize=(4,3))
        ax.scatter(_time, _rv, s=1, c='k', label='train data', zorder=3)
        ax.plot(t_long, rv_long, label='gp fit', lw=0.5, zorder=1)
        ax.legend()
        ax.set_xlabel('time')
        ax.set_ylabel('RV [m/s]')
        figpath = os.path.join(RESULTSDIR, 'HIP_67522',
                               f'training_data_long_Nobs{N_obs}.png')
        savefig(f, figpath, writepdf=0)


        ####################
        # plot #3: generated data

        plt.close('all')
        t0 = t_true.min()

        f, axs = plt.subplots(nrows=3, figsize=(4,3), sharex=True)
        axs[0].plot(t_true-t0, rv_true_spot, c='C0', lw=1)
        axs[1].plot(t_true-t0, rv_true_kep, c='C1', lw=1)

        axs[2].errorbar(t_obs-t0, rv_obs, yerr=sigma_obs, fmt='none', ecolor='k',
                        alpha=1, elinewidth=1, capsize=1, zorder=12)
        axs[2].plot(t_true-t0, rv_true, c='C2', lw=1)
        axs[2].plot(t_true-t0, rv_true_kep, c='C1', lw=1)

        axs[0].set_ylabel('RV spot')
        axs[1].set_ylabel('RV planet')
        axs[2].set_ylabel('Sum [m/s]')
        axs[2].set_xlabel('Days from 2021-Apr-28')

        figpath = os.path.join(RESULTSDIR, 'HIP_67522', f'sim_spot_Nobs{N_obs}.png')
        savefig(f, figpath, writepdf=0)

    obs_df = pd.DataFrame({
        't_obs':t_obs,
        'rv_obs':rv_obs,
        'rv_err_obs':rv_err_obs
    })

    true_df = pd.DataFrame({
        'rv_true_spot':rv_true_spot,
        'rv_true_kep':rv_true_kep,
        'rv_true':rv_true,
        't_true':t_true
    })

    return obs_df, true_df


def plot_korbpost(outpath, truth=None,
                  xlabel='K$_\mathrm{orb}$ [m$\,$s$^{-1}$]'):

    Nobs = [20, 30, 45, 90]

    Korb_d = {}

    from betty.paths import BETTYDIR

    for N in Nobs:
        pklpath = os.path.join(
            BETTYDIR, f'test_HIP67522_rvspotorbit_Nobs{N}.pkl'
        )

        with open(pklpath, 'rb') as f:
            d = pickle.load(f)
            trace = d['trace']

        samples = trace['K_orb']
        Korb_d[N] = samples

    # make plot

    set_style()
    plt.close('all')

    f, ax = plt.subplots(figsize=(4,3))

    for N in Nobs:
        ls = ':' if N in [30, 90] else '-'

        ax.hist(Korb_d[N], bins=np.arange(0,220,20),
                weights=np.ones_like(Korb_d[N])*1/len(Korb_d[N]),
                zorder=-1, histtype='step', alpha=1, ls=ls,
                label='N$_\mathrm{obs}$='+f'{N}')

    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    for ix, h in enumerate(handles):
        ls = ':' if ix in [1, 3] else '-'
        new_handles.append(Line2D([], [], c=h.get_edgecolor(), ls=ls))

    ax.legend(loc='best', fontsize='small', handles=new_handles, labels=labels)

    if isinstance(xlabel, str):
        ax.set_xlabel(xlabel)

    if truth is not None:
        ymin, ymax = ax.get_ylim()
        ax.vlines(
            truth, ymin, ymax, colors='C0', alpha=0.5,
            linestyles='--', zorder=-2, linewidths=0.5
        )
        ax.set_ylim((ymin, ymax))

    ax.set_ylabel('Relative probability')

    savefig(f, outpath, writepdf=0)
