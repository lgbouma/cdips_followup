"""
How many observations do we need over 21A to measure the mass of HIP 67522b?
"""

import os
import numpy as np, pandas as pd, pymc3 as pm

from numpy import array as nparr
from os.path import join

from astropy.time import Time
from astropy import units as u

from importlib.machinery import SourceFileLoader
from collections import OrderedDict

from betty.paths import BETTYDIR
import betty.plotting as bp
from betty.modelfitter import ModelFitter

from cdips_followup.paths import RESULTSDIR
from cdips_followup.rv_simulation import get_simulated_rvs, plot_korbpost

#
# parameters for the `exoplanet` planet model (mostly from Rizzuto+20)
# Table 5 and Table 1
#
ecc = 0.061
t0 = 2458604.02358
period = 6.95993 # days
b = 0.134
omega = -17*2*np.pi/360
m_star = 1.22
r_star = 1.38
K_orb = 93 # *u.m/u.s

#
# parameters for the `exoplanet` stellar rotation model
#
phot_rot_amp = 2e-2
P_rot = 1.418 # *u.day
rv_rot_amp = 1080 # *u.m/u.s
log_amp = np.log(51371) # log of variance of data
log_Q0 = 4.59
log_deltaQ = np.log(np.exp(4.59) - np.exp(2.00))
mix = 0.9905

xoparams = {
    # planet
    'ecc': ecc, 't0':t0, 'period':period, 'b':b, 'omega':omega,
    'm_star':m_star, 'r_star':r_star, 'K_orb':K_orb,
    # star
    'log_amp':log_amp, 'P_rot':P_rot,
    'log_Q0':log_Q0, 'log_deltaQ':log_deltaQ, 'mix':mix
}

sigma_obs = 20 # m/s
N_obs = 60 # really, asking for 45, but why not.

start_time = Time('2021-04-28 23:59:00')
end_time = Time('2021-05-26 23:59:00')

df = pd.read_csv('/Users/luke/Dropbox/proj/cdips_followup/data/lightcurves/hip_67522.csv')
df['time'] += 2400000
t_diff = start_time.jd - df.time.min()
df['time'] += t_diff

df['rv'] = (df['flux'] - np.nanmean(df['flux']))*rv_rot_amp**(3/2)
xo_df = df

##########################################
# generate / get the fake data
##########################################

plotdir = join(RESULTSDIR, 'HIP_67522')

obs_csvpath = join(plotdir, f'hip_67522_obs_df_Nobs{N_obs}.csv')
true_csvpath = join(plotdir, f'hip_67522_true_df_Nobs{N_obs}.csv')
if not os.path.exists(obs_csvpath):
    obs_df, true_df = get_simulated_rvs(
        sigma_obs, N_obs, start_time, end_time, xoparams, xo_df
    )
    obs_df.to_csv(obs_csvpath, index=False)
    true_df.to_csv(true_csvpath, index=False)
else:
    obs_df = pd.read_csv(obs_csvpath)
    true_df = pd.read_csv(true_csvpath)

##########################################
# fit the fake data
##########################################

starid = 'HIP67522'
modelid = 'rvspotorbit'

datasets = OrderedDict()
datasets['minerva'] = [nparr(obs_df.t_obs),
                       nparr(obs_df.rv_obs),
                       nparr(obs_df.rv_err_obs)]

pklpath = join(BETTYDIR, f'test_{starid}_{modelid}_Nobs{N_obs}.pkl')

priorpath = '/Users/luke/Dropbox/proj/cdips_followup/data/HIP_67522/HIP_67522_priors.py'
priormod = SourceFileLoader('prior', priorpath).load_module()
priordict = priormod.priordict

N_samples = 1000
m = ModelFitter(modelid, datasets, priordict, plotdir=plotdir,
                pklpath=pklpath, overwrite=0, N_samples=N_samples,
                N_cores=os.cpu_count(), target_accept=0.8)

print(pm.summary(m.trace, var_names=list(priordict)))

summdf = pm.summary(m.trace, var_names=list(priordict), round_to=10,
                    kind='stats', stat_funcs={'median':np.nanmedian},
                    extend=True)

cornerplot = 1
onedpostplot = 1
korbpostplot = 1

if cornerplot:
    outpath = join(plotdir, f'{starid}_{modelid}_cornerplot_Nobs{N_obs}.png')
    bp.plot_cornerplot(list(priordict), m, outpath)

if onedpostplot:
    outpath = join(plotdir, f'{starid}_{modelid}_1dpost_Nobs{N_obs}.png')
    bp.plot_1d_posterior(m.trace['K_orb'], outpath, truth=K_orb,
                         xlabel='K$_\mathrm{orb}$ [m$\,$s$^{-1}$]')

if korbpostplot:
    outpath = join(plotdir, f'{starid}_{modelid}_korbpost.png')
    plot_korbpost(outpath, truth=K_orb,
                  xlabel='K$_\mathrm{orb}$ [m$\,$s$^{-1}$]')
