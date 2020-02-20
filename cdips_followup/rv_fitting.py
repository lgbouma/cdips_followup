"""
Given RV vs time, does it phase up?
"""
import numpy as np, pandas as pd

from radvel import driver
import os
import emcee
if not emcee.__version__ == "2.2.1":
    raise AssertionError('radvel requires emcee v2')

from cdips.utils import today_YYYYMMDD
from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
DRIVERDIR = os.path.join(os.path.dirname(__path__[0]), 'drivers/radvel_drivers')



###########
# classes #
###########

class args_object(object):
    """
    a minimal version of the "parser" object that lets you work with the
    high-level radvel API from python. (without directly using the command line
    interface)
    """
    def __init__(self, setupfn, outputdir):
        # return args object with the following parameters set
        self.setupfn = setupfn
        self.outputdir = outputdir
        self.decorr = False
        self.plotkw = {}
        self.gp = False


def convert_vels_to_radvel_ready(ticid, is_bin):
    # "vels" is the Paul Butler pipeline format -> "radvel_vels" to be radvel
    # readable

    if is_bin:
        name = 'HD{}_PFSBIN.vels'.format(ticid)
    else:
        name = 'HD{}_PFS.vels'.format(ticid)

    datadir = os.path.join(DATADIR, 'PFS', 'butler_vels')
    velspath = os.path.join(datadir, name)

    if is_bin:
        names = ['time','mnvel','errvel']
    else:
        names = ['time','mnvel','errvel','idk0','idk1','idk2','exptime']

    df = pd.read_csv(velspath, delim_whitespace=True, names=names)

    df['tel'] = 'PFS'
    df['Name'] = 'TIC{}'.format(ticid)
    df['source'] = velspath

    outdir = os.path.join(DATADIR, 'PFS', 'radvel_vels', today_YYYYMMDD())
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outname = 'TIC{}.csv'.format(ticid)

    outpath = os.path.join(outdir, outname)

    df.to_csv(outpath, index=False, sep=',')
    print('made {}'.format(outpath))

    return outpath


def prepare_template(ticid, periodval, t0val, timebase,
                     rv_path, period_unc, t0_unc, k_prior_init=100):
    # make a radvel_driver for this target

    template_path = os.path.join(DRIVERDIR, 'template.py')
    with open(template_path, 'r') as f:
        lines = f.readlines()

    d = {
        'STARNAME': 'TIC{}'.format(ticid),
        'PERIODVAL': str(periodval),
        'T0VAL': str(t0val),
        'KPRIORINIT': str(k_prior_init),
        'MEDIANTIMEBASE': str(timebase),
        'RV_PATH': rv_path,
        'PERIOD_UNC': str(period_unc),
        'T0_UNC': str(t0_unc)
    }

    for ix, l in enumerate(lines):
        for k,v in d.items():
            if k in l:
                lines[ix] = l.replace(k, v)

    outlines = lines

    outdir = os.path.join(DRIVERDIR, today_YYYYMMDD())
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outpath = os.path.join(outdir, 'TIC{}.py'.format(ticid))

    with open(outpath, mode='w') as f:
        f.writelines(outlines)
    print('made {}'.format(outpath))

    return outpath


def run_radvel(driverpath, outdir):

    setupfn = driverpath
    outputdir = outdir # "/home/luke/Dropbox/proj/WASP-4b_anomaly/results/rv_fitting/LGB_20200202_fix_gammaddot"
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    args = args_object(setupfn, outputdir)

    # # perform max-likelihood fit. usually needed to be done first.
    # radvel fit -s $basepath
    driver.fit(args)

    # # plot the maxlikelihood fit
    # radvel plot -t rv -s $basepath
    args.type = ['rv']
    driver.plots(args)

    # # # perform mcmc to get uncertainties
    # # radvel mcmc -s $basepath
    # args.nsteps = 10000  # Number of steps per chain [10000]
    # args.nwalkers = 50   # Number of walkers. [50]
    # args.ensembles = 16   # Number of ensembles. Will be run in parallel on separate CPUs [8]
    # args.maxGR = 1.01    # Maximum G-R statistic for chains to be deemed well-mixed and halt the MCMC run [1.01]
    # args.burnGR = 1.03   # Maximum G-R statistic to stop burn-in period [1.03]
    # args.minTz = 1000    # Minimum Tz to consider well-mixed [1000]
    # args.minsteps = 1000 # Minimum number of steps per walker before convergence tests are performed [1000].
    #                      # Convergence checks will start after the minsteps threshold or the minpercent threshold has been hit.
    # args.minpercent = 5  # Minimum percentage of steps before convergence tests are performed [5]
    #                      # Convergence checks will start after the minsteps threshold or the minpercent threshold has been hit.
    # args.thin = 1        # Save one sample every N steps [default=1, save all samples]
    # args.serial = False  # If True, run MCMC in serial instead of parallel. [False]
    # driver.mcmc(args)

    # # # corner plot the samples
    # # radvel plot -t rv corner trend -s $basepath
    # args.type = ['rv','corner','trend']
    # driver.plots(args)

    # # # make a sick pdf report
    # # radvel report -s $basepath
    # args.comptype= 'ic' # Type of model comparison table to include. Default: ic
    # args.latex_compiler = 'pdflatex' # path to latex compiler
    # # driver.report(args)

    # # # optionally, include stellar parameters to derive physical parameters for the
    # # # planetary system
    # # radvel derive -s $basepath
    # driver.derive(args)

    # # # optionally, make corner plot for derived parameters
    # # radvel plot -t derived -s $basepath
    # args.type = ['derived']
    # driver.plots(args)

    # # # do model comparison. valid choices: ['nplanets', 'e', 'trend', 'jit', 'gp']
    # # radvel ic -t nplanets e trend -s $basepath
    # args.type = ['nplanets', 'e', 'trend', 'jit']
    # args.mixed = True      # flag to compare all models with the fixed parameters mixed and matched rather than treating each model comparison separately. This is the default.
    # args.unmixed = False   # flag to treat each model comparison separately (without mixing them) rather than comparing all models with the fixed parameters mixed and matched.
    # args.fixjitter = False # flag to fix the stellar jitters at the nominal model best-fit value
    # args.verbose = True    # get more details

    # driver.ic_compare(args)

    # # # make the final report
    # # radvel report -s $basepath
    # driver.report(args)
