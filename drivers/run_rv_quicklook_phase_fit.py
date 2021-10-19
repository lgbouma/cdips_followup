"""
Environment: py37
Given RV vs time, does it phase up?

Contents:
    "manualmain" and "automain" modes for producing radvel_drivers that then
    run the Keplerian fit.
"""

import os
import numpy as np, pandas as pd

from cdips_followup import __path__
from cdips_followup.manage_ephemerides import query_ephemeris

from cdips_followup.rv_fitting import (
    convert_vels_to_radvel_ready,
    prepare_template,
    run_radvel
)

##########
# config #
##########

def main():

    # ticid = '59859387'
    # instrument = 'CHIRON'
    # k_prior_init = 500
    # log_k_prior_high = np.log(1000)
    # log_k_prior_low = np.log(100)
    # is_bin = 0
    # phase_with_prot = 1

    # ticid = '59859387'
    # is_bin = 0
    # instrument = 'PFS'

    ticid = '268301217'
    is_bin = 1
    instrument = 'PFS'
    k_prior_init = 100
    log_k_prior_high = np.log(800)
    log_k_prior_low = np.log(10)
    phase_with_prot = 0

    do_mcmc = 0
    allowlineartrend = 1

    if not instrument in ['PFS', 'CHIRON']:
        errmsg = 'need to implement instrument specific read functions'
        raise NotImplementedError(errmsg)

    rv_path = convert_vels_to_radvel_ready(ticid, is_bin, instrument)

    rv_df = pd.read_csv(rv_path)
    timebase = np.nanmedian(rv_df.time)

    e = query_ephemeris(ticid=ticid)

    if phase_with_prot:
        period_mean = 3
        period_std = 1.5
        epoch_mean = 2458888.13511
        epoch_std = 2
        driver_path = prepare_template(ticid, period_mean, epoch_mean, timebase,
                                       rv_path, period_std, epoch_std,
                                       instrument, k_prior_init=k_prior_init,
                                       log_k_prior_high=log_k_prior_high,
                                       log_k_prior_low=log_k_prior_low,
                                       allowlineartrend=allowlineartrend)

    else:
        # phase with planet orbit
        driver_path = prepare_template(ticid, e['period'], e['epoch'], timebase,
                                       rv_path, e['period_unc'], e['epoch_unc'],
                                       instrument, k_prior_init=k_prior_init,
                                       log_k_prior_high=log_k_prior_high,
                                       log_k_prior_low=log_k_prior_low,
                                       allowlineartrend=allowlineartrend)

    RVRESULTDIR = os.path.join(os.path.dirname(__path__[0]), 'results',
                               'spec_analysis', instrument, 'phased_RVs')
    if not os.path.exists(RVRESULTDIR):
        os.mkdir(RVRESULTDIR)

    run_radvel(driver_path, RVRESULTDIR, do_mcmc=do_mcmc)


def manual_main():
    # If you manually wrote the driver file (rather than automatically
    # producing it). E.g., right now, this is for multi-instrument fitting.
    driver_path = os.path.join(
        os.path.dirname(__path__[0]),
        'drivers',
        'radvel_drivers',
        '20211014', 'HIP67522_CHIRON.py'
    )

    RVRESULTDIR = os.path.join(os.path.dirname(__path__[0]), 'results',
                               'spec_analysis', 'multi_instrument',
                               'phased_RVs_CHIRON')
    if not os.path.exists(RVRESULTDIR):
        os.mkdir(RVRESULTDIR)

    run_radvel(driver_path, RVRESULTDIR)


if __name__ == "__main__":
    automain = 0
    manualmain = 1

    if automain:
        main()
    if manualmain:
        manual_main()
