"""
Environment: py36_emcee2
Given RV vs time, does it phase up?
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

RVRESULTDIR = os.path.join(os.path.dirname(__path__[0]), 'results',
                           'spec_analysis', 'PFS', 'phased_RVs')
if not os.path.exists(RVRESULTDIR):
    os.mkdir(RVRESULTDIR)


def main():

    ticid = '59859387'
    is_bin = 0

    # ticid = '268301217'
    # is_bin = 1

    is_pfs = 1

    if not is_pfs:
        errmsg = 'need to implement instrument specific read functions'
        raise NotImplementedError(errmsg)

    rv_path = convert_vels_to_radvel_ready(ticid, is_bin)

    rv_df = pd.read_csv(rv_path)
    timebase = np.nanmedian(rv_df.time)

    e = query_ephemeris(ticid=ticid)

    driver_path = prepare_template(ticid, e['period'], e['epoch'], timebase,
                                   rv_path, e['period_unc'], e['epoch_unc'])

    run_radvel(driver_path, RVRESULTDIR)


if __name__ == "__main__":
    main()
