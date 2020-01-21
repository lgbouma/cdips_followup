"""
Given a list of source_ids, and a string defining the "run name", and start and
end times, make plots showing the number of hours targets are visible per
night.

Usage:
Copy-paste in source_ids for CANDIDATES below. Set START_TIME and END_TIME.
Then run via

    $ python run_candidate_observabilities.py
"""

###########
# imports #
###########

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
import os
from glob import glob

from cdips_followup.utils import (
    get_cdips_candidates, given_sourceid_get_radec
)

from get_observability_charts import (
    make_observability_chart_singlesite,
    make_observability_chart_multisite
)
from astropy.time import Time

##########
# config #
##########

#
# copy-paste these in from a google-sheet view of candidates.csv. for instance,
# you might take "specfu" candidates via:
# =QUERY(CANDIDATES_CURRENT!A:ZZ, "select * where K<2 and AD<7500 and not L
# contains 'SP2' and not L contains '--'", -1)
#
CANDIDATES = np.array([
    4844691297067063424, 5519619186857962112, 5525188767305211904,
    5251470948229949568, 5290721997195236480, 5557593814516968960,
    5838450865699668736, 2919143383943171200, 3340674976430098688,
    5290781443841554432, 5510676828723793920, 5516140233292943872,
    5596735638203997824, 5240531632175135616, 5246508676134913024,
    5254794943702272640, 5256717749007641344, 5304507540851982592,
    5312141072137904512, 5322083917107299712, 5329853517263005696,
    5338645555866749952, 5339389268061191040, 5489726768531119616,
    5514373043228306816, 5523449717870971776, 5544257185155788672,
    5548219137873425280, 5838183443852841216, 5765748511163751936
]).astype(str)
canddf = pd.DataFrame({'source_id': CANDIDATES})

START_TIME = Time('2019-09-13 20:00:00', format='iso')
END_TIME = Time('2020-09-13 20:00:00', format='iso')

RUN_NAME = '20200121_SPECFU_SP0_and_SP1'

###########
# drivers #
###########

def single_site(site='keck'):

    df = get_cdips_candidates()
    df.source_id = df.source_id.astype(str)
    df = canddf.merge(df, on='source_id', how='left')
    assert len(df) == len(canddf)

    for name, ra, dec in zip(
        nparr(df['targetid']), nparr(df['gaia_ra']), nparr(df['gaia_dec'])
    ):

        outdir = '../results/followup_planning/{}'.format(name)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        make_observability_chart_singlesite(
            name=name, site=site, ra=ra, dec=dec, outdir=outdir,
            start_time=START_TIME, end_time=END_TIME
        )


def multi_site(
    sites=['keck',
           'Las Campanas Observatory',
           'Kitt Peak National Observatory']
    ):
    """
    Campanas: PFS, PISCO
    Kitt Peak: TRES(?), NEID
    Keck: HIRES, Subaru
    La-Palma: HARPS-N
    """

    df = get_cdips_candidates()
    df.source_id = df.source_id.astype(str)
    df = canddf.merge(df, on='source_id', how='left')
    assert len(df) == len(canddf)

    for source_id, name, ra, dec in zip(
        CANDIDATES, nparr(df['targetid']),
        nparr(df['gaia_ra']), nparr(df['gaia_dec'])
    ):

        outdir = '../results/followup_planning/{}'.format(RUN_NAME)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if ra == dec == -1:
            ra, dec = given_sourceid_get_radec(source_id)

        print(name, ra, dec)
        make_observability_chart_multisite(
            name=name, sites=sites, ra=ra, dec=dec, outdir=outdir,
            start_time=START_TIME, end_time=END_TIME, save_csv=True
        )
        print('done with {}...'.format(name))


if __name__ == "__main__":
    multi_site()
