import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
import os
from glob import glob

from get_observability_charts import (
    make_observability_chart_singlesite,
    make_observability_chart_multisite
)
from astropy.time import Time

def single_site(site='keck'):

    df = pd.read_csv('../data/20190912_youngstar_cands_with_gaia.csv')

    for name, ra, dec in zip(
        nparr(df['source_id']), nparr(df['ra']), nparr(df['dec'])
    ):

        outdir = '../results/followup_planning/{}'.format(name)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        # 2020A: Feb 1 - Jul 31 2020
        make_observability_chart_singlesite(
            name=name, site=site,
            ra=ra, dec=dec,
            outdir=outdir,
            start_time=Time('2019-09-13 20:00:00', format='iso'),
            end_time=Time('2020-09-13 20:00:00', format='iso')
        )

def multi_site(
    sites=['keck',
           'Las Campanas Observatory',
           'Kitt Peak National Observatory']
):

    df = pd.read_csv('../data/20190912_youngstar_cands_with_gaia.csv')

    for name, ra, dec in zip(
        nparr(df['source_id']), nparr(df['ra']), nparr(df['dec'])
    ):

        outdir = '../results/followup_planning/{}'.format(name)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        # 2020A: Feb 1 - Jul 31 2020
        make_observability_chart_multisite(
            name=name, sites=sites,
            ra=ra, dec=dec,
            outdir=outdir,
            start_time=Time('2019-09-13 20:00:00', format='iso'),
            end_time=Time('2020-09-13 20:00:00', format='iso')
        )


if __name__ == "__main__":

    multi_site()

    # single_site(site='Las Campanas Observatory') # PFS, PISCO
    # single_site(site='Kitt Peak National Observatory') # TRES(?), NEID
    # single_site(site='keck') # HIRES, Subaru
    # single_site(site='lapalma') # HARPS-N
