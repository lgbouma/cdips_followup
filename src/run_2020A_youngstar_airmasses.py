import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
import os
from glob import glob

from get_airmass_charts import make_airmass_chart

from astropy.time import Time

def main(site='keck'):

    df = pd.read_csv('../data/20190912_youngstar_cands_with_gaia.csv')

    for name, ra, dec in zip(
        nparr(df['source_id']), nparr(df['ra']), nparr(df['dec'])
    ):

        outdir = '../results/followup_planning/{}'.format(name)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        make_airmass_chart(
            name=name, site=site,
            ra=ra, dec=dec,
            outdir=outdir,
            start_time=Time('2019-09-13 20:00:00', format='iso'),
            end_time=Time('2020-09-13 20:00:00', format='iso'),
            check_months_observable=False
        )

if __name__ == "__main__":
    # main(site='keck') # HIRES, Subaru
    main(site='Kitt Peak National Observatory') # TRES(?), NEID
    main(site='Las Campanas Observatory') # PFS, PISCO
    # main(site='lapalma') # HARPS-N
