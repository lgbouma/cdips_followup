"""
Make some quick-look light curves of SPOC 2-minute stars in Taurus.
"""
import os
from quicklooklc import quicklooklc
import numpy as np, pandas as pd
from astrobase.services.identifiers import gaiadr2_to_tic
from cdips_followup.paths import RESULTSDIR

df = pd.read_csv("/Users/luke/Dropbox/proj/cdips/data/cluster_data/"
                 "v07/Galli_2019_cut_confirmed_members_only.csv")

df['dr2_source_id'] = df['GaiaDR2'].str.replace('Gaia DR2 ','')

basedir = os.path.join(RESULTSDIR, 'quicklooklc', "TAURUS")

for ix, r in df.iterrows():

    try:
        ticid = gaiadr2_to_tic(r.dr2_source_id)
        print(42*'-')
        print(r.dr2_source_id, ticid)

        outdir = os.path.join(basedir, f'TIC{ticid}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        else:
            print(f'found {outdir} skip')
            continue

        quicklooklc(ticid,
            cdips = 0,
            outdir = outdir,
            spoc = 1,
            eleanor = 0,
            cdipspre = 0,
            kepler = 0,
            qlp = 0,
            detrend = 'biweight', # 'biweight' # None, 'biweight', or 'pspline'
            period = None,
            epoch = None,
            badtimewindows = None,
            do_mag_lcs = 0,
            do_eleanor_lcs = 0,
            do_flux_lcs = 1,
            do_periodogram = 1,
            do_pf = 0,
            require_quality_zero = 1,
            forceylim = None # [0.93, 1.07]# for the flux light curves
        )
    except Exception as e:
        print(42*'-')
        print(f'GAIA DR2 {r.dr2_source_id} FAILED!!!')
        print(42*'-')
