"""
For Rafael's CHAT db inserts, generate a comma separated text file with the
following rows for the new candidates so he can add them directly to the
database?

name, ra, dec, vmag, T0, P, dur, ecc, omega, band, focus, expt, posx, posy
TIC59859387.01,07:21:28.72,-45:34:03.85,10.735,2458492.61762789,1.12036920,0.0739,,,i,Focused,,,,
"""

##########
# config #
##########

import numpy as np, pandas as pd
from cdips_followup.manage_ephemerides import query_ephemeris
from cdips_followup.utils import get_cdips_candidates
from astroquery.gaia import Gaia
#
# 467522800.01 through 1030389830.01, taken from the candidate db
#
source_ids = [5240531632175135616, 5246508676134913024, 5252320149182351744,
              5254794943702272640, 5255504399239636992, 5256502966294604800,
              5256717749007641344, 5304507540851982592, 5308445029440522112,
              5312011948237411456, 5312141072137904512, 5321581959996214016,
              5322083917107299712, 5329853517263005696, 5338645555866749952,
              5339389268061191040, 5404579488593432576, 5405530566146401408,
              5406188928778889472, 5489726768531119616, 5514373043228306816,
              5523449717870971776, 5544257185155788672, 5548219137873425280,
              5790098845432414464, 5838183443852841216, 5869166341556812032 ]

##########################################
##########################################

df = get_cdips_candidates()
df = df.drop(['period'], axis=1)

ephems = []
for source_id in source_ids:
    d = query_ephemeris(source_id=str(source_id))
    ephems.append(d)

edf = pd.DataFrame(ephems)

assert len(edf) == len(source_ids)

edf.source_id = edf.source_id.astype(str)
df.source_id = df.source_id.astype(str)
bigdf = edf.merge(df, how='left', on=['source_id'], suffixes=('', '_canddb'))

assert len(bigdf) == len(edf)

outpath = (
    '../../results/SG1_candidate_sharing/20200121_s8_to_s11_cdips_cands.csv'
)
bigdf.to_csv(outpath, index=False, sep='|')
print('made {}'.format(outpath))

#
# rafael's format
#
scols = ['targetid', 'gaia_ra', 'gaia_dec', 'tic_Vmag', 'epoch', 'period',
         'duration']

rdf = bigdf[scols]

rafael_name_dict = {
    'name':'targetid', 'ra':'gaia_ra', 'dec':'gaia_dec', 'vmag':'tic_Vmag',
    'T0':'epoch', 'P':'period', 'dur':'duration'
}
rafael_names = {v: k for k, v in rafael_name_dict.items()}

rdf = rdf.rename(rafael_names, axis=1)

# fix ra/dec to be gaia precision
ras, decs = [], []
for ix, source_id in enumerate(source_ids):
    print('{}/{}'.format(ix, len(source_ids)))
    jobstr = (
        "select top 1 g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag from "
        "gaiadr2.gaia_source as g where g.source_id = {:d}".
        format(source_id)
    )

    job = Gaia.launch_job(jobstr)
    gaia_r = job.get_results()

    if len(gaia_r) != 1:
        raise AssertionError('gaia match failed')

    ra, dec = float(gaia_r['ra']), float(gaia_r['dec'])
    ras.append(ra)
    decs.append(dec)

rdf['ra'] = ras
rdf['dec'] = decs

# faint M dwarfs: make up the mag.
rdf.loc[rdf.vmag<0, 'vmag'] = 17.5

# duration should be in days for rafael
rdf['dur'] = rdf['dur']/24

nancols = ['ecc', 'omega', 'band', 'focus', 'expt', 'posx', 'posy']
for ncol in nancols:
    rdf[ncol] = np.nan

rdf['band'] = 'i'
outpath = (
    '../../results/SG1_candidate_sharing/20200121_s8_to_s11_cdips_cands_rafael.csv'
)
rdf.to_csv(outpath, index=False, sep=',')
print('made {}'.format(outpath))
