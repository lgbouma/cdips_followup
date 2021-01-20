"""
We want to quantify the output of the CDIPS project, for proposals.
Numbers will include (in latex format):

    numctois, numfalsepos, numneedobs, numconfirmable, numtoofaint
"""
import numpy as np, pandas as pd
from cdips.utils.catalogs import get_exofop_ctoi_catalog
from cdips_followup.utils import get_cdips_candidates

#
# Get CTOIs I've put on ExoFOP.  selection function is things that I've put up,
# with Gaia age provenance (my comment strings were inconsistent).
#
df_ctoi = get_exofop_ctoi_catalog()

sel_ctoi = (
    (df_ctoi.User == 'bouma')
    &
    (df_ctoi.Notes.str.contains('CG') |
     df_ctoi.Notes.str.contains('Zari') |
     df_ctoi.Notes.str.contains('KC') |
     df_ctoi.Notes.str.contains('Nice') |
     df_ctoi.Notes.str.contains('IC2602') |
     df_ctoi.Notes.str.contains('Collinder135')
    )
)

sdf_ctoi = df_ctoi[sel_ctoi]

#
# Get updates on them from candidate database.
#

df_cand = get_cdips_candidates()

mdf = sdf_ctoi.merge(df_cand, how='left', left_on='TIC ID', right_on='ticid')

scols = [
    'TIC ID', 'CTOI', 'Promoted to TOI', 'TESS Mag', 'TESS Mag err', 'RA',
    'Dec', 'Depth ppm', 'Depth ppm Error', 'Duration (hrs)', 'Duration (hrs) Error',
    'Radius (R_Earth)', 'Radius (R_Earth) Error', 'Stellar Distance (pc)',
    'Stellar Radius (R_Sun)', 'Stellar Radius (R_Sun) err', 'CTOI lastmod',
    'User', 'Group', 'Tag', 'Notes',
    'source_id', 'ticid', 'toi', 'targetid', 'iscdipstarget', 'reference',
    'name', 'age', 'nbhd_rating', 'init_priority', 'current_priority',
    'pending_spectroscopic_observations',
    'pending_photometry_observations', 'comment', 'rp', 'rp_unc', 'period',
    'gaia_ra', 'gaia_dec', 'gaia_plx', 'gaia_Gmag', 'gaia_Bmag', 'gaia_Rmag',
    'tic_Bmag', 'tic_Vmag', 'tic_Jmag', 'tic_Hmag', 'tic_Kmag', 'tic_Tmag',
    'tic_teff', 'tic_logg', 'tic_rstar', 'tic_mstar', 'candidate_provenance',
    'insert_time', 'last_update_time', 'isretired', 'disposition'
]

smdf = mdf[scols]

confirmableinV = (
    (smdf.isretired == 0)
    &
    (smdf.tic_Vmag < 14)
    &
    (smdf.tic_teff < 7000)
)

needobs = (
    (smdf.isretired == 0)
    &
    (smdf.pending_spectroscopic_observations.str.contains('SP0') |
     smdf.pending_spectroscopic_observations.str.contains('SP1') |
     smdf.pending_photometry_observations.str.contains('PP0') |
     smdf.pending_photometry_observations.str.contains('PP1')
    )
)

faintinT = (
    (smdf.isretired == 0)
    &
    (smdf.tic_Tmag > 15)
)

# the number of age-dated CTOIs produced by the CDIPS project, per ExoFOP-TESS
numctois = len(sdf_ctoi)
# the number of these that have significant evidence for not being planets
# (e.g., EBs, NEBs, pEBs, SB2s, dippers, etc)
numfalsepos = len(smdf[smdf.isretired == 1])
# the number of active candidates that are flagged as requiring either
# photometric or spectroscopic observations
numneedobs = len(smdf[needobs])
# the number that have V<14 and TIC Teff<7000K, suggesting masses might be
# measurable.
numconfirmableinV = len(smdf[confirmableinV])
# the number that have T>15.
numfaintinT = len(smdf[faintinT])

out_txt = (
r'\newcommand{\numctois}{'+str(numctois)+'\ }\n'+
r'\newcommand{\numfalsepos}{'+str(numfalsepos)+'\ }\n'+
r'\newcommand{\numneedobs}{'+str(numneedobs)+'\ }\n'+
r'\newcommand{\numconfirmableinV}{'+str(numconfirmableinV)+'\ }\n'+
r'\newcommand{\numfaintinT}{'+str(numfaintinT)+'\ }'
)

print(out_txt)
