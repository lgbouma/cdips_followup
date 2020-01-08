"""
After collecting basic source identifier information for candidate planets
around candidate young stars, we need to supplement the information for
observing proposals.

(Like: where are the objects? When are they observable? From where?)

For instance, in the 20190912_toi_youngstar_xmatching google sheet, the
cdips_KC19_xm_to_TOI_plus_my_CDIPS_cands sub-sheet contains candidates from:

    * 20190825 TOI list cross-matched to Kounkel & Covey 2019, and to the CDIPS
    target star list.

    * CDIPS planet candidates from paper 1.
"""

import pandas as pd
from astrobase.services.gaia import objectid_search

def match_google_sheet_to_gaia():

    df = pd.read_csv(r'../data/20190912_toi_youngstar_xmatching '
                     '- cdips_KC19_xm_to_TOI_plus_my_CDIPS_cands.csv')

    _results = []
    for ix, r in df.iterrows():

        res_dict = objectid_search(r['source_id'])
        _results.append(pd.read_csv(res_dict['result']))

    gaia_df = pd.concat(_results)

    mdf = df.merge(gaia_df, how='left', on='source_id')

    outpath = '../data/20190912_youngstar_cands_with_gaia.csv'
    mdf.to_csv(outpath, index=False)

import IPython; IPython.embed()
