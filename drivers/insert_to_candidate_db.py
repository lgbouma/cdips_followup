"""
One-time inserts of new young planet candidate to candidates.csv

insert_single: Manually insert 1 new candidate.
insert_many: Manually insert multiple new candidates.
"""
from cdips_followup.manage_candidates import insert_candidate
from cdips.utils import today_YYYYMMDD
import numpy as np, pandas as pd

from astrobase.services.identifiers import gaiadr2_to_tic

def insert_single_candidate(
    ticid = None,
    comment = '',
    nbhd_rating = 1,
    init_priority = 1,
    pending_photometry_observations = '--',
    pending_spectroscopic_observations = '--',
    raise_error_if_duplicate=True
):


    current_priority = init_priority

    d = {
        'nbhd_rating': 1,
        'init_priority': init_priority,
        'current_priority': current_priority,
        'pending_spectroscopic_observations': pending_spectroscopic_observations,
        'pending_photometry_observations': pending_photometry_observations,
        'comment': comment,
        'candidate_provenance': today_YYYYMMDD()+'_manual_insert',
        'isretired': 0
    }

    insert_candidate(
        source_id=None, ticid=ticid, manual_dict=d,
        raise_error_if_duplicate=raise_error_if_duplicate
    )


def manual_insert_many():

    # sourceids = [
    #     '5404435280772350848'
    # ]
    # comments = [
    #     "Probably multi. Prot 3d? is consistent with the youth, b/c Teff 6500K... but also consistent with a good deal older than 100 Myr.   Neighborhood plot looks ok.  Definitely worth getting Li, spectra, etc to check"
    # ]

    sourceids = ['184854434646126080', '212991482660701184',
                 '241728185423117184', '253885554213280896',
                 '332007951305099776', '431682700462002432',
                 '490451390799691392', '521454526326028928',
                 '526763312065766016', '542223510700499968',
                 '568619413331898240', '1056826907054257664',
                 '1203566251429799680', '1522821409742783616',
                 '1538570814298529664', '1686171213716517504',
                 '1728717091031847680', '1824446002246491136',
                 '1833519030401513984', '1935439776865267072',
                 '1942514069039535360', '1968083349388709120',
                 '1970742930580834176', '1972872684597722496',
                 '1984154395449640832', '1990082343676699136',
                 '2056621793787000704', '2067807061156062720',
                 '2122642764047814144', '2123377821930505472',
                 '2136224721950935168', '2152854796663244160',
                 '2155491051885121152', '2166937311541545856',
                 '2174182749574436096', '2190239295819005440',
                 '2204810126994432384', '2222159767637696896',
                 '2270404997834401664', '3447802181031092992',
                 '4519168132218628992', '4599686193337056256']
    comments = ['' for s in sourceids]

    for s, c in zip(sourceids, comments):

        ticid = gaiadr2_to_tic(s)

        insert_single_candidate(
            ticid = ticid,
            comment = c,
            nbhd_rating = 0,
            init_priority = 1,
            pending_photometry_observations = 'to_assess',
            pending_spectroscopic_observations = 'to_assess',
            raise_error_if_duplicate=False
        )


if __name__ == "__main__":

    insert_single = 0
    insert_many = 1

    if insert_single:

        insert_single_candidate(
            ticid = '56551765',
            comment = 'V1096 Tau = Herbig Anon 1.  Taurus 1-3 Myr.  Maybe HEB system, maybe 1-3Myr HJ..  Need to check for chromaticity, then careful RVs to untangle the SB2. ',
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'PP0.',
            pending_spectroscopic_observations = 'SP1.'
        )

    if insert_many:
        manual_insert_many()
