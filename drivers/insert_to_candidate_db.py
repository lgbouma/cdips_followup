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
    pending_spectroscopic_observations = '--'):


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
        raise_error_if_duplicate=True
    )


def manual_insert_many():

    sourceids = [
        '2103737241426734336',
        '2133589776693414016',
        '2132152474478407040',
        '772430527947893632',
        '2103714637018491136'
    ]

    comments = [
        'Kepler-1627 b, validated planet by Morton. PP0, SP0. Stephenson_1 is <~30 Myr old (KC19). The local KC19 members around this star show a clean main sequence. The G2V host has a ~2d rotation period, consistent with youth. With P=7.2d, Rp~=3.7Re, this is a high priority for SPECFU, incl RM.',
        'Kepler-52bcd, planets confirmed thru TTV. Compact multi in UBC_1 = Theia 520 (Curtis is finding ~350 Myr). Radii 2, 2, 2 Re. P 7.8, 16.4, 36.5.  Cluster also mentioned by Castro-Ginard+18. G=15.1 target star.',
        'Kepler-968bc. Confirmed (or validated?) system. 2 tranet in UBC_1 = Theia 520 (Curtis finding ~350 Myr). 2Re inner, 1.8Re outer. Teff 4500K host... G=14.3. Maybe can do things with Keck.',
        'GJ 1148bc. CARN (Carina-near; ~200Myr) membership not previously noted in literature. Inner ~Msaturn, e=0.3, M4V star HZ. Age~200Myr. Very metal rich ([Fe/H]=+0.32+/-0.05). Zuckerman+06 is the Gagne+18 reference.  I think this would be mainly cool to confirm because of the weird system architecture, and the high metallicity.  Also, if the metallicity of Carina-near in general is this high, that could be pretty useful for planet-hunting.',
        'PP2, SP1. 1.8Re, P=9.3d. Faint. Stephenson_1 is legit (per Kepler-1627b), but this star is faint and the planet is small. Also, the star rotation period (10ish days IIRC?) did not line up with the expectations for 30 Myr -- would be more like 1 Gyr.  If it were 30 Myr though, its size would place new limits on how SMALL a planet can be at 30 Myr.  Interesting, since most planets that young are large.'
    ]

    for s, c in zip(sourceids, comments):

        ticid = gaiadr2_to_tic(s)

        insert_single_candidate(
            ticid = ticid,
            comment = c,
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'to_assess',
            pending_spectroscopic_observations = 'to_assess'
        )


if __name__ == "__main__":

    insert_single = 0
    insert_many = 1

    if insert_single:

        insert_single_candidate(
            ticid = '235005571',
            comment = 'George Zhou found. V=12.41 late K / early-M dwarf in sector-28 with a 1.3d period hot Neptune at Rp=4Re. It is high in the HRD, potentially implying a photometric binary, or perhaps that it is on the PMS. Has archival UVES (ESO) spectrum. No sign of SB2, exhibits no Li absorption, but shows significant p-cygni structure in H-alpha, and a very complicated emission profile in Ca HK. The lack of Li is OK, given the speed of Li depletion in low mass stars. George suggests to keep on going with this target too. ',
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'PP0. try ip/gp.',
            pending_spectroscopic_observations = 'SP2.'
        )

    if insert_many:
        manual_insert_many()
