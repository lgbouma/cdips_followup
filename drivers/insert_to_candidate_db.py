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
        #'5805824988445654144',
        #'5838183443852841216',
        #'5838450865699668736',
        #'5878452576220299520',
        #'5885456190393040512',
        #'5889360246986113024',
        #'5940265981724168064',
        #'5940850376492524800',
        #'5952590785523816960',
        #'6024033615115982848',
        #'6037044827717626112',
        #'6042883578050870912',
        #'6637496339607744768',
        #'6658373007402886400',
        #'5239758155778687360',
        #'6539037542941988736',
        #'2974906868489280768',
        #'1476485996883837184',
        #'875071278432954240',
        #'6113920619134019456',
        '2003624341205485056',
        '2067161441672227968',
        '3033676761434888320',
        '5254512781523942912',
        '5290752787819629056',
        '5312912791851427840',
    ]

    comments = [
        #'Hot Jup or Nept, Theia 597 (KC19); cluster existence dubious.',
        #'Hot Nept, PMS M dwarf (Zari+18), Prot 2d, RotAmp 2%.',
        #'Warm subNept, Theia 64 (KC19). TOI 1097.',
        #'EB or Warm Jup, NGC 5617 (170Myr, CG18).',
        #'Likely EB, NGC 5925 (180Myr, CG18). Star properties unclear. ',
        #'Likely EB, Alessi 8 (CG18, 120Myr). Prot 6d; dubious bc Teff 6500K.',
        #'Crowded. HJ or EB. Collinder 307 (CG18, 30Myr). pmDec off.',
        #'10 Myr EB or HJ. PMS K dwarf host (NGC 6193, CG18).',
        #'PMS K dwarf HJ (Zari+18). Prot 2d, RotAmp 8%. Good target.',
        #'Likely EB. PMS (Zari+18), Teff 5000K, Prot 5d, ampltiude 2%.',
        #'Hot Jup. PMS M dwarf (Zari+18), Prot 2d?, Tdur a bit long.absHot',
        #'superNept. PMS M dwarf (Zari+18). Prot  0.8d, RotAmp 5%.',
        #'Host is a subgiant field star.',
        #'HATS-47. Theia 1098 (KC19); cluster existence dubious.',
        #'PATHOS-31: IC2602 Rp=4Re NON-GRAZING candidate. High quality.',
        #'TOI-251: Zhou+20 shows <200 Myr + validated. Rp=2.7, P=5d.',
        #'TOI-942: Zhou+20 shows <100 Myr + 2 planet validated. Both RM good.',
        #'TOI-1807 + TIC-27491137. 3 planets transiting 2749, 1 USP transiting 1807. Age 100-200Myr, in a Gaia comoving group. Hughes & C. Hedges, in prep.',
        #'TOI 1726, Mann+20 THYME3 confirmed. UMa MG member (400 Myr). V=6.9. Dai+20 submitted RM measurement.',
        #'HIP67522b, Rizzuto+20 THYME2 validated (b near 0). Claimed is HJ; but no mass. No RM yet either.'
        'NEB 34arcsec W (SG1).  KC19 Cluster might exist. HR diagram decently coeval.',
        'NEB 22arcsec NW (SG1). Group definitely real and KC19 membership seems good. MIT notes that Could be on neighbor',
        'SP2, PP2. Star is far from cluster core but could be part of fuzz-like population. Check Prot-Teff of the fuzz to see if the Kounkel members are real. FIES(1) shows Vrot=160km/s.',
        'SP0: check for Li. PP0: get multicolor zs, gp, ip. It seems to be grazing (boo). CDIPS vetting report looks good (not sure why missed). MIT odd-even is 4sigma --perhaps bc bad initial detrending. Prot is 2.5d for a 6000K star, correct for IC2602.  KC19 IC2602 structure seems clearly wrong because of one merged group at RA/dec 130-55. But this star seems like it could be in the outskirts of the main known cluster. MIT notes crowded field. No 2 min available. ',
        'NEB with eccentric secondary eclipse (SG2). Source is 17" SW. KC19 Membership seems good. MIT notes 9.9 day eccentric system with separations between events of ~4.72; 5.16',
        'SP2, PP2. MIT notes: probably an EB. Clearly not CG18 member. KC19s IC2602 structure seems wrong because of one merged group at RA/dec 130-55. But the HR diagram is pretty clean. So maybe the two populations are actually related. Pop Needs Prot sequence. Star is just at cluster turnoff, and is in the weird non-close population.'
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

    insert_single = 1
    insert_many = 0

    if insert_single:

        insert_single_candidate(
            ticid = '464646604',
            comment = 'ABDor MG member (50 Myr). Identified by George Zhou. Not yet a TOI. 500 or 600 ppm depth. V=8.3, perhaps good for RM.',
            nbhd_rating = 0,
            init_priority = 0,
            pending_photometry_observations = 'PP0.',
            pending_spectroscopic_observations = 'SP2.'
        )

    if insert_many:
        manual_insert_many()
