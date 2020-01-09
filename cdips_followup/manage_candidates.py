"""
Manage candidates.csv.

Each row in this CSV file is a planet candidate entry. This database is used to
track what needs observing, what has been and will be observed, and whether
more observations are needed.

Contents:
    insert_candidate
    query_candidate
    save_candidates_csv_file

Use cases:

* Add candidates given source_id or ticid: use `insert_candidate`.

* Update arbitrary columns given new information: use Libreoffice Calc
  (separated by "|", with "format quoted field as text"). Saving matches the
  pandas read/write formatting. The following columns are changed the most:

    current_priority
    pending_spectroscopic_observations
    pending_photometry_observations
    comment

* View spreadsheet to select viable candidates. Viewing is done in Libreoffice
Calc (separated by "|", with "format quoted field as text"), or else in Google
Spreadsheets. Google Sheets has nice auto-coloring logic, but cannot be used
for updates. (While Libreoffice calc can).

`candidates.csv` has columns:
    source_id, ticid, toi, targetid, reference, name, nbhd_rating,
    init_priority, current_priority, pending_spectroscopic_observations,
    pending_photometry_observations, comment, gaia_ra, gaia_dec, gaia_plx,
    gaia_Gmag, gaia_Bmag, gaia_Rmag, tic_Bmag, tic_Vmag, tic_Jmag, tic_Hmag,
    tic_Kmag, tic_Tmag, tic_teff, tic_logg, tic_rstar, tic_mstar,
    candidate_provenance, insert_time, last_update_time, isretired
"""

######################
# imports and config #
######################

import pandas as pd, numpy as np
import socket, os, csv
from datetime import datetime
from parse import search

from astrobase.services.identifiers import (
    tic_to_gaiadr2, gaiadr2_to_tic
)

from cdips.utils import today_YYYYMMDD
from cdips.utils.catalogs import (
    get_exofop_toi_catalog,
    get_toi_catalog,
    get_cdips_pub_catalog_entry,
    get_tic_star_information,
    get_exofop_toi_catalog_entry,
    get_exofop_ctoi_catalog_entry
)
from cdips_followup.utils import (
    ticid_to_toiid,
    ticid_and_toiid_to_targetid
)

HOMEDIR = os.path.expanduser('~')

if socket.gethostname() in ['ast1607-astro', 'brik']:
    CAND_PATH = os.path.join(
        HOMEDIR,
        'Dropbox/proj/cdips_followup/data/candidate_database/candidates.csv'
    )

if not os.path.exists(CAND_PATH):
    errmsg = 'Need to define CAND_PATH on this machine.'
    raise NotImplementedError(errmsg)


#############################
# general-purpose functions #
#############################

def insert_candidate(
    source_id=None, ticid=None, manual_dict=None
    ):
    """
    Insert a candidate to the candidates.csv database by passing either Gaia
    DR2 source_id or otherwise ticid (string).

    manual_dict (optional): with keys:

        nbhd_rating (0-2, or null),
        init_priority (0-2),
        current_priority (0-2),
        pending_spectroscopic_observations (str, '' if null),
        pending_photometry_observations (str, '' if null),
        comment (str, '' if null)
        candidate_provenance (str, '' if null)
        isretired (0 or 1)
    """

    #
    # Get identifiers (source_id, ticid, toi, targetid).
    #
    validate_source_id_ticid(source_id, ticid)

    if isinstance(source_id, str):
        ticid = gaiadr2_to_tic(source_id)
    elif isinstance(ticid, str):
        source_id = tic_to_gaiadr2(ticid)

    toiid = ticid_to_toiid(ticid)
    if isinstance(toiid, str):
        toiid = toiid.replace('.01','')

    targetid = ticid_and_toiid_to_targetid(ticid, toiid)

    #
    # Get CDIPS & GaiaDR2 catalog information, or else assign nans.
    #

    cdips_r = get_cdips_pub_catalog_entry(source_id)

    iscdipstarget = 1 if isinstance(cdips_r, pd.DataFrame) else 0

    if isinstance(cdips_r, pd.DataFrame):
        assert len(cdips_r) == 1
        cdips_r = cdips_r.iloc[0]

        for col in cdips_r.index:
            if pd.isnull(cdips_r[col]):
                cdips_r[col] = -1

    else:
        cdips_cols = [
            'cluster', 'reference', 'ext_catalog_name', 'ra', 'dec', 'pmra',
            'pmdec', 'parallax', 'phot_g_mean_mag', 'phot_bp_mean_mag',
            'phot_rp_mean_mag', 'k13_name_match', 'unique_cluster_name',
            'how_match', 'not_in_k13', 'comment', 'logt', 'e_logt',
            'logt_provenance'
        ]

        cdips_r = pd.Series({
            'source_id': source_id
        })

        for col in cdips_cols:
            cdips_r[col] = '-1'

    #
    # Get TIC information, or else assign nans
    #

    ticcols = ['ID', 'GAIA', 'Bmag', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'Tmag',
               'Teff', 'logg', 'rad', 'mass']

    tic_r = get_tic_star_information(ticid, desiredcols=ticcols)

    if isinstance(tic_r, pd.DataFrame):
        assert len(tic_r) == 1
        tic_r = tic_r.iloc[0]

        if not tic_r.GAIA == source_id:
            errmsg = (
                'expected tic GAIA ID ({}) to match my GAIA ID ({})'.
                format(tic_r.GAIA, source_id)
            )
            raise AssertionError(errmsg)

        for col in ticcols:
            if pd.isnull(tic_r[col]):
                tic_r[col] = -1

    else:
        tic_r = pd.Series({
            'source_id': source_id
        })

        for col in ticcols:
            tic_r[col] = -1

    #
    # Get the fit information as uploaded to ExoFOP-TESS. By default, use the
    # TOI table information. Otherwise, use the CTOI table.
    #

    plproperties_r = get_exofop_toi_catalog_entry(ticid)
    plkeyd = {
        'rp': 'Planet Radius (R_Earth)',
        'period': 'Period (days)',
        'depth': 'Depth (mmag)'
    }
    if plproperties_r is None:
        plproperties_r = get_exofop_ctoi_catalog_entry(ticid)
        plkeyd = {
            'rp': 'Radius (R_Earth)',
            'period': 'Period (days)',
            'depth': 'Depth mmag'
        }

    if isinstance(plproperties_r, pd.DataFrame):
        assert len(plproperties_r) == 1
        plproperties_r = plproperties_r.iloc[0]

        for col in plproperties_r.index:
            if pd.isnull(plproperties_r[col]):
                plproperties_r[col] = -1

    else:
        plproperties_r = pd.DataFrame(
            {v:-1 for v in plkeyd.values()}, index=[0]
        )

    #
    # Get: nbhd_rating, init_priority, current_priority,
    # pending_spectroscopic_observations, pending_photometry_observations,
    # comment, candidate_provenance, isretired.
    #

    if isinstance(manual_dict, dict):

        # Dictionary was passed containing the manually entries that SHOULD be
        # manually written.

        d = manual_dict
        nbhd_rating = d['nbhd_rating']
        init_priority = d['init_priority']
        current_priority = d['current_priority']
        pending_spectroscopic_observations = d['pending_spectroscopic_observations']
        pending_photometry_observations = d['pending_photometry_observations']
        comment = d['comment']
        candidate_provenance = d['candidate_provenance']
        isretired = d['isretired']

    else:

        # Set reasonable defaults. Raise warning.

        init_priority = 1
        nbhd_rating = init_priority
        current_priority = init_priority
        pending_spectroscopic_observations = ''
        pending_photometry_observations = ''
        comment = ''
        candidate_provenance = 'insert_candidate (w/out manual entries)'
        isretired = 0

        print(
            'WRN! For {}, did not get manual entries for PRIORITY, or COMMENT'.
            format(source_id)
        )

    #
    # Construct and insert the new row.
    #

    new_row = pd.DataFrame({
        'source_id': str(source_id),
        'ticid': str(ticid),
        'toi': str(toiid),
        'targetid': str(targetid),
        'iscdipstarget': iscdipstarget,
        'reference': cdips_r.reference,
        'name': cdips_r.cluster,
        'age': cdips_r.logt,
        'nbhd_rating': nbhd_rating if not pd.isnull(nbhd_rating) else '--',
        'init_priority': init_priority,
        'current_priority': current_priority,
        'pending_spectroscopic_observations': pending_spectroscopic_observations,
        'pending_photometry_observations': pending_photometry_observations,
        'comment': comment,
        'rp': plproperties_r[plkeyd['rp']],
        'period': plproperties_r[plkeyd['period']],
        'depth': plproperties_r[plkeyd['depth']],
        'gaia_ra': cdips_r.ra,
        'gaia_dec': cdips_r.dec,
        'gaia_plx': cdips_r.parallax,
        'gaia_Gmag': cdips_r.phot_g_mean_mag,
        'gaia_Bmag': cdips_r.phot_bp_mean_mag,
        'gaia_Rmag': cdips_r.phot_rp_mean_mag,
        'tic_Bmag': tic_r.Bmag,
        'tic_Vmag': tic_r.Vmag,
        'tic_Jmag': tic_r.Jmag,
        'tic_Hmag': tic_r.Hmag,
        'tic_Kmag': tic_r.Kmag,
        'tic_Tmag': tic_r.Tmag,
        'tic_teff': tic_r.Teff if not pd.isnull(tic_r.Teff) else -1,
        'tic_logg': tic_r.logg,
        'tic_rstar': tic_r.rad,
        'tic_mstar': tic_r.mass,
        'candidate_provenance': candidate_provenance,
        'insert_time': pd.Timestamp.now(),
        'last_update_time': pd.Timestamp.now(),
        'isretired': isretired
    }, index=[0])

    cand_df = pd.read_csv(CAND_PATH, sep='|')

    new_cand_df = pd.concat((cand_df, new_row), sort=False)

    #
    # convert default non-string columns to strings. ditto for float and
    # integer columns.
    #
    strcols = ['source_id', 'ticid', 'targetid', 'nbhd_rating',
               'pending_spectroscopic_observations',
               'pending_photometry_observations', 'comment',
               'candidate_provenance', 'insert_time', 'last_update_time']
    for strcol in strcols:
        new_cand_df[strcol] = new_cand_df[strcol].astype(str)

    floatcols = ['rp', 'period', 'gaia_ra', 'gaia_dec', 'gaia_plx',
                 'gaia_Gmag', 'gaia_Bmag', 'gaia_Rmag', 'tic_Bmag', 'tic_Vmag',
                 'tic_Jmag', 'tic_Hmag', 'tic_Kmag', 'tic_Tmag', 'tic_logg',
                 'tic_rstar', 'tic_mstar']
    for fcol in floatcols:
        new_cand_df[fcol] = new_cand_df[fcol].astype(float)


    intcols = ['iscdipstarget', 'tic_teff']
    for icol in intcols:
        new_cand_df[icol] = new_cand_df[icol].astype(int)

    #
    # replace empty and "nan" strings with "--" string, which LibreOffice calc
    # can handle.
    #
    new_cand_df = new_cand_df.replace(r'^\s*$', "--", regex=True)
    new_cand_df = new_cand_df.replace('nan', "--", regex=True)
    new_cand_df = new_cand_df.replace('None', "--", regex=True)

    save_candidates_csv_file(new_cand_df)


def validate_source_id_ticid(source_id, ticid):

    assert isinstance(source_id, str) or isinstance(ticid, str)

    if isinstance(source_id, str):
        assert not isinstance(ticid, str)

    if isinstance(ticid, str):
        assert not isinstance(source_id, str)


def query_candidate(source_id=None, ticid=None):
    """
    Query candidates.csv for the candidate entry. Use either Gaia DR2
    source_id or otherwise ticid. The NEWEST candidate is used (i.e., presumes
    that any new entries supersede old ones, for the same target).
    """

    validate_source_id_ticid(source_id, ticid)

    df = pd.read_csv(CAND_PATH, sep='|')

    # Get the newest entry for this source.
    if isinstance(source_id, str):
        seldf = df[df.source_id == source_id].iloc[-1]

    if isinstance(ticid, str):
        seldf = df[df.ticid == ticid].iloc[-1]

    return dict(seldf)


def save_candidates_csv_file(cand_df):
    """
    Save function does two things:
        1. If candidates_{YYYYMMDD}.csv doesn't exist, it will make it. Even
        though candidates.csv is version-controlled, this is a simple
        mechanism to allow periodic backups, rate-limited at one day.
        2. Overwrite candidates.csv
    """

    CAND_DAY_PATH = os.path.join(
        os.path.dirname(CAND_PATH),
        'candidates_{}.csv'.format(today_YYYYMMDD())
    )

    colorder = [
        'source_id', 'ticid', 'toi', 'targetid', 'iscdipstarget', 'reference',
        'name', 'age', 'nbhd_rating', 'init_priority', 'current_priority',
        'pending_spectroscopic_observations',
        'pending_photometry_observations', 'comment', 'rp', 'period',
        'gaia_ra', 'gaia_dec', 'gaia_plx', 'gaia_Gmag', 'gaia_Bmag',
        'gaia_Rmag', 'tic_Bmag', 'tic_Vmag', 'tic_Jmag', 'tic_Hmag',
        'tic_Kmag', 'tic_Tmag', 'tic_teff', 'tic_logg', 'tic_rstar',
        'tic_mstar', 'candidate_provenance', 'insert_time', 'last_update_time',
        'isretired'
    ]

    cand_df['insert_time'] = pd.to_datetime(cand_df.insert_time)
    out_df = cand_df[colorder].sort_values(by='insert_time')

    if not os.path.exists(CAND_DAY_PATH):
        out_df.to_csv(
            CAND_DAY_PATH, index=False, sep='|', quotechar='"',
            quoting=csv.QUOTE_NONNUMERIC, float_format='%.2f'
        )
        print('{}: backed up to {}'.
              format(datetime.utcnow(), CAND_DAY_PATH))

    out_df.to_csv(
        CAND_PATH, index=False, sep='|', quotechar='"',
        quoting=csv.QUOTE_NONNUMERIC, float_format='%.2f'
    )
    print('{}: overwrite {} with updated values'.
          format(datetime.utcnow(), CAND_PATH))
