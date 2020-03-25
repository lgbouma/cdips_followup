"""
Manage candidates.csv.

Each row in this CSV file is a planet candidate entry. This database is used to
track what needs observing, what has been and will be observed, and whether
more observations are needed.

Contents:
    insert_candidate
    query_candidate
    save_candidates_csv_file
    update_candidate_rot_params

Use cases:

* Add candidates given source_id or ticid: use `insert_candidate`.

* Update columns given new information: use Libreoffice Calc (separated by "|",
with "format quoted field as text"). Saving matches the pandas read/write
formatting. The following columns are changed the most:

    current_priority
    pending_spectroscopic_observations
    pending_photometry_observations
    comment

* Update manually gauged rotation parameters for particular candidates, with
update_candidate_rot_params. Note that this cannot be done through localc.

* View spreadsheet to select viable candidates. Viewing is done in Libreoffice
Calc (separated by "|", with "format quoted field as text"), or else in Google
Spreadsheets. Google Sheets has nice auto-coloring logic, but cannot be used
for updates. (While Libreoffice calc can).
"""

######################
# imports and config #
######################

import pandas as pd, numpy as np
import socket, os, csv
from datetime import datetime
from parse import search
from astropy import units as u

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
    source_id=None, ticid=None, manual_dict=None,
    raise_error_if_duplicate=True
    ):
    """
    Insert a candidate to the candidates.csv database by passing either Gaia
    DR2 source_id or otherwise ticid (string).

    Optional Arguments:
    ----------

    manual_dict: dict

        With keys:
        nbhd_rating (0-2, or -1 for null),
        init_priority (0-2),
        current_priority (0-2),
        pending_spectroscopic_observations (str, '--' if null),
        pending_photometry_observations (str, '--' if null),
        comment (str, '--' if null)
        candidate_provenance (str)
        isretired (0 or 1)

    raise_error_if_duplicate: boolean

        If attempting an insert on a source_id that already exists, an error
        will be raised, and the insert will not be performed. If false, a
        warning is raised, and the insert will not be performed.
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
        'rp_unc': 'Planet Radius (R_Earth) err',
        'period': 'Period (days)',
        'depth': 'Depth (mmag)'
    }
    if plproperties_r is None:
        plproperties_r = get_exofop_ctoi_catalog_entry(ticid)
        plkeyd = {
            'rp': 'Radius (R_Earth)',
            'rp_unc': 'Radius (R_Earth) Error',
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
        'rp_unc': plproperties_r[plkeyd['rp_unc']],
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
        'isretired': isretired,
        'disposition': 'PC',
        'rot_quality': '--',
        'Prot': '--',
        'vsini': '--',
        'rot_amp': '--',
        'sig_Prot': '--',
        'Tdur': '--',
        'sig_Tdur': '--',
        'Mp_pred': '--',
        'K_orb': '--',
        'K_RM': '--',
        'K_orb/sig_Prot': '--',
        'K_RM/sig_Tdur': '--'
    }, index=[0])

    cand_df = pd.read_csv(CAND_PATH, sep='|')

    if np.any(cand_df.source_id.astype(str).str.contains(str(source_id))):
        msg = (
            'Found existing candidates.csv entry for {}'.
            format(source_id)
        )
        if raise_error_if_duplicate:
            raise AssertionError('ERR! : '+msg)
        else:
            print('WRN! : '+msg)
            print('WRN! Not doing the insert.')
            return None

    new_cand_df = pd.concat((cand_df, new_row), sort=False)

    new_cand_df = format_candidates_file(new_cand_df)
    save_candidates_csv_file(new_cand_df)


def format_candidates_file(candidates_df):

    #
    # convert default non-string columns to strings. ditto for float and
    # integer columns.
    #
    strcols = ['source_id', 'ticid', 'targetid', 'nbhd_rating',
               'pending_spectroscopic_observations',
               'pending_photometry_observations', 'comment',
               'candidate_provenance', 'insert_time', 'last_update_time',
               'rot_quality', 'Prot', 'vsini', 'rot_amp', 'sig_Prot', 'Tdur',
               'sig_Tdur', 'Mp_pred', 'K_orb', 'K_RM', 'K_orb/sig_Prot',
               'K_RM/sig_Tdur' ]

    for strcol in strcols:
        candidates_df[strcol] = candidates_df[strcol].astype(str)

    floatcols = ['rp', 'rp_unc', 'period', 'gaia_ra', 'gaia_dec', 'gaia_plx',
                 'gaia_Gmag', 'gaia_Bmag', 'gaia_Rmag', 'tic_Bmag', 'tic_Vmag',
                 'tic_Jmag', 'tic_Hmag', 'tic_Kmag', 'tic_Tmag', 'tic_logg',
                 'tic_rstar', 'tic_mstar']
    for fcol in floatcols:
        candidates_df[fcol] = candidates_df[fcol].astype(float)

    intcols = ['iscdipstarget', 'tic_teff']
    for icol in intcols:
        candidates_df[icol] = candidates_df[icol].astype(int)

    #
    # replace empty and "nan" strings with "--" string, which LibreOffice calc
    # can handle.
    #
    candidates_df = candidates_df.replace(r'^\s*$', "--", regex=True)
    candidates_df = candidates_df.replace('nan', "--", regex=True)
    candidates_df = candidates_df.replace('None', "--", regex=True)

    return candidates_df


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
        seldf = df[df.ticid.astype(str) == ticid].iloc[-1]

    return dict(seldf)


def get_live_candidates(condition=None):
    """
    Return list of source_ids from candidates.csv that are "live".

    condition:
        None: filters on current_priority < 2
        'SP0SP1': & not pending_spectroscopic_observations contains 'SP2'
    """

    df = pd.read_csv(CAND_PATH, sep='|')

    sel = df.current_priority < 2

    if condition == 'SP0SP1':
        sel &= (~df.pending_spectroscopic_observations.str.contains('SP2'))

    print('Got {} live candidates meeting condition {}'.
          format(len(df[sel]), condition))

    return list(df[sel].source_id)


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
        'pending_photometry_observations', 'comment', 'rp', 'rp_unc', 'period',
        'gaia_ra', 'gaia_dec', 'gaia_plx', 'gaia_Gmag', 'gaia_Bmag',
        'gaia_Rmag', 'tic_Bmag', 'tic_Vmag', 'tic_Jmag', 'tic_Hmag',
        'tic_Kmag', 'tic_Tmag', 'tic_teff', 'tic_logg', 'tic_rstar',
        'tic_mstar', 'candidate_provenance', 'insert_time', 'last_update_time',
        'isretired', 'disposition', 'rot_quality', 'Prot', 'vsini', 'rot_amp',
        'sig_Prot', 'Tdur', 'sig_Tdur', 'Mp_pred', 'K_orb', 'K_RM',
        'K_orb/sig_Prot', 'K_RM/sig_Tdur'
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


def update_candidate_rot_params(ticid=None, source_id=None, rot_quality='--',
                                Prot='--', vsini='--', rot_amp='--',
                                Mp_pred='--', Tdur='--'):
    """
    Inputs:
        * rot_quality: 0 is good. 1 is so-so. 2 is very uncertain.
        * ticid or source_id: str.
        * At least one of Prot or vsini.
            (E.g., 3.14 and '--', or 42 and None, or 42 and 42)
        * Tdur (astropy quantity)
        * rot_amp: float
        * Mp_pred (astropy quantity).
    """

    if not isinstance(Prot, u.Quantity) and not isinstance(vsini, u.Quantity):
        msg = 'Need at least one of Prot or vsini.'
        raise ValueError(msg)

    if not isinstance(rot_amp, float):
        msg = 'Need rotation amplitude.'
        raise ValueError(msg)

    if not isinstance(Mp_pred, u.Quantity):
        msg = 'Need predicted mass.'
        raise ValueError(msg)

    if not isinstance(Tdur, u.Quantity):
        msg = 'Need transit duration.'
        raise ValueError(msg)

    cdict = query_candidate(source_id=source_id, ticid=ticid)

    # Require finite Rstar, Mstar, P_orb, and Rp for this calculation.
    keys = ['tic_rstar', 'tic_mstar', 'period', 'rp']
    for k in keys:
        if float(cdict[k]) < 0:
            raise ValueError('Got bad {} for {}'.
                             format(k, cdict['ticid']))

    # Calculate derived parameters
    # Potentially 'vsini',
    # definitely 'sig_Prot', 'K_orb', 'K_RM', 'K_orb/sig_Prot', 'K_RM/sig_Prot'

    sini = 1
    Rstar = float(cdict['tic_rstar'])*u.Rsun
    if not isinstance(vsini, u.Quantity):
        vsini = (
            2*np.pi*Rstar * sini / Prot
        ).to(u.km/u.s)

    Prot_is_photometric = True
    if not isinstance(Prot, u.Quantity):
        Prot = (
            2*np.pi*Rstar * sini / vsini
        ).to(u.day)
        Prot_is_photometric = False

    # Spot induced jitter at Prot [in m/s] ~= vsini*rotation amplitude.
    sig_Prot = (rot_amp * vsini).to(u.m/u.s)

    # Spot induced linear trend jitter at Tdur [in m/s]
    sig_Tdur = sig_Prot * (Tdur.to(u.day) / Prot.to(u.day)).cgs.value

    # K_orb: orbital semi-amplitude. Lovis & Fischer 2010, Eq 14.
    Mstar = float(cdict['tic_mstar'])*u.Msun
    P_orb = float(cdict['period'])*u.day
    Rplanet = float(cdict['rp'])*u.Rearth

    K_orb = (
        (28.4329*u.m/u.s)*
        (Mp_pred/(1*u.Mjup))*sini*
        ((Mp_pred + Mstar)/(1*u.Msun))**(-2/3)*
        (P_orb/(1*u.year))**(-1/3)
    ).to(u.m/u.s)

    # K_RM: Rossiter McLaughlin anomaly. Winn 2010, Eq 40.
    b = 0.7 # standard guess.
    depth = ((Rplanet/Rstar).cgs.value)**2
    K_RM = (
        depth * np.sqrt(1 - b**2) * vsini
    ).to(u.m/u.s)

    update_d = {}

    update_d['rot_quality'] = rot_quality,
    if Prot_is_photometric:
        update_d['Prot'] = Prot.value
    else:
        update_d['Prot'] = Prot
    update_d['vsini'] = vsini.value
    update_d['rot_amp'] = rot_amp
    update_d['sig_Prot'] = sig_Prot.value
    update_d['Tdur'] = Tdur.to(u.hour).value
    update_d['sig_Tdur'] = sig_Tdur.to(u.m/u.s).value
    update_d['Mp_pred'] Mp_pred.to(u.Mearth).value
    update_d['K_orb'] = K_orb.to(u.m/u.s).value
    update_d['K_RM'] = K_RM.to(u.m/u.s).value
    update_d['K_orb/sig_Prot'] = (K_orb.to(u.m/u.s)/sig_Prot.to(u.m/u.s)).value
    update_d['K_RM/sig_Tdur'] = (K_RM.to(u.m/u.s)/sig_Tdur.to(u.m/u.s)).value

    # Prepare the row to be updated.
    df = pd.read_csv(CAND_PATH, sep='|')
    inds = df.index.where(df.ticid.astype(str) == ticid)
    update_index = inds[~pd.isnull(inds)][0]

    for k,v in update_d.items():
        if isinstance(v, str):
            cdict[k] = v
        elif isinstance(v, float):
            if k == 'rot_amp':
                cdict[k] = '{:.4f}'.format(v)
            else:
                cdict[k] = '{:.2f}'.format(v)
        elif v is None:
            cdict[k] = '--'
        else:
            raise NotImplementedError

    updated_row = pd.DataFrame(cdict, index=[update_index])

    df.loc[update_index] = updated_row.values[0]

    df = format_candidates_file(df)

    save_candidates_csv_file(df)
