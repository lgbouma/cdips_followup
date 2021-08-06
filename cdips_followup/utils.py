"""
Contents:
    ticid_to_toiid
    ticid_and_toiid_to_targetid
    get_cdips_candidates
    given_sourceid_get_radec
    given_sourceid_get_gaiarow
"""
import os
import pandas as pd, numpy as np

from astroquery.gaia import Gaia
from cdips.utils.catalogs import ticid_to_toiid as cdips_ticid_to_toiid
from cdips_followup import __path__

def ticid_to_toiid(tic_id):

    return cdips_ticid_to_toiid(tic_id)


def ticid_and_toiid_to_targetid(tic_id, toi_id):
    """
    If TOI identifier exists, use that. Otherwise, try for TIC identifer. If
    neither was passed, return None.
    """

    if isinstance(tic_id, str) and isinstance(toi_id, str):

        return toi_id + '.01'

    elif isinstance(tic_id, str) and toi_id is None:

        return tic_id + '.01'

    else:

        return None


def get_cdips_candidates():

    datadir = os.path.join(os.path.dirname(__path__[0]), 'data')

    df = pd.read_csv(
        os.path.join(datadir, 'candidate_database/candidates.csv'),
        sep="|"
    )
    df.targetid = df.targetid.astype(str)

    return df


def given_sourceid_get_radec(source_id):

    jobstr = (
        "select top 1 g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag from "
        "gaiadr2.gaia_source as g where g.source_id = {:s}".
        format(source_id)
    )

    job = Gaia.launch_job(jobstr)
    gaia_r = job.get_results()

    if len(gaia_r) != 1:
        raise AssertionError('gaia match failed')

    ra, dec = float(gaia_r['ra']), float(gaia_r['dec'])

    return ra, dec


def given_sourceid_get_gaiarow(source_id, whichgaia='gaiadr2'):
    """
    args:
        source_id: from gaiadr2, gaiaedr3, etc.
    """

    jobstr = (
        "select top 1 g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag from "+
        f"{whichgaia}.gaia_source as g where g.source_id = {source_id:s}"
    )

    job = Gaia.launch_job(jobstr)
    gaia_r = job.get_results()

    if len(gaia_r) != 1:
        raise AssertionError('gaia match failed')

    return gaia_r
