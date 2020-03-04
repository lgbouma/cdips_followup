"""
Manage ephemerides.csv.

Each row in this CSV file is an ephemeris entry. Each target will typically
have multiple ephemeris entries, as knowledge improves with photometric
followup.

Tools:
* Insert new ephemerides.
    * From CDIPS pipeline.
    * From ExoFOP-TESS.
    * From Joel Hartman's .updateephem.txt files.
* Get latest ephemeris given target source_id, or ticid.
* Export to google sheets (?).
* Validate the spreadsheet, and save it.
"""

######################
# imports and config #
######################

import pandas as pd, numpy as np
import socket, os
from datetime import datetime
from parse import search

from astrobase.services.identifiers import tic_to_gaiadr2

from cdips.utils import today_YYYYMMDD
from cdips.utils.catalogs import (
    get_exofop_toi_catalog,
    get_exofop_ctoi_catalog
)

from cdips_followup.manage_candidates import validate_source_id_ticid

HOMEDIR = os.path.expanduser('~')
EPHEM_PATH = os.path.join(
    HOMEDIR,
    'Dropbox/proj/cdips_followup/data/ephemerides/ephemerides.csv'
)

if not os.path.exists(EPHEM_PATH):
    errmsg = 'Need to define EPHEM_PATH on this machine.'
    raise NotImplementedError(errmsg)


#############################
# general-purpose functions #
#############################

def insert_ephemeris(ephemeris_type=None, targetid=None, ephemsourcefile=None):
    """
    ephemeris_type: type of ephemeris source file. These could be Hartman's
    text files, cdips-pipeline, ExoFOP-TESS CSV files, or MIT QLP toi-plus
    catalog:  ['cdipspipeline', 'exofoptess_toi', 'hartmanupdate', 'toiplus']

    targetid: if type is MITQLP toi-plus, then e.g., 'TIC308538095.01' or
    '451.01'.

    ephemsourcefile: if type is Hartman's text files, then path to source file
    containing the ephemeris information.
    """

    valid_types = ['cdipspipeline', 'exofoptess_toi',
                   'exofoptess_ctoi','hartmanupdate', 'toiplus']

    if ephemeris_type not in valid_types:
        errmsg = (
            'got {} type. expected one of valid_types.'.
            format(ephemeris_type)
        )
        raise ValueError(errmsg)
    if ephemeris_type =='hartmanupdate':
        assert isinstance(ephemsourcefile, str)
    if ephemeris_type =='toiplus':
        assert isinstance(targetid, str)

    #
    # Retrieve dict with period, epoch, duration (and uncertainties).
    #
    if ephemeris_type == 'hartmanupdate':
        ephem_dict = read_hartman_updateephem_file(ephemsourcefile)

    elif ephemeris_type == 'exofoptess_toi':
        wrnmsg = (
            'WRN! EPHEMERIS FROM EXOFOPTESS TOI CATALOG. ENSURE LATEST '
            'VERSION HAS BEEN DOWNLOADED'
        )
        print(wrnmsg)
        ephem_dict = read_exofoptess_toi_ephem(targetid)

    elif ephemeris_type == 'exofoptess_ctoi':
        wrnmsg = (
            'WRN! EPHEMERIS FROM EXOFOPTESS CTOI CATALOG. ENSURE LATEST '
            'VERSION HAS BEEN DOWNLOADED'
        )
        print(wrnmsg)
        ephem_dict = read_exofoptess_ctoi_ephem(targetid)

    elif ephemeris_type == 'cdipspipeline':
        df = pd.read_csv(ephemsourcefile, sep="|")
        sdf = df[df.target == targetid]
        selcol = ['period', 'epoch', 'duration', 'period_unc', 'epoch_unc',
                  'duration_unc', 'depth', 'depth_unc']
        e_dict = sdf[selcol].to_dict('index')
        keys = ['period_val', 'epoch_val', 'duration_val', 'period_unc',
                'epoch_unc', 'duration_unc', 'depth_val', 'depth_unc']
        ephem_dict = {}
        for k, c in zip(keys, selcol):
            ephem_dict[k] = e_dict[list(e_dict.keys())[0]][c]

    else:
        raise NotImplementedError

    #
    # Construct and insert the new row.
    #

    if ephemeris_type == 'hartmanupdate':
        if 'TIC' not in ephemsourcefile:
            errmsg = (
                'assumed {} had TIC*.updatepehem.txt'.
                format(ephemsourcefile)
            )
            raise ValueError(errmsg)
        # assumes ../data/updated_ephemerides/20191030/TIC308538095.updateephem.txt
        ticid = os.path.basename(ephemsourcefile).split('.')[0].lstrip('TIC')
        source_id = tic_to_gaiadr2(ticid)
        targetid = 'TIC{}.01'.format(ticid)
        ephemeris_origin = os.path.abspath(ephemsourcefile)

    elif ephemeris_type == 'exofoptess_toi':

        toidf = get_exofop_toi_catalog()
        sel = toidf['TOI'].astype(str) == targetid
        targetrow = toidf[sel]

        ticid = str(targetrow['TIC ID'].iloc[0])
        source_id = tic_to_gaiadr2(ticid)
        targetid = targetid
        ephemeris_origin = get_exofop_toi_catalog(returnpath=True)

    elif ephemeris_type == 'exofoptess_ctoi':

        ctoidf = get_exofop_ctoi_catalog()
        sel = ctoidf['CTOI'].astype(str) == targetid.replace('TIC','')
        targetrow = ctoidf[sel]

        ticid = str(targetrow['TIC ID'].iloc[0])
        source_id = tic_to_gaiadr2(ticid)
        targetid = targetid
        ephemeris_origin = get_exofop_ctoi_catalog(returnpath=True)

    elif ephemeris_type == 'cdipspipeline':
        if not targetid.startswith('TIC'):
            raise NotImplementedError
        ticid = targetid.replace('TIC','').replace('.01','')
        source_id = str(sdf.source_id.iloc[0])
        targetid = targetid
        ephemeris_origin = ephemsourcefile

    else:
        raise NotImplementedError

    new_row = pd.DataFrame({
        'source_id': source_id,
        'ticid': ticid,
        'targetid': targetid,
        'insert_time': pd.Timestamp.now(),
        'period': ephem_dict['period_val'],
        'period_unc': ephem_dict['period_unc'],
        'epoch': ephem_dict['epoch_val'],
        'epoch_unc': ephem_dict['epoch_unc'],
        'duration': ephem_dict['duration_val'],
        'duration_unc': ephem_dict['duration_unc'],
        'depth': ephem_dict['depth_val'],
        'depth_unc': ephem_dict['depth_unc'],
        'ephemeris_origin': ephemeris_origin
    }, index=[0])

    ephem_df = pd.read_csv(EPHEM_PATH)

    new_ephem_df = pd.concat((ephem_df, new_row), sort=False)

    save_ephemerides_csv_file(new_ephem_df)


def query_ephemeris(source_id=None, ticid=None):
    """
    Query ephemerides.csv for the ephemeris entry. Use either Gaia DR2
    source_id or otherwise ticid. The NEWEST ephemeris is used (i.e., presumes
    that any new entries supersede old ones, for the same target).
    """

    validate_source_id_ticid(source_id, ticid)

    df = pd.read_csv(EPHEM_PATH)

    # Get the newest entry for this source.
    if isinstance(source_id, str):
        try:
            seldf = df[df.source_id == source_id].iloc[-1]
        except IndexError as e:
            msg = 'Failed to get ephemeris for {}'
            print(msg)
            raise IndexError(msg)

    if isinstance(ticid, str):
        df.ticid = df.ticid.astype(str)
        seldf = df[df.ticid == ticid].iloc[-1]

    return dict(seldf)


def save_ephemerides_csv_file(ephem_df):
    """
    Save function does two things:
        1. If ephemerides_{YYYYMMDD}.csv doesn't exist, it will make it. Even
        though ephemerides.csv is version-controlled, this is a simple
        mechanism to allow periodic backups, rate-limited at one day.
        2. Overwrite ephemerides.csv
    """

    EPHEM_DAY_PATH = os.path.join(
        os.path.dirname(EPHEM_PATH),
        'ephemerides_{}.csv'.format(today_YYYYMMDD())
    )
    colorder = ['source_id', 'ticid', 'targetid', 'insert_time', 'period',
                'period_unc', 'epoch', 'epoch_unc', 'depth', 'depth_unc',
                'duration', 'duration_unc', 'ephemeris_origin']

    ephem_df['insert_time'] = pd.to_datetime(ephem_df.insert_time)
    out_df = ephem_df[colorder].sort_values(by=['insert_time', 'source_id'])

    if not os.path.exists(EPHEM_DAY_PATH):
        out_df.to_csv(
            EPHEM_DAY_PATH, index=False, float_format='%.8f'
        )
        print('{}: backed up to {}'.
              format(datetime.utcnow(), EPHEM_DAY_PATH))

    out_df.to_csv(
        EPHEM_PATH, index=False, float_format='%.8f'
    )
    print('{}: overwrite {} with updated values'.
          format(datetime.utcnow(), EPHEM_PATH))


##############################
# dedicated format functions #
##############################

def read_hartman_updateephem_file(updateephempath, verbose=0):
    # read in the new ephemeris provided by Joel Hartman

    with open(updateephempath, 'r') as f:
        lines = f.readlines()

    epoch = [l for l in lines if '- Epoch: ' in l][0]
    period = [l for l in lines if '- Period: ' in l][0]
    dur = [l for l in lines if '- Transit duration: ' in l][0]
    rprs = [l for l in lines if '- Rp/R*: ' in l][0]

    if verbose:
        print(epoch, period, dur)

    epoch_val = float(search('{} - Epoch: {} +/- {} (20', epoch)[1].strip())
    epoch_unc = float(search('{} - Epoch: {} +/- {} (20', epoch)[2].strip())

    period_val = float(search('{} - Period: {} +/- {}\n', period)[1].strip())
    period_unc = float(search('{} - Period: {} +/- {}\n', period)[2].strip())

    dur_val = float(search('{} - Transit duration: {} +/- {} hours', dur)[1].strip())
    dur_unc = float(search('{} - Transit duration: {} +/- {} hours', dur)[2].strip())

    #  - Rp/R*:  0.18033 +-  0.00274
    rprs_val = float(search('{} - Rp/R*: {} +- {}\n', rprs)[1].strip())
    rprs_unc = float(search('{} - Rp/R*: {} +- {}\n', rprs)[2].strip())

    depth_val = rprs_val**2
    depth_unc = np.abs(
        2 * rprs_unc / rprs_val
    )

    if verbose:
        print(epoch_val, epoch_unc, period_val, period_unc, dur_val, dur_unc)

    ephem_dict = {
        'period_val':period_val,
        'epoch_val':epoch_val,
        'duration_val':dur_val,
        'period_unc':period_unc,
        'epoch_unc':epoch_unc,
        'duration_unc':dur_unc,
        'depth_val':depth_val*1e6,
        'depth_unc':depth_unc*1e6
    }

    return ephem_dict


def read_exofoptess_toi_ephem(targetid):
    # Get TOI ephemerides

    toidf = get_exofop_toi_catalog()

    if targetid.startswith('TIC'):
        raise NotImplementedError

    else:
        # targetid is Full TOI ID
        sel = toidf['TOI'].astype(str) == targetid
        targetrow = toidf[sel]

        if not len(targetrow) == 1:
            raise AssertionError(
                'failed to get ID match for {}'.format(targetid)
            )

    keymatchdict = {
        'period_val': 'Period (days)',
        'epoch_val': 'Epoch (BJD)',
        'duration_val': 'Duration (hours)',
        'period_unc': 'Period (days) err',
        'epoch_unc': 'Epoch (BJD) err',
        'duration_unc': 'Duration (hours) err',
        'depth_val': 'Depth (ppm)',
        'depth_unc': 'Depth (ppm) err'
    }

    ephem_dict  = {}
    for k,v in keymatchdict.items():
        result = np.float64(targetrow[v])
        ephem_dict[k] = result

    return ephem_dict


def read_exofoptess_ctoi_ephem(targetid):
    # Get TOI ephemerides

    ctoidf = get_exofop_ctoi_catalog()

    if not targetid.startswith('TIC'):
        raise NotImplementedError

    else:

        ticid = targetid.replace('TIC','')

        # stripped targetid is CTOI ID
        sel = ctoidf['CTOI'].astype(str) == ticid
        targetrow = ctoidf[sel]

        if not len(targetrow) == 1:
            raise AssertionError(
                'failed to get ID match for {}'.format(targetid)
            )

    keymatchdict = {
        'period_val': 'Period (days)',
        'epoch_val': 'Midpoint (BJD)',
        'duration_val': 'Duration (hrs)',
        'period_unc': 'Period (days) Error',
        'epoch_unc': 'Midpoint err',
        'duration_unc': 'Duration (hrs) Error',
        'depth_val': 'Depth ppm',
        'depth_unc': 'Depth ppm Error'
    }

    ephem_dict  = {}
    for k,v in keymatchdict.items():
        result = np.float64(targetrow[v])
        ephem_dict[k] = result

    return ephem_dict





def read_mitqlp_ephem(targetid):
    # NOTE: WIP
    raise NotImplementedError

    toidf = get_exofop_toi_catalog()

    if targetid.startswith('TIC'):
        raise NotImplementedError

    else:
        # targetid is Full TOI ID
        sel = toidf['Full TOI ID'].astype(str) == targetid
        targetrow = toidf[sel]

        if not len(targetrow) == 1:
            raise AssertionError(
                'failed to get ID match for {}'.format(targetid)
            )

    keymatchdict = {
        'period_val': 'Orbital Period Value',
        'epoch_val': 'Orbital Epoch Value',
        'duration_val': 'Transit Duration Value',
        'period_unc': 'Orbital Period Error',
        'epoch_unc': 'Orbital Epoch Error',
        'duration_unc': 'Transit Duration Error',
        'depth_val': 'Transit Depth Value',
        'depth_unc': 'Transit Depth Error'
    }

    ephem_dict  = {}
    for k,v in keymatchdict.items():

        if k == 'epoch_val':
            # BTJD := BJD - 2457000
            result = 2457000 + np.float64(targetrow[v])
        else:
            result = np.float64(targetrow[v])

        ephem_dict[k] = result

    import IPython; IPython.embed()

    return ephem_dict
