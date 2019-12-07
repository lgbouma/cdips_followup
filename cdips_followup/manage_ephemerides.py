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

import pandas as pd
import socket
from datetime import datetime
from parse import search

if socket.gethostname() == 'ast1607-astro':
    EPHEM_PATH = '/Users/luke/Dropbox/proj/cdips_followup/data/ephemerides/ephemerides.csv'
elif socket.gethostname() == 'brik':
    EPHEM_PATH = '/home/luke/Dropbox/proj/cdips_followup/data/ephemerides/ephemerides.csv'
else:
    errmsg = 'Need to define EPHEM_PATH on this machine.'
    raise NotImplementedError(errmsg)


#############################
# general-purpose functions #
#############################

def insert_ephemeris(ephemsourcefile, ephemeris_type=None):
    """
    ephemsourcefile: path to source file containing the ephemeris information.

    ephemeris_type: type of ephemeris source file. These could be Hartman's
    text files, cdips-pipeline's pickles, ExoFOP-TESS CSV files, or similar.
    """

    valid_types = ['cdips-pipeline', 'ExoFOP-TESS', 'hartman_updatephem']

    if ephemeris_type not in valid_types:
        errmsg = (
            'got {} type. expected one of valid_types.'.
            format(ephemeris_type)
        )
        raise ValueError(errmsg)

    #
    # Retrieve dict with period, epoch, duration (and uncertainties).
    #
    if ephemeris_type == 'hartman_updatephem':
        ephem_dict = read_hartman_updateephem_file(ephemsourcefile)

    elif ephemeris_type == 'cdips-pipeline':
        raise NotImplementedError

    elif ephemeris_type == 'ExoFOP-TESS':
        raise NotImplementedError

    #
    # Construct and insert the new row.
    #
    ephem_df = pd.read_csv(EPHEM_PATH)

    import IPython; IPython.embed()
    assert 0
    save_ephemerides_csv_file(ephem_df)


def get_ephemeris(source_id=None, ticid=None):
    """
    Query ephemerides.csv for the ephemeris entry. Use either Gaia DR2
    source_id or otherwise ticid.
    """

    assert isinstance(source_id, str) or isinstance(ticid, str)

    if isinstance(source_id, str):
        assert not isinstance(ticid, str)

    if isinstance(ticid, str):
        assert not isinstance(source_id, str)

    pass


def save_ephemerides_csv_file(ephem_df):

    ephem_df.sort_values(by='insert_time').to_csv(
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

    if verbose:
        print(epoch, period, dur)

    epoch = float(search('{} - Epoch: {} +/- {}', epoch)[1].strip())
    period = float(search('{} - Period: {} +/- {}', period)[1].strip())
    dur = float(search('{} - Transit duration: {} +/- {}', dur)[1].strip())

    #TODO get unc
    import IPython; IPython.embed()

    if verbose:
        print(epoch, period, dur)

    ephem_dict = {
        'period':period,
        'epoch':epoch,
        'duration':dur
    }

    return ephem_dict
