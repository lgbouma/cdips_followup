"""
To insert CDIPS candidates into the "CHIRON planning" spreadsheet,
a particular format is required...

TICID Rp P Teff TIC RA DEC Vmag N_exp
"""

import numpy as np, pandas as pd
from astrobase.services.identifiers import tic_to_gaiadr2
from astrobase.services.mast import tic_objectsearch
from numpy import array as nparr
from astropy.coordinates import ICRS
from astropy import units as u
import json

def given_ticid_get_quinn_row(ticid, disposition='STAR'):

    dr2_source_id = tic_to_gaiadr2(ticid)

    tic_res = tic_objectsearch(ticid)
    res_path = tic_res['cachefname']
    with open(res_path) as f:
        _d = json.load(f)
    assert len(_d['data']) == 1

    tic8_data = _d['data'][0]

    ra = tic8_data['ra']
    dec = tic8_data['dec']
    vmagstr = str(tic8_data['Vmag'])

    c = ICRS(nparr(ra)*u.deg, nparr(dec)*u.deg)
    rahmsstr = c.ra.to_string(u.hour, sep=' ', pad=True)
    decdmsstr = c.dec.to_string(u.degree, sep=' ', pad=True)

    quinnrow = (
        f'{ticid}\t{disposition}\t' +
        f'T{str(ticid).zfill(10)}\t'+
        f'{rahmsstr}\t{decdmsstr}\t'+
        f'{vmagstr}'
    )

    return quinnrow


def main():

    runid = 'toi1937_compstars_ngc2516'
    disposition = 'STAR'

    ticids = [
        "364398101",
        "340455080",
        "267893340",
        "267892931",
        "264873451",
        "300510057",
        "382625643",
        "382507393"
    ]

    for t in ticids:
        r = given_ticid_get_quinn_row(t, disposition=disposition)
        print(r)


if __name__ == "__main__":
    main()
