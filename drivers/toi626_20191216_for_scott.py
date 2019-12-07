import os
import numpy as np, matplotlib.pyplot as plt, pandas as pd

from astropy.io.votable import from_table, writeto, parse
from astropy.coordinates import SkyCoord
from astropy import units as u

datadir = '../data/'

def given_votable_get_df(votablepath, assert_equal='source_id'):

    vot = parse(votablepath)
    tab = vot.get_first_table().to_table()
    df = tab.to_pandas()

    if isinstance(assert_equal, str):
        np.testing.assert_array_equal(tab[assert_equal], df[assert_equal])

    return df

def main(extra_overplot=0):

    datapath = os.path.join(datadir,
                            'group449_xmatch-result.vot.gz')

    df = given_votable_get_df(datapath, assert_equal='source_id')

    import IPython; IPython.embed()

if __name__ == "__main__":
    main()
