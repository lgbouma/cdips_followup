import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from tessmaps.get_time_on_silicon import (
        given_cameras_get_stars_on_silicon as gcgss
)
from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from numpy import array as nparr

import cdips as cd

import os, json

datadir = os.path.join(os.path.dirname(cd.__path__[0]),
                       'data/skymap_data')

from cdips.utils.catalogs import (
    get_cdips_pub_catalog
)

def get_n_observations(dirnfile, outpath, merged=False,
                       withgaps=True,
                       aligncelestial=False):

    #targetdf = pd.read_csv('../../data/cg18_cdips_table1_subset.csv', sep=';')
    targetdf = get_cdips_pub_catalog()

    # targetdf = targetdf.sample(n=2300)

    ras = nparr(targetdf.ra)*u.deg
    decs = nparr(targetdf.dec)*u.deg

    coords = SkyCoord(ra=ras, dec=decs, frame='icrs')

    if merged:
        df_pri = pd.read_csv(os.path.join(
            datadir,'primary_mission_truenorth.csv'), sep=';')
        df_ext = pd.read_csv(dirnfile, sep=';')
        df = pd.concat([df_pri, df_ext])

    else:
        df = pd.read_csv(dirnfile, sep=';')

    lats = nparr([
        nparr(df['cam1_elat']),
        nparr(df['cam2_elat']),
        nparr(df['cam3_elat']),
        nparr(df['cam4_elat'])]).T
    lons = nparr([
        nparr(df['cam1_elon']),
        nparr(df['cam2_elon']),
        nparr(df['cam3_elon']),
        nparr(df['cam4_elon'])]).T

    cam_directions = []
    for lat, lon in zip(lats, lons):

        c1lat,c2lat,c3lat,c4lat = lat[0],lat[1],lat[2],lat[3]
        c1lon,c2lon,c3lon,c4lon = lon[0],lon[1],lon[2],lon[3]

        this_cam_dirn = [(c1lat, c1lon),
                         (c2lat, c2lon),
                         (c3lat, c3lon),
                         (c4lat, c4lon)]

        cam_directions.append(this_cam_dirn)

    df['camdirection'] = cam_directions

    n_observations = np.zeros_like(coords)

    for ix, row in df.iterrows():

        print(row['start'])
        cam_direction = row['camdirection']

        onchip = gcgss(coords, cam_direction, verbose=False, withgaps=withgaps,
                       aligncelestial=aligncelestial)

        n_observations += onchip

    outdf = pd.DataFrame({'ra':coords.ra.value,
                          'dec':coords.dec.value,
                          'elat':coords.barycentrictrueecliptic.lat.value,
                          'elon':coords.barycentrictrueecliptic.lon.value,
                          'n_observations': n_observations })
    outdf[['ra','dec','elon','elat','n_observations']].to_csv(
        outpath, index=False, sep=';')
    print('saved {}'.format(outpath))


if __name__=="__main__":

    dirnfile = os.path.join( datadir,'primary_mission_truenorth.csv')
    outpath = '../../results/yield_calculation/cdips_v0.4_n_lc_expected.csv'

    if not os.path.exists(outpath):
        get_n_observations(dirnfile, outpath, merged=False, withgaps=True,
                           aligncelestial=False)

    else:
        df = pd.read_csv(outpath, sep=';')

    print('N cdips targets: {}'.format(len(df)))

    print('N southern cdips targets: {}'.format(len(df[df.elat < 0])))

    print('N southern cdips targets on silicon: {}'.
          format(len(df[(df.elat < 0) & (df.n_observations > 0)])))

