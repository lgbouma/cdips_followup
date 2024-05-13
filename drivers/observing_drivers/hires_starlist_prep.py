import os, json
import numpy as np, pandas as pd
from astropy.coordinates import ICRS, SkyCoord
from astropy import units as u

from astrobase.services.identifiers import tic_to_gaiadr2
from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data
from astrobase.services.mast import tic_xmatch
from exposure import exposure_time

def get_log_row(star_id, coord, vmag, exptime, counts):

    ra_hms = coord.split(' ')[0]
    dec_dms = coord.split(' ')[1]
    try:
        exptime = int(min([exptime, 1800]))
    except ValueError:
        exptime = 1800
    max_exptime = int(min([np.round(exptime+300,-1), 1800]))
    ra_str = ra_hms.replace('h',' ').replace('m',' ').replace('s','')
    dec_str = dec_dms.replace('d',' ').replace('m',' ').replace('s','')

    s_str = str(int(np.round(float(dec_str.split(" ")[-1])))).zfill(2)
    _dec_str = f'{dec_str.split(" ")[0]} {dec_str.split(" ")[1]} {s_str}'

    pstr = 'p1'

    print(f"{star_id} {ra_str} {_dec_str}")

    row_oldfmt = (
        f"{star_id} "
        f"{ra_str} "
        f"{_dec_str} "
        "2000 "
        ""
        f"vmag={vmag:.1f}  "
        f"{exptime}/{max_exptime} "
        f"{int(counts)}k "
        "C2 "
        "1x "
        "out "
        f"{pstr} "
        "LB  "
        "\n"
    )

    N_obs_full_semester = 1
    row_newfmt = (
        f"{N_obs_full_semester},"
        f"{star_id},"
        f"{ra_hms.split('h')[0]},"
        f"{ra_hms.split('h')[1].split('m')[0]},"
        f"{ra_hms.split('h')[1].split('m')[1].replace('s','')},"
        f"{dec_dms.split('d')[0]},"
        f"{dec_dms.split('d')[1].split('m')[0]},"
        f"{dec_dms.split('d')[1].split('m')[1].replace('s','')},"
        "2000,"
        "vmag=,"
        f"{vmag:.1f},"
        f"{exptime},"
        f"{max_exptime},"
        f"{int(counts)},"
        "C2,"
        "1,"
        "out,"
        f"{pstr},"
        "LB,"
        "C245\n"
        #"Single-shot recon\n"
    )

    return row_oldfmt, row_newfmt

def main():
    # given TIC ID's, get HIRES info

    datestr = '20240516'
    jobstr = '24A_kepler'

    targetdir = './targetlists/'
    targetpath = os.path.join(
        targetdir,
        f'{datestr}_kepler_targets.csv'
    )
    df = pd.read_csv(targetpath)

    # NOTE: THIS IS THE MANUAL SELECTION STEP.
    #df = df[df.hires_priority == 2]

    star_ids, coords, vmags, exptimes, counts = [], [], [], [], []

    for ticid in df.ticid:

        print(42*'-')
        print(f'Starting TIC {ticid}...')

        star_ids.append(f"TIC{ticid}")

        dr2_source_id = tic_to_gaiadr2(str(ticid))

        gdf = given_source_ids_get_gaia_data(
            np.array([dr2_source_id]).astype(np.int64),
            str(dr2_source_id), n_max=1, overwrite=False
        )

        assert len(gdf) == 1

        radius = 4
        if dr2_source_id == '1450067137649449728':
            radius = 10

        xmatches = tic_xmatch(
            np.array(gdf.ra), np.array(gdf.dec), radius_arcsec=radius
        )
        if xmatches is None:
            import IPython; IPython.embed()
        with open(xmatches['cachefname'], 'r') as f:
            d = json.load(f)
        tdf = pd.DataFrame(d['data'])

        if len(tdf) > 1:
            print(42*'!')
            print(f'WRN! TIC{ticid} has multiple matches.  Brightest selected...')
            print(tdf[['MatchID','GAIAmag']])
            temp = tdf.sort_values(by='GAIAmag', ascending=True)
            tdf = temp.head(n=1)

        assert len(tdf) == 1

        c = SkyCoord(ra=np.array(gdf.ra)*u.deg,
                     dec=np.array(gdf.dec)*u.deg,
                     frame='icrs')

        coords.append(c.to_string('hmsdms', precision=1)[0])

        if tdf.Vmag.iloc[0] is not None:
            vmag = float(tdf.Vmag)
        else:
            # Eq 3 Stassun+2019, which is from Evans et al., 2018. DR2.
            # scatter 0.046mag
            assert len(gdf) == 1
            vmag = float(
                gdf.phot_g_mean_mag +
                0.01760 + 0.006860*gdf.bp_rp +
                0.1732*(gdf.bp_rp)**2
            )
        vmags.append(vmag)

        if vmag <= 13:
            count = 20
        elif vmag <= 15:
            count = 10
        else:
            count = 5
        counts.append(count)

        exptime_sec = exposure_time(vmag, count, iod=False)
        exptimes.append(np.round(exptime_sec, -1))

        print(f'Ending TIC {ticid}...')
        print(42*'-')


    oldlines, newlines = [], []
    for star_id, coord, vmag, exptime, count in zip(
        star_ids, coords, vmags, exptimes, counts
    ):

        print(star_id)

        oldl, newl = get_log_row(star_id, coord, vmag, exptime, count)

        oldlines.append(oldl)
        newlines.append(newl)

    outpath = f'scripts/{datestr}_{jobstr}_script_newfmt.txt'
    with open(outpath, 'w') as f:
        f.writelines(newlines)

    outpath = f'scripts/{datestr}_{jobstr}_script_oldfmt.txt'
    with open(outpath, 'w') as f:
        f.writelines(oldlines)

    print(f"Made {outpath}")

if __name__ == "__main__":
    main()
