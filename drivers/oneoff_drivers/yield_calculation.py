"""
Estimate how many planets can be detected by CDIPS. Assume the target star
sample is just Cantat-Gaudin+2018.
"""

###########
# imports #
###########

import os, subprocess, shlex
import multiprocessing as mp
import numpy as np, pandas as pd
from glob import glob

from astropy.coordinates import SkyCoord
import astropy.units as u, astropy.constants as const

from numpy import array as nparr

from astrobase.services.identifiers import (
    tic_to_gaiadr2, gaiadr2_to_tic
)

from cdips.utils.catalogs import (
    get_cdips_pub_catalog,
    get_tic_star_information
)

##########
# config #
##########

TICCOLS = ['ID', 'version', 'HIP', 'TYC', 'UCAC', 'TWOMASS', 'SDSS', 'ALLWISE',
           'GAIA', 'APASS', 'KIC', 'objType', 'typeSrc', 'ra', 'dec',
           'POSflag', 'pmRA', 'e_pmRA', 'pmDEC', 'e_pmDEC', 'PMflag', 'plx',
           'e_plx', 'PARflag', 'gallong', 'gallat', 'eclong', 'eclat', 'Bmag',
           'e_Bmag', 'Vmag', 'e_Vmag', 'umag', 'e_umag', 'gmag', 'e_gmag',
           'rmag', 'e_rmag', 'imag', 'e_imag', 'zmag', 'e_zmag', 'Jmag',
           'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'TWOMflag', 'prox',
           'w1mag', 'e_w1mag', 'w2mag', 'e_w2mag', 'w3mag', 'e_w3mag', 'w4mag',
           'e_w4mag', 'GAIAmag', 'e_GAIAmag', 'Tmag', 'e_Tmag', 'TESSflag',
           'SPFlag', 'Teff', 'e_Teff', 'logg', 'e_logg', 'MH', 'e_MH', 'rad',
           'e_rad', 'mass', 'e_mass', 'rho', 'e_rho', 'lumclass', 'lum',
           'e_lum', 'd', 'e_d', 'ebv', 'e_ebv', 'numcont', 'contratio',
           'disposition', 'duplicate_id', 'priority', 'eneg_EBV', 'epos_EBV',
           'EBVflag', 'eneg_Mass', 'epos_Mass', 'eneg_Rad', 'epos_Rad',
           'eneg_rho', 'epos_rho', 'eneg_logg', 'epos_logg', 'eneg_lum',
           'epos_lum', 'eneg_dist', 'epos_dist', 'distflag', 'eneg_Teff',
           'epos_Teff', 'TeffFlag', 'gaiabp', 'e_gaiabp', 'gaiarp', 'e_gaiarp',
           'gaiaqflag', 'starchareFlag', 'VmagFlag', 'BmagFlag', 'splists',
           'e_RA', 'e_Dec', 'RA_orig', 'Dec_orig', 'e_RA_orig', 'e_Dec_orig',
           'raddflag', 'wdflag', 'objID']

TICDTYPE = { 'ID': str, 'version': str, 'HIP': str, 'TYC': str, 'UCAC': str,
            'TWOMASS': str, 'SDSS': str, 'ALLWISE': str, 'GAIA': str, 'APASS':
            str, 'KIC': str, 'objType': str, 'typeSrc': str, 'ra': float,
            'dec': float, 'POSflag': str, 'pmRA': float, 'e_pmRA': float,
            'pmDEC': float, 'e_pmDEC': float, 'PMflag': str, 'plx': float,
            'e_plx': float, 'PARflag': str, 'gallong': float, 'gallat': float,
            'eclong': float, 'eclat': float, 'Bmag': float, 'e_Bmag': float,
            'Vmag': float, 'e_Vmag': float, 'umag': float, 'e_umag': float,
            'gmag': float, 'e_gmag': float, 'rmag': float, 'e_rmag': float,
            'imag': float, 'e_imag': float, 'zmag': float, 'e_zmag': float,
            'Jmag': float, 'e_Jmag': float, 'Hmag': float, 'e_Hmag': float,
            'Kmag': float, 'e_Kmag': float, 'TWOMflag': str, 'prox': float,
            'w1mag': float, 'e_w1mag': float, 'w2mag': float, 'e_w2mag': float,
            'w3mag': float, 'e_w3mag': float, 'w4mag': float, 'e_w4mag': float,
            'GAIAmag': float, 'e_GAIAmag': float, 'Tmag': float, 'e_Tmag':
            float, 'TESSflag': str, 'SPFlag': str, 'Teff': float, 'e_Teff':
            float, 'logg': float, 'e_logg': float, 'MH': float, 'e_MH': float,
            'rad': float, 'e_rad': float, 'mass': float, 'e_mass': float,
            'rho': float, 'e_rho': float, 'lumclass': str, 'lum': float,
            'e_lum': float, 'd': float, 'e_d': float, 'ebv': float, 'e_ebv':
            float, 'numcont': str, 'contratio': float, 'disposition': str,
            'duplicate_id': str, 'priority': float, 'eneg_EBV': float,
            'epos_EBV': float, 'EBVflag': str, 'eneg_Mass': float, 'epos_Mass':
            float, 'eneg_Rad': float, 'epos_Rad': float, 'eneg_rho': float,
            'epos_rho': float, 'eneg_logg': float, 'epos_logg': float,
            'eneg_lum': float, 'epos_lum': float, 'eneg_dist': float,
            'epos_dist': float, 'distflag': str, 'eneg_Teff': float,
            'epos_Teff': float, 'TeffFlag': str, 'gaiabp': float, 'e_gaiabp':
            float, 'gaiarp': float, 'e_gaiarp': float, 'gaiaqflag': str,
            'starchareFlag': str, 'VmagFlag': str, 'BmagFlag': str, 'splists':
            str, 'e_RA': float, 'e_Dec': float, 'RA_orig': float, 'Dec_orig':
            float, 'e_RA_orig': float, 'e_Dec_orig': float, 'raddflag': str,
            'wdflag': str, 'objID': str, }

#############
# functions #
#############

def get_cg18_df():

    outpath = '../../data/cg18_cdips_table1_subset.csv'

    if not os.path.exists(outpath):

        cdips_df = get_cdips_pub_catalog()
        cg18_df = cdips_df[cdips_df.reference.str.contains('CantatGaudin_2018')]

        # Stassun+2019, Eq 1
        Tmag_pred = (cg18_df['phot_g_mean_mag']
                    - 0.00522555 * (cg18_df['phot_bp_mean_mag'] - cg18_df['phot_rp_mean_mag'])**3
                    + 0.0891337 * (cg18_df['phot_bp_mean_mag'] - cg18_df['phot_rp_mean_mag'])**2
                    - 0.633923 * (cg18_df['phot_bp_mean_mag'] - cg18_df['phot_rp_mean_mag'])
                    + 0.0324473)

        cg18_df['Tmag_pred'] = Tmag_pred

        cg18_df.to_csv(outpath, sep=';')

    else:

        cg18_df = pd.read_csv(outpath, sep=';')

    return cg18_df


def get_cg18_stars_above_cutoff_T_mag(cutoff=16):

    cg18_df = get_cg18_df()

    print('{} CG18 stars with T<{}'.format(
        len(cg18_df[cg18_df.Tmag_pred < cutoff]), cutoff
    ))

    return cg18_df[cg18_df.Tmag_pred < cutoff]


def compute_dilution_fractions(source_ids, do_it_correctly=False,
                               naive_serial=False, naive_parallel=True):
    """
    args:
        source_ids : list of Gaia DR2 source_id strings.

    --------------------
    The correct way to estimate flux dilution is to make synthetic images. For
    any target star (G_Rp<16), collect all neighbors (down to say G_Rp=19). You
    need their (ra,dec,G,Bp,Rp). Combine the G,Bp,Rp to get the total counts
    (Tmag) for each star. Use a pointing model to convert (ra,dec)->(x,y).

    To make the image, you need a PSF model. Stassun+2018,2019 used a 2d gaussian,
    because it can be analytically integrated (and even though it can be ~30%
    wrong at the individual target level).

    A * int_-inf^+inf exp(-[x-x0]^2/2sigma^2 ) dx = A sigma sqrt(pi).

    The idea of the calculation is solve an analytic equation like:

    N_counts (from Tmag) = "A sigma sqrt(pi) " = integral of the 2d gaussian PSF

    Then solve for A, the normalization of the 2d gaussian, in terms of sigma
    (the FWHM of the assumed gaussian PSF) and N_counts.

    Then you have a vector of (x,y,A,sigma) for each target (sigma is
    constant). Discretize the scene at the sub-pixel level (e.g., for a 10x10
    TESS pixel grid, maybe make the image as 1000x1000). Add in the stars. Put
    down apertures at the center of the target star. Calculate the fraction!

    --------------------
    The way this function calculates the dilution fraction is a crude
    approximation. It queries neighbors for (ra,dec,G,Bp,Rp). And it posits the
    PSF is a delta function. For a *population-level* guess at the dilution
    fractions, this is probably OK.

    With naive_serial, it computes ~600 per minute, or takes 6 hours for the ~2e5
    CG18 stars.

    With naive_parallel, it computes ~10k per minute (over 50 cores), so it's
    done in like 20 minutes.
    """

    if do_it_correctly:
        errmsg = compute_dilution_fractions.__doc__
        raise NotImplementedError(errmsg)

    if naive_serial:
        for source_id in source_ids:
            compute_dilution_fraction(source_id)

    if naive_parallel:
        nworkers = 52
        maxworkertasks = 1000
        pool = mp.Pool(nworkers, maxtasksperchild=maxworkertasks)

        tasks = list(source_ids)

        # fire up the pool of workers
        results = pool.map(compute_dilution_fraction, tasks)

        # wait for the processes to complete work
        pool.close()
        pool.join()


def get_ra_dec_Tmag_given_sourceid(source_id,
                                   outdir='../../data/dilution_fractions/idqueries/'):

    outfile = os.path.join(outdir, '{}.csv'.format(source_id))

    if not os.path.exists(outfile):

        gaia2readcmd = (
            'gaia2read --id {source_id} --header --extra -o {outfile}'.
            format(source_id=source_id, outfile=outfile)
        )

        gaia2readproc = subprocess.Popen(shlex.split(gaia2readcmd),
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
        g_stdout, g_stderr = gaia2readproc.communicate()

    if not os.path.exists(outfile):
        print(g_stderr)
        print(g_stdout)
        raise ValueError('no gaia cone match for {}'.format(source_id))

    df = pd.read_csv(outfile, delim_whitespace=True)

    Tmag = (
        df['phot_g_mean_mag[20]']
        - 0.00522555 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])**3
        + 0.0891337 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])**2
        - 0.633923 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])
        + 0.0324473
    )

    return float(df['RA[deg][2]']), float(df['Dec[deg][3]']), float(Tmag)



def conequery_given_radec_sourceid(ra, dec, source_id, boxrad=0.03333,
                                   faintrmag=19,
                                   outdir='../../data/dilution_fractions/conequeries'):

    outfile = os.path.join(outdir, '{}.csv'.format(source_id))

    if not os.path.exists(outfile):

        gaia2readcmd = (
            "gaia2read -r {ra:f} -d {dec:f} -s {boxrad:f} --MR {faintrmag:f} "
            "--circ --header --extra -o {outfile}"
        ).format(
            ra=ra,
            dec=dec,
            boxrad=boxrad,
            faintrmag=faintrmag,
            outfile=outfile
        )

        gaia2readproc = subprocess.Popen(shlex.split(gaia2readcmd),
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
        g_stdout, g_stderr = gaia2readproc.communicate()

        if 'stars found' not in g_stdout.decode('utf-8'):
            print(g_stderr)
            print(g_stdout)
            raise ValueError('no gaia cone match for {}'.format(source_id))

    return pd.read_csv(outfile, delim_whitespace=True)


def compute_dilution_fraction(
    source_id,
    aperture_radii=[0.75,1.,1.25,1.5,1.75,2.,2.25,2.5],
    ):

    assert isinstance(source_id, str)

    outpath = (
        '../../data/dilution_fractions/dilutionvalues/{}.csv'.
        format(source_id)
    )

    if os.path.exists(outpath):
        print('found {}, skip'.format(outpath))
        return
    #
    # Cone query down to G_Rp=19, within say 6 pixels = 120 arcseconds =
    # 0.03333 degrees.  Do it w/ gaia2read, b/c it's faster than internet-based
    # queries.
    #

    ra, dec, target_Tmag = get_ra_dec_Tmag_given_sourceid(source_id)
    df = conequery_given_radec_sourceid(ra, dec, source_id)

    #
    # Calculate separations and Tmags.
    #
    c_obj = SkyCoord(ra, dec, unit=(u.deg), frame='icrs')

    c_cone = SkyCoord(nparr(df['RA[deg][2]']), nparr(df['Dec[deg][3]']),
                      unit=(u.deg), frame='icrs')

    cone_seps = c_cone.separation(c_obj).to(u.arcsec).value

    df['sep_arcsec'] = cone_seps

    px_scale = 20.25 # arcsec/px
    df['sep_px'] = cone_seps / px_scale

    Tmag_pred = (df['phot_g_mean_mag[20]']
                - 0.00522555 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])**3
                + 0.0891337 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])**2
                - 0.633923 * (df['phot_bp_mean_mag[25]'] - df['phot_rp_mean_mag[30]'])
                + 0.0324473)

    df['Tmag_pred'] = Tmag_pred

    #
    # Compute dilution fraction for each aperture radius.
    #

    dilutions = []
    nstars = []

    for ap_radius in aperture_radii:

        sdf = df[df.sep_px < ap_radius]

        numerator = 10**(-0.4 * target_Tmag)
        denominator = np.sum( 10**(-0.4 * nparr(sdf[~pd.isnull(sdf.Tmag_pred)].Tmag_pred) ) )
        dilution = numerator/denominator

        dilutions.append(dilution)
        nstars.append(len(sdf))

    #
    # Save dataframe of dilutions and nstars for each target star.
    #
    outdf = pd.DataFrame({
        'ap_radius': aperture_radii,
        'dilution': dilutions,
        'nstars': nstars
    })

    outdf.to_csv(outpath, index=False)
    print('made {}'.format(outpath))


def get_merge_dilution_fractions(
    cg18_df, source_ids, aperture_radii=[0.75,1.,1.25,1.5,1.75,2.,2.25,2.5]
    ):
    """
    add columns of
       dilution_apX.XX
    and
       nstar_apX.XX
    for whatever apertures were used.
    """

    outpath = '../../data/cg18_cdips_table1_subset_with_dilution.csv'

    if not os.path.exists(outpath):

        dilutiond = {}
        nstard = {}

        for ap in aperture_radii:
            dilutiond['dilution_ap{:.2f}'.format(ap)] = []
            nstard['nstar_ap{:.2f}'.format(ap)] = []

        for ix, source_id in enumerate(source_ids):
            print('{}/{}'.format(ix, len(source_ids)))

            inpath = (
                '../../data/dilution_fractions/dilutionvalues/{}.csv'.
                format(source_id)
            )

            dil_df = pd.read_csv(inpath)

            for ap in aperture_radii:
                dilutiond['dilution_ap{:.2f}'.format(ap)] = (
                    float(dil_df[dil_df.ap_radius==ap].dilution)
                )

                nstard['nstar_ap{:.2f}'.format(ap)] = (
                    float(dil_df[dil_df.ap_radius==ap].nstars)
                )

        for ap in aperture_radii:
            dilkey = 'dilution_ap{:.2f}'.format(ap)
            starkey = 'nstar_ap{:.2f}'.format(ap)

            cg18_df[dilkey] = dilutiond[dilkey]
            cg18_df[starkey] = nstard[starkey]

        cg18_df.to_csv(outpath, index=False)

    return pd.read_csv(outpath)


def _get_tic_star_information(ticid, source_id):

    ticcols = ['ID', 'GAIA', 'Vmag', 'Tmag', 'Teff', 'logg', 'rad', 'e_rad',
               'mass', 'ebv', 'e_ebv']
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

    return tic_r

def get_star_info(cg18_df):

    cg18_df.source_id = cg18_df.source_id.astype(str)

    ticdir = '/nfs/phtess2/ar0/TESS/CAT/TIC8'
    ticpaths = glob(os.path.join(ticdir, 'tic_*00*.csv.gz'))

    selcols = ['ID', 'GAIA', 'Vmag', 'Tmag', 'Teff', 'logg', 'rad', 'e_rad',
               'mass', 'ebv', 'e_ebv', 'source_id']

    for ticchunk in ticpaths:

        # set low_memory=False to not chunk the file. force correct data types.
        tic_df = pd.read_csv(ticchunk, names=TICCOLS, low_memory=False,
                             dtype=TICDTYPE)

        tic_df = tic_df[tic_df.Tmag < 16.5]

        mdf = cg18_df.merge(tic_df, how='inner', left_on='source_id',
                            right_on='GAIA')

        outdf = mdf[selcols]

        outstr = os.path.basename(ticchunk).replace('.csv.gz','_xmatch.csv')
        outpath = '../../data/cg18_ticv8_xmatch/{}'.format(outstr)
        outdf.to_csv(outpath, index=False)
        print('made {}'.format(outpath))


if __name__ == "__main__":

    calc_dilution = 0
    get_merge_dilution = 1
    merge_tic = 1

    cg18_df = get_cg18_stars_above_cutoff_T_mag()

    if calc_dilution:
        compute_dilution_fractions(
            nparr(cg18_df.source_id.astype(str))
        )

    if get_merge_dilution:
        cg18_df = get_merge_dilution_fractions(
            cg18_df, nparr(cg18_df.source_id.astype(str))
        )

    if merge_tic:
        df = get_star_info(cg18_df)
