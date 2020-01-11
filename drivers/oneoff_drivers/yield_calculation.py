"""
Estimate how many planets can be detected by CDIPS. Assume the target star
sample is just Cantat-Gaudin+2018.
"""

import os, subprocess, shlex
import multiprocessing as mp
import numpy as np, pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units as u, astropy.constants as const

from numpy import array as nparr

from cdips.utils.catalogs import get_cdips_pub_catalog

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
        denominator = np.sum( 10**(-0.4 * nparr(sdf.Tmag_pred) ) )
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


if __name__ == "__main__":

    run = 1
    test = 0

    if test:
        compute_dilution_fraction('5525188767305211904')

    if run:
        cg18_df = get_cg18_stars_above_cutoff_T_mag()

        compute_dilution_fractions(nparr(cg18_df.source_id.astype(str)))
