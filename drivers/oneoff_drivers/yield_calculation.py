"""
Estimate how many planets can be detected by CDIPS. Assume the target star
sample is just Cantat-Gaudin+2018.
"""

import os
import numpy as np, pandas as pd

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


def compute_dilution_fractions(source_ids, do_it_correctly=False):
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
    """

    if do_it_correctly:
        errmsg = compute_dilution_fractions.__doc__
        raise NotImplementedError(errmsg)


def compute_dilution_fraction(
    source_id,
    aperture_radii=[0.75,1.,1.25,1.5,1.75,2.,2.25,2.5]
    ):

    #
    # Cone query down to G_Rp=19, within say 8 pixels = 160 arcseconds.
    #
    pass





if __name__ == "__main__":

    cg18_df = get_cg18_stars_above_cutoff_T_mag()
