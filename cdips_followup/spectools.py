"""
READ:
    read_nextgen: Read NextGen Model Atmosphere Grid from Hauschildt et al 1999.
    read_veloce: Read Veloce FITS file.
    read_fies: Read FIES FITS file.
    read_tres: Read TRES FITS file.
    read_hires: Read HIRES FITS file.
    read_pfs: Read PFS IDL SAV file.
    read_feros: Read FEROS FITS file.
    read_galah: Read GALAH (DR3) spectrum FITS file.
        read_galah_given_sobject_id
    read_gaiaeso: Read GAIA-ESO (DR4) spectrum FITS file. (GIRAFFE and UVES)
    read_winered: Read WINERED 1d FITS file.
    read_dbsp: Read DBSP 1d FITS file, extracted by PypeIt.

SPECMATCH WRAPPERS:
    specmatch_analyze: shift+cross-correlate to get vsini, Rstar, FeH w/ SME.
    specmatchsyn_analyze: ditto, but using specmatch-synth.

CALCULATE:
    fit_continuum: fits a quadratic polynomial across the order.
    get_Li_6708_EW: measure Li EW.
    get_line_EW: measure EW for an arbitrary line.
    measure_veloce_vsini: wraps below.
        measure_vsini: compares target spectrum to broadend NextGen spectra
    given_deltawvlen_get_vsys: convert Δλ to velocity.
    air_to_vac: given air wavelengths, get vacuum wavelengths
    get_naive_rv: order-by-order CCF synthetic template match RV

VISUALIZE:
    viz_1d_spectrum: flux vs wavelength (with select lins underplotted).
        plot_orders: wrapper to viz_1d_spectrum
        inspect_pfs: DEPRECATED
    plot_stack_comparison: visual comparison of spectra in a stack.
    specmatch_viz_compare: SME comparison of target and HIRES library spectra.
    plot_spec_vs_dwarf_library: SME ditto.

WIP:
    get_Ca_HK_emission: in theory, for CaHK emission line widths.
"""
import os, re
from datetime import datetime
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
from matplotlib.transforms import blended_transform_factory
from os.path import join
from typing import Tuple

from astropy.io import fits
from astropy import units as u, constants as const
from astropy.modeling.polynomial import Chebyshev1D
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import StdDevUncertainty
from astropy.modeling.fitting import NonFiniteValueError

from scipy.io import readsav
from scipy.interpolate import interp1d
import scipy.ndimage as nd

from copy import deepcopy

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum, fit_lines
from specutils.analysis import equivalent_width, centroid
from specutils.manipulation.utils import excise_regions
from specutils.analysis import correlation

from astropy.modeling import models

from aesthetic.plot import savefig, format_ax, set_style

from cdips_followup import __path__
from cdips_followup.paths import SPECDIR

from cdips.utils.lcutils import p2p_rms

from specutils.manipulation import (
    FluxConservingResampler, LinearInterpolatedResampler,
    SplineInterpolatedResampler
)
fluxcon = FluxConservingResampler()
linear = LinearInterpolatedResampler()

def air_to_vac(wavelength):
    """
    Implements the air to vacuum wavelength conversion described in eqn 65 of
    Griesen 2006

    Requires dimensionful (air) wavelength input.
    """
    wlum = wavelength.to(u.um).value

    return (1+1e-6*(287.6155+1.62887/wlum**2+0.01360/wlum**4)) * wavelength

def vac_to_air(wavelength):
    """
    Griesen 2006 reports that the error in naively inverting Eqn 65 is less
    than 10^-9 and therefore acceptable.  This is therefore eqn 67
    """
    wlum = wavelength.to(u.um).value
    nl = (1+1e-6*(287.6155+1.62887/wlum**2+0.01360/wlum**4))
    return wavelength/nl

# usable orders in Veloce spectra. 0-based count.
VELOCE_ORDERS = [
    6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
    27, 28, 29, 30, 31, 32, 33, 34, 35, 36
]

# usable orers in Veloce spectra for vsini measurements against the synthetic
# "NextGen" model atmospheres. (The ones that are 20 years old).
VELOCE_VSINI_ORDERS_VS_NEXTGEN = [
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36
]

# a few common lines, in air...
LINE_D = [
    ['BalmerLimit', 3645],
    ['H30', 3662.22],
    ['H29', 3663.41],
    ['H28', 3664.65],
    ['H27', 3666.08],
    ['H26', 3667.73],
    ['H25', 3669.45],
    ['H24->2', 3671.32],
    ['H23->2', 3673.81],
    ['H22->2', 3676.376],
    ['H21->2', 3679.370],
    ['H20->2', 3682.823],
    ['H19->2', 3686.831],
    ['H18->2', 3691.551],
    ['H17->2', 3697.157],
    ['H16->2', 3703.859],
    ['H15->2', 3711.978],
    ['H14->2', 3721.946],
    ['H13->2', 3734.369],
    ['H12->2', 3750.151], # 12->2
    ['H11->2', 3770.633], # 11->2
    ['H10->2', 3797.909], # 10->2, NIST
    ['Fe_L', 3820.44],
    ['H$\eta$', 3835.397], # 9->2
    ['H$\zeta$', 3889.096], # 8->2
    ['He', 3926.7],
    ['Ca_K', 3933.66],
    ['Ca_H', 3968.47],
    ['H$\epsilon$', 3970.075], # 7->2
    ['He[I]', 4026.4],
    ['H$\delta$', 4101.75],
    ['CaI', 4226],
    ['Fe-CR', 4233.3], # Curtiss+1916 (?!)
    ['Fe[I]', 4383], # K5 in absorption
    ['H$\gamma$', 4340.47],
    ['He[I]', 4471],
    [r'H$\beta$', 4861.35],
    ['FeI', 4920.51],
    ['FeI', 4957.61],
    #['FeI', 5169.00], #rough
    ['Feb3', 5168.91],
    ['Mgb4', 5167.33],
    ['Mgb1', 5183.62],
    ['Mgb2', 5172.70],
    ['CaI+FeI', 5269],
    ['Fe_E2', 5270.39],
    ['TiO', 5446],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['Hg_e', 5460.73],
    ['DIB', 5780.6],
    ['He', 5875.618],
    ['NaI_D2', 5889.95],
    ['NaI_D1', 5895.92],
    ['Halpha', 6562.8], # 3->2
    ['DIB', 6613.7],
    ['He (&Fe)', 6678],
    ['FeI_a', 6703.58],
    ['FeI_b', 6705.1],
    ['FeI_c', 6707.44],
    ['LiI_a', 6707.76],
    ['LiI_b', 6707.91],
    ['CaI$\lambda$', 6718],
    ['CaH', 6750],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['O [b]', 6860], # fraunhoffer oxygen
    ['K I', 7665],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['TiO', 7666],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['VO', 7851],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['K I', 7699],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['Rb I', 7800], # Rubidium I resonance line
    ['Rb I', 7947], # Rubidium I resonance line
    ['Na', 8183],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['Na', 8195],# m-dwarfs https://www.stsci.edu/~inr/ldwarf.html
    ['FeI', 8470.7],
    ['CaII(a)', 8500],
    ['SiI', 8504.5],
    ['FeI', 8516.3],
    ['FeI', 8529.0],
    ['CaII(b)', 8544],
    ['SiI', 8559],
    ['FeI', 8614.2],
    ['CaII(c)', 8664.5],
    ['FeI', 8677.1],
    ['NI', 8685.4],
    ['FeI' , 8691.0],
    ['S[I]', 8696.7],
    ['Pa η', 9017.4], # 10->3
    ['Pa ζ', 9231.5], # 9->3
    ['Pa ε', 9548.6], # 8->3
    ['Pa δ', 10052.1], # 7->3
    ['Pa γ', 10941.1], # 6->3
    ['Pa β', 12821.6], # 5->3
    ['Pa α', 18756.1], # 4->3
    #['Ca[I]', air_to_vac(10343.8194*u.angstrom).value], # Muirhead+2020,m from VALD3 database
    ['Ti[I]', air_to_vac(10396.802*u.angstrom).value],  # ''
    ['Ti[I]', air_to_vac(10496.113*u.angstrom).value],  # ''
    ['Ti[I]', air_to_vac(10584.633*u.angstrom).value],  # ''
    ['Ti[I]', air_to_vac(10677.047*u.angstrom).value],  # ''
    ['Ti[I]', air_to_vac(10726.391*u.angstrom).value],  # ''
    ['Ti[I]', air_to_vac(10774.866*u.angstrom).value],  # ''
    ['[S III]', air_to_vac(9531.10052*u.angstrom).value], # http://astronomy.nmsu.edu/drewski/tableofemissionlines.html
    ['[C I]', air_to_vac(9824.130*u.angstrom).value], #''
    ['[C I]', air_to_vac(9850.260*u.angstrom).value], #''
    ['[S VIII]', air_to_vac(9913.000*u.angstrom).value], #''
    ['[Fe XIII]', air_to_vac(10746.800*u.angstrom).value], #''

]

LINELISTDIR = join(os.path.dirname(__path__[0]), 'data/linelists')
# from https://physics.nist.gov/PhysRefData/ASD/lines_form.html, 9000A to 12000A
nir_csvpath = join(LINELISTDIR, "nist_linelist_9000A_to_12000A.csv")

def retrieve_integer_digits(string):
    digits = re.findall(r'\d+', string)
    if digits:
        return int(digits[0])
    else:
        return None
def retrieve_floats(string):
    floats = re.findall(r'\d+\.\d+', string)
    if floats:
        return float(floats[0])
    else:
        return None

df = pd.read_csv(nir_csvpath)
df['intensint'] = df['intens'].apply(retrieve_integer_digits)
df['obs_wl_air'] = df['obs_wl_air(A)'].apply(retrieve_floats)

# add hydrogen lines btwn 9000-12000A
sdf = df[df.element == 'H'].drop_duplicates('obs_wl_air')
for element, wl in zip(sdf.element, sdf.obs_wl_air):
    LINE_D.append([element, air_to_vac(wl*u.angstrom).value])

# add sodium, calcium, potassium, magnesium, helium, silicon lines btwn 9000-12000A
#for el in ['Na', 'Ca', 'K', 'Mg', 'He', 'Si']:
for el in ['Ca', 'He', 'Si']:
    sdf = df[df.element == el].drop_duplicates('obs_wl_air')
    for element, wl in zip(sdf.element, sdf.obs_wl_air):
        LINE_D.append([element, air_to_vac(wl*u.angstrom).value])

# take the top 5% strongest iron lines
sdf = df[df.element == 'Fe'].drop_duplicates('obs_wl_air')
fe_95 = np.nanpercentile(sdf.intensint, 95)
sdf = sdf[sdf.intensint > fe_95]
for element, wl in zip(sdf.element, sdf.obs_wl_air):
    LINE_D.append([element, air_to_vac(wl*u.angstrom).value])

# take the top 10% strongest Ti lines
sdf = df[df.element == 'Ti'].drop_duplicates('obs_wl_air')
ti_90 = np.nanpercentile(sdf.intensint, 90)
sdf = sdf[sdf.intensint > ti_90]
for element, wl in zip(sdf.element, sdf.obs_wl_air):
    LINE_D.append([element, air_to_vac(wl*u.angstrom).value])

# directories
DATADIR = join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = join(os.path.dirname(__path__[0]), 'results/spec_analysis')
TESTOUTDIR = join(OUTDIR, 'tests')
for d in [OUTDIR, TESTOUTDIR]:
    if not os.path.exists(d): os.mkdir(d)

########
# read #
########
def read_nextgen(teff=5500, vsini=5, logg=4.5, vturb=2):
    """
    Read the NextGen Model Atmosphere Grid from Hauschildt et al 1999.

    teff: 4000, 4250, 4500, ... up to 7000
    vsini: 0, 2, 5, 10, 30, 60.
    logg: 4.5
    vturb: 2 or 4.

    returns:
        wavelength, flux
    """


    NEXTGENDIR = join(
        os.path.expanduser('~'),
        'local',
        'HighResGrid'
    )

    if logg==4.5:
        logg = 45
    else:
        raise NotImplementedError

    spectrum_name = '{}_{}_{}_{}'.format(
        str(teff).zfill(5),
        str(logg),
        str(vturb).zfill(3),
        str(vsini).zfill(3)
    )

    spectrum_path = join(
        NEXTGENDIR, spectrum_name
    )

    df = pd.read_csv(spectrum_path, comment='#', delim_whitespace=True,
                     names=['wav', 'physical_flux', 'continuum_flux',
                            'residual_flux'])

    # f,ax = plt.subplots(figsize=(10,3))
    # ax.plot(df.wav*10, df.residual_flux)
    # f.savefig('temp_nextgen.png')

    return 10*nparr(df.wav), nparr(df.residual_flux)


def read_veloce(spectrum_path, start=0, end=None, return_err=False):
    """
    Read Veloce FITS file.
    Return (39x4112) arrays of flux and wavelength.
    39 orders, 4112 pixels in the cross-dispersion direction.

    start and end can be e.g., 200 and -200 to preclean each order.
    """

    hdul = fits.open(spectrum_path)

    # optionally "pre-clean" the data. the edge-pixels are, as a general rule,
    # wonky.
    flux = hdul[0].data[:, start:end]
    flux_err = hdul[1].data[:, start:end]
    wav = hdul[2].data[:, start:end]

    hdul.close()

    if return_err:
        return flux, wav, flux_err
    else:
        return flux, wav


def read_fies(spectrum_path, start=0, end=None, return_err=False):
    """
    Read FIES FITS file.
    Return (90x2062) arrays of flux and wavelength.
    90 orders, 2062 pixels in the cross-dispersion direction.
    """

    hdul = fits.open(spectrum_path)

    # optionally "pre-clean" the data. edge-pixels are, generally, wonky
    flux = hdul[0].data[0, :, start:end]
    wav = hdul[0].data[1, :, start:end]
    flux_err = hdul[0].data[2, :, start:end]

    hdul.close()

    if return_err:
        return flux, wav, flux_err
    else:
        return flux, wav


def read_tres(spectrum_path, start=0, end=None, return_err=False):
    """
    Read TRES FITS file containing the individual orders, extracted and
    rectified (intensity corrected for the blaze and with a wavelength
    solution).

    Return (51x2304) arrays of flux and wavelength.
    51 orders, 2304 pixels in the cross-dispersion direction.
    """

    from cdips_followup.readmultispec import readmultispec

    d = readmultispec(spectrum_path)

    # optionally "pre-clean" the data. edge-pixels are, generally, wonky
    flux = d['flux'][:, start:end]
    wav = d['wavelen'][:, start:end]

    if return_err:
        raise NotImplementedError
    else:
        return flux, wav


def read_hires(spectrum_path, start=0, end=None, is_registered=1, return_err=1):
    """
    Read HIRES FITS file. Assume was "registered" (blaze-corrected, and then
    shifted).
    """
    hdul = fits.open(spectrum_path)

    if is_registered:
        flux = hdul[3].data[:, start:end]
        flux_err = hdul[4].data[:, start:end]
        wav = hdul[5].data[:, start:end]
    else:
        flux = hdul[0].data[:, start:end]
        flux_err = hdul[1].data[:, start:end]
        wav = hdul[2].data[:, start:end]

    hdul.close()

    if return_err:
        return flux, wav, flux_err
    else:
        return flux, wav


def _get_full_hires_spectrum(spectrum_path):

    k = os.path.basename(spectrum_path).replace('.fits','')
    flxs, norm_flxs, wavs= [], [], []
    print('WRN! Assuming HIRES spectrum is already deblazed')
    flx_2d, wav_2d = read_hires(spectrum_path, is_registered=0, return_err=0)
    start = 10
    end = -10
    for order in range(wav_2d.shape[0]):
        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
        flxs.append(flx)
        norm_flxs.append(flx/np.nanmedian(flx))
        wavs.append(wav)

    flx_arr = np.array(flxs)
    norm_flx_arr = np.array(norm_flxs)
    wav_arr = np.array(wavs)

    # get date time RA DEC for barycorr.
    hdul = fits.open(spectrum_path)
    hdr = hdul[0].header
    dateobs = hdr['DATE-OBS']
    utc_time = hdr['UTC']
    mjd = hdr['MJD']
    dec = hdr['DEC']
    ra = hdr['RA']
    hdul.close()

    return flx_arr, norm_flx_arr, wav_arr, mjd, ra, dec


def read_pfs(spectrum_path, wvlen_soln, verbose=False, is_template=False):
    """
    Read PFS IDL SAV file.
    Return (73, 3520) arrays of flux and wavelength.
    73 orders, 3520 pixels in the cross-dispersion direction.

    Args:
        spectrum_path: path to spectrum
        wvlen_soln: path to wavelength solution (might not exist, if you're
        given a template)
        is_template: if True, uses spectrum_path to get the wavelength
        solution.
    """

    if not is_template:
        sp = readsav(spectrum_path, python_dict=True, verbose=verbose)
        wvlen = readsav(wvlen_soln, python_dict=True, verbose=verbose)
        return sp['sp'], wvlen['w']

    else:
        s = readsav(spectrum_path, python_dict=True, verbose=verbose)
        return s['star'], s['w']


def read_feros(spectrum_path):

    hdul = fits.open(spectrum_path)
    d = hdul[0].data
    wav = d[0,0]
    flx = d[3,0]
    hdul.close()
    return flx, wav


def read_gaiaeso(spectrum_path):

    hdul = fits.open(spectrum_path)
    d = hdul[1].data
    wav = 10*d.WAVE # convert to angstrom
    flx = d.FLUX
    hdul.close()
    return flx.flatten(), wav.flatten()

def read_winered(spectrum_path):
    # per warp/Spec1Dtools.py, in the WARP WINERED pipeline

    if spectrum_path.find("fits") == -1:
        spectrum_path += ".fits"

    spfits = fits.open(spectrum_path)
    splength = spfits[0].header["NAXIS1"]
    spdata = spfits[0].data

    rcrval1 = float(spfits[0].header["CRVAL1"])
    rcdelt1 = float(spfits[0].header["CDELT1"])
    rcrpix1 = float(spfits[0].header["CRPIX1"])

    lamx = np.array([rcrval1 + rcdelt1 * (l - rcrpix1 + 1.) for l in range(splength)])
    spfits.close()

    return spdata, lamx


def read_dbsp(spectrum_path):

    hdul = fits.open(spectrum_path)

    d = hdul[1].data

    wav = hdul[1].data['OPT_WAVE']
    flx = hdul[1].data['OPT_COUNTS']

    hdul.close()

    return flx.flatten(), wav.flatten()


def read_neid(filename, read_ccf=True,
              skip_orders=(slice(0, 3, None), slice(118, 123, None))):
    '''
    Read a NEID spectrum.

    NEID Data Format can be found here:
    https://neid.ipac.caltech.edu/docs/NEID-DRP/dataformat.html

    Skips orders 0-2, 118-122, setting them to nan

    Args:
    -----
    filename : str
        Spectrum filename (NEID L1 or L2 .fits file)
    '''
    hdu_list = fits.open(filename)
    header = hdu_list[0].header
    level = header['DATALVL']
    if level < 2:
        read_ccf = False

    w = np.copy(hdu_list[7].data)
    s = np.copy(hdu_list[1].data)
    s_err = np.sqrt(hdu_list[4].data)

    for order in skip_orders:
        s[order] = np.nan
        s_err[order] = np.nan

    if read_ccf:
        ccf_hdu = hdu_list[12]
        ccf_header = ccf_hdu.header
        ccf_start, ccf_step = ccf_header['CCFSTART'], ccf_header['CCFSTEP']
        n_pixels = ccf_header['NAXIS1']
        #  gamma_rv = float(header['QRV'])
        ccf_vels = np.arange(n_pixels) * ccf_step + ccf_start
        ccf_data = np.copy(ccf_hdu.data)
        for order in skip_orders:
            ccf_data[order] = np.nan

        hdu_list.close()
        return w, s, s_err, ccf_vels, ccf_data
    else:
        hdu_list.close()
        return w, s, s_err


def read_galah_given_sobject_id(sobject_id, working_directory, verbose=True,
                                requirefind=True, single_ccd=None):
    """
    Read in all available GALAH DR3 CCDs and return a dictionary.

    Cite Buder+2021 if you use this, since it was based on Sven Buder's GALAH
    DR3 tutorial notebooks:
        https://github.com/svenbuder/GALAH_DR3/blob/master/tutorials/tutorial3_plotting_reduced_spectra.ipynb

    Args:
        sobject_id (np.int64): GALAH DR3 identifier. The most useful way to
        find and work with these is using the GALAH_DR3_main_allstar_v1
        catalog.

        working_directory (str): directory to which the GALAH DR3 spectra have
        already been downloaded.

        requirefind (bool): if True, will raise ValueError if the spectra are
        not found.

        single_ccd (int or None): if passed, will return a single CCD. E.g.,
        for Li 6708, we only need CCD 3.

    Returns:
        specdict (dict): dictionary with the spectra read in, if single_ccd is
        None. Else standard `wav, flx` tuple.
    """

    from glob import glob

    # Check if FITS files already available in working directory
    fits_files = [[], [], [], []]
    ccdlist = [1,2,3,4] if single_ccd is None else [single_ccd]
    for each_ccd in ccdlist:
        globstr = join(
                working_directory, str(sobject_id)+str(each_ccd)+'.fits'
        )
        fits_files[each_ccd-1] = glob(globstr)

    # If not already available, try to download
    for each_ccd in ccdlist:
        if fits_files[each_ccd-1] == []:
            globstr = join(
                    working_directory, str(sobject_id)+str(each_ccd)+'.fits'
            )
            msg = f'Did not find {globstr}'
            raise ValueError(msg)

    spectrum = dict()
    for each_ccd in ccdlist:
        if fits_files[each_ccd-1]!=[]:
            hdul = fits.open(fits_files[each_ccd-1][0])

            # Extension 0: Reduced spectrum
            # Extension 1: Relative error spectrum
            # Extension 4: Normalised spectrum, NB: cut for CCD4

            # Extract wavelength grid for the reduced spectrum
            start_wavelength = hdul[0].header["CRVAL1"]
            dispersion       = hdul[0].header["CDELT1"]
            nr_pixels        = hdul[0].header["NAXIS1"]
            reference_pixel  = hdul[0].header["CRPIX1"]
            if reference_pixel == 0:
                reference_pixel = 1
            spectrum['wave_red_'+str(each_ccd)] = ((np.arange(0,nr_pixels)--reference_pixel+1)*dispersion+start_wavelength)

            # Extract wavelength grid for the normalised spectrum
            start_wavelength = hdul[4].header["CRVAL1"]
            dispersion       = hdul[4].header["CDELT1"]
            nr_pixels        = hdul[4].header["NAXIS1"]
            reference_pixel  = hdul[4].header["CRPIX1"]
            if reference_pixel == 0:
                reference_pixel=1
            spectrum['wave_norm_'+str(each_ccd)] = ((np.arange(0,nr_pixels)--reference_pixel+1)*dispersion+start_wavelength)

            # Extract flux and flux error of reduced spectrum
            spectrum['sob_red_'+str(each_ccd)]  = np.array(hdul[0].data)
            spectrum['uob_red_'+str(each_ccd)]  = np.array(hdul[0].data * hdul[1].data)

            # Extract flux and flux error of reduced spectrum
            spectrum['sob_norm_'+str(each_ccd)] = np.array(hdul[4].data)
            if each_ccd != 4:
                spectrum['uob_norm_'+str(each_ccd)] = np.array(hdul[4].data * hdul[1].data)
            else:
                # for normalised error of CCD4, only used appropriate parts of error spectrum
                spectrum['uob_norm_4'] = np.array(hdul[4].data * (hdul[1].data)[-len(spectrum['sob_norm_4']):])

            hdul.close()
        else:
            spectrum['wave_red_'+str(each_ccd)] = []
            spectrum['wave_norm_'+str(each_ccd)] = []
            spectrum['sob_red_'+str(each_ccd)] = []
            spectrum['sob_norm_'+str(each_ccd)] = []
            spectrum['uob_red_'+str(each_ccd)] = []
            spectrum['uob_norm_'+str(each_ccd)] = []

    spectrum['wave_red'] = np.concatenate(
        ([spectrum['wave_red_'+str(each_ccd)] for each_ccd in ccdlist])
    )
    spectrum['wave_norm'] = np.concatenate(
        ([spectrum['wave_norm_'+str(each_ccd)] for each_ccd in ccdlist])
    )
    spectrum['sob_red'] = np.concatenate(
        ([spectrum['sob_red_'+str(each_ccd)] for each_ccd in ccdlist])
    )
    spectrum['sob_norm'] = np.concatenate(
        ([spectrum['sob_norm_'+str(each_ccd)] for each_ccd in ccdlist])
    )
    spectrum['uob_red'] = np.concatenate(
        ([spectrum['uob_red_'+str(each_ccd)] for each_ccd in ccdlist])
    )
    spectrum['uob_norm'] = np.concatenate(
        ([spectrum['uob_norm_'+str(each_ccd)] for each_ccd in ccdlist])
    )

    if verbose:
        outstr = """
        VERBOSE:
        The spectra of the 4 CCDs are now read in and saved in a dictionary.
        For convenience, the wavelength (wave_*), observed signal (sob_*), and
        uncertainties of the observed signal (uob_*) for both reduced (*red*)
        and normalised (*norm*) are saved both for each CCD (*_1/*_2/*_3/*_4)
        as well as concatenated.
        The dictionary therefore has the keywords
        """
        outstr = re.sub("(?m)^\s+", "", outstr) # fancier than textwrap.dedent
        print(outstr)
        print(spectrum.keys())

    if single_ccd is None:
        return spectrum

    else:
        wav = spectrum[f'wave_red_{single_ccd}']
        flx = spectrum[f'sob_norm_{single_ccd}']
        return flx, wav


def read_galah(spectrum_path, single_ccd, verbose=True):
    """
    Read in GALAH DR3 spectrum from a single CCD.

    Args:
        spectrum_path (str): GALAH DR3 spectrum apth.

        single_ccd (int): GALAH DR3 CCD number. CCD 4 has some specifics.
        Default for Li measurements is 3.

    Returns:
        flx, wav
    """

    hdul = fits.open(spectrum_path)

    # Extension 0: Reduced spectrum
    # Extension 1: Relative error spectrum
    # Extension 4: Normalised spectrum, NB: cut for CCD4

    spectrum = dict()
    each_ccd = single_ccd

    # Extract wavelength grid for the reduced spectrum
    start_wavelength = hdul[0].header["CRVAL1"]
    dispersion       = hdul[0].header["CDELT1"]
    nr_pixels        = hdul[0].header["NAXIS1"]
    reference_pixel  = hdul[0].header["CRPIX1"]
    if reference_pixel == 0:
        reference_pixel = 1
    spectrum['wave_red_'+str(each_ccd)] = ((np.arange(0,nr_pixels)--reference_pixel+1)*dispersion+start_wavelength)

    # Extract wavelength grid for the normalised spectrum
    start_wavelength = hdul[4].header["CRVAL1"]
    dispersion       = hdul[4].header["CDELT1"]
    nr_pixels        = hdul[4].header["NAXIS1"]
    reference_pixel  = hdul[4].header["CRPIX1"]
    if reference_pixel == 0:
        reference_pixel=1
    spectrum['wave_norm_'+str(each_ccd)] = ((np.arange(0,nr_pixels)--reference_pixel+1)*dispersion+start_wavelength)

    # Extract flux and flux error of reduced spectrum
    spectrum['sob_red_'+str(each_ccd)]  = np.array(hdul[0].data)
    spectrum['uob_red_'+str(each_ccd)]  = np.array(hdul[0].data * hdul[1].data)

    # Extract flux and flux error of reduced spectrum
    spectrum['sob_norm_'+str(each_ccd)] = np.array(hdul[4].data)
    if each_ccd != 4:
        spectrum['uob_norm_'+str(each_ccd)] = np.array(hdul[4].data * hdul[1].data)
    else:
        # for normalised error of CCD4, only used appropriate parts of error spectrum
        spectrum['uob_norm_4'] = np.array(hdul[4].data * (hdul[1].data)[-len(spectrum['sob_norm_4']):])

    wav = spectrum[f'wave_red_{single_ccd}']
    flx = spectrum[f'sob_norm_{single_ccd}']
    return flx, wav


##########
# shifts #
##########
def given_vsys_get_li_target_wv(vsys=27*u.km/u.s):

    target_wav = 6708*u.AA

    # deltaf = deltav / c * f0
    # v = f*lambda
    f0 = const.c/target_wav
    deltaf = vsys/ const.c * f0

    new_f = f0+deltaf

    new_lambda = const.c/new_f

    return new_lambda.to(u.AA)


def given_deltawvlen_get_vsys(deltawvlen=2.3*u.AA, wvlen_0=5180*u.AA,
                              time=None, verbose=True):

    # NOTE: currently manual, thru Wright & Eastman's code
    # from astropy.time import Time
    # t = Time(['2020-02-04T00:00:00'], format='isot', scale='utc')
    # print(t.jd)

    # note: varies by ~500m/s due to earth's rotation. (need the exact time to
    # account for this better)
    if time is None:
        barycorr = 0 # m/s from BARYCORR
        print('WRN! barycorr = 0')

    # deltawvlen probably good to ~20-30%, so delta_v no better.
    delta_v = const.c * (deltawvlen / wvlen_0)

    v_star = delta_v + barycorr*(u.m/u.s)
    if verbose:
        print(v_star.to(u.km/u.s))

    return v_star



#################
# visualization #
#################

def viz_1d_spectrum(flx, wav, outpath, xlim=None, vlines=None, names=None,
                    ylim=None, norm_median=False, ylabel=None, xlabel=None,
                    fig=None, ax=None, axtitle=None, show_vel=0, wav0=None):
    """
    Only nonintuitive kwargs:

        show_vel (bool): whether to transform to velocity on x-axis.  If so,
        requires wav0, the reference wavelength, in angstrom.
    """

    if show_vel:
        assert isinstance(wav0, (float, int))

    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        wav = wav[sel]
        flx = flx[sel]

    if norm_median:
        median_flx = np.nanmedian(flx)
        flx /= median_flx

    FNSAVEPLOT = 1
    if fig is None and ax is None:
        plt.close('all')
        fig,ax = plt.subplots(figsize=(10,3))
    else:
        # if you are plotting onto some axis grid, pass the figure too
        assert not (fig is None) and not (ax is None)
        FNSAVEPLOT = 0

    if show_vel:
        def get_vel(wav, wav0):
            deltawvlen = ( wav - wav0 )
            delta_v = const.c * (deltawvlen / wav0)
            delta_v_kms = delta_v.to(u.km/u.s)
            return delta_v_kms.value
        xvals = get_vel(wav, wav0)
    else:
        xvals = wav

    ax.plot(xvals, flx, c='k', zorder=3, lw=0.2)

    y_90 = np.nanpercentile(flx, 90)
    y_10 = np.nanpercentile(flx, 10)
    y_median = np.nanmedian(flx)
    y_diff = y_90 - y_10

    # ax.set_ylim( (y_median-1.1*y_diff, y_median+1.1*y_diff) )

    if not show_vel:
        ax.set_xlabel('λvac [Angstrom]')
    else:
        ax.set_xlabel('Δv [km/s]')

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if xlabel == '':
        ax.set_xticklabels([])

    ax.set_ylabel('flux [e-]')
    if norm_median:
        ax.set_ylabel('flux (median normalized)')
    if ylabel:
        ax.set_ylabel(ylabel)

    if isinstance(xlim, (list, tuple)):
        xmin = xlim[0]
        xmax = xlim[1]
    else:
        xmin = min(xvals)
        xmax = max(xvals)

    this_d = []
    for k, v in LINE_D:
        if v > xmin and v<xmax:
            this_d.append([k, v])

    if isinstance(ylim, (list, tuple)):
        ax.set_ylim(ylim)

    if isinstance(axtitle, str):
        ax.set_title(axtitle)

    if norm_median:
        txt = f'med(f)={median_flx:.1f}e$^-$'
        #trans = blended_transform_factory(ax.transAxes, ax.transData)
        bbox = dict(facecolor='white', edgecolor='none', alpha=0.7)
        ax.text(0.98, 0.98, txt, va='top', ha='right', bbox=bbox,
                transform=ax.transAxes, fontsize='small')

    if len(this_d) > 0:
        for k, v in this_d:
            ylim = ax.get_ylim()
            delta_y = 0.9*(max(ylim) - min(ylim))
            ax.set_ylim(ylim)

            tform = blended_transform_factory(ax.transData, ax.transAxes)
            if show_vel:
                ax.vlines(get_vel(v,wav0), min(ylim)+delta_y, max(ylim), zorder=-3,
                          linestyles=':', color='k', lw=0.3)
                ax.text(get_vel(v,wav0), 0.95, k, ha='center', va='top',
                        transform=tform, fontsize=4)
            else:
                ax.vlines(v, min(ylim)+delta_y, max(ylim), zorder=-3,
                          linestyles=':', color='k', lw=0.3)
                ax.text(v, 0.95, k, ha='center', va='top', transform=tform,
                        fontsize=4)

    if isinstance(xlim, (list, tuple)):
        if not show_vel:
            ax.set_xlim(xlim)

    if isinstance(vlines, list):
        sel = (nparr(vlines)>min(xlim)) & (nparr(vlines)<max(xlim))
        vlines, names = nparr(vlines)[sel], nparr(names)[sel]
        ylim = ax.get_ylim()
        delta_y = 0.9*(max(ylim) - min(ylim))
        if show_vel:
            ax.vlines(get_vel(vlines,wav0), min(ylim)+delta_y, max(ylim),
                      zorder=-3, linestyles=':', color='k', lw=0.3)
        else:
            ax.vlines(vlines, min(ylim)+delta_y, max(ylim), zorder=-3,
                      linestyles=':', color='k', lw=0.3)
        ax.set_ylim(ylim)

        tform = blended_transform_factory(ax.transData, ax.transAxes)
        for x, n in zip(vlines, names):
            if show_vel:
                ax.text(get_vel(x,wav0), 0.95, n, ha='center', va='top',
                        transform=tform, fontsize=4)
            else:
                ax.text(x, 0.95, n, ha='center', va='top', transform=tform,
                        fontsize=4)

    ax.grid(which='both', axis='x', zorder=-3, lw=0.3, ls='--')

    if FNSAVEPLOT:
        format_ax(ax)
        savefig(fig, outpath, writepdf=False)
        plt.close()


def plot_orders(spectrum_path, wvsol_path=None, outdir=None, idstring=None,
                is_template=False, xshift=0, flat_path=None):

    if not isinstance(outdir, str):
        raise ValueError
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path,
                                  is_template=is_template)

        if isinstance(flat_path, str):
            _f = readsav(flat_path, python_dict=True, verbose=True)
            flat_2d = _f['nf']

    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path)

    elif "FIES" in spectrum_path:
        flx_2d, wav_2d = read_fies(spectrum_path)

    elif 'hires' in spectrum_path.lower():
        if 'registered' in spectrum_path:
            flx_2d, wav_2d = read_hires(spectrum_path, is_registered=1,
                                        return_err=0)
        elif 'deblazed' in spectrum_path:
            flx_2d, wav_2d = read_hires(spectrum_path, is_registered=0,
                                        return_err=0)
        else:
            print('WRN! Assuming HIRES spectrum is already deblazed')
            flx_2d, wav_2d = read_hires(spectrum_path, is_registered=0,
                                        return_err=0)

    elif 'TRES' in spectrum_path:
        flx_2d, wav_2d = read_tres(spectrum_path)

    elif 'RVS' in spectrum_path:
        # Gaia RVS
        df = pd.read_csv(spectrum_path)
        flx_2d = np.array(df.flux).reshape((1, len(df)))
        wav_2d = 10*np.array(df.wavelength).reshape((1, len(df)))

    elif 'WINERED' in spectrum_path:
        flx, wav = read_winered(spectrum_path)
        raise NotImplementedError

    else:
        raise NotImplementedError

    for order in range(wav_2d.shape[0]):

        if 'Veloce' in spectrum_path:
            start = 200
            end = -200
        elif "RVS" in spectrum_path:
            start = 0
            end = None
        else:
            start = 10
            end = -10

        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
        if isinstance(flat_path, str):
            flat = flat_2d[order, start:end]

        if 'PFS' in spectrum_path:
            sel = (flx > 0)
            flx, wav = flx[sel], wav[sel]
            if isinstance(flat_path, str):
                flat = flat[sel]

        if isinstance(xshift, (float, int)):
            wav = deepcopy(wav) - xshift

        outname = f'{idstring}_order{str(order).zfill(2)}.png'
        outpath = join(outdir, outname)
        viz_1d_spectrum(flx, wav, outpath)

        if isinstance(flat_path, str):
            outname = f'{idstring}_order{str(order).zfill(2)}_flat.png'
            outpath = join(outdir, outname)
            viz_1d_spectrum(flat, wav, outpath)


def plot_stack_comparison(spectrum_paths, wvsol_path=None, outdir=None,
                          idstring=None, savstr='', is_template=False, xshift=0,
                          flat_path=None, subtractmodel=True, viz_1d=0,
                          labeldict=None, synth_paths=None,
                          rot_broadening_vsini=None):
    """
    Do a visual comparison of spectra in a stack.
    (This might be the same star, with many time-series spectra.  If so,
    subtractmodel=True)
    (Or it could be of different stars, for which you want to compare the
    spectra. Then subtractmodel=False).

    spectrum_paths: list of HIRES spectra to compare against

    synth_paths: list of synthetic spectra (PHOENIX) to compare against, if
        subtractmodel=False
    """

    if not isinstance(outdir, str):
        raise ValueError
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    assert 'hires' in spectrum_paths[0].lower()

    flx_arrs, norm_flx_arrs, wav_arrs = [], [], []

    # collect all spectra; trim and shift if requested
    for sp in spectrum_paths:
        k = os.path.basename(sp).replace('.fits','')
        flxs, norm_flxs, wavs= [], [], []
        print('WRN! Assuming HIRES spectrum is already deblazed')
        flx_2d, wav_2d = read_hires(sp, is_registered=0, return_err=0)
        start = 10
        end = -10
        for order in range(wav_2d.shape[0]):
            flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
            if isinstance(xshift, (float, int)):
                wav = deepcopy(wav) - xshift
            flxs.append(flx)
            norm_flxs.append(flx/np.nanmedian(flx))
            wavs.append(wav)

            if viz_1d:
                outname = f'{idstring}_{k}_order{str(order).zfill(2)}.png'
                outpath = join(outdir, outname)
                if not os.path.exists(outpath):
                    ylim = None
                    if order == 7 and 'bj' in sp.lower():
                        # Ca-HK
                        _flx = deepcopy(flx)
                        viz_1d_spectrum(
                            flx, wav, outpath.replace('.png','_normmedian.png'),
                            ylim=[0,3], norm_median=True
                        )
                        flx = _flx
                        ylim = [0,600]
                    viz_1d_spectrum(flx, wav, outpath, ylim=ylim)
                else:
                    print(f'Found {outpath}')

        flx_arr = np.array(flxs)
        norm_flx_arr = np.array(norm_flxs)
        wav_arr = np.array(wavs)
        flx_arrs.append(flx_arr)
        norm_flx_arrs.append(norm_flx_arr)
        wav_arrs.append(wav_arr)

    # stack, take median (Nspectra x Norders x Nwavelengths)
    flx_arrs = np.array(flx_arrs)
    norm_flx_arrs = np.array(norm_flx_arrs)
    wav_arrs = np.array(wav_arrs)

    # model is the median over all orders.
    if subtractmodel:
        model_flx = np.nanmedian(norm_flx_arrs, axis=0)

    # same viz, on the differences...
    if subtractmodel:
        for sp in spectrum_paths:
            k = os.path.basename(sp).replace('.fits','')
            flx_2d, wav_2d = read_hires(sp, is_registered=0, return_err=0)
            start, end = 10, -10
            for order in range(wav_2d.shape[0]):
                flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
                _model_f = model_flx[order, :]
                assert _model_f.shape == flx.shape
                if isinstance(xshift, (float, int)):
                    wav = deepcopy(wav) - xshift
                norm_flx = flx/np.nanmedian(flx)

                diff_f = norm_flx - _model_f

                outname = f'{idstring}_{k}_order{str(order).zfill(2)}_modelsubtracted.png'
                outpath = join(outdir, outname)
                csvpath = outpath.replace('.png','.csv')
                if not os.path.exists(outpath):
                    ylim = [-0.5,0.5]
                    if subtractmodel:
                        ylabel = 'f/(med_λ(f)) - med_t(f/med_λ(f))'
                    else:
                        ylabel = 'f/(med_λ(f))'
                    viz_1d_spectrum(diff_f, wav, outpath, ylim=ylim,
                                    ylabel=ylabel)
                else:
                    print(f'Found {outpath}')

                if not os.path.exists(csvpath):
                    outdf = pd.DataFrame({
                        'wav': wav,
                        'norm_flx': norm_flx,
                        'model_flx': _model_f,
                        'diff_flx': diff_f
                    })
                    outdf.to_csv(csvpath, index=False)
                    print(f'Made {csvpath}')

    else:

        if isinstance(synth_paths, list):

            syn_flxs, syn_norm_flxs, syn_wavs= [], [], []
            for sp in synth_paths:

                hl = fits.open(sp)
                dirname = os.path.dirname(sp)
                if "HiRes" in dirname:
                    _hl = fits.open(
                        join(dirname, "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")
                    )
                    syn_wav = _hl[0].data
                    _hl.close()
                elif "MedRes" in dirname:
                    # downloaded the δλ = 1A model, fe/h=0, alpha/M=0
                    syn_wav = np.arange(int(3e3), int(1e4), 0.1)
                else:
                    NotImplementedError
                syn_flx = hl[0].data
                # NOTE: might want to actually continuum normalize...
                syn_norm_flx = syn_flx/np.nanmedian(syn_flx)
                hl.close()

                if rot_broadening_vsini is not None:
                    #FIXME this convolution turned out to be much too expensive
                    #over synthetic model grids
                    raise NotImplementedError
                    # python setup.py develop from hd_34382
                    from specmatchtools.wavsol import wav_to_dvel
                    from specmatchtools.kernels import rot

                    dvel = wav_to_dvel(syn_wav)
                    dvel0 = dvel[0]

                    if np.allclose(dvel, dvel[0], rtol=1e-3, atol=1) is False:
                        print("wav not uniform in loglambda, using mean dvel")
                        dvel0 = np.mean(dvel)

                    n = 151 # Correct for VsinI up to ~50 km/s
                    vsini = rot_broadening_vsini*1.
                    varr, M = rot(n, dvel, vsini)
                    syn_flx = nd.convolve1d(syn_flx, M)
                    syn_norm_flx = nd.convolve1d(syn_norm_flx, M)

                syn_flxs.append(syn_flx)
                syn_norm_flxs.append(syn_norm_flx)
                syn_wavs.append(syn_wav)

            # stack, take median (Nspectra x Nwavelengths)
            syn_flx_arrs = np.array(syn_flxs)
            syn_norm_flx_arrs = np.array(syn_norm_flxs)
            syn_wav_arrs = np.array(syn_wav)

        n_spec = flx_arrs.shape[0]
        if isinstance(synth_paths, list):
            n_syn_spec = syn_flx_arrs.shape[0]
        else:
            n_syn_spec = 0
        n_orders = flx_arrs.shape[1]
        n_rows = n_spec + n_syn_spec

        if n_syn_spec == 0:
            star_ids = [
                os.path.basename(os.path.dirname(sp)) for sp in
                spectrum_paths
            ]
        else:
            synth_star_ids = [
                os.path.basename(sp).replace(
                    "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits", ""
                ) for sp in synth_paths
            ]
            obs_star_ids = [
                os.path.basename(os.path.dirname(sp)) for sp in
                spectrum_paths
            ]
            star_ids = synth_star_ids + obs_star_ids

        obs_ids = [
            os.path.basename(sp).replace('.fits','') for sp in
            spectrum_paths
        ]

        from aesthetic.plot import set_style, savefig
        from scipy.ndimage import gaussian_filter1d
        fn = lambda x: gaussian_filter1d(x, sigma=1)

        for ij in range(n_orders):

            ostr = str(ij).zfill(2)
            outname = idstring + "_".join(star_ids) + f"order{ostr}{savstr}.png"
            outpath = join(outdir, outname)

            plt.close("all")
            set_style("clean")

            fig, axs = plt.subplots(nrows=n_rows, figsize=(12,3*n_rows))
            axs = axs.flatten()

            # only empirically observed spectra
            if n_syn_spec == 0:

                for ix in range(n_spec):

                    star_id = star_ids[ix]
                    flx, wav = flx_arrs[ix, ij, :], wav_arrs[ix, ij, :]
                    norm_flx = flx/np.nanmedian(flx)

                    axs[ix].plot(wav, fn(np.array(norm_flx)), c='k', zorder=3, lw=0.2)
                    if ix == 0:
                        wav_min, wav_max = np.nanmin(wav), np.nanmax(wav)

                    # Star name
                    label = labeldict[star_id]
                    props = dict(boxstyle='square', facecolor='white', alpha=0.95,
                                 pad=0.15, linewidth=0)
                    axs[ix].text(0.98, 0.04, label, transform=axs[ix].transAxes,
                                 ha='right',va='bottom', color='k', zorder=43,
                                 fontsize='small', bbox=props)

            # mixing empirical spectra and synthetic PHOENIX spectra
            else:

                for ix in range(n_rows):

                    star_id = star_ids[ix]
                    if ix < n_syn_spec:
                        flx, wav = syn_flx_arrs[ix, :], syn_wav_arrs[:]
                        obs_wav = wav_arrs[-1, ij, :]
                        sel = (wav > min(obs_wav)) & (wav < max(obs_wav))
                        flx, wav = flx[sel], wav[sel]
                        if ix == 0:
                            wav_min, wav_max = np.nanmin(wav), np.nanmax(wav)
                    else:
                        flx, wav = flx_arrs[ix-n_rows, ij, :], wav_arrs[ix-n_rows, ij, :]

                    norm_flx = flx/np.nanmedian(flx)

                    if ix < n_syn_spec:
                        # lte12000-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
                        axs[ix].plot(wav, norm_flx, c='k', zorder=3, lw=0.2)
                        label = (
                            synth_star_ids[ix].
                            replace("lte", "Teff ").replace("-"," logg ")
                        )
                    else:
                        # data!
                        axs[ix].plot(wav, fn(np.array(norm_flx)), c='k', zorder=3, lw=0.2)
                        label = labeldict[star_id]

                    # Star name
                    props = dict(boxstyle='square', facecolor='white', alpha=0.95,
                                 pad=0.15, linewidth=0)
                    axs[ix].text(0.98, 0.04, label, transform=axs[ix].transAxes,
                                 ha='right',va='bottom', color='k', zorder=43,
                                 fontsize='small', bbox=props)

            for ax in axs:
                ax.set_xlim((wav_min, wav_max))

            savefig(fig, outpath, dpi=300, writepdf=0)




def inspect_pfs(nightstr, targetline, xlim=None):
    # NIST has the Li I resonance doublet listed with one transition at 6707.76
    # and the other at 6707.91 A.

    raise DeprecationWarning

    targetid = 'TIC268301217.01'
    datestr = '20200203_{}'.format(nightstr)
    xshift = 2.30
    # note: wavelength soln might change by night...
    spectrum_path = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.{}'.
        format(nightstr)
    )
    wvsol_path = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
    )

    flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)

    if xlim == 'fullorder':
        xlim = None
        vlines = None
        names = None
        if targetline == 'Halpha':
            target_wav = 6562.8
        elif targetline == 'LiI':
            target_wav = 6707.85
        elif targetline == 'Mgb1':
            target_wav = 5183.62

    elif xlim == 'assign':
        if targetline == 'Halpha':
            target_wav = 6562.8
            vlines = [target_wav]
            names = ['Halpha']
        elif targetline == 'LiI':
            target_wav = 6707.85
            # first two are guesstimates from Berger+18 fig3.
            vlines = [6703.58, 6705.1, 6707.44, 6707.76, 6707.91, 6718]
            names = ['FeI', 'FeI', 'FeI', 'Li', '', 'CaI$\lambda$']
        elif targetline == 'Mgb1':
            target_wav = 5183.62
            vlines = [target_wav]
            names = ['Mgb1']

        delta_wav = 5
        xlim = [target_wav-delta_wav, target_wav+delta_wav]

    # retrieve the order corresponding to target wavelength
    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    order = int(np.argwhere((_preorder != wav_2d.shape[1]-1) & (_preorder != 0)))

    flx, wav = flx_2d[order, :], wav_2d[order, :]
    shiftstr = ''
    if isinstance(xshift, (float, int)):
        wav = deepcopy(wav) - xshift
        shiftstr = '_shift{:.2f}'.format(float(xshift))

    xstr = ''
    if isinstance(xlim, list) or isinstance(xlim, tuple):
        xstr = '_{}_xlim{:.1f}-{:.1f}'.format(
            targetline, xlim[0], xlim[1]
        )

    outname = '{}_{}_order{}{}{}.png'.format(
        datestr, targetid, str(order).zfill(2), xstr, shiftstr
    )

    outpath = '../results/spec_analysis/PFS/spec_viz/{}'.format(outname)

    viz_1d_spectrum(flx, wav, outpath, xlim=xlim, vlines=vlines, names=names)


def specmatch_viz_compare(wavlim=[5160,5210]):
    """
    Visualize HIRES spectra for dwarf stars within the wavelengths.
    """

    import specmatchemp.plots as smplot
    import specmatchemp.library

    lib = specmatchemp.library.read_hdf(wavlim=wavlim)

    cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
    g = cut.groupby(pd.cut(cut.Teff,bins=np.arange(5000,7000,250)))
    cut = g.first()

    fig = plt.figure()
    plt.plot(lib.library_params.Teff, lib.library_params.radius,'b.',
             label='_nolegend_')
    plt.plot(cut.Teff, cut.radius,'ro', label='Selected Stars')
    plt.legend()
    smplot.label_axes('Teff','radius')
    wvstr = '{}-{}'.format(wavlim[0], wavlim[1])
    outpath = (
        '../results/spec_analysis/specmatch/quickstart-library-selected-stars_{}.png'.
        format(wvstr)
    )
    savefig(fig, outpath, writepdf=False)

    plt.close('all')

    fig,ax = plt.subplots(figsize=(8,4))
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    bbox = dict(facecolor='white', edgecolor='none',alpha=0.8)
    step = 1
    shift = 0
    for _,row in cut.iterrows():
        spec = lib.library_spectra[row.lib_index,0,:]
        plt.plot(lib.wav, spec.T + shift,color='RoyalBlue',lw=0.5)
        s = "{cps_name:s}, Teff={Teff:.0f}".format(**row)
        plt.text(0.01, 1+shift, s, bbox=bbox, transform=trans)
        shift+=step

    plt.grid()
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Normalized Flux (Arbitrary Offset)')

    outpath = (
        '../results/spec_analysis/specmatch/quickstart-library-selected-stars-spectra_{}.png'.
        format(wvstr)
    )
    savefig(fig, outpath, writepdf=False)


def plot_spec_vs_dwarf_library(wavlim, teff, outdir, idstring, sm_res=None,
                               spectrum_path=None, wvsol_path=None,
                               is_template=False):

    import specmatchemp.library

    lib = specmatchemp.library.read_hdf(wavlim=wavlim)
    if teff < 6500:
        cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
        g = cut.groupby(pd.cut(
            cut.Teff,
            bins=np.arange(np.round(teff-300,-2), np.round(teff+300,-2), 100)
        ))
    else:
        cut = lib.library_params.query('radius < 2.3 and -0.25 < feh < 0.25')
        g = cut.groupby(pd.cut(
            cut.Teff,
            bins=np.arange(6300, 6800, 100)
        ))

    cut = g.first()

    plt.close('all')

    fig,ax = plt.subplots(figsize=(8,4))
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    bbox = dict(facecolor='white', edgecolor='none',alpha=0.8)
    step = 1
    shift = 0
    for _,row in cut.iterrows():
        spec = lib.library_spectra[row.lib_index,0,:]
        plt.plot(lib.wav, spec.T + shift,color='RoyalBlue',lw=0.5)
        s = "{cps_name:s}, Teff={Teff:.0f}".format(**row)
        plt.text(0.01, 1.3+shift, s, bbox=bbox, transform=trans, fontsize=6)
        shift+=step

    if not sm_res is None:
        sm_res.target.plot(offset=shift,
                           plt_kw={'color':'k', 'lw':0.5})
        plt.text(0.01, 1+shift, 'shifted',
                 bbox=bbox, transform=trans, fontsize=6)
        shift += step
        sm_res.target_unshifted.plot(
            offset=shift, normalize=True, plt_kw={'color':'k', 'lw':0.5}
        )
        plt.text(0.01, 1.3+shift, 'Target (unshifted) {}'.format(idstring),
                 bbox=bbox, transform=trans, fontsize=6)
    if not spectrum_path is None:
        ###########################
        # copy in a bunch of code #
        ###########################
        if 'Veloce' in spectrum_path:
            flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
        elif 'PFS' in spectrum_path:
            flx_2d, wav_2d = flx_2d, wav_2d = read_pfs(
                spectrum_path, wvsol_path, is_template=is_template
            )
        else:
            raise NotImplementedError
        target_wav = np.mean(wavlim)
        _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
        viable_orders = np.argwhere((_preorder != wav_2d.shape[1]-1) &
                                    (_preorder != 0))
        order = int(viable_orders[np.argmin(np.abs(
            _preorder[viable_orders] - wav_2d.shape[1]/2))])
        flx, wav = flx_2d[order, :], wav_2d[order, :]
        sel = (wav > wavlim[0]) & (wav < wavlim[1])
        flx, wav = flx[sel], wav[sel]
        if 'Veloce' in spectrum_path:
            wav, flx = wav[::-1], flx[::-1]
        spec = Spectrum1D(spectral_axis=wav*u.AA,
                          flux=flx*u.dimensionless_unscaled)
        cont_flx = fit_generic_continuum(spec)(spec.spectral_axis)
        cont_norm_spec = spec / cont_flx
        flx = cont_norm_spec.flux

        ################
        # done copying #
        ################

        plt.plot(wav, flx+shift, color='k', lw=0.5)
        s = 'target (unshifted), teff={}'.format(teff)
        plt.text(0.01, 1.3+shift, s, bbox=bbox, transform=trans, fontsize=6)
        shift+=step

    plt.grid(True)
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Normalized Flux (Arbitrary Offset)')

    outpath = join(outdir, '{}_compare_speclib.png'.format(idstring))
    savefig(fig, outpath, writepdf=False)


#####################
# continuum fitting #
#####################

@custom_model
def sum_of_Chebyshev1D_deg3(x, c0=0, c1=0, c2=0, d0=0, d1=0, d2=0):
    # NOTE: does not work.

    m1 = Chebyshev1D(3)
    m2 = Chebyshev1D(3)

    m1.c0 = c0
    m1.c1 = c1
    m1.c2 = c2

    m2.c0 = d0
    m2.c1 = d1
    m2.c2 = d2

    return (
        m1.evaluate(x, c0, c1, c2)
        +
        m2.evaluate(x, d0, d1, d2)
    )

def Chebyshev1D_deg3(x, c0=0, c1=0, c2=0):
    # NOTE: also does not work.

    m1 = Chebyshev1D(3)

    return (
        m1.evaluate(x, c0, c1, c2)
    )



def fit_continuum(flx, wav, instrument=None):
    """
    ----------
    Args:

        flx, wav: 1d flux [arbitrary units] and wavelength [angstrom].

    ----------
    Returns:

        cont_flx, cont_norm_spec: fitted continuum flux, and
        continuum-normalized Spectrum1D instance.

    ----------
    Description:

    By default, fit a quadratric polynomial (in fancy speak, a 3rd order
    Chebyshev polynomial series).

    For Veloce, the continuum is better described as a linear combination of
    two quadratics, because there is a notch rougly in the middle of each
    order.
    """

    spec = Spectrum1D(spectral_axis=wav*u.AA,
                      flux=flx*u.dimensionless_unscaled)

    if instrument is None:
        deg = 3
        model = [Chebyshev1D(deg)]
        cont_flx = fit_generic_continuum(
            spec, model=model
        )(spec.spectral_axis)

    if instrument == 'Veloce':
        #model = [sum_of_Chebyshev1D_deg3()]
        #model = [Chebyshev1D_deg3()]
        model = [Chebyshev1D(4)]
        fit_m = fit_generic_continuum(
            spec, model=model
        )
        cont_flx = fit_m[0](spec.spectral_axis)

    cont_norm_spec = spec / cont_flx

    return cont_flx, cont_norm_spec



###########################################
# measure quantities (Li EW, vsini, ....) #
###########################################
def get_Ca_HK_emission(spectrum_path, wvsol_path=None, xshift=None,
                       delta_wav=5, outpath=None, is_template=False):
    """
    spectrum_path: path to PFS or Veloce spectrum

    wvsol_path: path to PFS wavelength solution (optional)

    xshift: angstrom shift required to get into source frame (not vacuum frame).

    delta_wav: window to do the measurement over (angstrom)

    outpath: summary figure is written here.

    is_template: special PFS arg for template path I/O.
    """

    if not isinstance(outpath, str):
        raise ValueError

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path,
                                  is_template=is_template)
        instrument = 'PFS'
    elif "Veloce" in spectrum_path:
        raise ValueError('Veloce does not cover Ca HK')
    else:
        raise NotImplementedError

    # target_wav = 3951.065
    target_wavs = [3933.66, 3968.47]
    vlines = [3933.66, 3968.47]
    names = ['Ca_K', 'Ca_H']

    for target_wav, vline, name in zip(target_wavs, vlines, names):

        xlim = [target_wav-delta_wav, target_wav+delta_wav]

        #
        # retrieve the order corresponding to target wavelength.
        # then shift the wavelength solution to source frame, if needed.
        #
        if instrument == 'PFS':
            if target_wav == 3933.66:
                order = 0
            elif target_wav == 3968.47:
                order = 2
        else:
            raise NotImplementedError

        flx, wav = flx_2d[order, :], wav_2d[order, :]

        thispath = outpath.replace('.png', '{}_1d_check.png'.format(name))
        viz_1d_spectrum(flx, wav, thispath, xlim=xlim, vlines=vlines,
                        names=names)

        shiftstr = ''
        if isinstance(xshift, (float, int)):
            wav = deepcopy(wav) - xshift
            shiftstr = '_shift{:.2f}'.format(float(xshift))

        #
        # cut spectrum to region of interest
        #
        if isinstance(xlim, list):
            xmin = xlim[0]
            xmax = xlim[1]
            sel = (wav > xmin) & (wav < xmax)

            wav = wav[sel]
            flx = flx[sel]

        spec = Spectrum1D(spectral_axis=wav*u.AA,
                          flux=flx*u.dimensionless_unscaled)

        #
        # fit continuum. when doing so, exclude absorption lines.
        #
        if isinstance(xlim, list):
            exclude_regions = []
            for _wv in vlines:
                if xmin < _wv-0.5 and xmax > _wv+0.5:
                    exclude_regions.append(
                        SpectralRegion((_wv-0.5)*u.AA, (_wv+0.5)*u.AA)
                    )

        cont_flx = (
            fit_generic_continuum(spec,
                                  exclude_regions=exclude_regions
            )(spec.spectral_axis)
        )

        cont_norm_spec = spec / cont_flx

        #
        # to fit gaussians, look at 1-flux.
        #
        full_spec = Spectrum1D(spectral_axis=cont_norm_spec.wavelength,
                               flux=(1-cont_norm_spec.flux))

        # #
        # # get the Ca H and Ca K equivalent widths
        # #
        # region = SpectralRegion((target_wav-0.5)*u.AA, (target_wav+0.5)*u.AA)
        # li_equiv_width = equivalent_width(cont_norm_spec, regions=region)
        # li_centroid = centroid(full_spec, region)

        # #
        # # fit a gaussian too, and get ITS equiv width
        # # https://specutils.readthedocs.io/en/stable/fitting.html
        # #
        # g_init = models.Gaussian1D(amplitude=0.2*u.dimensionless_unscaled,
        #                            mean=target_wav*u.AA, stddev=0.5*u.AA)
        # g_fit = fit_lines(full_spec, g_init, window=(region.lower, region.upper))
        # y_fit = g_fit(full_spec.wavelength)

        # fitted_spec = Spectrum1D(spectral_axis=full_spec.wavelength,
        #                          flux=(1-y_fit)*u.dimensionless_unscaled)
        # fitted_li_equiv_width = equivalent_width(fitted_spec, regions=region)

        # #
        # # print bestfit params
        # #
        # print(42*'=')
        # print('got Li equiv width of {}'.format(li_equiv_width))
        # print('got fitted Li equiv width of {}'.format(fitted_li_equiv_width))
        # print('got Li centroid of {}'.format(li_centroid))
        # print('fit gaussian1d params are\n{}'.format(repr(g_fit)))
        # print(42*'=')

        #
        # plot the results
        #
        f,axs = plt.subplots(nrows=2, ncols=1, figsize=(6,8))

        axs[0].plot(wav, flx, c='k', zorder=3)
        axs[0].plot(wav, cont_flx, c='r', zorder=2)

        axs[1].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

        # axs[2].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

        # axs[3].plot(full_spec.wavelength, full_spec.flux, c='k')
        # axs[3].plot(full_spec.wavelength, y_fit, c='g')

        # txt = (
        #     'gaussian1d\namplitude:{:.3f}\nmean:{:.3f}\nstd:{:.3f}\nEW:{:.3f}'.
        #     format(g_fit.amplitude.value,
        #            g_fit.mean.value,
        #            g_fit.stddev.value,
        #            fitted_li_equiv_width)
        # )
        # axs[3].text(
        #     0.95, 0.95, txt, ha='right', va='top', transform=axs[3].transAxes,
        #     fontsize='xx-small'
        # )

        axs[0].set_ylabel('flux')
        axs[1].set_ylabel('contnorm flux')
        # axs[2].set_ylabel('contnorm flux [zoom]')
        # axs[3].set_ylabel('1 - (contnorm flux)')

        if isinstance(xlim, list):
            for ax in axs:
                ax.set_xlim(xlim)

        # axs[2].set_xlim([target_wav-1, target_wav+1])
        # axs[3].set_xlim([target_wav-1, target_wav+1])
        axs[-1].set_xlabel('wavelength [angstrom]')

        for ax in axs:
            format_ax(ax)

        savefig(f, outpath.replace('.png', '_{}.png'.format(name)),
                writepdf=False)



def _given_2d_get_1d_order(wav_2d, flx_2d, target_wav):

    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    viable_orders = np.argwhere(
        (_preorder != wav_2d.shape[1]-1) & (_preorder != 0)
    )
    order = int(
        viable_orders[np.argmin(
            np.abs(_preorder[viable_orders] - wav_2d.shape[1]/2)
        )]
    )

    flx, wav = flx_2d[order, :], wav_2d[order, :]
    return flx, wav


def get_Li_6708_EW(spectrum_path, wvsol_path=None, xshift=None, delta_wav=7.5,
                   outpath=None, is_template=False, writecsvresults=True,
                   verbose=True, montecarlo_errors=True):
    """
    Shift the spectrum, fit a gaussian, and numerically integrate it over a
    window (+/-1 angstrom, centered on 6707.835A) to get the Li 6708A
    equivalent width.  Follows methods from
    https://specutils.readthedocs.io/en/stable/fitting.html

    For uncertainties, if `montecarlo_errors` is True, will perform the
    procedure 20 times, adding continuum noise to the spectrum and repeating
    the fit.  This takes ~40-60 seconds per spectrum, because the
    implementation is not at all optimized.

    Args:

        spectrum_path: path to PFS, Veloce, TRES, or GALAH spectrum

        wvsol_path: path to PFS wavelength solution (optional)

        xshift: angstrom shift required to get into source frame (not vacuum
        frame).  Sign such that if `xshift = 1`, you blueshift the spectrm by 1
        angstrom.  Can also be "find", in which case xshift will be determined
        automatically.  This option currently only works for HIRES.

        delta_wav: window to do the measurement over (angstrom)

        montecarlo_errors: whether to fit the EW many times, adding random
        noise to the spectrum over each iteration.

        outpath: summary figure is written here.

    Returns:
        pandas DataFrame with keys:

            'LiEW_mA': EW of the spectrum sans interpolation or gaussian
            fitting.  At typical R~=20k resolution, underestimates by ~50mA.

            'Fitted_Li_EW_mA': EW of the gaussian interpolated fit. Calibrated
            against Randich+18 for GALAH EW measurements, and agrees at >~50 mA.

            'Fitted_Li_EW_mA_perr', 'Fitted_Li_EW_mA_merr': 84-50 percentile
            and 50-16 percentile of the Monte-Carlo'ed `Fitted_Li_EW_mA`.

            Other keys include ['Li_centroid_A', 'gaussian_fit_amplitude'
            'gaussian_fit_mean', 'gaussian_fit_stddev'].
    """

    if not isinstance(outpath, str):
        raise ValueError

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path,
                                  is_template=is_template)
        instrument = 'PFS'
    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
        instrument = 'Veloce'
    elif "TRES" in spectrum_path:
        flx_2d, wav_2d = read_tres(spectrum_path)
        instrument = 'TRES'
    elif "HIRES" in spectrum_path:
        print('WRN! Assuming HIRES spectrum is already deblazed')
        flx_2d, wav_2d = read_hires(spectrum_path, is_registered=0,
                                    return_err=0)
        instrument = 'HIRES'
    elif "HARPS" in spectrum_path:
        print('WRN! Assuming Louise special format')
        df = pd.read_csv(
                spectrum_path, delim_whitespace=True, names=['wav', 'flx']
        )
        instrument = 'HARPS'
        flx, wav = np.array(df.flx), np.array(df.wav)
    elif "CORALIE" in spectrum_path:
        print('WRN! Assuming Louise special format')
        df = pd.read_csv(
                spectrum_path, delim_whitespace=True, names=['wav', 'flx']
        )
        instrument = 'CORALIE'
        flx, wav = np.array(df.flx), np.array(df.wav)
    elif "galah" in spectrum_path.lower():
        single_ccd = 3 # since we're doing Li measurements
        flx, wav = read_galah(spectrum_path, single_ccd)
        instrument = 'HERMES'
    elif 'gaiaeso' in spectrum_path.lower():
        flx, wav = read_gaiaeso(spectrum_path)
        instrument = 'UVES/GIRAFFE'
    elif 'NEID' in spectrum_path:
        wav_2d, flx_2d, _ = read_neid(spectrum_path, read_ccf=False)
        instrument = 'NEID'
    else:
        raise NotImplementedError

    # FeI at 6703.58 and 6705.1 (Berger+18 Fig3).
    # FeI at 6707.44
    # Li I doublet at 6707.76, 6707.91 -> avg 6707.835
    # CaI lambda at ~6718.
    target_wav = 6707.835
    # "~" here is not a line -- it's a spacer to ensure we don't cut the
    # continuum too closer.
    vlines = [6703.58, 6705.1, 6707.44, 6707.76, 6707.91, 6708.1,
              6710.2, 6713.1, 6718, 6711, 6713.65, 6715.9]
    names = ['FeI', 'FeI', 'FeI', 'Li', 'Li', '~', 'FeI', '?',
             'CaI$\lambda$', '?', "?", "?"]
    xlim = [target_wav-delta_wav, target_wav+delta_wav]

    if instrument in ['Veloce', 'TRES', 'PFS', 'HIRES', 'NEID']:
        #
        # retrieve the order corresponding to target wavelength.
        # then shift the wavelength solution to source frame, if needed.
        #
        flx, wav = _given_2d_get_1d_order(wav_2d, flx_2d, target_wav)

    if instrument == 'Veloce':
        wav = wav[::-1]
        flx = flx[::-1]
        thispath = outpath.replace('.png', '_1d_check.png')
        viz_1d_spectrum(flx, wav, thispath, xlim=(6700, 6725), vlines=vlines,
                        names=names)

    origxshift = deepcopy(xshift)
    if xshift == 'find':
        # automate finding the xshift, by cross-correlating against CK03864,
        # which has RV~=0, no lithium, and nice neighboring iron lines.
        # NOTE: this procedure will probably fail on M-dwarfs.  you might want
        # to do more careful template selection for such cases.

        assert instrument == 'HIRES'
        template_path = join(SPECDIR, "HIRES", "TEMPLATES", "ij172.1101.fits")
        tflx_2d, twav_2d = read_hires(template_path, is_registered=0,
                                      return_err=0)
        tflx, twav = _given_2d_get_1d_order(twav_2d, tflx_2d, target_wav)

        # clip wild positive outliers, ie cosmic rays or emission lines
        nflx = flx/np.nanmedian(flx)
        ntflx = tflx/np.nanmedian(tflx)

        sel = nflx < 1.1
        nflx, nwav = nflx[sel], wav[sel]

        uncertainty = StdDevUncertainty(np.ones(len(nflx)))
        tuncertainty = StdDevUncertainty(np.ones(len(tflx)))

        spec = Spectrum1D(spectral_axis=nwav*u.AA,
                          flux=nflx*u.dimensionless_unscaled,
                          uncertainty=uncertainty)

        tspec = Spectrum1D(spectral_axis=twav*u.AA,
                          flux=ntflx*u.dimensionless_unscaled,
                          uncertainty=tuncertainty)

        from specutils.analysis import template_correlate

        corr, lag = template_correlate(spec, tspec)

        # TODO: add a smoother?  corr and lag can have spikes
        shift_lag = lag[np.argmax(corr)]

        shift_wvlen_arr = (lag * (target_wav*u.AA) / (const.c)).to(u.AA)

        shift_wvlen = (shift_lag * (target_wav*u.AA) / (const.c)).to(u.AA)

        xshift = np.round(shift_wvlen.value, 2)

        MAKE_SHIFTPLOT = 0
        if MAKE_SHIFTPLOT:
            set_style("clean")
            fig = plt.figure(figsize=(10,3))
            axd = fig.subplot_mosaic(
                """
                AAAB
                """
            )
            ax = axd['A']

            ax.plot(nwav, nflx, lw=0.5, label='target')
            ax.plot(twav, 0.5+ntflx, lw=0.5, label='template')
            ax.plot(nwav - xshift, 1+nflx, lw=0.5, label='target (shifted)')
            ax.legend(loc='best', fontsize='xx-small')
            ax.update({'xlabel': 'λ [Angstrom]', 'ylabel': 'f',
                       'ylim':[0.5, 2.1], 'xlim':[6690, 6730]})

            ax = axd['B']
            ax.plot(shift_wvlen_arr.value, corr)
            _sel = (shift_wvlen_arr.value > -3) & (shift_wvlen_arr.value < 3)
            ylim = [min(corr[_sel]), max(corr[_sel])]
            ax.vlines(
                xshift, ylim[0], ylim[1], color='C0', zorder=-1, alpha=0.7, ls=":"
            )
            ax.update({'xlabel': 'shift [Angstrom]', 'ylabel': 'CCF',
                       'title': f'xshift: {xshift}', 'xlim':[-3,3],
                       'ylim':ylim})

            _outpath = outpath.replace(".png", "_auto-xshift.png")
            savefig(fig, _outpath, writepdf=False)
            plt.close("all")

    else:
        ntflx, nflx, nwav, twav = flx*1., flx*1., wav*1., wav*1.


    shiftstr = ''
    if isinstance(xshift, (float, int)):
        print(f"Shifting by {xshift} Angstrom")
        wav = deepcopy(wav) - xshift
        shiftstr = '_shift{:.2f}'.format(float(xshift))
    elif isinstance(xshift, astropy.units.quantity.Quantity):
        raise NotImplementedError('convert xshift to float or int')

    #
    # cut spectrum to region of interest
    #
    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        if verbose:
            print(f'Cutting spectrum to {xmin:.3f} - {xmax:.3f} Angstrom')

        wav = wav[sel]
        flx = flx[sel]

        sel = np.isfinite(flx) & np.isfinite(wav)
        wav, flx = wav[sel], flx[sel]

        assert len(wav) >= 5
        assert len(flx) >= 5

        # exclude cosmic rays, by requiring flux < median(flux) + 5*STD_MAD
        stdev_hat = 1.483 * np.nanmedian(np.abs(flx - np.nanmedian(flx)))
        sel = ~(flx>np.nanmedian(flx)+5*stdev_hat)

        wav, flx = wav[sel], flx[sel]

        if verbose:
            thispath = outpath.replace('.png', '_fittedslice.csv')
            outdf = pd.DataFrame({
                'wav': wav,
                'flx': flx
            })
            outdf.to_csv(thispath, index=False)
            print(f'Cacheing to {thispath}')

    #
    # fit continuum. when doing so, exclude absorption lines.
    #
    if isinstance(xlim, list):
        exclude_regions = []
        for _wv in vlines:
            if xmin < _wv-0.5 and xmax > _wv+0.5:
                exclude_regions.append(
                    SpectralRegion((_wv-0.5)*u.AA, (_wv+0.5)*u.AA)
                )

    region = SpectralRegion((target_wav-1.0)*u.AA, (target_wav+1.0)*u.AA)

    spec = Spectrum1D(spectral_axis=wav*u.AA,
                      flux=flx*u.dimensionless_unscaled)

    cont_flx = fit_generic_continuum(
        spec, exclude_regions=exclude_regions
    )(spec.spectral_axis)

    cont_norm_spec = spec / cont_flx

    #
    # to fit gaussians, look at 1-flux.
    #
    full_spec = Spectrum1D(spectral_axis=cont_norm_spec.wavelength,
                           flux=(1-cont_norm_spec.flux))

    # get the Li EW, assuming that the spectrum has been correctly shifted to
    # within a 1A region.
    li_equiv_width = equivalent_width(cont_norm_spec, regions=region)
    li_centroid = centroid(full_spec, region)

    # fit a gaussian, and integrate to get ITS equivalent width
    g_init = models.Gaussian1D(
        amplitude=0.2*u.dimensionless_unscaled,
        mean=target_wav*u.AA, stddev=0.5*u.AA
    )
    g_fit = fit_lines(
        full_spec, g_init, window=(region.lower, region.upper)
    )

    min_x, max_x, N_x = (
        min(full_spec.wavelength), max(full_spec.wavelength), int(1e4)
    )
    x_fit = np.linspace(min_x, max_x, N_x)
    y_fit = g_fit(x_fit)

    fitted_spec = Spectrum1D(
        spectral_axis=x_fit, flux=(1-y_fit)*u.dimensionless_unscaled
    )
    fitted_li_equiv_width = equivalent_width(
        fitted_spec, regions=region
    )


    if montecarlo_errors:
        N_montecarlo = int(2e1)

        if exclude_regions is not None:
            excised_spectrum = excise_regions(spec, exclude_regions)
            continuum_rms = np.nanstd(excised_spectrum.flux).value

        fitted_li_equiv_widths = []

        for i in range(N_montecarlo):
            if verbose:
                print(f'{datetime.utcnow().isoformat()}: {i}/{N_montecarlo}')

            np.random.seed(i)
            err = np.random.normal(loc=0, scale=continuum_rms, size=len(wav))

            spec = Spectrum1D(spectral_axis=wav*u.AA,
                              flux=(flx+err)*u.dimensionless_unscaled)

            cont_flx = fit_generic_continuum(
                spec, exclude_regions=exclude_regions
            )(spec.spectral_axis)

            cont_norm_spec = spec / cont_flx

            # to fit gaussians, look at 1-flux.
            # fit a gaussian, and integrate to get ITS equivalent width
            full_spec = Spectrum1D(
                spectral_axis=cont_norm_spec.wavelength,
                flux=(1-cont_norm_spec.flux)
            )

            g_init = models.Gaussian1D(
                amplitude=0.2*u.dimensionless_unscaled,
                mean=target_wav*u.AA, stddev=0.5*u.AA
            )
            try:
                g_fit = fit_lines(
                    full_spec, g_init, window=(region.lower, region.upper)
                )
            except NonFiniteValueError as e:
                txt = (
                    f'{datetime.utcnow().isoformat()}: {spectrum_path} failed!'
                    'Got error:\n'
                    f'{repr(e)}'
                )
                print(txt)
                fitted_li_equiv_widths.append(
                    np.nan
                )
                continue

            y_fit = g_fit(x_fit)

            fitted_spec = Spectrum1D(
                spectral_axis=x_fit,
                flux=(1-y_fit)*u.dimensionless_unscaled
            )
            try:
                _fitted_li_equiv_width = equivalent_width(
                    fitted_spec, regions=region
                )
                fitted_li_equiv_widths.append(
                    _fitted_li_equiv_width.to(u.angstrom).value
                )
            except IndexError as e:
                txt = (
                    f'{datetime.utcnow().isoformat()}: {spectrum_path} failed!'
                    'Got error:\n'
                    f'{repr(e)}'
                )
                print(txt)
                fitted_li_equiv_widths.append(
                    np.nan
                )

        fitted_li_equiv_widths = nparr(fitted_li_equiv_widths)
        # drop any nans
        fitted_li_equiv_widths = (
            fitted_li_equiv_widths[~pd.isnull(fitted_li_equiv_widths)]
        )

        fitted_li_equiv_widths_perr = (
            np.percentile(fitted_li_equiv_widths, 84)
            -
            np.percentile(fitted_li_equiv_widths, 50)
        )
        fitted_li_equiv_widths_merr = (
            np.percentile(fitted_li_equiv_widths, 50)
            -
            np.percentile(fitted_li_equiv_widths, 16)
        )

    else:
        fitted_li_equiv_widths_perr = np.nan
        fitted_li_equiv_widths_merr = np.nan


    #
    # print bestfit params
    #
    if verbose:
        print(42*'=')
        print(f'got Li equiv width of {li_equiv_width:.3f}')
        print(f'got fitted Li equiv width of {fitted_li_equiv_width:.3f}, '
              f'from gaussian over {min_x:.1f} - {max_x:.1f}, with '
              f'{N_x} points. ')
        print(f'got Li centroid of {li_centroid:.3f}')
        print(f'fit gaussian1d params are\n{repr(g_fit)}')
        print(42*'=')

    #
    # plot the results
    #
    set_style("science")

    f = plt.figure(figsize=(6, 10))
    if origxshift == 'find':
        axd = f.subplot_mosaic(
            """
            AAAB
            CCCC
            DDDD
            EEEE
            FFFF
            """
        )
    else:
        axd = f.subplot_mosaic(
            """
            CCCC
            DDDD
            EEEE
            FFFF
            """
        )

    # only when origxshift == 'find'... show CCF vs RV
    if origxshift == 'find':
        # raw flux vs wavelength
        ax = axd['A']
        ax.plot(nwav, nflx, lw=0.5, label='target')
        ax.plot(twav, 0.5+ntflx, lw=0.5, label='template')
        ax.plot(nwav - xshift, 1+nflx, lw=0.5, label='target (shifted)')
        ax.legend(loc='best', fontsize=4)
        ax.update({'ylabel': '$f$',
                   'ylim':[0.4, 2.1], 'xlim':[6691, 6729]})

        ax = axd['B']
        _corr = corr/np.nanmax(corr)
        ax.plot(shift_wvlen_arr.value, _corr)
        _sel = (shift_wvlen_arr.value > -3) & (shift_wvlen_arr.value < 3)
        ylim = [min(_corr[_sel]), max(_corr[_sel])]
        ax.vlines(
            xshift, ylim[0], ylim[1], color='C0', zorder=-1, alpha=0.7, ls=":"
        )
        txt = f'Δx: {xshift}'
        bbox = dict(facecolor='white', alpha=0.9, pad=0, edgecolor='white')
        ax.text(
            0.95, 0.95, txt, ha='right', va='top', transform=axd['B'].transAxes,
            fontsize='xx-small', zorder=1, bbox=bbox
        )
        ax.update({'ylabel': 'CCF', 'xlim':[-3,3], 'ylim':ylim})
        ax.set_yticklabels([])

    ax = axd['C']
    ax.plot(wav, flx, c='k', zorder=3)

    if montecarlo_errors:
        ax.scatter(excised_spectrum.wavelength, excised_spectrum.flux,
                    c='gray', zorder=4, s=5)

    ax.plot(wav, cont_flx, c='r', zorder=2)

    ax = axd['D']
    ax.plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    ax = axd['E']
    ax.plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    ylim = ax.get_ylim()
    ax.vlines([region.lower.value, region.upper.value], min(ylim), max(ylim),
              color='orangered', linestyle='--', zorder=-2, lw=0.5, alpha=0.3)

    ax = axd['F']
    ax.plot(full_spec.wavelength, full_spec.flux, c='k')
    ax.plot(x_fit, y_fit, c='g')
    ylim = ax.get_ylim()
    ax.vlines([region.lower.value, region.upper.value], min(ylim), max(ylim),
              color='g', linestyle='--', zorder=-2, lw=0.5, alpha=0.3)

    if not montecarlo_errors:
        txt = (
            f'gaussian1d\namplitude:{g_fit.amplitude.value:.3f}'
            f'\nmean:{g_fit.mean.value:.3f}'
            f'\nstd:{g_fit.stddev.value:.3f}'
            f'\nEW:{(li_equiv_width*1e3).value:.1f}mA'
            f'\nGaussFitEW:{(fitted_li_equiv_width*1e3).value:.1f}mA'
        )
    else:
        txt = (
            f'gaussian1d\namplitude:{g_fit.amplitude.value:.3f}'
            f'\nmean:{g_fit.mean.value:.3f}'
            f'\nstd:{g_fit.stddev.value:.3f}'
            f'\nEW:{(li_equiv_width*1e3).value:.1f}mA'
            f'\nGaussFitEW:{(fitted_li_equiv_width*1e3).value:.1f}'
            f'+{fitted_li_equiv_widths_perr*1e3:.1f}'
            f'-{fitted_li_equiv_widths_merr*1e3:.1f}mA'
        )

    ax.text(
        0.95, 0.95, txt, ha='right', va='top', transform=axd['F'].transAxes,
        fontsize='xx-small'
    )

    axd['C'].set_ylabel('$f$')
    axd['D'].set_ylabel('$f_\mathrm{contnorm}$')
    axd['E'].set_ylabel('$f_\mathrm{contnorm}$ [zoom]')
    axd['F'].set_ylabel('$1-f_\mathrm{contnorm}$')

    assert isinstance(xlim, list)
    for ax in [axd['C'], axd['D'], axd['E'], axd['F']]:
        ax.set_xlim(xlim)

    # Plot known lines
    if isinstance(vlines, list):
        sel = (nparr(vlines)>min(xlim)) & (nparr(vlines)<max(xlim))
        vlines, names = nparr(vlines)[sel], nparr(names)[sel]
        ylim = axd['C'].get_ylim()
        delta_y = 0.9*(max(ylim) - min(ylim))
        axd['C'].vlines(vlines, min(ylim)+delta_y, max(ylim), zorder=-3,
                      linestyles=':', color='k', lw=0.3)
        axd['C'].set_ylim(ylim)

        tform = blended_transform_factory(axd['C'].transData, axd['C'].transAxes)
        for x, n in zip(vlines, names):
            axd['C'].text(x, 0.95, n, ha='center', va='top', transform=tform,
                        fontsize=4)

    # Plot xshift diagnostic bars
    target_wav = 6707.835
    for xbar in np.arange(0, 2.0, 0.3333):

        for ax in [axd['C'], axd['D'], axd['E']]:
            ylim = ax.get_ylim()
            delta_y = 0.7*(max(ylim) - min(ylim))
            ax.vlines([target_wav], min(ylim), max(ylim)-delta_y, zorder=-3,
                      linestyles=':', color='C0', lw=0.6)
            ax.vlines([target_wav-0.5], min(ylim), max(ylim)-delta_y, zorder=-3,
                      linestyles=':', color='C1', lw=0.6)
            ax.vlines([target_wav+0.5], min(ylim), max(ylim)-delta_y, zorder=-3,
                      linestyles=':', color='C1', lw=0.6)
            ymid = np.mean([min(ylim), max(ylim)-delta_y])
            ax.vlines([target_wav-1.0], min(ylim), max(ylim)-delta_y, zorder=-3,
                      linestyles=':', color='C2', lw=0.3)
            ax.vlines([target_wav+1.0], min(ylim), max(ylim)-delta_y, zorder=-3,
                      linestyles=':', color='C2', lw=0.3)
            ax.plot(
                [target_wav - xbar, target_wav + xbar],
                [ymid, ymid],
                c='C0',
                lw=5,
                alpha=0.2
            )
            ax.set_ylim(ylim)


    axd['E'].set_xlim([target_wav-1.5, target_wav+1.5])
    axd['F'].set_xlim([target_wav-1.5, target_wav+1.5])
    axd['F'].set_xlabel('$\lambda$ [angstrom] '+ f'(xshift={xshift:.4f})')

    f.tight_layout()

    savefig(f, outpath, writepdf=False)

    if writecsvresults:
        outpath = outpath.replace('.png','_results.csv')
        outdict = {
            'xshift': xshift,
            # LiEW_mA: result from adding up the flux
            'LiEW_mA': np.round(li_equiv_width*1e3, 3),
            # Fitted_Li_EW_mA: result from Gaussian fit
            'Fitted_Li_EW_mA': np.round(fitted_li_equiv_width*1e3, 3),
            'Fitted_Li_EW_mA_perr': np.round(fitted_li_equiv_widths_perr*1e3, 3),
            'Fitted_Li_EW_mA_merr': np.round(fitted_li_equiv_widths_merr*1e3, 3),
            'Li_centroid_A': np.round(li_centroid, 3),
            'gaussian_fit_amplitude': np.round(g_fit.amplitude.value, 4),
            'gaussian_fit_mean': np.round(g_fit.mean.value, 4),
            'gaussian_fit_stddev': np.round(g_fit.stddev.value, 4)
        }
        outdf = pd.DataFrame(outdict, index=[0])
        outdf.to_csv(outpath, index=False)
        print(f'Made {outpath}')


def get_line_EW(spectrum_path, xshift=None, delta_wav=10, outpath=None,
                writecsvresults=True, verbose=True, montecarlo_errors=True,
                linewindowwidth_angstrom=1, target_wav=4861.35,
                isemissionline=True, dintwav=4, dblock=6, linename='hbeta'):
    """
    Shift the spectrum, fit a gaussian, and numerically integrate it over a
    window (+/- `linewindowwidth_angstrom` angstrom, centered on A) to get the
    line EW.  See e.g. get_Li_6708_EW

    Args:

        linecenter: in angstroms, e.g. 6562.8 for halpha.

        isemissionline: sets default guess for fit line amplitude

        ...rest is as in get_Li_6708_EW

    Returns:
        pandas DataFrame with keys:

            'EW_mA': EW of the spectrum sans interpolation or gaussian
            fitting.  At typical R~=20k resolution, underestimates by ~50mA.

            'Fitted_EW_mA': EW of the gaussian interpolated fit. Calibrated
            against Randich+18 for GALAH EW measurements, and agrees at >~50 mA.

            'Fitted_EW_mA_perr', 'Fitted_EW_mA_merr': 84-50 percentile
            and 50-16 percentile of the Monte-Carlo'ed `Fitted_EW_mA`.

            Other keys include ['centroid_A', 'gaussian_fit_amplitude'
            'gaussian_fit_mean', 'gaussian_fit_stddev'].
    """

    if not isinstance(outpath, str):
        raise ValueError

    if "DBSP" in spectrum_path:
        flx, wav = read_dbsp(spectrum_path)
        instrument = 'DBSP'
    elif "HIRES" in spectrum_path:
        print('WRN! Assuming HIRES spectrum is already deblazed')
        flx_2d, wav_2d = read_hires(spectrum_path, is_registered=0,
                                    return_err=0)
        flx, wav = _given_2d_get_1d_order(wav_2d, flx_2d, target_wav)
        instrument = 'HIRES'
    else:
        raise NotImplementedError

    xlim = [target_wav-delta_wav, target_wav+delta_wav]

    origxshift = deepcopy(xshift)

    shiftstr = ''
    if isinstance(xshift, (float, int)):
        print(f"Shifting by {xshift} Angstrom")
        wav = deepcopy(wav) - xshift
        shiftstr = '_shift{:.2f}'.format(float(xshift))
    elif isinstance(xshift, astropy.units.quantity.Quantity):
        raise NotImplementedError('convert xshift to float or int')

    #
    # cut spectrum to region of interest
    #
    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        if verbose:
            print(f'Cutting spectrum to {xmin:.3f} - {xmax:.3f} Angstrom')

        wav = wav[sel]
        flx = flx[sel]

        sel = np.isfinite(flx) & np.isfinite(wav)
        wav, flx = wav[sel], flx[sel]

        ## exclude cosmic rays, by requiring flux < median(flux) + 5*STD_MAD
        #stdev_hat = 1.483 * np.nanmedian(np.abs(flx - np.nanmedian(flx)))
        #sel = ~(flx>np.nanmedian(flx)+5*stdev_hat)

        #wav, flx = wav[sel], flx[sel]

        if verbose:
            thispath = outpath.replace('.png', '_fittedslice.csv')
            outdf = pd.DataFrame({
                'wav': wav,
                'flx': flx
            })
            outdf.to_csv(thispath, index=False)
            print(f'Cacheing to {thispath}')

    #
    # fit continuum. when doing so, exclude absorption lines.
    #
    vlines = [target_wav]
    names = [linename]
    if isinstance(xlim, list):
        exclude_regions = []
        for _wv in vlines:
            if xmin < _wv-dblock and xmax > _wv+dblock:
                exclude_regions.append(
                    SpectralRegion((_wv-dblock)*u.AA, (_wv+dblock)*u.AA)
                )

    region = SpectralRegion((target_wav-dintwav)*u.AA, (target_wav+dintwav)*u.AA)

    spec = Spectrum1D(spectral_axis=wav*u.AA,
                      flux=flx*u.dimensionless_unscaled)

    cont_flx = fit_generic_continuum(
        spec, exclude_regions=exclude_regions
    )(spec.spectral_axis)

    cont_norm_spec = spec / cont_flx

    if isemissionline:
        full_spec = Spectrum1D(spectral_axis=cont_norm_spec.wavelength,
                               flux=(cont_norm_spec.flux-1))
    else:
        # to fit gaussians to absorption lines... look at 1-flux.
        full_spec = Spectrum1D(spectral_axis=cont_norm_spec.wavelength,
                               flux=(1-cont_norm_spec.flux))


    # get the Li EW, assuming that the spectrum has been correctly shifted to
    # within a `2*dintwav` Angstrom region.
    equiv_width = equivalent_width(cont_norm_spec, regions=region)
    _centroid = centroid(full_spec, region)

    # fit a gaussian, and integrate to get ITS equivalent width
    g_init = models.Gaussian1D(
        amplitude=0.2*u.dimensionless_unscaled,
        mean=target_wav*u.AA, stddev=0.5*u.AA
    )
    g_fit = fit_lines(
        full_spec, g_init, window=(region.lower, region.upper)
    )

    min_x, max_x, N_x = (
        min(full_spec.wavelength), max(full_spec.wavelength), int(1e4)
    )
    x_fit = np.linspace(min_x, max_x, N_x)
    y_fit = g_fit(x_fit)

    if isemissionline:
        fitted_spec = Spectrum1D(
            spectral_axis=x_fit, flux=(y_fit-1)*u.dimensionless_unscaled
        )
    else:
        fitted_spec = Spectrum1D(
            spectral_axis=x_fit, flux=(1-y_fit)*u.dimensionless_unscaled
        )
    fitted_equiv_width = equivalent_width(
        fitted_spec, regions=region
    )


    if montecarlo_errors:
        N_montecarlo = int(2e1)

        if exclude_regions is not None:
            excised_spectrum = excise_regions(spec, exclude_regions)
            continuum_rms = np.nanstd(excised_spectrum.flux).value

        fitted_equiv_widths = []

        for i in range(N_montecarlo):
            if verbose:
                print(f'{datetime.utcnow().isoformat()}: {i}/{N_montecarlo}')

            np.random.seed(i)
            err = np.random.normal(loc=0, scale=continuum_rms, size=len(wav))

            spec = Spectrum1D(spectral_axis=wav*u.AA,
                              flux=(flx+err)*u.dimensionless_unscaled)

            cont_flx = fit_generic_continuum(
                spec, exclude_regions=exclude_regions
            )(spec.spectral_axis)

            cont_norm_spec = spec / cont_flx

            # to fit gaussians, look at 1-flux.
            # fit a gaussian, and integrate to get ITS equivalent width
            if isemissionline:
                full_spec = Spectrum1D(
                    spectral_axis=cont_norm_spec.wavelength,
                    flux=(cont_norm_spec.flux-1)
                )
            else:
                full_spec = Spectrum1D(
                    spectral_axis=cont_norm_spec.wavelength,
                    flux=(1-cont_norm_spec.flux)
                )

            g_init = models.Gaussian1D(
                amplitude=0.2*u.dimensionless_unscaled,
                mean=target_wav*u.AA, stddev=0.5*u.AA
            )
            g_fit = fit_lines(
                full_spec, g_init, window=(region.lower, region.upper)
            )

            y_fit = g_fit(x_fit)

            if isemissionline:
                fitted_spec = Spectrum1D(
                    spectral_axis=x_fit,
                    flux=(y_fit-1)*u.dimensionless_unscaled
                )
            else:
                fitted_spec = Spectrum1D(
                    spectral_axis=x_fit,
                    flux=(1-y_fit)*u.dimensionless_unscaled
                )
            try:
                _fitted_equiv_width = equivalent_width(
                    fitted_spec, regions=region
                )
                fitted_equiv_widths.append(
                    _fitted_equiv_width.to(u.angstrom).value
                )
            except IndexError as e:
                txt = (
                    f'{datetime.utcnow().isoformat()}: {spectrum_path} failed!'
                    'Got error:\n'
                    f'{repr(e)}'
                )
                print(txt)
                fitted_equiv_widths.append(
                    np.nan
                )

        fitted_equiv_widths = nparr(fitted_equiv_widths)
        # drop any nans
        fitted_equiv_widths = (
            fitted_equiv_widths[~pd.isnull(fitted_equiv_widths)]
        )

        fitted_equiv_widths_perr = (
            np.percentile(fitted_equiv_widths, 84)
            -
            np.percentile(fitted_equiv_widths, 50)
        )
        fitted_equiv_widths_merr = (
            np.percentile(fitted_equiv_widths, 50)
            -
            np.percentile(fitted_equiv_widths, 16)
        )

    else:
        fitted_equiv_widths_perr = np.nan
        fitted_equiv_widths_merr = np.nan


    #
    # print bestfit params
    #
    if verbose:
        print(42*'=')
        print(f'got line equiv width of {equiv_width:.3f}')
        print(f'got fitted equiv width of {fitted_equiv_width:.3f}, '
              f'from gaussian over {min_x:.1f} - {max_x:.1f}, with '
              f'{N_x} points. ')
        print(f'got centroid of {_centroid:.3f}')
        print(f'fit gaussian1d params are\n{repr(g_fit)}')
        print(42*'=')

    #
    # plot the results
    #
    set_style("science")

    f = plt.figure(figsize=(6, 10))
    if origxshift == 'find':
        axd = f.subplot_mosaic(
            """
            AAAB
            CCCC
            DDDD
            EEEE
            FFFF
            """
        )
    else:
        axd = f.subplot_mosaic(
            """
            CCCC
            DDDD
            EEEE
            FFFF
            """
        )

    ax = axd['C']
    ax.plot(wav, flx, c='k', zorder=3)

    if montecarlo_errors:
        ax.scatter(excised_spectrum.wavelength, excised_spectrum.flux,
                    c='gray', zorder=4, s=5)

    ax.plot(wav, cont_flx, c='r', zorder=2)

    ax = axd['D']
    ax.plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    ax = axd['E']
    ax.plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    ylim = ax.get_ylim()
    ax.vlines([region.lower.value, region.upper.value], min(ylim), max(ylim),
              color='orangered', linestyle='--', zorder=-2, lw=0.5, alpha=0.3)

    ax = axd['F']
    ax.plot(full_spec.wavelength, full_spec.flux, c='k')
    ax.plot(x_fit, y_fit, c='g')
    ylim = ax.get_ylim()
    ax.vlines([region.lower.value, region.upper.value], min(ylim), max(ylim),
              color='g', linestyle='--', zorder=-2, lw=0.5, alpha=0.3)

    if not montecarlo_errors:
        txt = (
            f'gaussian1d\namplitude:{g_fit.amplitude.value:.3f}'
            f'\nmean:{g_fit.mean.value:.3f}'
            f'\nstd:{g_fit.stddev.value:.3f}'
            f'\nEW:{(equiv_width*1e3).value:.1f}mA'
            f'\nGaussFitEW:{(fitted_equiv_width*1e3).value:.1f}mA'
        )
    else:
        txt = (
            f'gaussian1d\namplitude:{g_fit.amplitude.value:.3f}'
            f'\nmean:{g_fit.mean.value:.3f}'
            f'\nstd:{g_fit.stddev.value:.3f}'
            f'\nEW:{(equiv_width*1e3).value:.1f}mA'
            f'\nGaussFitEW:{(fitted_equiv_width*1e3).value:.1f}'
            f'+{fitted_equiv_widths_perr*1e3:.1f}'
            f'-{fitted_equiv_widths_merr*1e3:.1f}mA'
        )

    ax.text(
        0.95, 0.95, txt, ha='right', va='top', transform=axd['F'].transAxes,
        fontsize='xx-small'
    )

    axd['C'].set_ylabel('$f$')
    axd['D'].set_ylabel('$f_\mathrm{contnorm}$')
    axd['E'].set_ylabel('$f_\mathrm{contnorm}$ [zoom]')
    axd['F'].set_ylabel('$1-f_\mathrm{contnorm}$')

    assert isinstance(xlim, list)
    for ax in [axd['C'], axd['D'], axd['E'], axd['F']]:
        ax.set_xlim(xlim)

    # Plot known lines
    if isinstance(vlines, list):
        sel = (nparr(vlines)>min(xlim)) & (nparr(vlines)<max(xlim))
        vlines, names = nparr(vlines)[sel], nparr(names)[sel]
        ylim = axd['C'].get_ylim()
        delta_y = 0.9*(max(ylim) - min(ylim))
        axd['C'].vlines(vlines, min(ylim)+delta_y, max(ylim), zorder=-3,
                      linestyles=':', color='k', lw=0.3)
        axd['C'].set_ylim(ylim)

        tform = blended_transform_factory(axd['C'].transData, axd['C'].transAxes)
        for x, n in zip(vlines, names):
            axd['C'].text(x, 0.95, n, ha='center', va='top', transform=tform,
                        fontsize=4)

    ## Plot xshift diagnostic bars
    ##target_wav = 6707.835
    #for xbar in np.arange(0, 2.0, 0.3333):

    #    for ax in [axd['C'], axd['D'], axd['E']]:
    #        ylim = ax.get_ylim()
    #        delta_y = 0.7*(max(ylim) - min(ylim))
    #        ax.vlines([target_wav], min(ylim), max(ylim)-delta_y, zorder=-3,
    #                  linestyles=':', color='C0', lw=0.6)
    #        ax.vlines([target_wav-0.5], min(ylim), max(ylim)-delta_y, zorder=-3,
    #                  linestyles=':', color='C1', lw=0.6)
    #        ax.vlines([target_wav+0.5], min(ylim), max(ylim)-delta_y, zorder=-3,
    #                  linestyles=':', color='C1', lw=0.6)
    #        ymid = np.mean([min(ylim), max(ylim)-delta_y])
    #        ax.vlines([target_wav-1.0], min(ylim), max(ylim)-delta_y, zorder=-3,
    #                  linestyles=':', color='C2', lw=0.3)
    #        ax.vlines([target_wav+1.0], min(ylim), max(ylim)-delta_y, zorder=-3,
    #                  linestyles=':', color='C2', lw=0.3)
    #        ax.plot(
    #            [target_wav - xbar, target_wav + xbar],
    #            [ymid, ymid],
    #            c='C0',
    #            lw=5,
    #            alpha=0.2
    #        )
    #        ax.set_ylim(ylim)


    axd['E'].set_xlim([target_wav-dblock, target_wav+dblock])
    axd['F'].set_xlim([target_wav-dblock, target_wav+dblock])
    axd['F'].set_xlabel('$\lambda$ [angstrom] '+ f'(xshift={xshift:.4f})')

    f.tight_layout()

    savefig(f, outpath, writepdf=False)

    if writecsvresults:
        outpath = outpath.replace('.png','_results.csv')
        outdict = {
            'xshift': xshift,
            'EW_mA': np.round(equiv_width*1e3, 3),
            'Fitted_EW_mA': np.round(fitted_equiv_width*1e3, 3),
            'Fitted_EW_mA_perr': np.round(fitted_equiv_widths_perr*1e3, 3),
            'Fitted_EW_mA_merr': np.round(fitted_equiv_widths_merr*1e3, 3),
            'centroid_A': np.round(_centroid, 3),
            'gaussian_fit_amplitude': np.round(g_fit.amplitude.value, 4),
            'gaussian_fit_mean': np.round(g_fit.mean.value, 4),
            'gaussian_fit_stddev': np.round(g_fit.stddev.value, 4)
        }
        outdf = pd.DataFrame(outdict, index=[0])
        outdf.to_csv(outpath, index=False)
        print(f'Made {outpath}')



def measure_vsini(wav, flx, flxerr=None, teff=6000, logg=4.5, vturb=2, outdir=None,
                  targetname=None):

    if not isinstance(outdir, str):
        raise ValueError
    if not isinstance(targetname, str):
        raise ValueError

    model_teff_grid = np.arange(4000, 7000+250, 250)
    use_teff = model_teff_grid[np.argmin(np.abs(teff - model_teff_grid))]

    # vsini: 0, 2, 5, 10, 30, 60.
    model_d = {}

    vsinis = [60, 30, 10, 5, 2, 0]
    for vsini in vsinis:
         m_wav, m_flx = read_nextgen(teff=use_teff, vsini=vsini,
                                     logg=logg, vturb=vturb)
         sel = (m_wav > np.nanmin(wav.value)) & (m_wav < np.nanmax(wav.value))
         if len(m_wav[sel])==0:
             return None
         model_d['wav_vsini{}'.format(vsini)] = m_wav[sel]
         model_d['flx_vsini{}'.format(vsini)] = m_flx[sel]

    # calculate chi^2 using the models as model spectra. will probably also
    # want to have exclude regions (per order?) where the fit is obviously
    # terrible.
    chi2_d = {}
    shifts = np.arange(-1, 1+0.05, 0.05)
    for vsini in vsinis:
        for shift in shifts:

            # interpolate the high-resolution model grid onto data grid
            m_wav_highres = model_d['wav_vsini{}'.format(vsini)]
            m_flx_highres = model_d['flx_vsini{}'.format(vsini)]
            interp_fn = interp1d(m_wav_highres+shift, m_flx_highres, kind='quadratic',
                                 bounds_error=False, fill_value=np.nan)
            m_wav_lowres = wav
            m_flx_lowres = interp_fn(m_wav_lowres)

            k = 'vsini{:d}_shift{:.2f}'.format(vsini, shift)
            if not flxerr is None:
                chi2_d[k] = np.nansum(
                    (flx - m_flx_lowres)**2 / (flxerr**2)
                ).value
            else:
                flxerr = 1e-2
                chi2_d[k] = np.nansum(
                    (flx - m_flx_lowres)**2 / (flxerr**2)
                ).value

    # now plot the preferred shift and template
    min_key = min(chi2_d, key=chi2_d.get)
    best_vsini = float(min_key.split('_')[0].split('vsini')[1])
    best_shift = float(min_key.split('shift')[1])
    best_chi2 = chi2_d[min_key]

    # make a stacked spectrum plot showing target spectrum on top, and the
    # candidate model spectra below. Nice for "chi by eye"
    outpath = join(
        outdir, targetname+'_vsini-stack-shift{:.2f}.png'.format(best_shift)
    )
    if not os.path.exists(outpath):
        plt.close('all')
        f, ax = plt.subplots(figsize=(10,7))
        trans = blended_transform_factory(ax.transAxes, ax.transData)
        bbox = dict(facecolor='white', edgecolor='none', alpha=0.7)
        n_spectra = len(vsinis) + 1

        offset = 0
        for i, vsini in enumerate(vsinis):
            ax.plot(model_d['wav_vsini{}'.format(vsini)]+best_shift,
                    offset+model_d['flx_vsini{}'.format(vsini)], lw=0.5, c='k')
            s = "Teff {} vsini {}\nshift {:.2f}".format(use_teff, vsini, best_shift)
            ax.text(0.01, 0.9+offset, s, bbox=bbox, transform=trans,
                    fontsize='xx-small')
            offset +=1

        ax.plot(wav, offset+flx, lw=0.5, c='C0')
        s = "{} (Teff {})".format(targetname, teff)
        ax.text(0.01, 0.9+offset, s, bbox=bbox, transform=trans,
                fontsize='xx-small')

        ax.text(0.97, 0.9+offset, '$\chi^2$ vsini: {}'.format(best_vsini),
                bbox=bbox, transform=trans, ha='right', va='top', fontsize='small')

        ax.set_xlabel('wavelength [angstrom]')
        ax.set_ylabel('relative flux')
        format_ax(ax)
        savefig(f, outpath, writepdf=False)
    else:
        print('found {}'.format(outpath))

    return best_vsini, best_shift, best_chi2


def get_synth_spectrum(synth_path):

    hl = fits.open(synth_path)
    dirname = os.path.dirname(synth_path)
    if "HiRes" in dirname:
        _hl = fits.open(
            os.path.join(dirname,
                         "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")
        )
        syn_wav = _hl[0].data
        _hl.close()
    elif "MedRes" in dirname:
        # downloaded the δλ = 1A model, fe/h=0, alpha/M=0
        syn_wav = np.arange(int(3e3), int(1e4), 0.1)
    else:
        NotImplementedError
    syn_flx = hl[0].data
    # NOTE: might want to actually continuum normalize...
    syn_norm_flx = syn_flx/np.nanmedian(syn_flx)
    hl.close()

    return syn_flx, syn_wav


def _compute_chi_sq(
    args: Tuple[int, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    ) -> Tuple[int, float, np.ndarray, np.ndarray]:
    """
    Computes chi squared between a shifted observed spectrum and
    the template spectrum, returning only primitive types/arrays.

    Args:
        args: (index, wav_shift, rv_shift, obs_wav, obs_flx, tpl_wav, tpl_flx)
            index (int): index of this shift
            wav_shift (float): wavelength shift in Angstrom
            rv_shift (float): radial velocity shift (not used inside, but included
                for clarity)
            obs_wav (np.ndarray): 1D array of the object wavelength
            obs_flx (np.ndarray): 1D array of the object flux
            tpl_wav (np.ndarray): 1D array of the template wavelength
            tpl_flx (np.ndarray): 1D array of the template flux

    Returns:
        (index, chi_sq, shifted_wav, shifted_flx)
    """
    (_ix, wav_shift, _rv_shift,
     obs_wav, obs_flx, tpl_wav, tpl_flx) = args

    # Shift the observation
    _wav = obs_wav + wav_shift

    # Make a selection so we only compare in region overlapping the template
    sel = (_wav > np.min(tpl_wav)) & (_wav < np.max(tpl_wav))
    _wav = _wav[sel]
    _flx = obs_flx[sel]

    # For demonstration: simple “uncertainty” as the scatter of obs
    _unc = np.ones_like(_flx) * np.std(_flx)

    # Interpolate the shifted obs onto the same grid as the template
    # NOTE: do a quick-and-dirty linear interpolation
    #       (you could do a more elaborate resampler if you want).
    interp_flx = np.interp(tpl_wav, _wav, _flx, left=np.nan, right=np.nan)
    interp_unc = np.interp(tpl_wav, _wav, _unc, left=np.nan, right=np.nan)

    # Compute chi^2
    # Only compare finite overlap region
    finite_sel = ~np.isnan(interp_flx) & ~np.isnan(tpl_flx)
    data = interp_flx[finite_sel]
    model = tpl_flx[finite_sel]
    unc = interp_unc[finite_sel]
    chi_sq = np.sum((data - model)**2 / (unc**2))

    # Return arrays and floats only
    return _ix, chi_sq, tpl_wav[finite_sel], data


def get_naive_rv(spectrum_path, synth_path, outdir, make_plot=1,
                 run_in_parallel=0):
    """
    Given a target HIRES spectrum `spectrum_path` and a PHOENIX model spectrum
    `synth_path` guess the RV using an order-by-order CCF.

    Example PHOENIX spectra at
    /Users/luke/Dropbox/proj/hd_34382/data/synthetic_spectra/PHOENIX_MedRes
    # downloaded the δλ = 1A model, fe/h=0, alpha/M=0
    """

    tname = os.path.basename(spectrum_path).replace(".fits", "")
    sname = os.path.basename(synth_path).replace(".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits", "")

    outdir = os.path.join(outdir, tname+"_v_"+sname)
    if not os.path.exists(outdir): os.mkdir(outdir)

    # get the relevant spectra
    flx_2d, norm_flx_2d, wav_2d, mjd, ra, dec = _get_full_hires_spectrum(spectrum_path)
    from astropy.coordinates import SkyCoord
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    ra, dec = coord.ra.deg, coord.dec.deg
    print(ra, dec)

    syn_flx, syn_wav = get_synth_spectrum(synth_path)

    # barycentric correction:
    # account for position and velocity of the Earth geocenter
    # with respect to the Solar System barycenter, and Earth's rotation, and
    # drop higher order terms.
    from barycorrpy import get_BC_vel
    jd = float(mjd) + 2400000.5
    epoch = 2451545.0
    rv = 0
    pmra = 0
    pmdec = 0
    px = 1
    bcresult = get_BC_vel(JDUTC=jd, obsname='Keck', ephemeris='de430',
                          zmeas=0.0, ra=ra, dec=dec, pmra=pmra, pmdec=pmdec,
                          px=px, rv=rv, predictive=True)
    bc_ms = bcresult[0] # barycentric correction in m/s
    bc_kms = bc_ms / 1e3


    # iterate over orders and fit for the RV
    n_order = flx_2d.shape[0]
    orders, rvs_ccf, rvs_chisq = [], [], []

    outcsvpath = os.path.join(outdir, f'rv_{tname}_{sname}.csv')
    if os.path.exists(outcsvpath):
        print(f"Found {outcsvpath}, return.")
        return pd.read_csv(outcsvpath)

    print(f"naive_rv: Will write to {outdir}")

    #good_orders = [1,2,6,19]
    all_orders = range(n_order)
    iterable = all_orders # NOTE can tweak
    for ix in iterable:

        wav = wav_2d[ix, :]
        flx = norm_flx_2d[ix, :] / np.nanmedian(norm_flx_2d[ix, :])
        obs_sel = np.isfinite(wav) & np.isfinite(flx)
        wav, flx = wav[obs_sel], flx[obs_sel]
        min_wav, max_wav = np.nanmin(wav), np.nanmax(wav)
        wav_0 = np.nanmean(wav)*u.AA

        syn_sel = (syn_wav > min_wav) & (syn_wav < max_wav)
        sflx = syn_flx[syn_sel] / np.nanmedian(syn_flx[syn_sel])
        swav = syn_wav[syn_sel]

        # define object spectrum, smooth, and continuum normalize

        # smooth
        from scipy.ndimage import gaussian_filter1d
        fn = lambda x: gaussian_filter1d(x, sigma=1)
        flx = fn(flx)

        o_unc = StdDevUncertainty(
            np.ones_like(flx)*p2p_rms(flx)*u.dimensionless_unscaled
        )

        o_spec = Spectrum1D(spectral_axis = wav*u.AA,
                            flux = fn(flx)*u.dimensionless_unscaled,
                            uncertainty = o_unc,
                            velocity_convention = 'doppler_optical',
                            rest_value = wav_0)

        # continuum normalize
        cont_flx = fit_generic_continuum(o_spec)(o_spec.spectral_axis)
        cont_norm_spec = o_spec / cont_flx
        std_data = np.nanstd(cont_norm_spec.flux)
        o_spec = Spectrum1D(spectral_axis = wav*u.AA,
                            flux = cont_norm_spec.flux,
                            uncertainty = o_unc,
                            velocity_convention = 'doppler_optical',
                            rest_value = wav_0)

        # template spectrum (tiny uncertainty)
        t_unc = StdDevUncertainty(
            1e-3*np.ones_like(sflx)*p2p_rms(sflx)*u.dimensionless_unscaled
        )
        t_spec = Spectrum1D(spectral_axis = swav*u.AA,
                            flux = sflx*u.dimensionless_unscaled,
                            uncertainty = t_unc)

        # continuum normalize, and match stdevn of data
        cont_flx = fit_generic_continuum(t_spec)(t_spec.spectral_axis)
        cont_norm_spec = t_spec / cont_flx
        std_model = np.nanstd(cont_norm_spec.flux)
        fudge = std_data/std_model
        t_spec = Spectrum1D(spectral_axis = swav*u.AA,
                            flux = fudge*cont_norm_spec.flux,
                            uncertainty = t_unc)

        # object spectrum resampled onto the template wavelength scale
        spl = SplineInterpolatedResampler()
        ospl_spec = spl(o_spec, swav*u.AA)

        # degrade object spectrum resolution to ~1Å (FWHM=1Å)
        dw = np.median(np.diff(swav))
        sigma_pix = 1.0 / (2.355 * dw)  # convert 1Å FWHM into pixels
        print(f'Degrading target by {sigma_pix:.1f} px to 1angstrom resolution')
        degraded_flux = gaussian_filter1d(ospl_spec.flux.value, sigma_pix)
        ospl_spec = Spectrum1D(
            spectral_axis=ospl_spec.spectral_axis,
            flux=degraded_flux*u.dimensionless_unscaled,
            uncertainty=ospl_spec.uncertainty,
            velocity_convention=ospl_spec.velocity_convention,
            rest_value=ospl_spec.rest_value
        )

        # template vs object cross-correlation
        # NOTE default 0.5 apodization_window.. what if 1? should this number
        # matter?
        corr, lag = correlation.template_correlate(
            ospl_spec, t_spec, apodization_window=0.5
        )

        # negative RV: blueshifted, star coming towards you
        # ccf-based RV
        rv_ccf = lag[np.argmax(corr)]

        equiv = u.doppler_optical(wav_0)
        ccf_shifted_mean = rv_ccf.to(u.AA, equivalencies=equiv)
        ccf_shift = wav_0 - ccf_shifted_mean # positive = redshift

        ospl_shift_spec = Spectrum1D(
            spectral_axis = ospl_spec.wavelength + ccf_shift,
            flux = ospl_spec.flux
        )

        # inverse approach.  try a grid of shifts... and chi^2 minimize.
        # +/- ~240km/s, or +/- 3 angstrom

        # Serial implementation has been tested to be good to ~5km/s, but SLOW
        #
        # 100 points ->  5km/s grid resolution...
        # 1000 points ->  0.5km/s grid resolution...
        # 10000 points ->  0.05km/s grid resolution...
        #
        # Parallel implementation running...

        N_grid = int(1e4)

        if not run_in_parallel:

            wav_shift_grid = np.linspace(-4, 4, N_grid)
            rv_shift_grid = (wav_0 + wav_shift_grid*u.AA).to(u.km/u.s, equivalencies=equiv)
            chi_sqs = []

            N_tot = len(wav_shift_grid)

            for _ix, wav_shift, rv_shift in (
                zip(range(N_tot), wav_shift_grid, rv_shift_grid)
            ):

                # redefine data to ensure correct trimming
                _wav = ospl_spec.wavelength + wav_shift*u.AA
                _flx = ospl_spec.flux

                sel = (
                    (_wav > min(t_spec.wavelength)) &
                    (_wav < max(t_spec.wavelength))
                )
                _wav, _flx = _wav[sel], _flx[sel]
                _unc = np.ones_like(_flx) * p2p_rms(_flx)

                # resample, since the default "shift" produces odd pixelation
                # effects
                o_unc = StdDevUncertainty(
                     _unc * u.dimensionless_unscaled
                )
                ospl_shift_spec = linear(ospl_spec, _wav)

                # define model
                # sum((data-model)^2/unc^2)
                data = ospl_shift_spec.flux
                unc = ospl_shift_spec.uncertainty.array
                _sel = (
                    (t_spec.wavelength > min(ospl_shift_spec.wavelength))
                    &
                    (t_spec.wavelength > max(ospl_shift_spec.wavelength))
                )
                model = t_spec.flux[sel]

                chi_sq = np.sum((data-model)**2/(unc**2))

                chi_sqs.append(chi_sq)

                if _ix == 0:
                    best_shift_spec = ospl_shift_spec
                else:
                    if chi_sq < np.nanmin(chi_sqs):
                        best_shift_spec = ospl_shift_spec

                if _ix % 1000 == 0:
                    print(f"{_ix}/{N_tot}, {wav_shift:.4f}, {rv_shift:.2f}, χ^2={chi_sq:.1f}...")

            rv_chisq = rv_shift_grid[np.argmin(np.array(chi_sqs))]
            shift_chisq = wav_shift_grid[np.argmin(np.array(chi_sqs))]

            orders.append(ix)
            rvs_chisq.append(rv_chisq.value)
            rvs_ccf.append(rv_ccf.value)

        else:

            # Explanation: Why we pass only NumPy arrays into multiprocessing
            #
            # Python’s multiprocessing uses “pickling” to send data from the main
            # process to worker processes. NumPy arrays are straightforward to pickle
            # because NumPy provides a well-defined interface for serializing array
            # metadata (shape, strides, dtype) and the raw memory buffer. Hence, it’s
            # easy for Python to re-create the same array in each worker process.
            #
            # On the other hand, complex objects like Spectrum1D may contain internal
            # references to C/C++ code, custom data structures, or Python-callable
            # attributes that are not pickle-friendly. These objects cannot be
            # automatically decomposed and reassembled by the pickle module. Thus,
            # the workaround is to pass only primitive data (arrays, floats, etc.)
            # to multiprocessing, and then rebuild or reconstruct any specialized
            # objects in the parent process as needed.

            ###################################################
            # 1) Extract your arrays from the Spectrum1D objects
            ###################################################
            #   E.g. “ospl_spec” might be your object spectrum,
            #   while “t_spec” is the template.  Convert them
            #   into arrays that can be passed to multiprocessing:
            ###################################################
            obs_wav_arr = ospl_spec.wavelength.value  # 1D numpy array
            obs_flx_arr = ospl_spec.flux.value
            tpl_wav_arr = t_spec.wavelength.value
            tpl_flx_arr = t_spec.flux.value

            ###################################################
            # 2) Define your shift grids as normal
            ###################################################
            wav_shift_grid = np.linspace(-4, 4, N_grid)
            rv_shift_grid = (
                (wav_0 + wav_shift_grid*u.AA)
                .to(u.km/u.s, equivalencies=equiv)
            )
            N_tot = len(wav_shift_grid)

            from multiprocessing import Pool

            args_list = [
                (i,
                 wav_shift_grid[i],
                 rv_shift_grid[i].value,     # pass as float
                 obs_wav_arr,
                 obs_flx_arr,
                 tpl_wav_arr,
                 tpl_flx_arr)
                for i in range(N_tot)
            ]

            chi_sqs = []
            min_chi = np.inf
            best_result = None

            # Start pool
            with Pool() as pool:
                results = pool.map(_compute_chi_sq, args_list)

            # results is a list of tuples: (_ix, chi_sq, shifted_wav, shifted_flx)
            for _ix, chi_sq, shifted_wav, shifted_flx in results:
                chi_sqs.append(chi_sq)
                if chi_sq < min_chi:
                    min_chi = chi_sq
                    best_result = (shifted_wav, shifted_flx)

                if _ix % 1000 == 0:
                    print(f"{_ix}/{N_tot}, {wav_shift_grid[_ix]:.4f}, "
                          f"{rv_shift_grid[_ix]:.2f} km/s, chi^2={chi_sq:.1f}...")

            # pick best solution
            best_ix = np.argmin(chi_sqs)
            rv_chisq = rv_shift_grid[best_ix]
            shift_chisq = wav_shift_grid[best_ix]

            # Rebuild a Spectrum1D only for the best shift
            best_wav, best_flx = best_result
            best_unc = StdDevUncertainty(np.std(best_flx)*np.ones_like(best_flx))
            best_shift_spec = Spectrum1D(
                spectral_axis=best_wav*u.AA,
                flux=best_flx*u.dimensionless_unscaled,
                uncertainty=best_unc
            )

            orders.append(ix)
            rvs_chisq.append(-1. * rv_chisq.value) # not sure why, but needed.
            rvs_ccf.append(rv_ccf.value)


        if make_plot:
            set_style('science')
            plt.close("all")
            fig = plt.figure(figsize=(16, 4))
            axd = fig.subplot_mosaic(
                """
                AAAAAABB
                AAAAAACC
                """
            )

            ax = axd['A']
            ax.plot(best_shift_spec.spectral_axis,#-shift_chisq*u.AA,
                    best_shift_spec.flux,
                    lw=0.5, label=tname + " Δλ=1A best-shift")

            dflux = 0.3*(
                np.nanmax(best_shift_spec.flux) -
                np.nanmin(best_shift_spec.flux)
            )

            ax.plot(ospl_spec.spectral_axis,
                    ospl_spec.flux - dflux,
                    lw=0.5, label=tname + " Δλ=1A observed")

            ax.plot(t_spec.spectral_axis, t_spec.flux - 2*dflux, lw=0.5,
                    label=sname)

            tstr = (
                f'order {ix}: RVccf={rv_ccf:.2f} RVχ2={rv_chisq:.2f}, '
                f'datashift (χ2) λ={shift_chisq:.3f}'
            )
            ax.set_title(tstr)
            ax.grid(True, which='both', ls='--', lw=0.3)
            ax.set_xlabel('Wavelength [Angstrom]')
            ax.set_ylabel('Relative flux')

            ax.legend(loc='lower right', fontsize='x-small')

            ax = axd['B']
            ax.plot(lag, corr)
            ax.set_xlim([-250,250])
            ax.set_xlabel('Velocity [km/s]')
            ax.set_ylabel('CCF')
            ax.grid(True, which='both', ls='--', lw=0.3)

            ax = axd['C']
            ax.plot(rv_shift_grid, chi_sqs)
            ax.set_xlim([-250,250])
            ax.set_xlabel('Velocity [km/s]')
            ax.set_ylabel('χ^2')
            ax.grid(True, which='both', ls='--', lw=0.3)

            fig.tight_layout()
            ostr = str(ix).zfill(2)
            outpath = os.path.join(
                outdir, f"shiftcompare_ord{ostr}_{tname}_{sname}.png"
            )
            savefig(fig, outpath, writepdf=False, dpi=300)

        DEBUG = 0

        if DEBUG:
            plt.close("all")
            fig, ax = plt.subplots()
            ax.plot(o_spec.spectral_axis, o_spec.flux)
            ax.plot(ospl_spec.spectral_axis, ospl_spec.flux)
            plt.savefig("temp_resample.png", dpi=300)

        if DEBUG:
            plt.close("all")
            fig, ax = plt.subplots()
            ax.plot(lag, corr)
            ax.set_xlim([-300,300])
            plt.savefig("temp_correlate.png", dpi=300)

    outdf = pd.DataFrame({
        'order': orders,
        'rv_ccf_kms': np.round(rvs_ccf, 4),
        'rv_chisq_kms': np.round(rvs_chisq, 4),
        'bc_kms': np.round(np.ones(len(orders))*bc_kms, 4),
    })

    # note: the correction is defined s.t. you subtract it off the measured
    # redshift.
    outdf['rv_chisq_minus_bc_kms'] = np.round(
        outdf['rv_chisq_kms'] - outdf['bc_kms'], 4
    )

    rvs = np.array(outdf['rv_chisq_minus_bc_kms'])

    rchip_good_orders = [0,2,3,4,5,6,11] # previous: [0,2,3,4,5,6,10,11,13]
    sel_rvs = rvs[np.array(rchip_good_orders)]

    outdf['meangoodorder_rv_chisq_minus_bc_kms'] = np.round(np.nanmean(sel_rvs), 4)

    outdf.to_csv(outcsvpath, index=False)
    print(f"wrote {outcsvpath}")

    return pd.read_csv(outcsvpath)



def measure_veloce_vsini(specname, targetname, teff, outdir):
    """
    wrapper to measure_vsini, which uses a chi-squared grid. returns mean vsini
    over selected orders, and mean wavelength shift over the same.
    """

    spectrum_path = join(
        DATADIR, 'Veloce', specname
    )

    flx_2d, wav_2d, flxerr_2d = read_veloce(spectrum_path, start=400, end=-600,
                                            return_err=True)

    o_vsini, o_shift, o_chi2 = [], [], []
    for o in VELOCE_VSINI_ORDERS_VS_NEXTGEN:

        flx, wav, flxerr = flx_2d[o, :], wav_2d[o, :], flxerr_2d[o, :]
        wav, flx, flxerr = wav[::-1], flx[::-1], flxerr[::-1]
        cont_flx, cont_norm_spec = fit_continuum(flx, wav, instrument='Veloce')

        sel = (
            (cont_flx > 0) &
            (cont_norm_spec.flux < 2) &
            (cont_norm_spec.flux > -1)
        )

        flxerr = flxerr[sel] / cont_flx[sel]
        wav, flx = cont_norm_spec.wavelength[sel], cont_norm_spec.flux[sel]

        res = measure_vsini(wav, flx, flxerr=flxerr, teff=teff, logg=4.5,
                            vturb=2, outdir=outdir,
                            targetname=targetname+'_order{}'.format(o))

        if res is None:
            msg = (
                'WRN! Order {} failed to get template and target overlap. '
                'Skipping!'.format(o)
            )
            print(msg)

        fit_vsini, fit_shift, fit_chi2 = res
        o_vsini.append(fit_vsini)
        o_shift.append(fit_shift)
        o_chi2.append(fit_chi2)

    out_df = pd.DataFrame({
        'order': VELOCE_VSINI_ORDERS_VS_NEXTGEN,
        'vsini': o_vsini,
        'shift': o_shift,
        'chi2': o_chi2
    })
    outpath = join(outdir, specname.replace('.fits','_vsini_fit.csv'))
    out_df.to_csv(outpath, index=False)
    print('made {}'.format(outpath))

    gamma = given_deltawvlen_get_vsys(deltawvlen=np.mean(o_shift)*u.AA,
                                      wvlen_0=6500*u.AA, verbose=False)

    print('Average vsini: {}km/s'.format(np.mean(o_vsini)))
    print('Average shift: {:.2f}A'.format(np.mean(o_shift)))
    print('Average shift: {}'.format(gamma.to(u.km/u.s)))

    return np.mean(o_vsini), np.mean(o_shift), gamma


######################
# specmatch analysis #
######################

def specmatch_analyze(spectrum_path, wvsol_path=None, region=None, outdir=None,
                      idstring=None, is_template=False):

    if 'PFS' in spectrum_path:
        instrument = 'PFS'
    elif 'Veloce' in spectrum_path:
        instrument = 'Veloce'
    else:
        raise NotImplementedError('check wavlim works for non-pfs spectrum')

    for s in [outdir, region, idstring]:
        if not isinstance(s, str):
            raise ValueError
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path,
                                  is_template=is_template)
    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
    else:
        raise NotImplementedError

    #
    # define target region
    #

    if region == 'Mgb1' and instrument=='PFS':
        target_wav = 5183.62
        wavlim = [5160,5210]
    elif region == 'Mgb1' and instrument!='PFS':
        raise NotImplementedError('veloce doesnt cover Mg b1')

    if 'order' in region and instrument=='PFS':
        ix = int(region.split('order')[-1])
        target_wav = np.nanmedian(wav_2d[ix, :])
        wavlim = [target_wav-21, target_wav+21]
    elif 'order' in region and instrument!='PFS':
        raise NotImplementedError()

    if region == '6300' and instrument=='Veloce':
        target_wav = 6300
        wavlim = [6280,6350]
    elif region == '6300' and instrument!='Veloce':
        raise NotImplementedError(
            'does this region work on {}?'.format(instrument)
        )

    #
    # first, find orders that contain the target wavelength. select the best
    # order (the one for which the target wavelength is closest to the center
    # of the order). then trim the 1d spectrum from that order
    #
    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    viable_orders = np.argwhere(
        (_preorder != wav_2d.shape[1]-1) & (_preorder != 0)
    )
    order = int(
        viable_orders[np.argmin(
            np.abs(_preorder[viable_orders] - wav_2d.shape[1]/2)
        )]
    )

    flx, wav = flx_2d[order, :], wav_2d[order, :]
    sel = (wav > wavlim[0]) & (wav < wavlim[1])
    flx, wav = flx[sel], wav[sel]

    if instrument == 'Veloce':
        # the australians drive on the other side of the road?
        wav = wav[::-1]
        flx = flx[::-1]

    #
    # continuum normalize, and then check you did it ok.
    #
    spec = Spectrum1D(spectral_axis=wav*u.AA,
                      flux=flx*u.dimensionless_unscaled)

    # avoid manual exclude regions if possible
    exclude_regions = []
    if len(exclude_regions)>=1:
        cont_flx = (
            fit_generic_continuum(spec,
                                  exclude_regions=exclude_regions
            )(spec.spectral_axis)
        )
    else:
        cont_flx = (
            fit_generic_continuum(spec)(spec.spectral_axis)
        )

    cont_norm_spec = spec / cont_flx

    outpath = join(outdir, f'{idstring}_{region}_cont_norm_check.png')

    f,axs = plt.subplots(nrows=2, ncols=1, figsize=(6,4))
    axs[0].plot(wav, flx, c='k', zorder=3, lw=0.5)
    axs[0].plot(wav, cont_flx, c='r', zorder=2, lw=0.5)
    axs[1].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k', lw=0.5)

    axs[0].set_ylabel('flux')
    axs[1].set_ylabel('contnorm flux')

    axs[-1].set_xlabel('wavelength [angstrom]')
    for ax in axs:
        format_ax(ax)
    savefig(f, outpath, writepdf=0)
    plt.close('all')

    flx = cont_norm_spec.flux

    #
    # shift and cross-correlate w/ specmatch
    #
    import specmatchemp.library
    lib = specmatchemp.library.read_hdf(wavlim=wavlim)

    from specmatchemp.spectrum import Spectrum
    s_spectrum = Spectrum(wav, flx)
    s_spectrum.name = idstring

    from specmatchemp.specmatch import SpecMatch
    sm_res = SpecMatch(s_spectrum, lib)
    sm_res.shift()

    try:
        match_row = lib.library_params[
            (lib.library_params.cps_name == sm_res.shift_ref.name)
        ]
        print(match_row)
        match_name = match_row.source_name.iloc[0]
    except IndexError:
        match_name = sm_res.shift_ref.name

    outpath = join(outdir, f'{idstring}_{region}_shift_check.png')
    fig = plt.figure(figsize=(10,5))
    sm_res.target_unshifted.plot(normalize=True, plt_kw={'color':'forestgreen'}, text='Target (unshifted)')
    sm_res.target.plot(offset=1.0, plt_kw={'color':'royalblue'}, text='Target (shifted): {}'.format(idstring))
    sm_res.shift_ref.plot(offset=2.0, plt_kw={'color':'firebrick'}, text='Reference: '+match_name)
    plt.xlim((wavlim[0],wavlim[1]))
    plt.ylim(0,3.0)
    ax = plt.gca()
    format_ax(ax)
    savefig(fig, outpath, writepdf=0)
    plt.close('all')

    wavlimshifted = (min(sm_res.target.w), max(sm_res.target.w))

    #
    # cross-correlate against the templates to fit for vsini, Rstar, FeH.
    #
    sm_res.match(wavlim=wavlimshifted)

    # Plot chi-squared surfaces
    outpath =  join(outdir, f'{idstring}_{region}_chisq.png')

    fig = plt.figure(figsize=(12, 8))
    sm_res.plot_chi_squared_surface()
    ax = plt.gca()
    format_ax(ax)
    savefig(fig, outpath, writepdf=0)
    plt.close('all')

    sm_res.lincomb()

    print('Derived Parameters: ')
    print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
        sm_res.results['Teff'], sm_res.results['radius'], sm_res.results['feh'])
    )

    #
    # make a plot comparing you spectrum to dwarf star spectra of comparable
    # Teff
    #
    plot_spec_vs_dwarf_library(
        wavlim,
        sm_res.results['Teff'],
        outdir,
        f'{idstring}_{region}',
        sm_res=sm_res
    )

