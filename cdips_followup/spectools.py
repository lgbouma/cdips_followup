"""
READ:
    read_nextgen: Read NextGen Model Atmosphere Grid from Hauschildt et al 1999.
    read_veloce: Read Veloce FITS file.
    read_fies: Read FIES FITS file.
    read_tres: Read TRES FITS file.
    read_hires: Read HIRES FITS file.
    read_pfs: Read PFS IDL SAV file.
    read_feros: Read FEROS FITS file.

SPECMATCH WRAPPERS:
    specmatch_analyze: shift+cross-correlate to get vsini, Rstar, FeH w/ SME.
    specmatchsyn_analyze: ditto, but using specmatch-synth.

CALCULATE:
    fit_continuum: fits a quadratic polynomial across the order.
    get_Li_6708_EW: measure Li EW.
    measure_veloce_vsini: wraps below.
        measure_vsini: compares target spectrum to broadend NextGen spectra
    given_deltawvlen_get_vsys: convert Δλ to velocity.

VISUALIZE:
    viz_1d_spectrum: flux vs wavelength (with select lins underplotted).
    plot_orders
    inspect_pfs
    specmatch_viz_compare: SME comparison of target and HIRES library spectra.
    plot_spec_vs_dwarf_library: SME ditto.

WIP:
    get_Ca_HK_emission: in theory, for CaHK emission line widths.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
from matplotlib.transforms import blended_transform_factory

from astropy.io import fits
from astropy import units as u, constants as const
from astropy.modeling.polynomial import Chebyshev1D
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter

from scipy.io import readsav
from scipy.interpolate import interp1d

import os
from copy import deepcopy

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum, fit_lines
from specutils.analysis import equivalent_width, centroid
from astropy.modeling import models

from specmatchemp.spectrum import Spectrum
from specmatchemp.specmatch import SpecMatch
import specmatchemp.library
import specmatchemp.plots as smplot

from aesthetic.plot import savefig, format_ax

from cdips_followup import __path__

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

# a few common lines
line_d = {
    'Fe_L': 3820.44,
    'Ca_K': 3933.66,
    'Ca_H': 3968.47,
    'H$\delta$': 4101.75,
    'H$\gamma$': 4340.47,
    r'H$\beta$': 4861.35,
    'FeI': 4920.51,
    'FeI': 4957.61,
    'FeI': 5169.00, #rough
    'Mgb1': 5183.62,
    'Mgb2': 5172.70,
    'Feb3': 5168.91,
    'Mgb4': 5167.33,
    'CaI+FeI': 5269,
    'Fe_E2': 5270.39,
    'Hg_e': 5460.73,
    'He': 5875.618,
    'NaI_D2': 5889.95,
    'NaI_D1': 5895.92,
    'FeI_a': 6703.58,
    'FeI_b': 6705.1,
    'FeI_c': 6707.44,
    'LiI_a': 6707.76,
    'LiI_b': 6707.91,
    'CaI$\lambda$': 6718,
    'Halpha': 6562.8,
    'CaI_IRT_a': 8498,
    'CaI_IRT_b': 8542,
    'CaI_IRT_c': 8662
}

# directories
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = os.path.join(os.path.dirname(__path__[0]), 'results/spec_analysis')
TESTOUTDIR = os.path.join(OUTDIR, 'tests')
if not os.path.exists(TESTOUTDIR):
    os.mkdir(TESTOUTDIR)

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


    NEXTGENDIR = os.path.join(
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

    spectrum_path = os.path.join(
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

    if return_err:
        return flux, wav, flux_err
    else:
        return flux, wav


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
    return wav, flx


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

def viz_1d_spectrum(flx, wav, outpath, xlim=None, vlines=None, names=None):

    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        wav = wav[sel]
        flx = flx[sel]

    plt.close('all')
    f,ax = plt.subplots(figsize=(10,3))
    ax.plot(wav, flx, c='k', zorder=3, lw=0.2)

    y_90 = np.nanpercentile(flx, 90)
    y_10 = np.nanpercentile(flx, 10)
    y_median = np.nanmedian(flx)
    y_diff = y_90 - y_10

    # ax.set_ylim( (y_median-1.1*y_diff, y_median+1.1*y_diff) )

    ax.set_xlabel('wavelength [angstrom]')
    ax.set_ylabel('flux [e-]')

    xmin = min(wav)
    xmax = max(wav)

    this_d = {}
    for k, v in line_d.items():
        if v > xmin and v<xmax:
            this_d[k] = v

    if len(this_d) > 0:
        for k, v in this_d.items():
            ylim = ax.get_ylim()
            delta_y = 0.9*(max(ylim) - min(ylim))
            ax.vlines(v, min(ylim)+delta_y, max(ylim), zorder=-3,
                      linestyles=':', color='k', lw=0.3)
            ax.set_ylim(ylim)

            tform = blended_transform_factory(ax.transData, ax.transAxes)
            ax.text(v, 0.95, k, ha='center', va='top', transform=tform,
                    fontsize=4)

    if isinstance(xlim, (list, tuple)):
        ax.set_xlim(xlim)

    if isinstance(vlines, list):
        sel = (nparr(vlines)>min(xlim)) & (nparr(vlines)<max(xlim))
        vlines, names = nparr(vlines)[sel], nparr(names)[sel]
        ylim = ax.get_ylim()
        delta_y = 0.9*(max(ylim) - min(ylim))
        ax.vlines(vlines, min(ylim)+delta_y, max(ylim), zorder=-3,
                  linestyles=':', color='k', lw=0.3)
        ax.set_ylim(ylim)

        tform = blended_transform_factory(ax.transData, ax.transAxes)
        for x, n in zip(vlines, names):
            ax.text(x, 0.95, n, ha='center', va='top', transform=tform,
                    fontsize=4)

    ax.grid(b=True, which='both', axis='x', zorder=-3, lw=0.3, ls='--')

    format_ax(ax)
    savefig(f, outpath, writepdf=False)
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
            raise NotImplementedError
    elif 'TRES' in spectrum_path:
        flx_2d, wav_2d = read_tres(spectrum_path)
    else:
        raise NotImplementedError

    for order in range(wav_2d.shape[0]):

        if 'Veloce' in spectrum_path:
            start = 200
            end = -200
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
        outpath = os.path.join(outdir, outname)
        viz_1d_spectrum(flx, wav, outpath)

        if isinstance(flat_path, str):
            outname = f'{idstring}_order{str(order).zfill(2)}_flat.png'
            outpath = os.path.join(outdir, outname)
            viz_1d_spectrum(flat, wav, outpath)


def inspect_pfs(nightstr, targetline, xlim=None):
    # NIST has the Li I resonance doublet listed with one transition at 6707.76
    # and the other at 6707.91 A.

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

    outpath = os.path.join(outdir, '{}_compare_speclib.png'.format(idstring))
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




def get_Li_6708_EW(spectrum_path, wvsol_path=None, xshift=None, delta_wav=5,
                   outpath=None, is_template=False):
    """
    spectrum_path: path to PFS, Veloce, or TRES spectrum

    wvsol_path: path to PFS wavelength solution (optional)

    xshift: angstrom shift required to get into source frame (not vacuum frame).

    delta_wav: window to do the measurement over (angstrom)

    outpath: summary figure is written here.
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
    else:
        raise NotImplementedError

    # FeI at 6703.58 and 6705.1 (Berger+18 Fig3).
    # FeI at 6707.44
    # Li I doublet at 6707.76, 6707.91
    # CaI lambda at ~6718.
    target_wav = 6707.85
    vlines = [6703.58, 6705.1, 6707.44, 6707.76, 6707.91, 6718]
    names = ['FeI', 'FeI', 'FeI', 'Li', '', 'CaI$\lambda$']
    xlim = [target_wav-delta_wav, target_wav+delta_wav]

    #
    # retrieve the order corresponding to target wavelength.
    # then shift the wavelength solution to source frame, if needed.
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

    if instrument == 'Veloce':
        wav = wav[::-1]
        flx = flx[::-1]
        thispath = outpath.replace('.png', '_1d_check.png')
        viz_1d_spectrum(flx, wav, thispath, xlim=(6700, 6725), vlines=vlines,
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

    #
    # get the Li EW
    #
    region = SpectralRegion((target_wav-0.5)*u.AA, (target_wav+0.5)*u.AA)
    li_equiv_width = equivalent_width(cont_norm_spec, regions=region)
    li_centroid = centroid(full_spec, region)

    #
    # fit a gaussian too, and get ITS equiv width
    # https://specutils.readthedocs.io/en/stable/fitting.html
    #
    g_init = models.Gaussian1D(amplitude=0.2*u.dimensionless_unscaled,
                               mean=target_wav*u.AA, stddev=0.5*u.AA)
    g_fit = fit_lines(full_spec, g_init, window=(region.lower, region.upper))
    y_fit = g_fit(full_spec.wavelength)

    fitted_spec = Spectrum1D(spectral_axis=full_spec.wavelength,
                             flux=(1-y_fit)*u.dimensionless_unscaled)
    fitted_li_equiv_width = equivalent_width(fitted_spec, regions=region)

    #
    # print bestfit params
    #
    print(42*'=')
    print('got Li equiv width of {}'.format(li_equiv_width))
    print('got fitted Li equiv width of {}'.format(fitted_li_equiv_width))
    print('got Li centroid of {}'.format(li_centroid))
    print('fit gaussian1d params are\n{}'.format(repr(g_fit)))
    print(42*'=')

    #
    # plot the results
    #
    f,axs = plt.subplots(nrows=4, ncols=1, figsize=(6,8))

    axs[0].plot(wav, flx, c='k', zorder=3)
    axs[0].plot(wav, cont_flx, c='r', zorder=2)

    axs[1].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    axs[2].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k')

    axs[3].plot(full_spec.wavelength, full_spec.flux, c='k')
    axs[3].plot(full_spec.wavelength, y_fit, c='g')

    txt = (
        'gaussian1d\namplitude:{:.3f}\nmean:{:.3f}\nstd:{:.3f}\nEW:{:.3f}'.
        format(g_fit.amplitude.value,
               g_fit.mean.value,
               g_fit.stddev.value,
               fitted_li_equiv_width)
    )
    axs[3].text(
        0.95, 0.95, txt, ha='right', va='top', transform=axs[3].transAxes,
        fontsize='xx-small'
    )

    axs[0].set_ylabel('flux')
    axs[1].set_ylabel('contnorm flux')
    axs[2].set_ylabel('contnorm flux [zoom]')
    axs[3].set_ylabel('1 - (contnorm flux)')

    if isinstance(xlim, list):
        for ax in axs:
            ax.set_xlim(xlim)

    axs[2].set_xlim([target_wav-1, target_wav+1])
    axs[3].set_xlim([target_wav-1, target_wav+1])
    axs[-1].set_xlabel('wavelength [angstrom]')

    for ax in axs:
        format_ax(ax)

    savefig(f, outpath)


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
    outpath = os.path.join(
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


def measure_veloce_vsini(specname, targetname, teff, outdir):
    """
    wrapper to measure_vsini, which uses a chi-squared grid. returns mean vsini
    over selected orders, and mean wavelength shift over the same.
    """

    spectrum_path = os.path.join(
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
    outpath = os.path.join(outdir, specname.replace('.fits','_vsini_fit.csv'))
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

    outpath = os.path.join(outdir, f'{idstring}_{region}_cont_norm_check.png')

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
    lib = specmatchemp.library.read_hdf(wavlim=wavlim)

    s_spectrum = Spectrum(wav, flx)
    s_spectrum.name = idstring

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

    outpath = os.path.join(outdir, f'{idstring}_{region}_shift_check.png')
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
    outpath =  os.path.join(outdir, f'{idstring}_{region}_chisq.png')

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
