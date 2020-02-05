from numpy import array as nparr
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u, constants as const
from scipy.io import readsav
from matplotlib.transforms import blended_transform_factory
from copy import deepcopy

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum, fit_lines
from specutils.analysis import equivalent_width, centroid
from astropy.modeling import models

from stringcheese.plotutils import savefig, format_ax

########
# read #
########
def read_veloce(spectrum_path, start=0, end=None):
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

    return flux, wav


def read_pfs(spectrum_path, wvlen_soln, verbose=False):
    """
    Read PFS IDL SAV file.
    Return (73, 3520) arrays of flux and wavelength.
    73 orders, 3520 pixels in the cross-dispersion direction.
    """

    sp = readsav(spectrum_path, python_dict=True, verbose=verbose)
    wvlen = readsav(wvlen_soln, python_dict=True, verbose=verbose)

    return sp['sp'], wvlen['w']


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
    f,ax = plt.subplots(figsize=(4,3))
    ax.plot(wav, flx, c='k', zorder=3, lw=0.2)

    ax.set_xlabel('wavelength [angstrom]')
    ax.set_ylabel('flux [e-]')

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

    format_ax(ax)
    savefig(f, outpath, writepdf=False)
    plt.close()


def viz_compare(wavlim=[6700,6725]):
    # NOTE: doesnt work, b/c the install fails

    import specmatchemp.library
    import specmatchemp.plots as smplot

    lib = specmatchemp.library.read_hdf(wavlim=wavlim)

    cut = lib.library_params.query('radius < 1.5 and -0.25 < feh < 0.25')
    g = cut.groupby(pd.cut(cut.Teff,bins=arange(5000,7000,250)))
    cut = g.first()

    fig = plt.figure()
    plt.plot(lib.library_params.Teff, lib.library_params.radius,'b.',
             label='_nolegend_')
    plt.plot(cut.Teff, cut.radius,'ro', label='Selected Stars')
    plt.legend()
    smplot.label_axes('Teff','radius')
    wvstr = '{}-{}'.format(wavlim[0], wavlim[1])
    fig.savefig(
        '../results/spec_viz/quickstart-library-selected-stars_{}.png'.  wvstr
    )


##############
# measure EW #
##############
def get_Li_6708_EW(spectrum_path, wvsol_path=None, xshift=None, delta_wav=5,
                   outpath=None):
    """
    xshift: angstrom shift required to get into source frame (not vacuum frame).

    delta_wav: window to do the measurement over (angstrom)
    """

    if not isinstance(outpath, str):
        raise ValueError

    flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)

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
    order = int(np.argwhere((_preorder != wav_2d.shape[1]-1) & (_preorder != 0)))

    flx, wav = flx_2d[order, :], wav_2d[order, :]
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
