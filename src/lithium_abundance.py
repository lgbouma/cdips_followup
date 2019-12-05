"""
Measure lithium abundances.
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum, fit_lines
from specutils.analysis import equivalent_width, centroid
from astropy.modeling import models

from stringcheese.plotutils import savefig, format_ax

def plot_feros_spectra(spectrum_path, outpath=None, xlim=None):
    hdul = fits.open(spectrum_path)

    d = hdul[0].data

    wav = d[0,0]
    flx = d[3,0]

    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        wav = wav[sel]
        flx = flx[sel]

    spec = Spectrum1D(spectral_axis=wav*u.AA, flux=flx*u.dimensionless_unscaled)

    f,ax = plt.subplots()
    ax.plot(wav, flx, c='k', zorder=3)

    if isinstance(xlim, list):
        exclude_regions = []
        if xmin < 6709.2 and xmax > 6710.2:
            exclude_regions.append(SpectralRegion(6709.2*u.AA, 6710.2*u.AA))
        if xmin < 6679 and xmax > 6681:
            exclude_regions.append(SpectralRegion(6679*u.AA, 6681*u.AA))

        cont_flx = fit_generic_continuum(spec, exclude_regions=exclude_regions)(spec.spectral_axis)
        ax.plot(wav, cont_flx, c='r', zorder=2)

    ax.set_xlabel('wavelength [angstrom]')
    ax.set_ylabel('relative flux')

    if isinstance(xlim, list):
        ax.set_xlim(xlim)

    format_ax(ax)
    savefig(f, outpath)


def make_toi837_feros_plots():

    spectrum_path = '../data/spectra/TOI-837_FEROS.fits'

    plotpath = '../results/TOI_837/feros_spectrum_quicklook.png'
    plot_feros_spectra(spectrum_path, outpath=plotpath)

    xlim = [6702, 6715]
    plotpath = (
        '../results/TOI_837/feros_spectrum_xlim{}.png'.
        format(repr(xlim).replace(' ','_').replace(',',''))
    )
    plot_feros_spectra(spectrum_path, outpath=plotpath, xlim=xlim)

    xlim = [6708, 6712]
    plotpath = (
        '../results/TOI_837/feros_spectrum_xlim{}.png'.
        format(repr(xlim).replace(' ','_').replace(',',''))
    )
    plot_feros_spectra(spectrum_path, outpath=plotpath, xlim=xlim)

    xlim = [6670, 6713]
    plotpath = (
        '../results/TOI_837/feros_spectrum_xlim{}.png'.
        format(repr(xlim).replace(' ','_').replace(',',''))
    )
    plot_feros_spectra(spectrum_path, outpath=plotpath, xlim=xlim)


def get_toi837_li_equivalent_width():

    spectrum_path = '../data/spectra/TOI-837_FEROS.fits'
    plotpath = '../results/TOI_837/feros_spectrum_get_li_equivalent_width.png'

    #
    # fit out the continuum to get the continuum normalized flux, over the
    # window of 6670 angstrom to 6713 angstrom.
    #
    xlim = [6670, 6713]

    hdul = fits.open(spectrum_path)
    d = hdul[0].data
    wav = d[0,0]
    flx = d[3,0]

    if isinstance(xlim, list):
        xmin = xlim[0]
        xmax = xlim[1]
        sel = (wav > xmin) & (wav < xmax)

        wav = wav[sel]
        flx = flx[sel]

    spec = Spectrum1D(spectral_axis=wav*u.AA, flux=flx*u.dimensionless_unscaled)

    if isinstance(xlim, list):
        exclude_regions = []
        if xmin < 6709.2 and xmax > 6710.2:
            exclude_regions.append(SpectralRegion(6709.2*u.AA, 6710.2*u.AA))
        if xmin < 6679 and xmax > 6681:
            exclude_regions.append(SpectralRegion(6679*u.AA, 6681*u.AA))

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
    region = SpectralRegion(6708.5*u.AA, 6711.5*u.AA)
    li_equiv_width = equivalent_width(cont_norm_spec, regions=region)
    li_centroid = centroid(full_spec, region)

    #
    # fit a gaussian too, and get ITS equiv width
    # https://specutils.readthedocs.io/en/stable/fitting.html
    #
    g_init = models.Gaussian1D(amplitude=0.2*u.dimensionless_unscaled,
                               mean=6709.7*u.AA, stddev=0.5*u.AA)
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

    axs[2].set_xlim([6708.5, 6711.5])
    axs[3].set_xlim([6708.5, 6711.5])
    axs[-1].set_xlabel('wavelength [angstrom]')

    for ax in axs:
        format_ax(ax)
    outpath = '../results/TOI_837/toi837_li_equivalent_width_routine.png'
    savefig(f, outpath)






if __name__ == "__main__":

    get_toi837_li_equivalent_width()
    assert 0

    make_toi837_feros_plots()

