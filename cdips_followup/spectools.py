from numpy import array as nparr
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

from astropy.io import fits
from astropy import units as u, constants as const
from scipy.io import readsav

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

from stringcheese.plotutils import savefig, format_ax

line_d = {
    'Mgb1': 5183.62,
    'Mgb2': 5172.70,
    'Feb3': 5168.91,
    'Mgb4': 5167.33
}

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


def given_deltawvlen_get_vsys(deltawvlen=2.3*u.AA, wvlen_0=5180*u.AA):

    # NOTE: currently manual, thru Wright & Eastman's code
    from astropy.time import Time
    t = Time(['2020-02-04T00:00:00'], format='isot', scale='utc')
    print(t.jd)

    # note: varies by ~500m/s due to earth's rotation. (need the exact time to
    # account for this better)
    barycorr = 2203.497975544 # m/s from BARYCORR

    # deltawvlen probably good to ~20-30%, so delta_v no better.
    delta_v = const.c * (deltawvlen / wvlen_0)

    v_star = delta_v + barycorr*(u.m/u.s)
    print(v_star.to(u.km/u.s))



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


def plot_orders(spectrum_path, wvsol_path=None, outdir=None, idstring=None):

    if not isinstance(outdir, str):
        raise ValueError
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)
    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path)
    else:
        raise NotImplementedError

    for order in range(wav_2d.shape[0]):

        if 'PFS' in spectrum_path:
            start = 10
            end = -10
        elif 'Veloce' in spectrum_path:
            start = 200
            end = -200

        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]

        outname = '{}_order{}.png'.format(
            idstring, str(order).zfill(2),
        )

        outpath = os.path.join(outdir, outname)

        viz_1d_spectrum(flx, wav, outpath)




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
                               spectrum_path=None, wvsol_path=None):

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
            flx_2d, wav_2d = flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)
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


##############
# measure EW #
##############
def get_Li_6708_EW(spectrum_path, wvsol_path=None, xshift=None, delta_wav=5,
                   outpath=None):
    """
    spectrum_path: path to PFS or Veloce spectrum

    wvsol_path: path to PFS wavelength solution (optional)

    xshift: angstrom shift required to get into source frame (not vacuum frame).

    delta_wav: window to do the measurement over (angstrom)

    outpath: summary figure is written here.
    """

    if not isinstance(outpath, str):
        raise ValueError

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)
        instrument = 'PFS'
    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
        instrument = 'Veloce'
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


######################
# specmatch analysis #
######################

def specmatch_analyze(spectrum_path, wvsol_path=None, region=None, outdir=None,
                      idstring=None):

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

    if region == 'Mgb1' and instrument=='PFS':
        target_wav = 5183.62
        wavlim = [5160,5210]
    elif region == 'Mgb1' and instrument!='PFS':
        raise NotImplementedError('veloce doesnt cover Mg b1')

    if region == '6300' and instrument=='Veloce':
        target_wav = 6300
        wavlim = [6280,6350]
    elif region == '6300' and instrument!='Veloce':
        raise NotImplementedError(
            'does this region work on {}?'.format(instrument)
        )

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)
    elif "Veloce" in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
    else:
        raise NotImplementedError

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

    outpath = os.path.join(outdir, '{}_cont_norm_check.png'.format(idstring))

    f,axs = plt.subplots(nrows=2, ncols=1, figsize=(6,4))
    axs[0].plot(wav, flx, c='k', zorder=3, lw=0.5)
    axs[0].plot(wav, cont_flx, c='r', zorder=2, lw=0.5)
    axs[1].plot(cont_norm_spec.wavelength, cont_norm_spec.flux, c='k', lw=0.5)

    axs[0].set_ylabel('flux')
    axs[1].set_ylabel('contnorm flux')

    axs[-1].set_xlabel('wavelength [angstrom]')
    for ax in axs:
        format_ax(ax)
    savefig(f, outpath)
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

    outpath = os.path.join(outdir, '{}_shift_check.png'.format(idstring))
    fig = plt.figure(figsize=(10,5))
    sm_res.target_unshifted.plot(normalize=True, plt_kw={'color':'forestgreen'}, text='Target (unshifted)')
    sm_res.target.plot(offset=1.0, plt_kw={'color':'royalblue'}, text='Target (shifted): {}'.format(idstring))
    sm_res.shift_ref.plot(offset=2.0, plt_kw={'color':'firebrick'}, text='Reference: '+match_name)
    plt.xlim((wavlim[0],wavlim[1]))
    plt.ylim(0,3.0)
    ax = plt.gca()
    format_ax(ax)
    savefig(fig, outpath)
    plt.close('all')


    #
    # cross-correlate against the templates to fit for vsini, Rstar, FeH.
    #
    sm_res.match(wavlim=(wavlim[0],wavlim[1]))

    # Plot chi-squared surfaces
    outpath =  os.path.join(outdir, '{}_chisq.png'.format(idstring))
    fig = plt.figure(figsize=(12, 8))
    sm_res.plot_chi_squared_surface()
    ax = plt.gca()
    format_ax(ax)
    savefig(fig, outpath)
    plt.close('all')

    sm_res.lincomb()

    print('Derived Parameters: ')
    print('Teff: {0:.0f}, Radius: {1:.2f}, [Fe/H]: {2:.2f}'.format(
        sm_res.results['Teff'], sm_res.results['radius'], sm_res.results['feh']))

    #
    # make a plot comparing you spectrum to dwarf star spectra of comparable
    # Teff
    #
    plot_spec_vs_dwarf_library(
        wavlim,
        sm_res.results['Teff'],
        outdir,
        idstring,
        sm_res=sm_res
    )

