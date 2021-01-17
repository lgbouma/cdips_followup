"""
Tools for working with specmatch-syn (note: you need to be using a py27
environment for it to work)

    specmatchsyn_analyze: main driver

plotting:
    plot_data_model_comparison
    plot_chisq
"""
import os, sys
if not sys.version_info[0] == 2:
    raise AssertionError(
        'specmatch-syn depends on python2.7, and deprecated scipy routines.'
    )

import numpy as np
from scipy.io import readsav

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
from datetime import datetime

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'specmatch-syn')
COELHO05_PATH = os.path.join(LOCALDIR, 'coelho.h5')

def format_ax(ax):
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize('small')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize('small')

def savefig(fig, figpath, writepdf=True, dpi=450):
    fig.savefig(figpath, dpi=dpi, bbox_inches='tight')
    print('{}: made {}'.format(datetime.utcnow().isoformat(), figpath))

    if writepdf:
        pdffigpath = figpath.replace('.png','.pdf')
        fig.savefig(pdffigpath, bbox_inches='tight', dpi=dpi)
        print('{}: made {}'.format(datetime.utcnow().isoformat(), pdffigpath))

    plt.close('all')


def read_pfs(spectrum_path, wvlen_soln, verbose=False, is_template=False):
    """
    Duplicate from cdips_followup.spectools because of py27 compatability.
    """

    if not is_template:
        sp = readsav(spectrum_path, python_dict=True, verbose=verbose)
        wvlen = readsav(wvlen_soln, python_dict=True, verbose=verbose)
        return sp['sp'], wvlen['w']

    else:
        s = readsav(spectrum_path, python_dict=True, verbose=verbose)
        return s['star'], s['w']


def specmatchsyn_analyze(
    spectrum_path, wvsol_path=None, regions=None, outdir=None, idstring=None,
    is_template=False, flat_path=None, START=200, END=-200,
    init_starparam_dict = {'teff': 5850, 'logg': 4.4, 'fe': 0.1, 'vsini': 9,
                           'psf': None, 'rot_method': 'rotmacro'}):

    if 'PFS' in spectrum_path:
        instrument = 'PFS'
    else:
        raise NotImplementedError('check wavlim works for non-pfs spectrum')

    for s in [outdir, idstring]:
        if not isinstance(s, str):
            raise ValueError
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if "PFS" in spectrum_path:
        flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path,
                                  is_template=is_template)
    else:
        raise NotImplementedError

    #
    # continuum normalize using the flat
    #
    if not isinstance(flat_path, str):
        errmsg = (
            'specmatchsynth needs some form of continuum normalization. '
            'at the moment, i only have implemented the flat spectrum division.'
        )
        raise NotImplementedError(errmsg)

    if isinstance(flat_path, str):
        _f = readsav(flat_path, python_dict=True, verbose=True)
        flat_2d = _f['nf']

    # regions are like `[f'order{ix}' for ix in range(35, 53)]`
    _f2d, _w2d, _n2d, _u2d = [], [], [], []
    for region in regions:

        order = int(region.split('order')[-1])

        start = START
        end = END

        flx, wav = flx_2d[order, start:end], wav_2d[order, start:end]
        if isinstance(flat_path, str):
            flat = flat_2d[order, start:end]

        sel = (flx > 0)
        flx, wav = flx[sel], wav[sel]
        if isinstance(flat_path, str):
            flat = flat[sel]
        #u_flx = 1/np.sqrt(flx)

        flx_norm = np.nanpercentile(flx, 95)
        flat_norm = np.nanpercentile(flat, 95)
        flx /= flx_norm
        flat /= flat_norm

        unc_prefactor = 1e-2
        u_flx = np.ones_like(flx)*unc_prefactor

        cont_norm_flx = flx / flat

        outpath = os.path.join(
            outdir, '{}_{}_cont_norm_check.png'.format(idstring, region)
        )

        if not os.path.exists(outpath):

            f,axs = plt.subplots(nrows=2, ncols=1, figsize=(6,4))
            axs[0].plot(wav, flx, c='k', zorder=3, lw=0.5)
            axs[0].plot(wav, flat, c='r', zorder=2, lw=0.5)
            axs[1].plot(wav, cont_norm_flx, c='k', lw=0.5)

            axs[0].set_ylabel('flux')
            axs[1].set_ylabel('contnorm flux')

            axs[-1].set_xlabel('wavelength [angstrom]')
            for ax in axs:
                format_ax(ax)
            savefig(f, outpath, writepdf=0)
            plt.close('all')

        flx = cont_norm_flx

        _f2d.append(flx)
        _w2d.append(wav)
        _n2d.append(flat)
        _u2d.append(u_flx)

    s_flx_2d = np.vstack(_f2d)
    s_wav_2d = np.vstack(_w2d)
    s_flat_2d = np.vstack(_n2d)
    s_u_flx_2d = np.vstack(_u2d)

    #
    # shift spectrum onto model wavelength scale
    #
    import smsyn.library
    import smsyn.io.spectrum
    import smsyn.echelle

    lib = smsyn.library.read_hdf(COELHO05_PATH, wavlim=[4000,7000])

    ref_wav = np.logspace(np.log10(s_wav_2d[0,0]),np.log10(s_wav_2d[-1,-1]),64000)
    _i = init_starparam_dict

    if _i['psf'] is None:
        print(42*'-')
        print('WRN! Ignoring PSF broadening, b/c I assume its less than vsini')
        print(42*'-')

    ref_flux = lib.synth(ref_wav, _i['teff'], _i['logg'], _i['fe'],
                         _i['vsini'], _i['psf'], _i['rot_method'])

    # now SHIFT. key outputs are flux_shift, uflux_shift.
    # save the output spectrum to a fits file.
    NSEG = 5 # N segments to break each order into when doing the cross-correlation.
    PIXVELSHIFT_METHOD = 'spline'
    ech = smsyn.echelle.Echelle(s_wav_2d, s_flx_2d, s_u_flx_2d)
    pixmid, vshift = smsyn.echelle.vshift(ech, ref_wav, ref_flux, nseg=NSEG)
    pvs = smsyn.echelle.PixelVelocityShift(ech.nord, ech.npix, pixmid, vshift)
    dvel = pvs.caculate_dvel(method=PIXVELSHIFT_METHOD)
    ech_shift = smsyn.echelle.shift_orders(ech, ref_wav, dvel)
    flux_shift, uflux_shift = smsyn.echelle.flatten(ech_shift)

    name = idstring
    obs = regions[0] + '_thru_' + regions[-1]
    spec = smsyn.io.spectrum.Spectrum(ref_wav, flux_shift, uflux_shift,
                                      header=dict(name=idstring, obs=obs))
    smspecpath = os.path.join(outdir, '{}_{}.sm.fits'.format(name, obs))
    spec.to_fits(smspecpath)

    # now run the SpecMatch algorithm
    from smsyn.inst.hires.pipeline import (
        Pipeline, grid_search, lincomb, polish, read_pickle
    )
    segfile = os.path.join(smsyn.__path__[0], 'inst', 'pfs', 'segments.csv')
    pipe = Pipeline(smspecpath, COELHO05_PATH, segfile)

    pklpath = pipe.smfile.replace('.fits', '.pkl')

    out = None
    if not os.path.exists(pklpath):
        grid_search(pipe, debug=False)
        lincomb(pipe)
        polish(pipe)
        pipe.to_pickle()
    else:
        pipe, out = read_pickle(pklpath, name=name, obs=obs)

    if out is None:
        out = pd.concat(
            [pd.DataFrame(pipe.polish_output[k]['params_out'], index=[k])
             for k in pipe.polish_output.keys()]
        )

    print(out.describe())

    # chisq plot
    plt.close('all')
    fig = plot_chisq(out, columns=['teff','logg','fe','vsini'])
    outpath = os.path.join(
        outdir, '{}_chisq.png'.format(idstring)
    )
    savefig(fig, outpath)

    # data vs model comparison

    # ignore botched orders for best-fit params
    sel = (
        (out.logg != 3.5)
        &
        (out.logg != 5.0)
    )
    cols = ['teff','logg','fe','vsini','psf','rchisq1','rchisq0']
    s_out = out[sel][cols]
    bestpars = s_out.mean()
    bestparuncs = s_out.std()

    plt.close('all')
    plot_data_model_comparison(
        pipe, bestpars, out, regions, idstring, outdir, fulldata=1,
        feh_references=1
    )
    plt.close('all')
    plot_data_model_comparison(
        pipe, bestpars, out, regions, idstring, outdir
    )
    plt.close('all')
    plot_data_model_comparison(
        pipe, bestpars, out, regions, idstring, outdir, fulldata=1
    )

    # print params
    for n, p, u in zip(cols, bestpars, bestparuncs):
        print('{}: {:.2f} +/- {:.2f}'.format(n, p, u))

    return


def plot_data_model_comparison(pipe, bestpars, out, regions, idstring, outdir,
                               fulldata=0, feh_references=0):
    """
    pieces poached from smsyn.plotting.output
    """

    teff = bestpars['teff']
    logg = bestpars['logg']
    feh = bestpars['fe']
    vsini = bestpars['vsini']
    psf = bestpars['psf']

    segs = pipe.segments

    import smsyn.library, smsyn.io
    lib = smsyn.library.read_hdf(pipe.libfile, wavlim=(segs[0][0], segs[-1][-1]))
    spec = smsyn.io.spectrum.read_fits(pipe.smfile)

    fitwav = []
    fitflux = []
    fitres = []
    fitmod = []
    for s in pipe.segments:
        segment0 = s[0]
        output = pipe.polish_output[segment0]
        wav = output['wav']
        flux = output['flux']
        resid = output['resid']
        model = flux - resid

        fitwav.append(wav)
        fitflux.append(flux)
        fitres.append(resid)
        fitmod.append(model)

    wav = np.hstack(fitwav).flatten()
    flux = np.hstack(fitflux).flatten()
    resid = np.hstack(fitres).flatten()
    model = np.hstack(fitmod).flatten()

    fullwav = spec.wav
    # fullmod = lib.synth(lib.wav, teff, logg, feh, vsini, psf, rot_method='rotmacro')
    fullmod = lib.synth(spec.wav, teff, logg, feh, vsini, psf, rot_method='rotmacro')
    if feh_references:
        feh_metlo, feh_methi = 0.0, 0.3
        fullmod_metlo = lib.synth(spec.wav, teff, logg, feh_metlo, vsini, psf, rot_method='rotmacro')
        fullmod_methi = lib.synth(spec.wav, teff, logg, feh_methi, vsini, psf, rot_method='rotmacro')

    fullspec = spec.flux
    allres = spec.flux - fullmod
    if feh_references:
        allres_metlo = spec.flux - fullmod_metlo
        allres_methi = spec.flux - fullmod_methi
    wavmask = pipe.wav_exclude

    #
    # make the by-order plot!
    #
    for r, seg in zip(regions, segs):

        _s = '' if not fulldata else '_fulldata'
        if feh_references:
            _s += '_fehreferences'
        outpath = os.path.join(outdir, '{}_{}_modelresid{}.png'.
                               format(idstring, r, _s))

        plt.close('all')
        f, axs = plt.subplots(nrows=2, figsize=(8,4), sharex=True)

        crop = (wav >= seg[0]) & (wav <= seg[1])
        pltflux = flux[crop]
        pltwav = wav[crop]
        pltmod = model[crop]

        if fulldata:
            crop = (fullwav >= seg[0]) & (fullwav <= seg[1])
            pltwav = fullwav[crop]
            pltflux = fullspec[crop]
            pltmod = fullmod[crop]
            pltflux *= (np.mean(pltmod)/np.mean(pltflux))

        if feh_references:
            pltmod_metlo = fullmod_metlo[crop]
            pltmod_methi = fullmod_methi[crop]

        o = np.argsort(pltwav)
        mstr = (
            'Model (teff {:d}, logg {:.2f}, feh {:.2f} vsini {:.2f} psf {:.2f})'
            .format(int(teff), logg, feh, vsini, psf)
        )

        axs[0].plot(pltwav[o], pltflux[o], color='k', linewidth=0.7,
                    label='Data', zorder=2)
        axs[0].plot(pltwav[o], pltmod[o], color='C0', linewidth=0.7,
                    label=mstr, zorder=1)
        if feh_references:
            mstr = (
                'Model (teff {:d}, logg {:.2f}, feh {:.2f} vsini {:.2f} psf {:.2f})'
                .format(int(teff), logg, feh_metlo, vsini, psf)
            )
            axs[0].plot(pltwav[o], pltmod_metlo[o], color='C1', linewidth=0.7,
                        label=mstr, zorder=1)
            # mstr = (
            #     'Model (teff {:d}, logg {:.2f}, feh {:.2f} vsini {:.2f} psf {:.2f})'
            #     .format(int(teff), logg, feh_methi, vsini, psf)
            # )
            # axs[0].plot(pltwav[o], pltmod_methi[o], color='C2', linewidth=0.7,
            #             label=mstr, zorder=1)

        axs[0].axhline(1, color='gray', linewidth=0.7, zorder=-2)
        axs[1].axhline(0, color='gray', linewidth=0.7, zorder=-2)
        if feh_references:
            axs[1].plot(pltwav[o], pltflux[o] - pltmod[o], color='C0',
                        linewidth=0.7)
            axs[1].plot(pltwav[o], pltflux[o] - pltmod_metlo[o], color='C1',
                        linewidth=0.7)
            #axs[1].plot(pltwav[o], pltflux[o] - pltmod_methi[o], color='C2',
            #            linewidth=0.7)

        else:
            axs[1].plot(pltwav[o], pltflux[o] - pltmod[o], color='k',
                        linewidth=0.7)


        for w0,w1 in wavmask:
            if w0 > seg[0] or w1 < seg[1]:
                for a in axs:
                    a.axvspan(w0,w1, color='LightGray', zorder=0)

        axs[0].set_ylim(-0.2, 1.2)
        axs[1].set_ylim(-0.2, 0.2)
        for a in axs:
            a.set_xlim(seg[0], seg[1])
            format_ax(a)

        axs[0].legend(loc='lower left', fontsize='xx-small')

        axs[1].set_xlabel('Wavelength [$\\AA$]')
        axs[1].set_ylabel('Residual')
        axs[0].set_ylabel('Flux')
        f.tight_layout()

        savefig(f, outpath, writepdf=0)

    return


def plot_chisq(df, fig=None, columns=['teff','logg','fe'], **kwargs):
    """
    Make a multi-panel plot of chisq
    Poached from smsyn.plotting.output
    """
    ncols = len(columns)
    if fig is None:
        fig,axL = plt.subplots(ncols=ncols, figsize=(4,3), sharey=True)
    else:
        axL = fig.get_axes()

    i = 0
    for col in columns:
        plt.sca(axL[i])
        plt.semilogy()
        plt.scatter(df[col],df['rchisq1'],**kwargs)
        plt.xlabel(col)
        i+=1
    plt.ylabel('rchisq1')
    return fig





