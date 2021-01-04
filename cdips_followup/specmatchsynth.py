"""
Tools for working with specmatch-syn (note: you need to be using a py27
environment for it to work)
"""
import os
import numpy as np
from scipy.io import readsav

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from numpy import array as nparr
from matplotlib.transforms import blended_transform_factory

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'specmatch-syn')
COELHO05_PATH = os.path.join(LOCALDIR, 'coelho.h5')

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
                           'psf': None, 'rot_method': 'rotmacro'}
):

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
        u_flx = np.sqrt(flx)

        flx_norm = np.nanpercentile(flx, 95)
        flat_norm = np.nanpercentile(flat, 95)
        flx /= flx_norm
        u_flx /= flx_norm
        flat /= flat_norm

        # spec = Spectrum1D(spectral_axis=wav*u.AA,
        #                   flux=flx*u.dimensionless_unscaled)
        # cont_norm_spec = spec / flat

        cont_norm_flx = flx / flat

        outpath = os.path.join(
            outdir, '{}_{}_cont_norm_check.png'.format(idstring, region)
        )

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

        _f2d.append(flx.value)
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
    from smsyn.inst.hires.pipeline import Pipeline, grid_search, lincomb
    segfile = os.path.join(smsyn.__path__[0], 'inst', 'pfs', 'segments.csv')
    pipe = Pipeline(smspecpath, COELHO05_PATH, segfile)
    grid_search(pipe, debug=False)
    lincomb(pipe)



    import IPython; IPython.embed()
    assert 0

    if return_velocities:
        return flux_shift, uflux_shift, dvel


    #FIXME: ok great. now from here SHIFT.

    flux_shift, uflux_shift = smsyn.inst.hires.shift.shift(
        wav, flux, uflux, ref_wav, ref_flux
    )

    spec = smsyn.io.spectrum.Spectrum(
        ref_wav, flux_shift, uflux_shift, header=dict(name=name,obs=obs)
    )
    spec.to_fits(outfile) 

    import IPython; IPython.embed()

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

    outpath = os.path.join(outdir, '{}_{}_shift_check.png'.
                           format(idstring, region))
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
    outpath =  os.path.join(outdir,
                            '{}_{}_chisq.png'.format(idstring, region))
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


