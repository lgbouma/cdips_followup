"""
Given spectrum from PFS / Veloce / CHIRON / FEROS...
* Measure Li EW at 6708A.
* Vis Halpha at 6562.8A
* Measure Teff, vsini, logg.
"""
import os
from cdips_followup import __path__
from cdips_followup.spectools import (
    get_Li_6708_EW, inspect_pfs, specmatch_viz_compare, specmatch_analyze,
    plot_orders, plot_spec_vs_dwarf_library, measure_veloce_vsini
)

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = os.path.join(os.path.dirname(__path__[0]), 'results/spec_analysis')

class argclass(object):
    pass


def main():
    args = argclass()

    args.do_orders = 0          # plot all orders
    args.do_sm_analysis = 0     # get Teff, Rstar + compare spectra
    args.do_sm_viz = 0          # specmatch-emp check
    args.do_inspect = 0         # inspect to figure out require rest-frame shift
    args.do_li_ew = 0           # once rest-frame shift is known
    args.do_vsini = 1           # measure vsini

    args.is_pfs = 0
    args.is_veloce = 1

    if args.is_pfs:
        nightind =  3556 # 4257
        args.spectrum_name = 'rn56.{}'.format(nightind)
        args.wvsol_name = 'w_n56.dat'
        args.idstring = 'TIC268301217_20200206_{}'.format(nightind)
        if args.do_sm_analysis:
            args.teff = 5700
        if args.do_li_ew:
            args.xshift = 2.30 # set 
        main_pfs(args)

    elif args.is_veloce:
        main_veloce(args)


def main_pfs(args):

    spectrum_path = os.path.join(
        DATADIR, 'PFS', '{}'.format(args.spectrum_name)
    )
    wvsol_path = os.path.join(
        DATADIR, 'PFS', '{}'.format(args.wvsol_name)
    )

    if args.do_orders:
        outdir = os.path.join(OUTDIR, 'PFS', 'spec_viz_orders')
        plot_orders(spectrum_path, wvsol_path=wvsol_path,
                    outdir=outdir, idstring=args.idstring)

    if args.do_sm_viz:
        # 4990 through 6410 in the SpecMatch library. b/c HIRES needs different
        # echelle settings >~6200A (https://www2.keck.hawaii.edu/inst/hires/)
        specmatch_viz_compare(wavlim=[5160,5210])

    if args.do_inspect:
        for xlim in ['assign', 'fullorder']:
            night = args.spectrum_name.split('.')[-1]
            inspect_pfs(night, 'Halpha', xlim)
            inspect_pfs(night, 'LiI', xlim)
            inspect_pfs(night, 'Mgb1', xlim)

    if args.do_li_ew:
        outdir = os.path.join(OUTDIR, 'PFS', 'Li_EW')
        outpath = os.path.join(
            outdir,
            '{}_Li_EW_shift{:.2f}.png'.format(args.idstring, args.xshift)
        )
        get_Li_6708_EW(spectrum_path, wvsol_path=wvsol_path,
                       xshift=args.xshift, outpath=outpath)

    if args.do_sm_analysis:
        outdir = os.path.join(OUTDIR, 'PFS', 'specmatch')
        plot_spec_vs_dwarf_library([5160, 5210], args.teff, outdir,
                                   args.idstring, spectrum_path=spectrum_path,
                                   wvsol_path=wvsol_path)
        specmatch_analyze(spectrum_path, wvsol_path=wvsol_path, region='Mgb1',
                          outdir=outdir, idstring=args.idstring)


def main_veloce(args):

    # idstring = '20200131_837.01'
    # teff = 6100
    # shift = -0.10

    # idstring = '20200131_TIC308538095'
    # teff = 6140
    # shift = -0.20

    idstring = '20200131_TIC146129309'
    teff = 6633
    shift = -0.02

    specname = '{}_Bouma_final_combined.fits'.format(idstring)
    spectrum_path = os.path.join(
        DATADIR, 'Veloce', specname
    )

    if args.do_orders:
        outdir = '../results/spec_analysis/Veloce/spec_viz_orders/'
        plot_orders(spectrum_path, outdir=outdir, idstring=idstring)

    if args.do_sm_analysis:
        outdir = '../results/spec_analysis/Veloce/specmatch/'
        plot_spec_vs_dwarf_library([6280,6350], teff, outdir, idstring,
                                   spectrum_path=spectrum_path)
        specmatch_analyze(spectrum_path, region='6300', outdir=outdir,
                          idstring=idstring)
    if args.do_sm_viz:
        #specmatch_viz_compare(wavlim=[5160,5210])
        specmatch_viz_compare(wavlim=[6270,6350])

    if args.do_inspect:
        for xlim in ['assign', 'fullorder']:
            for night in ['3555', '3556']:
                inspect_pfs(night, 'Halpha', xlim)
                inspect_pfs(night, 'LiI', xlim)
                inspect_pfs(night, 'Mgb1', xlim)

    if args.do_li_ew:
        outdir = '../results/spec_analysis/Veloce/Li_EW/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        outpath = os.path.join(
            outdir, '{}_Li_EW_shift{:.2f}.png'.format(idstring, shift)
        )

        get_Li_6708_EW(spectrum_path, xshift=shift, outpath=outpath)

    if args.do_vsini:
        outdir = '../results/spec_analysis/Veloce/vsini/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        vsini, shift, gamma = measure_veloce_vsini(
            specname, idstring, teff, outdir)


if __name__ == "__main__":
    main()
