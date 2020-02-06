"""
Given spectrum from PFS / Veloce / CHIRON / FEROS...
* Measure Li EW at 6708A.
* Vis Halpha at 6562.8A
* Measure Teff, vsini, logg.
"""
import os
from cdips_followup.spectools import (
    get_Li_6708_EW, inspect_pfs, specmatch_viz_compare, specmatch_analyze,
    plot_orders, plot_spec_vs_dwarf_library
)

class argclass(object):
    pass

def main_pfs(args):

    if args.do_orders:

        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outdir = '../results/spec_analysis/PFS/spec_viz_orders/'
        idstring = 'TIC268301217_20200203_3556'

        plot_orders(spectrum_path, wvsol_path=wvsol_path, outdir=outdir,
                    idstring=idstring)

    if args.do_sm_viz:
        # 4990 through 6410 in the SpecMatch library. b/c HIRES needs different
        # echelle settings >~6200A (https://www2.keck.hawaii.edu/inst/hires/)
        specmatch_viz_compare(wavlim=[5160,5210])

    if args.do_inspect:
        for xlim in ['assign', 'fullorder']:
            for night in ['3555', '3556']:
                inspect_pfs(night, 'Halpha', xlim)
                inspect_pfs(night, 'LiI', xlim)
                inspect_pfs(night, 'Mgb1', xlim)

    if args.do_li_ew:
        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outdir = '../results/spec_analysis/PFS/Li_EW/'
        outpath = os.path.join(
            outdir, '{}_Li_EW_shift{:.2f}.png'.format(idstring, shift)
        )


        get_Li_6708_EW(spectrum_path, wvsol_path=wvsol_path, xshift=2.30,
                       outpath=outpath)

    if args.do_sm_analysis:
        # spectrum_path = (
        #     '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3555'
        # )
        # wvsol_path = (
        #     '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        # )
        # outdir = '../results/spec_analysis/PFS/specmatch/'
        # idstring = 'TIC268301217_20200203_3555'

        # specmatch_analyze(spectrum_path, wvsol_path=wvsol_path, region='Mgb1',
        #                   outdir=outdir, idstring=idstring)

        spectrum_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/rn56.3556'
        )
        wvsol_path = (
            '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/PFS/w_n56.dat'
        )
        outdir = '../results/spec_analysis/PFS/specmatch/'
        idstring = 'TIC268301217_20200203_3556'

        specmatch_analyze(spectrum_path, wvsol_path=wvsol_path, region='Mgb1',
                          outdir=outdir, idstring=idstring)


def main_veloce(args):

    idstring = '20200131_837.01'
    teff = 6100
    shift = -0.10

    # idstring = '20200131_TIC308538095'
    # teff = 6140
    # shift = -0.20

    # idstring = '20200131_TIC146129309'
    # teff =6633
    # shift = -0.02

    spectrum_path = (
        '/Users/luke/Dropbox/proj/cdips_followup/data/spectra/Veloce/'
        '{}_Bouma_final_combined.fits'.format(idstring)
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
        #
        # 4990 through 6410 in the SpecMatch library. b/c HIRES needs different
        # echelle settings >~6200A (https://www2.keck.hawaii.edu/inst/hires/).
        # fine for PFS, but for Veloce, we go from like 6000 to 9450.
        #
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




if __name__ == "__main__":

    args = argclass()

    args.do_orders = 0       # plot all orders
    args.do_sm_analysis = 0  # for Teff, Rstar, comparing spectra
    args.do_sm_viz = 0       # specmatch-emp check
    args.do_inspect = 0      # inspect to figure out require rest-frame shift
    args.do_li_ew = 1        # once rest-frame shift is known

    main_veloce(args)
    # main_pfs(args)
