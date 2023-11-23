"""
Given spectrum from PFS / Veloce / CHIRON / FEROS / TRES:

* Plot all orders
* Measure Li EW at 6708A.
* Visualize Halpha at 6562.8A
* Measure Teff, vsini, logg.
* Measure Ca HK
* Run specmatch-emp.
* Run specmatch-synth (nb. this functionality requires py27 environment)
"""
import os
from os.path import join
import numpy as np
from glob import glob
from cdips_followup import __path__
from cdips_followup.spectools import (
    get_Li_6708_EW, get_Halpha_EW, inspect_pfs, specmatch_viz_compare,
    specmatch_analyze, plot_orders, plot_spec_vs_dwarf_library,
    plot_stack_comparison, measure_veloce_vsini, get_Ca_HK_emission
)

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = os.path.join(os.path.dirname(__path__[0]), 'results/spec_analysis')

class argclass(object):
    pass


def main():
    args = argclass()

    args.do_orders = 1           # plot all orders
    args.do_sms_analysis = 0     # run specmatch-syn analysis
    args.do_sme_analysis = 0     # specmatch-emp for Teff, Rstar, spec compare
    args.do_sme_viz = 0          # specmatch-emp check
    args.do_inspect = 0          # inspect to figure out require rest-frame shift
    args.do_li_ew = 0            # once rest-frame shift is known
    args.do_halpha_ew = 0        # emission, absorption, either way!
    args.do_vsini = 0            # measure vsini
    args.do_ca_hk = 0            # get Ca HK emission properties
    args.do_stack_comparison = 0 # compare versus stack

    args.is_pfs = 0
    args.is_veloce = 0
    args.is_fies = 0
    args.is_tres = 0
    args.is_hires = 0
    args.is_neid = 0
    args.is_harps = 0
    args.is_coralie = 0
    args.is_rvs = 0
    args.is_winered = 0
    args.is_dbsp = 1

    if args.is_pfs:

        # Night-specific quicklooks
        nightind =  3556 # 3555 # # 4258 # 4257 # # 
        args.spectrum_name = 'rn56.{}'.format(nightind)
        args.wvsol_name = 'w_n56.dat'
        args.idstring = 'TIC268301217_20200206_{}'.format(nightind)
        args.is_template = False

        # Template quicklooks
        #args.spectrum_name = 'TIC268301217_template_spectrum_v20201111.dat'
        #args.flat_name = 'nf_n58_10.dat' # or none
        args.spectrum_name = 'TOI3364.dat'
        args.flat_name = None
        args.wvsol_name = None
        args.idstring = 'TOI3364'
        args.is_template = True

        if args.do_sme_analysis:
            args.teff = 5700
        if args.do_li_ew or args.do_ca_hk or args.do_orders:
            args.xshift = 1.80 # set 
            if args.is_template:
                args.idstring += '_shift{:.2f}'.format(args.xshift)

        main_pfs(args)

    elif args.is_veloce:
        main_veloce(args)

    elif args.is_fies:
        main_fies(args)

    elif args.is_tres:
        main_tres(args)

    elif args.is_hires:
        main_hires(args)

    elif args.is_neid:
        main_neid(args)

    elif args.is_harps:
        main_harps(args)

    elif args.is_coralie:
        main_coralie(args)

    elif args.is_rvs:
        main_rvs(args)

    elif args.is_winered:
        main_winered(args)

    elif args.is_dbsp:
        main_dbsp(args)


def main_pfs(args):

    spectrum_path = os.path.join(
        DATADIR, 'PFS', '{}'.format(args.spectrum_name)
    )
    wvsol_path = os.path.join(
        DATADIR, 'PFS', '{}'.format(args.wvsol_name)
    )
    if isinstance(args.flat_name, str):
        flat_path = os.path.join(
        DATADIR, 'PFS', '{}'.format(args.flat_name)
    )

    if args.do_orders:
        outdir = os.path.join(OUTDIR, 'PFS', 'spec_viz_orders')
        plot_orders(spectrum_path, wvsol_path=wvsol_path,
                    outdir=outdir, idstring=args.idstring,
                    is_template=args.is_template, xshift=args.xshift,
                    flat_path=flat_path)

    if args.do_sme_viz:
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
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for deltawav in [7.5, 10]:
            outpath = os.path.join(
                outdir,
                f'{args.idstring}_Li_EW_shift{args.xshift:.2f}_deltawav{deltawav}.png'
            )
            get_Li_6708_EW(spectrum_path, wvsol_path=wvsol_path,
                           xshift=args.xshift, outpath=outpath,
                           is_template=args.is_template, delta_wav=deltawav)

    if args.do_ca_hk:
        outdir = os.path.join(OUTDIR, 'PFS', 'Ca_HK')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        outpath = os.path.join(
            outdir,
            '{}_Ca_HK_shift{:.2f}.png'.format(args.idstring, args.xshift)
        )
        get_Ca_HK_emission(spectrum_path, wvsol_path=wvsol_path,
                           xshift=args.xshift, outpath=outpath,
                           is_template=args.is_template)


    if args.do_sme_analysis:
        outdir = os.path.join(OUTDIR, 'PFS', 'specmatch')
        plot_spec_vs_dwarf_library([5160, 5210], args.teff, outdir,
                                   args.idstring, spectrum_path=spectrum_path,
                                   wvsol_path=wvsol_path,
                                   is_template=args.is_template)
        # r'H$\beta$' outside of HIRES region, 4990 - 6409 A.
        regions = ['Mgb1']
        regions = ['order{}'.format(ix) for ix in range(35, 53)]
        for r in regions:
            specmatch_analyze(spectrum_path, wvsol_path=wvsol_path, region=r,
                              outdir=outdir, idstring=args.idstring,
                              is_template=args.is_template)

    if args.do_sms_analysis:
        raise NotImplementedError(
            'use drivers.run_specmatchsynth (py27 compatability)'
        )


def main_fies(args):

    specname = 'F4_TOI-1000_FICj170127_2019-10-18T05-51-50.779.spec.fits'
    spectrum_path = os.path.join(
        DATADIR, 'FIES', 'toi_1000', specname
    )
    idstring = 'toi1000'

    if args.do_orders:
        outdir = '../results/spec_analysis/FIES/spec_viz_orders/'
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        plot_orders(spectrum_path, outdir=outdir, idstring=idstring)

    else:
        raise NotImplementedError


def main_tres(args):

    # nb. any of these are good
    # specname = '269-003518_2016-01-03_03h36m58s_cb.spec.fits'
    # specname = '269-003518_2016-12-12_05h30m47s_cb.spec.fits'
    # specname = '269-003518_2015-10-01_08h35m20s_cb.spec.fits'
    # specname = '269-003518_2015-10-02_08h09m59s_cb.spec.fits'
    # specname = '269-003518_2017-01-19_04h03m24s_cb.spec.fits'
    specname = 'TRES_spectrum_Kepler1627.fits'
    specname = 'KOI_7368_ab20150601.fits'

    spectrum_path = os.path.join(
        DATADIR, 'TRES', specname
    )
    idstring = 'KOI_7368'

    if args.do_orders:
        outdir = os.path.join(OUTDIR, 'TRES', 'spec_viz_orders')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        plot_orders(spectrum_path, outdir=outdir, idstring=idstring)

    if args.do_li_ew:
        outdir = os.path.join(OUTDIR, 'TRES', 'Li_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        args.xshift = -0.3
        args.idstring = idstring
        outpath = os.path.join(
            outdir,
            '{}_Li_EW_shift{:.2f}.png'.format(args.idstring, args.xshift)
        )
        get_Li_6708_EW(spectrum_path, wvsol_path=None,
                       xshift=args.xshift, outpath=outpath)

    else:
        raise NotImplementedError


def main_hires(args):

    # single star; single spectrum
    idstring = 'HD_34382'
    specname = 'bj488.337.fits'
    specname = 'bj487.135.fits'
    specname = 'bj486.531.fits'

    idstring = 'TIC376116071'
    specname = 'ij490.326.fits'

    idstring = 'TIC146539195'
    #specname = 'bj501.88.fits'
    #specname = 'ij501.88.fits'
    #specname = 'rj501.88.fits'

    # CPVs
    idstring = 'TIC59129133'
    idstring = 'TIC264599508'
    idstring = 'TIC408188366'
    idstring = 'TIC141146667'
    idstring = 'TIC402980664'

    spectrum_paths = glob(os.path.join(
        DATADIR, 'HIRES', idstring, '?j*fits'
    ))

    if args.do_orders:
        outdir = os.path.join(OUTDIR, 'HIRES', 'spec_viz_orders')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for spectrum_path in spectrum_paths:
            specname = os.path.basename(spectrum_path)
            _idstr = idstring + "_" + specname.replace(".fits","")
            plot_orders(spectrum_path, outdir=outdir, idstring=_idstr)

    if args.do_stack_comparison:
        outdir = os.path.join(OUTDIR, 'HIRES', 'stack_comparisons')
        #NOTE: with iodine, Kep1627 2021/08/09.
        # bj: CaHK, ij: Halpha, Li. rj: Na D and Mg b.
        specglob = 'ij423*.fits'
        spectrum_paths = glob(os.path.join(
            '/Users/luke/Dropbox/proj/rudolf/data/spec/20210807_HIRES/CKS_REDUC',
            specglob
        ))
        idstring = 'Kepler1627'
        assert len(spectrum_paths) > 1
        plot_stack_comparison(spectrum_paths, outdir=outdir, idstring=idstring)

    if args.do_li_ew:

        outdir = os.path.join(OUTDIR, 'HIRES', 'Li_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        args.xshift = -0.5
        args.idstring = idstring
        for delta_wav in [2.5,5,7.5]:
            outname = (
                f"{args.idstring}_{specname.replace('.fits','')}_"
                f"Li_EW_shift{args.xshift:.2f}_deltawav{delta_wav:.1f}.png"
            )
            outpath = os.path.join(outdir, outname)
            get_Li_6708_EW(spectrum_path, wvsol_path=None, delta_wav=delta_wav,
                           xshift=args.xshift, outpath=outpath)

    if args.do_halpha_ew:
        outdir = os.path.join(OUTDIR, 'HIRES', 'Halpha_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        args.xshift = -0.05
        args.idstring = idstring
        for delta_wav in [2.5,5,7.5]:
            outname = (
                f"{args.idstring}_{specname.replace('.fits','')}_"
                f"Halpha_shift{args.xshift:.2f}_deltawav{delta_wav:.1f}.png"
            )
            outpath = os.path.join(outdir, outname)
            get_Halpha_EW(spectrum_path, wvsol_path=None, delta_wav=delta_wav,
                          xshift=args.xshift, outpath=outpath)



def main_harps(args):

    # single star; single spectrum
    idstring = 'TOI-858'
    specname = 'TOI-858_spec.dat'
    spectrum_path = os.path.join(
        DATADIR, 'HARPS', idstring, specname
    )

    if args.do_li_ew:

        outdir = os.path.join(OUTDIR, 'HARPS', 'Li_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        args.xshift = 0
        args.idstring = idstring
        for delta_wav in [2.5,5,7.5]:
            outname = (
                f"{args.idstring}_{specname.replace('.fits','')}_"
                f"Li_EW_shift{args.xshift:.2f}_deltawav{delta_wav:.1f}.png"
            )
            outpath = os.path.join(outdir, outname)
            get_Li_6708_EW(spectrum_path, wvsol_path=None, delta_wav=delta_wav,
                           xshift=args.xshift, outpath=outpath)


def main_coralie(args):

    # single star; single spectrum
    idstring = 'TOI-858A'
    specname = 'TIC198008002_spec.dat'
    spectrum_path = os.path.join(
        DATADIR, 'CORALIE', idstring, specname
    )

    if args.do_li_ew:

        outdir = os.path.join(OUTDIR, 'CORALIE', 'Li_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        args.xshift = 0
        args.idstring = idstring
        for delta_wav in [2.5,5,7.5]:
            outname = (
                f"{args.idstring}_{specname.replace('.fits','')}_"
                f"Li_EW_shift{args.xshift:.2f}_deltawav{delta_wav:.1f}.png"
            )
            outpath = os.path.join(outdir, outname)
            get_Li_6708_EW(spectrum_path, wvsol_path=None, delta_wav=delta_wav,
                           xshift=args.xshift, outpath=outpath)


def main_neid(args):

    # single star; single spectrum
    idstring = 'TOI4145'
    specname = 'TOI4145S-sy20211128-neidL2.fits'
    spectrum_path = os.path.join(
        DATADIR, 'NEID', idstring, specname
    )

    if args.do_orders:
        outdir = os.path.join(OUTDIR, 'NEID', 'spec_viz_orders')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        plot_orders(spectrum_path, outdir=outdir, idstring=idstring)

    if args.do_li_ew:
        outdir = os.path.join(OUTDIR, 'NEID', 'Li_EW')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        args.xshift = 1.80
        args.idstring = idstring
        outpath = os.path.join(
            outdir,
            '{}_{}_Li_EW_shift{:.2f}.png'.format(
                args.idstring, specname.replace('.fits',''), args.xshift)
        )
        get_Li_6708_EW(spectrum_path, wvsol_path=None,
                       xshift=args.xshift, outpath=outpath)


def main_veloce(args):

    # idstring = '20200131_837.01'
    # teff = 6100
    # shift = -0.10

    # idstring = '20200131_TIC308538095'
    # teff = 6140
    # shift = -0.20

    # idstring = '20200131_TIC146129309'
    # teff = 6633
    # shift = -0.02

    specname = '{}_Bouma_final_combined.fits'.format(idstring)
    spectrum_path = os.path.join(
        DATADIR, 'Veloce', specname
    )

    if args.do_orders:
        outdir = '../results/spec_analysis/Veloce/spec_viz_orders/'
        plot_orders(spectrum_path, outdir=outdir, idstring=idstring)

    if args.do_sme_analysis:
        outdir = '../results/spec_analysis/Veloce/specmatch/'
        plot_spec_vs_dwarf_library([6280,6350], teff, outdir, idstring,
                                   spectrum_path=spectrum_path,
                                   is_template=args.is_template)
        specmatch_analyze(spectrum_path, region='6300', outdir=outdir,
                          idstring=idstring, is_template=args.is_template)
    if args.do_sme_viz:
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


def main_rvs(args):

    ####################
    # change these #
    DOWNLOAD_SPECTRA = 0
    #cache_id = 'bprp_gt2_dlt20pc'
    cache_id = '20230411_good_CPV_dlt150pc_hasRVS'
    ####################

    cachedir = f'/Users/luke/.gaia_cache/RVS/{cache_id}'

    from cdips.utils.gaiaqueries import given_source_ids_get_gaia_data

    if DOWNLOAD_SPECTRA:
        from cdips.utils.gaiaqueries import run_query_to_get_rvs_spectra
        csvpaths = run_query_to_get_rvs_spectra()
    else:
        csvpaths = glob(join(cachedir, "*.csv.gz"))

    if args.do_orders:
        outdir = cachedir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for spectrum_path in csvpaths:

            dr3_source_id = os.path.basename(spectrum_path).split("_")[1]

            gdf = given_source_ids_get_gaia_data(
                np.array([dr3_source_id]).astype(np.int64),
                f"gdr3_{dr3_source_id}", n_max=5, overwrite=False,
                enforce_all_sourceids_viable=True, savstr='',
                which_columns='*', table_name='gaia_source_lite',
                gaia_datarelease='gaiadr3'
            )

            bp_rp = float(gdf["bp_rp"])
            g = float(gdf["phot_g_mean_mag"])

            idstring = f"bprp{bp_rp:.3f}_G{g:.3f}_gdr3_{dr3_source_id}"

            plot_orders(spectrum_path, outdir=outdir, idstring=idstring)


def main_winered(args):

    if args.do_orders:
        # just plot the orders
        ##########################################
        do0 = 0
        do1 = 1
        if do0:
            outputdir =(
                '/Users/luke/Dropbox/proj/cpv/data/spectra/WINERED/'
                'bouma_june10_hiresy/'
                'bd-185550_WINA00039003_output_v0/'
            )
            object_id = 'bd-185550'
        elif do1:
            outputdir =(
                '/Users/luke/Dropbox/proj/cpv/data/spectra/WINERED/'
                'bouma_june03_hiresy/'
                'cd-38245_WINA00037499_output_v0/'
            )
            object_id = 'cd-38245'
        ##########################################
        airvac = "VAC"
        datadir = join(outputdir, f'{object_id}_sum/{airvac}_norm/fsr1.05')
        spectrum_paths = glob(join(datadir, f'{object_id}*fits'))

        outdir = join(outputdir, f'order_plots_{airvac}_norm_fsr1.05')
        if not os.path.exists(outdir): os.mkdir(outdir)

        from cdips_followup.spectools import read_winered, viz_1d_spectrum

        for spectrum_path in spectrum_paths:
            flx, wav = read_winered(spectrum_path)
            outname = os.path.basename(spectrum_path.replace('.fits','.png'))
            outpath = join(outdir, outname)
            viz_1d_spectrum(flx, wav, outpath, ylabel='flx [cont norm]')

    if args.do_stack_comparison:
        # compare science spectra against telluric refernce stars
        ##########################################
        do0 = 0
        do1 = 1
        if do0:
            basedir = (
                '/Users/luke/Dropbox/proj/cpv/data/spectra/WINERED/'
                'bouma_june10_hiresy/'
            )
            ref_objectid = 'bd-185550'
            ref_frameid = 'WINA00039003'
            refdir = join(basedir, f'{ref_objectid}_{ref_frameid}_output_v0/')
            sci_objectid = 'TIC_167664935'
            sci_frameid = 'WINA00038979'
            scidir = join(basedir, f'{sci_objectid}_{sci_frameid}_output_v0/')

        elif do1:
            basedir = (
                '/Users/luke/Dropbox/proj/cpv/data/spectra/WINERED/'
                'bouma_june10_hiresy/'
            )

            #ref_objectid = 'cd-38245'
            #ref_frameid = 'WINA00037499'
            #refdir = join(basedir, f'{ref_objectid}_{ref_frameid}_output_v0/')
            #ref_objectid = 'TIC_89026133'
            #ref_frameid = 'WINA00037475'
            #refdir = join(basedir, f'{ref_objectid}_{ref_frameid}_output_v0/')

            ref_objectid = 'bd-185550'
            ref_frameid = 'WINA00039003'
            refdir = join(basedir, f'{ref_objectid}_{ref_frameid}_output_v0/')

            sci_objectid = 'TIC_89026133'
            sci_frameid = 'WINA00038992'
            scidir = join(basedir, f'{sci_objectid}_{sci_frameid}_output_v1/')

        ##########################################

        # make directories
        airvac = "VAC"
        outdir = join(basedir, f'stackcomparison_plots_{airvac}_norm_fsr1.05')
        if not os.path.exists(outdir): os.mkdir(outdir)
        outdir = join( outdir,
            f'{sci_objectid}_{sci_frameid}_vs_{ref_objectid}_{ref_frameid}'
        )
        if not os.path.exists(outdir): os.mkdir(outdir)

        # collect the paths
        sci_datadir = join(scidir, f'{sci_objectid}_sum/{airvac}_norm/fsr1.05')
        ref_datadir = join(refdir, f'{ref_objectid}_sum/{airvac}_norm/fsr1.05')

        sci_sppaths = np.sort(glob(join(sci_datadir, f'{sci_objectid}*fits')))
        ref_sppaths = np.sort(glob(join(ref_datadir, f'{ref_objectid}*fits')))

        assert len(sci_sppaths) == len(ref_sppaths)
        assert len(sci_sppaths) > 0

        from cdips_followup.spectools import read_winered, viz_1d_spectrum
        import matplotlib.pyplot as plt
        from aesthetic.plot import set_style, savefig

        for sci_sppath, ref_sppath in zip(sci_sppaths, ref_sppaths):

            plt.close('all')
            set_style("science")
            fig,axs = plt.subplots(nrows=2, ncols=1, figsize=(10,6))

            flx, wav = read_winered(ref_sppath)
            axtitle0 = os.path.basename(ref_sppath).rstrip(".fits").replace("_fsr1.05_VAC_norm","")
            viz_1d_spectrum(flx, wav, None, ylabel='flx [cont norm]',
                            xlabel='', fig=fig, ax=axs[0], axtitle=axtitle0,
                            ylim=[0.7,1.3])

            flx, wav = read_winered(sci_sppath)
            axtitle1 = (
                os.path.basename(sci_sppath).rstrip(".fits").replace("_fsr1.05_VAC_norm","")
                +
                "_"
                +
                sci_frameid
            )
            viz_1d_spectrum(flx, wav, None, ylabel='flx [cont norm]',
                            fig=fig, ax=axs[1], axtitle=axtitle1, ylim=[0.7,1.3])

            outname = f"sci-{axtitle1}_vs_ref-{axtitle0}.png"
            outpath = join(outdir, outname)

            fig.tight_layout()
            savefig(fig, outpath, writepdf=0)


def main_dbsp(args):

    if args.do_orders:
        # just plot the orders
        ##########################################
        outputdir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/20231111/p200_dbsp_blue_A/Viz'
        object_id = 'LP_12-502'
        ##########################################

        datadir = '/Users/luke/Dropbox/proj/cpv/data/spectra/DBSP_REDUX/20231111/p200_dbsp_blue_A/Science'

        spectrum_paths = glob(join(datadir, f'spec1d_blue*{object_id}*fits'))

        outdir = outputdir

        if not os.path.exists(outdir): os.mkdir(outdir)

        from cdips_followup.spectools import read_dbsp, viz_1d_spectrum

        for spectrum_path in spectrum_paths:
            flx, wav = read_dbsp(spectrum_path)
            outname = os.path.basename(spectrum_path.replace('.fits','.png'))
            outpath = join(outdir, outname)
            viz_1d_spectrum(flx, wav, outpath, ylabel='counts')


if __name__ == "__main__":
    main()
