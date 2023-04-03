"""
Given a TICID, make diagnostic plots to help do an initial analysis of any
light-curve of interest. (E.g., phase-folding, river plots, different wide
figure sizes).
"""

from glob import glob
import os, multiprocessing
import numpy as np, pandas as pd, matplotlib.pyplot as plt

from cdips.utils.lcutils import _given_mag_get_flux

from cdips_followup.quicklooktools import (
    get_tess_data, explore_flux_lightcurves, explore_eleanor_lightcurves,
    explore_mag_lightcurves, make_periodogram, get_kepler_data
)

def quicklooklc(
    ticid,
    outdir = None,
    cdips = 0,
    spoc = 0,
    eleanor = 1,
    cdipspre = 0,
    kepler = 0,
    qlp = 0,
    detrend = 'minimal',#'best', None, 'biweight', 'locor', 'notch', 'minimal'
    do_mag_lcs = 0,
    do_eleanor_lcs = 0,
    do_flux_lcs = 0,
    do_periodogram = 1,
    do_pf = 0,
    require_quality_zero = 0,
    forceylim = None, # [0.93, 1.07]# for the flux light curves
    period = None,
    epoch = None,
    badtimewindows = None,
    slideclipdict = None, #{'window_length':5, 'high':4, 'low':4},
    mask_orbit_edges = True
):

    if do_pf:
        assert (
            isinstance(period, (float,int)) and isinstance(epoch, (float,int))
        )

    pipedict = {'cdips': cdips, 'spoc':spoc, 'eleanor':eleanor,
                'cdipspre': cdipspre, 'kepler':kepler, 'qlp':qlp}
    for k,v in pipedict.items():
        if v:
            pipeline = k

    data, hdrs = get_tess_data(ticid, outdir=outdir, cdips=cdips, spoc=spoc,
                               cdipspre=cdipspre, eleanor=eleanor, qlp=qlp)
    if data is None and kepler:
        data = get_kepler_data(ticid, outdir=outdir)

    if do_eleanor_lcs:
        explore_eleanor_lightcurves(data, hdrs, ticid, period=period, epoch=epoch,
                                    require_quality_zero=require_quality_zero)
        if detrend:
            explore_eleanor_lightcurves(
                data, hdrs, ticid, period=period, epoch=epoch,
                require_quality_zero=require_quality_zero, detrend=detrend,
                do_phasefold=do_pf
            )

    if do_mag_lcs:
        explore_mag_lightcurves(data, hdrs, ticid, period=period, epoch=epoch,
                                do_phasefold=do_pf)

    if do_flux_lcs:

        explore_flux_lightcurves(data, hdrs, ticid, pipeline=pipeline, period=period,
                                 outdir=outdir,
                                 epoch=epoch,
                                 require_quality_zero=require_quality_zero,
                                 forceylim=forceylim,
                                 slideclipdict=slideclipdict,
                                 do_phasefold=do_pf,
                                 mask_orbit_edges=mask_orbit_edges)
        if detrend:
            t, f = explore_flux_lightcurves(
                data, hdrs, ticid, pipeline=pipeline, period=period, epoch=epoch,
                outdir=outdir,
                detrend=detrend, do_phasefold=do_pf,
                badtimewindows=badtimewindows,
                require_quality_zero=require_quality_zero, get_lc=1,
                forceylim=forceylim, slideclipdict=slideclipdict,
                mask_orbit_edges=mask_orbit_edges
            )

    if do_periodogram:
        if pipeline == 'cdips':
            time = data[0]['TMID_BJD']
            flux, err = _given_mag_get_flux(data[0]['IRM1'], data[0]['IRE1'])
            _data = {'time': time, 'flux': flux, 'err': err}
            _data = [_data]
        elif pipeline in ['spoc', 'kepler']:
            _data = data
        else:
            raise NotImplementedError

        for ix, _d in enumerate(_data):
            make_periodogram(_d, ticid, pipeline, id_str=f"{ticid}_{ix}",
                             period_min=0.1, period_max=20, nterms_0=1,
                             outdir=outdir)


if __name__ == "__main__":

    #delLyr_ticids = ['68833286', '28477591', '336409172', '425705688',
    #                 '390058041', '425914815', '390134160', '28769219',
    #                 '377724027', '377724518', '28765338', '28772706',
    #                 '258279124', '258349043', '258349041', '237163448',
    #                 '237163474', '237163213', '237185823', '237184017',
    #                 '237183630', '237183962', '237184441', '29064620',
    #                 '237195399', '1715456535', '237193016', '237195537',
    #                 '237195915', '350989117', '350992827', '350992203',
    #                 '120048930', '120097480', '120046404', '120044575',
    #                 '120044726', '20529532', '20530184', '120098840',
    #                 '120099346', '20534565', '120182387', '120251661',
    #                 '120181767', '120261671', '20697988', '20818771',
    #                 '120423937', '20988031', '399822428', '20995907',
    #                 '399791681', '120498272', '21150335', '21149747',
    #                 '21150307', '41867756', '41874592', '120757393',
    #                 '120688854', '120688481', '42198961', '120896000',
    #                 '120900020', '120972481', '121085990']
    #rsg5_ticids = ['185461568', '185668188', '185667262', '185848343',
    #               '185847064', '185848569', '185948737', '186129360',
    #               '378170387', '378173315', '186178469', '186140431',
    #               '186237466', '186255176', '186237579', '295799540',
    #               '186238215', '186253914', '186238820', '416970247',
    #               '416972085', '406076288', '256066070', '406087385',
    #               '256184740', '256346500', '351939704', '351937745',
    #               '352162814', '193335018', '193539544', '278355321',
    #               '265250815']


    ticid =  '34488204' #
    ticid = '389423271' # speedy mic
    # ticid = '245821931' # hyades 40pc 2Re maybe from s5
    # ticid = '245833065' # another hyad with snr=7.5
    # ticid = '268397995' # V830 Tau (claimed RV planet, rot period 2.741d)
    # ticid = '460950389' # IC 2602, PATHOS-31
    # ticid = '220322660'
    # ticid = '11286209' # sigma Ori E
    ticid = '411614400' # weirdo from s11
    ticid = '11893637' # GJ 1151
    ticid = '281582156'
    ticid = '12359079'
    ticid = '283410775' # V374 Peg -- 0.44d M dwarf, cited in Stauffer+17
    ticid = '302773669' # HD 17156 b
    ticid = '158385467' # KIC 7740983
    ticid = '434567277' # weirdo
    ticid = '296206104' # dip 0.04mag
    ticid = '440113053' # qatar-4
    ticid = '123755508' # pathos-38
    ticid = '374732772' # pathos-42
    ticid = '411614400' # by urs truly
    ticid = '264593815' # CVSO 17
    ticid = '264593828' # CVSO 17 neighbor
    ticid = '264593828' # CVSO 17 neighbor
    ticid = '324101484' # SV Centauri
    ticid = '286347420' # looks like planet?
    ticid = '340486476' # EB + dips?
    ticid = '268301217' # TOI 1937b
    ticid = '360156606' # TOI 1227
    ticid = '96533063' # Hillenbrand weirdo
    ticid = '53682439' # CTOI from CDIPS, need S33
    ticid = '268431671' # manual
    ticid = '339308290' # similar to CD-48
    ticid = '438790187' # CD-48 8201
    ticid = '424470353' # 2MASS J15460752-6258042 50 Myr TT
    ticid = '409525054' # BD+45 598, edge on disk
    ticid = '120105470' # Kepler 1627b
    ticid = '179367009' # J1407, V1400 Cen, Mamajek's object
    ticid = '464405850' # Schneiderman's object
    ticid = '4294779' # toi 2451
    ticid = '238597707' # trojan candidate
    ticid = '438790187' # from Montet and Elsa, 10 Myr LCC
    ticid = '56551765' # Tau1 =V1096 Tau = Anon1
    ticid = '56655841' # Tau2
    ticid = '332207549'
    ticid = '302773669' # HD 17156b
    ticid = '198456933' # complex rotator
    ticid = '364075855' # CR steph-1
    ticid = '440686535' # CR in pleiades
    ticid = '94478915'
    ticid = '329335127' # V1716 Cyg
    ticid = '56655841'
    ticid = '219234987' # gamma dor
    ticid = '27009706'
    ticid = '96680681' # AB Aur
    ticid = '142276270' # TOI-1136
    ticid = '166089209' # toi-1136 friend1
    ticid = '219467520' # toi-1136 friend2
    ticid = '99904467' # toi-1136 friend3
    ticid = '427685831' # WASP-52
    ticid = '198008005' # TOI-858
    ticid = '268301217' # TOI 1937b
    ticid = '1232360'
    ticid = '2234723' # B-star complex rotator
    ticid = '188109809' # CR in alpha-Per?
    ticid = '26412438' # CDIPS CTOI
    ticid = '117689799' # TOI-3504 (probably not in Alessi12)
    ticid = '12996623'
    ticid = '403085618' # the plx/pmra/pmdec missing alpha Per star
    ticid = '427394876' # orion EB1
    ticid = '427395094' # orion EB2
    ticid = '427393298' # orion EB3
    ticid = '160329609' # AP Col, closest PMS star (and second youngest)
    ticid = '277888583' # GJ 393, candidate <10pc PMS star
    ticid = '247985094' # GJ 896A, candidate <10pc PMS star
    ticid = '167756732' # debris disk host
    ticid = '54147487' # debris disk host TWA7 / CE Ant
    ticid = '134873712' # debris disk host TWA25 / V1249 Cen
    ticid = '73540072'
    ticid = '185336364' # tabby's star
    ticid = '368129164' # brightest CPV known?
    ticid = '243921117'

    # # optional #
    # period = 1.395733 # None
    # epoch = 2458597.22243 - 2457000 + 0.23*1.395733 #None
    # badtimewindows = None

    # # # CD-48 8201 / TIC 438790187
    # period = 2.099291 # None
    # #epoch = 2458622.89036 # default
    # epoch = 2458622.8904 #None
    # badtimewindows = [(1616.75,1617.0),(1617.6,1617.8)]

    # # HD 34382 (TIC 2234723)
    # period = 2.46156932
    # epoch = 2459000

    # # Kepler1627
    # period, epoch, badtimewindows = 7.20280608, 2454953.790531, None

    period, epoch = None, None
    #period, epoch = 1.035844, 2458684.473340 # a little long in EM
    #period, epoch = 1.037, 2458684.473340
    period, epoch = 6.44/24, 1640.

    #for ticid in delLyr_ticids:
    quicklooklc(ticid, period=period, epoch=epoch)
