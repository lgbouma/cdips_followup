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

def main():

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
    ticid = '438790187' # from Montet and Elsa, 10 Myr LCC
    ticid = '464405850' # Schneiderman's object
    ticid = '4294779' # toi 2451

    # # optional #
    # period = 1.395733 # None
    # epoch = 2458597.22243 - 2457000 + 0.23*1.395733 #None
    # badtimewindows = None

    # # # CD-48 8201 / TIC 438790187
    # period = 2.099291 # None
    # #epoch = 2458622.89036 # default
    # epoch = 2458622.8904 #None
    # badtimewindows = [(1616.75,1617.0),(1617.6,1617.8)]

    # # Kepler1627
    # period, epoch, badtimewindows = 7.20280608, 2454953.790531, None

    cdips = 0
    spoc = 1
    eleanor = 0
    cdipspre = 0
    kepler = 0
    qlp = 0

    detrend = 'biweight' # None, 'biweight', or 'pspline'
    period, epoch, badtimewindows = None, None, None

    do_mag_lcs = 0
    do_eleanor_lcs = 0
    do_flux_lcs = 1

    do_periodogram = 0
    do_pf = 0
    do_riverplot = 0

    require_quality_zero = 0

    ####################

    pipedict = {'cdips': cdips, 'spoc':spoc, 'eleanor':eleanor,
                'cdipspre': cdipspre, 'kepler':kepler, 'qlp':qlp}
    for k,v in pipedict.items():
        if v:
            pipeline = k

    data = get_tess_data(ticid, outdir=None, cdips=cdips, spoc=spoc,
                         cdipspre=cdipspre, eleanor=eleanor, qlp=qlp)
    if data is None and kepler:
        data = get_kepler_data(ticid, outdir=None)

    if do_eleanor_lcs:
        explore_eleanor_lightcurves(data, ticid, period=period, epoch=epoch,
                                    require_quality_zero=require_quality_zero)
        if detrend:
            explore_eleanor_lightcurves(
                data, ticid, period=period, epoch=epoch,
                require_quality_zero=require_quality_zero, detrend=detrend,
                do_phasefold=do_pf
            )

    if do_mag_lcs:
        explore_mag_lightcurves(data, ticid, period=period, epoch=epoch)

    if do_flux_lcs:

        explore_flux_lightcurves(data, ticid, pipeline=pipeline, period=period,
                                 epoch=epoch,
                                 require_quality_zero=require_quality_zero)
        if detrend:
            t, f = explore_flux_lightcurves(
                data, ticid, pipeline=pipeline, period=period, epoch=epoch,
                detrend=detrend, do_phasefold=do_pf,
                badtimewindows=badtimewindows,
                require_quality_zero=require_quality_zero, get_lc=1
            )

    if do_periodogram:
        if pipeline == 'cdips':
            time = data[0]['TMID_BJD']
            flux, err = _given_mag_get_flux(data[0]['IRM1'], data[0]['IRE1'])
            _data = {'time': time, 'flux': flux, 'err': err}
        elif pipeline in ['spoc', 'kepler']:
            # just look at a single sector
            _data = data[0]
        else:
            raise NotImplementedError
        make_periodogram(_data, ticid, pipeline, period_min=1, period_max=20,
                         nterms_0=1)

    if do_riverplot:
        make_riverplot(data)


if __name__ == "__main__":
    main()
