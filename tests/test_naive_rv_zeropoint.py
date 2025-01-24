import os
from os.path import join
import numpy as np, matplotlib.pyplot as plt
from cdips_followup.spectools import get_naive_rv
from cdips_followup.paths import RESULTSDIR, DATADIR

####################
# change these
run_in_parallel = 1
chip = 'i'
####################

outdir = join(RESULTSDIR, 'spec_analysis/naive_rv_test_results')
if not os.path.exists(outdir): os.mkdir(outdir)
outdir = join(outdir, f"{chip}_chip")
if not os.path.exists(outdir): os.mkdir(outdir)

TESTDICT = {
    # file name, grid teff, %2f grid logg, expected RV
    #'HD22879': [f'{chip}j59.2010.fits', 5700, 4.50, 120.34], # DROP bc anomalously blue F9 good: all except 14
    #'HD39881': [f'{chip}j63.450.fits', 5700, 4.50, 0.32], # G5 DROP - struggles
    #'WASP-4': [f'{chip}j161.506.fits', 5400, 4.50, 58.84], # technically nonstandard, G7 good: 0,2,3,4,5,6, 11. (nb cf did fine?)
    'HD101259': [f'{chip}j152.901.fits', 5000, 3.50, 96.812], # G6 good: 0-6,8, 9, 10,11,13
    'HD97658': [f'{chip}j124.477.fits', 5100, 4.50, -1.71], # K1 good: 0-6, 10,11,13
    'HD139323': [f'{chip}j122.698.fits', 5100, 4.50, -67.105], # K3 good: 0-6, 10,11,13
    'HD32147': [f'{chip}j246.514.fits', 4600, 4.50, 21.54], # K3 good: 0-6, 10,11,13
    'HD36003': [f'{chip}j06.1051.fits', 4500, 4.50, -55.56], # K5 good: 0-7, 10,11,13
    'HD245409': [f'{chip}j222.308.fits', 4000, 4.50, 22.19], # K7 good: 0-7, 10,11,13
    'HD84035': [f'{chip}j69.131.fits', 4500, 4.50, -12.27], # K7 good: 0-6, 10,11,13
    'HD97101': [f'{chip}j185.1188.fits', 4100, 4.50, -16.16], #K8
    'HD199305': [f'{chip}j59.1033.fits', 3700, 5.00, -17.144], # M0.5 low vsini
    'GJ686': [f'{chip}j59.1274.fits', 3600, 5.00, -9.499], # M1 vsini=2, good: 0-7, 11, 13 (not 10!)
    'HD36395': [f'{chip}j318.78.fits', 3700, 5.00, 8.68], # M1.5 good: 0-7, 10,11,13
    'GJ908': [f'{chip}j131.1109.fits', 3600, 5.00, -71.084], # M2
    'HD265866': [f'{chip}j149.110.fits', 3500, 5.00, 22.942], # M3
    'GJ687': [f'{chip}j172.1088.fits', 3400, 5.00, -28.720], # M3.5 vsini=79, 0-8, 11, 13
    'GJ388': [f'{chip}j229.80.fits', 3400, 5.00, 12.45], # M4.5 vsini=66, 0-7, 11, 13
    'GJ699': [f'{chip}j130.2075.fits', 3200, 5.00, -110.416], # M4 vsini 113.2 Barnard's star!, 0-8, 11, 13
}

def main():

    rv_list = []
    drvs = []
    std_rvs = []

    for starname,v in TESTDICT.items():

        fitsname, teff, logg, rv_expected = v

        spectrum_path = join(DATADIR, f'spectra/HIRES/{starname}/{fitsname}')
        localdir = join(os.path.expanduser('~'), 'local')
        synth_path = join(localdir,
                          f'synthetic_spectra/PHOENIX_MedRes/'
                          f'lte0{teff:d}-{logg:.2f}-0.0'
                          f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

        df = get_naive_rv(spectrum_path, synth_path, outdir, chip, make_plot=1,
                          run_in_parallel=run_in_parallel)

        # recalculate!
        rvs = np.array(df['rv_chisq_kms'])
        bc = np.array(df['bc_kms'])
        ZEROPOINTS = {
            'r': 81.527, # HIRES instrumental velocity ZP on r-chip is 81.527 +/- 0.664 km/s
            'i': 81.527+5.443, # shift based on M dwarf TiO bands; makes scales consistent
        }
        ZP = ZEROPOINTS[chip]

        rvs = rvs - bc + ZP  # this is the systemic RV!!! (at the order level)
        if chip == 'r':
            chip_good_orders = [0,2,3,4,5,6,11] # previous: [0,2,3,4,5,6,10,11,13]
        elif chip == 'i':
            if teff >= 4000:
                chip_good_orders = [1]
            elif teff < 4000:
                # TiO and K 7700 order 8
                chip_good_orders = [4] # previous: [0,2,3,4,5,6,10,11,13]

        sel_rvs = rvs[np.array(chip_good_orders)]
        df['meangoodorder_rv_chisq_minus_bc_kms'] = np.round(np.nanmean(sel_rvs), 4)

        rv_std = np.round(np.nanstd(sel_rvs), 4)
        df['stdgoodorder_rv_chisq_minus_bc_kms'] = rv_std

        rv = df['meangoodorder_rv_chisq_minus_bc_kms'].iloc[0]

        drv = rv_expected - rv

        print(42*'*')
        print(f"{starname} {fitsname}: drv = {drv:.1f}")
        print(42*'*')

        rv_list.append(rv)
        drvs.append(drv)
        std_rvs.append(rv_std)

    drvs = np.array(drvs)
    std_rvs = np.array(std_rvs)

    weights = 1.0 / (std_rvs**2)
    # calculate zero-point, and its uncertainty, using a weighted sum based on
    # each uncertainty (measured via order to order scatter)
    wmean = np.sum(drvs * weights) / np.sum(weights)
    wstd = np.sqrt(1.0 / np.sum(weights))

    zeropoint = 1. * wmean  + ZP
    zeropoint_unc = 1. * wstd

    # single order case
    if np.isnan(zeropoint):
        zeropoint = ZP
        zeropoint_unc = 0.664

    plt.close("all")
    fig, ax = plt.subplots()
    ax.errorbar(np.arange(len(drvs)), drvs, yerr=std_rvs, ls=":", c='k',
                markersize=2, marker='o')
    ax.set_xlabel("star index (G8 left, M4.5 right)")
    ax.set_ylabel("expected RV - my RV - zeropoint")
    title = f'zeropoint ({chip} chip) = {zeropoint:.3f} +/- {zeropoint_unc:.3f} km/s'
    ax.set_title(title)
    outpath = os.path.join(outdir, "zeropoint.png")
    fig.savefig(outpath, bbox_inches='tight', dpi=300)

    print(title)
    print(f"saved {outpath}")

    import IPython; IPython.embed()



if __name__ == "__main__":
    main()
