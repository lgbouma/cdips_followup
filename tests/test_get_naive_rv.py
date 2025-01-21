import numpy as np
from cdips_followup.spectools import get_naive_rv

outdir = '/Users/luke/Dropbox/proj/cdips_followup/results/spec_analysis/test_results'

TESTDICT = {
    # file name, grid teff, %2f grid logg, expected RV
    'HD22879': ['rj59.2010.fits', 5700, 4.50, 120.34], # F9 good: all except 14
    'HD39881': ['rj63.450.fits', 5700, 4.50, 0.32], # G5
    'HD101259': ['rj152.901.fits', 5000, 3.50, 96.812], # G6 good: 0-6,8, 9, 10,11,13
    'WASP-4': ['rj161.506.fits', 5400, 4.50, 58.84], # G7 good: 0,2,3,4,5,6, 11. (nb cf did fine?)
    'HD97658': ['rj124.477.fits', 5100, 4.50, -1.71], # K1 good: 0-6, 10,11,13
    'HD139323': ['rj122.698.fits', 5100, 4.50, -67.105], # K3 good: 0-6, 10,11,13
    'HD32147': ['rj246.514.fits', 4600, 4.50, 21.54], # K3 good: 0-6, 10,11,13
    'HD36003': ['rj06.1051.fits', 4500, 4.50, -55.56], # K5 good: 0-7, 10,11,13
    'HD245409': ['rj222.308.fits', 4000, 4.50, 22.19], # K7 good: 0-7, 10,11,13
    'HD84035': ['rj69.131.fits', 4500, 4.50, -12.27], # K7 good: 0-6, 10,11,13
    'HD36395': ['rj318.78.fits', 3700, 5.00, 8.68], # M1.5 good: 0-7, 10,11,13
    'HD265866': ['rj149.110.fits', 3500, 5.00, 22.942], #M3
    'GJ388': ['rj229.80.fits', 3400, 5.00, 12.45], # M4.5 vsini=66
    #'GJ687': []
}

drvs = []

for starname,v in TESTDICT.items():

    fitsname, teff, logg, rv_expected = v

    spectrum_path = (f'/Users/luke/Dropbox/proj/cdips_followup/data/'
                     f'spectra/HIRES/{starname}/{fitsname}')
    synth_path = (f'/Users/luke/local/synthetic_spectra/PHOENIX_MedRes/'
                  f'lte0{teff:d}-{logg:.2f}-0.0'
                  f'.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

    df = get_naive_rv(spectrum_path, synth_path, outdir, make_plot=1)

    # recalculate!
    rvs = np.array(df['rv_chisq_minus_bc_kms'])
    rchip_good_orders = [0,2,3,4,5,6,11] # previous: [0,2,3,4,5,6,10,11,13]
    sel_rvs = rvs[np.array(rchip_good_orders)]
    df['meangoodorder_rv_chisq_minus_bc_kms'] = np.round(np.nanmean(sel_rvs), 4)

    rv = df['meangoodorder_rv_chisq_minus_bc_kms'].iloc[0]

    drv = rv_expected - rv

    print(42*'*')
    print(f"{starname}: drv = {drv:.1f}")
    print(42*'*')


    drvs.append(drv)

print(drvs)
import IPython; IPython.embed()

