"""
Count how many CDIPS light-curves have been put out.
Also though, count how many unique stars there are!
"""

from glob import glob
import os
import numpy as np

lcdir = '/nfs/phtess2/ar0/TESS/PROJ/lbouma/CDIPS_LCS'

lcpaths1 = glob(os.path.join(lcdir, 'sector-?', 'cam?_ccd?', 'hlsp*.fits'))
lcpaths2 = glob(os.path.join(lcdir, 'sector-??', 'cam?_ccd?', 'hlsp*.fits'))

lcpaths = lcpaths1 + lcpaths2
lcnames = list(map(os.path.basename, lcpaths))

gaiaids = list(map(
    lambda x: x.split('gaiatwo')[1].split('-')[0].lstrip('0'), lcnames
))

print('N cdips lightcurves: {}'.format(len(lcpaths)))

print('N unique cdips targets: {}'.format(len(np.unique(gaiaids))))
