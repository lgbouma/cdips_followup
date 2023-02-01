from cdips_followup.spectools import get_Li_6708_EW

delta_wav = 7.5
outpath = "./li_test.png"

# random G star
spectrum_path = '/Users/luke/local/young-KOIs_HIRES_lithium/ij127.710.fits'

# Kepler-1627 / KOI-5245
spectrum_path = '/Users/luke/local/young-KOIs_HIRES_lithium/ij405.85.fits'

# KOI-7368, Teff 5240 K
spectrum_path = '/Users/luke/local/young-KOIs_HIRES_lithium/ij481.291.fits'

# KOI-7913 N = the primary, Teff 4324 K
spectrum_path = '/Users/luke/local/young-KOIs_HIRES_lithium/ij440.75.fits'

# KOI-7913 B = the secondary, Teff 4000 K
spectrum_path = '/Users/luke/local/young-KOIs_HIRES_lithium/ij438.81.fits'

get_Li_6708_EW(spectrum_path, wvsol_path=None, delta_wav=delta_wav,
               outpath=outpath, xshift='find')

assert 0
