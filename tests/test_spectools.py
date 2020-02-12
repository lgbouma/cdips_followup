from cdips_followup.spectools import given_deltawvlen_get_vsys

from cdips_followup.spectools import (
    read_veloce, read_pfs, fit_continuum
)

debug = 0
test = 1

datadir = '/Users/luke/Dropbox/proj/cdips_followup/data/spectra'

#FIXME
def test_veloce_continuum_fit():


    if 'Veloce' in spectrum_path:
        flx_2d, wav_2d = read_veloce(spectrum_path, start=200, end=-200)
    elif 'PFS' in spectrum_path:
        flx_2d, wav_2d = flx_2d, wav_2d = read_pfs(spectrum_path, wvsol_path)
    else:
        raise NotImplementedError
    target_wav = np.mean(wavlim)
    _preorder = np.argmin(np.abs(wav_2d - target_wav), axis=1)
    viable_orders = np.argwhere((_preorder != wav_2d.shape[1]-1) &
                                (_preorder != 0))
    order = int(viable_orders[np.argmin(np.abs(
        _preorder[viable_orders] - wav_2d.shape[1]/2))])
    flx, wav = flx_2d[order, :], wav_2d[order, :]
    sel = (wav > wavlim[0]) & (wav < wavlim[1])
    flx, wav = flx[sel], wav[sel]
    if 'Veloce' in spectrum_path:
        wav, flx = wav[::-1], flx[::-1]


    fit_continuum(flx, wav, instrument='Veloce')


if debug:
    given_deltawvlen_get_vsys()

if test:

