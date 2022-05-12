import os
from astrobase.services.tesslightcurves import get_unpopular_lightcurve
from cdips_followup.paths import PHOTDIR

tic_id = "120105470" # Kepler-1627
download_dir = os.path.join(PHOTDIR, f"TIC{tic_id}")
if not os.path.exists(download_dir):
    os.mkdir(download_dir)

lcfiles = get_unpopular_lightcurve(tic_id, download_dir=download_dir)
