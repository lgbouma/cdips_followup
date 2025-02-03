import os
from astrobase.services.tesslightcurves import get_unpopular_lightcurve
from cdips_followup.paths import PHOTDIR

tic_id = "120105470" # Kepler-1627
tic_id = "107012050" # DENIS J1048
tic_id = "396740648" # Luke Handley's TOI
tic_id = "169138338" # iras lynne's object
tic_id = "167913198" # Gliese710

download_dir = os.path.join(PHOTDIR, f"TIC{tic_id}")
if not os.path.exists(download_dir):
    os.mkdir(download_dir)

lcfiles = get_unpopular_lightcurve(tic_id, ffi_dir=download_dir)
