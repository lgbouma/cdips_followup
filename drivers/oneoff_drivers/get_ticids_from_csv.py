import os
import pandas as pd
from astrobase.services.identifiers import gaiadr2_to_tic

csvpath = os.path.join(
    '~/Dropbox/proj/earhart/results/glue_ngc2516_manual_rotator_viz',
    'close_in_prot_and_color_chiron_five.csv'
)

df = pd.read_csv(csvpath)

ticids = []

for s in df.source_id:
    ticids.append(gaiadr2_to_tic(s))

outdf = pd.DataFrame({
    'source_id': df.source_id,
    'ticid': ticids
})

import IPython; IPython.embed()
