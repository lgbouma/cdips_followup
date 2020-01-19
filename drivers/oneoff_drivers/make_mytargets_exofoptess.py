import os
import numpy as np, pandas as pd
from cdips_followup.utils import get_cdips_candidates
from cdips.utils import today_YYYYMMDD
from cdips_followup import __path__ as cfpath

datadir = os.path.join(os.path.dirname(cfpath[0]), 'data')

df = get_cdips_candidates()

outdf = pd.DataFrame()
outdf['tic'] = df.ticid
outdf['data'] = 'A'
outdf['days'] = 30

outpath = os.path.join(datadir, 'exofoptess_mytargets_uploads',
                       'mytargets_{}.txt'.format(today_YYYYMMDD()))
outdf.to_csv(outpath, index=False, sep='|')
print('made {}'.format(outpath))
