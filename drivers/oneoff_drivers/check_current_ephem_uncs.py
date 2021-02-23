from cdips_followup.utils import get_ephemeris_uncertainty
import pandas as pd

df = pd.read_csv('/Users/luke/Dropbox/proj/cdips_followup/data/ephemerides/ephemerides.csv')

for ix, r in df.iterrows():
    eph_unc = get_ephemeris_uncertainty(r['epoch'], r['epoch_unc'],
                                        r['period'], r['period_unc'],
                                        epoch_obs='today')

    print(f"{ix}: {r['source_id']}. {24*eph_unc:.1f} hr.")
