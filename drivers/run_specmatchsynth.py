import os
from cdips_followup import __path__
from cdips_followup.specmatchsynth import specmatchsyn_analyze

DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data/spectra')
OUTDIR = os.path.join(os.path.dirname(__path__[0]), 'results/spec_analysis')

class argclass(object):
    pass

args = argclass()

args.spectrum_name = 'TIC268301217_template_spectrum_v20201111.dat'
args.flat_name = 'nf_n58_10.dat' # or none
args.wvsol_name = None
args.idstring = 'TIC268301217_20201111_template'
args.is_template = True

spectrum_path = os.path.join(
    DATADIR, 'PFS', '{}'.format(args.spectrum_name)
)
wvsol_path = os.path.join(
    DATADIR, 'PFS', '{}'.format(args.wvsol_name)
)
if isinstance(args.flat_name, str):
    flat_path = os.path.join(
    DATADIR, 'PFS', '{}'.format(args.flat_name)
)

outdir = os.path.join(OUTDIR, 'PFS', 'synthetic_specmatch')
# pick regions...
regions = ['order{}'.format(ix) for ix in range(35, 53)]

specmatchsyn_analyze(spectrum_path, wvsol_path=wvsol_path, regions=regions,
                     outdir=outdir, idstring=args.idstring,
                     is_template=args.is_template, flat_path=flat_path)
