import os
from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = os.path.join(os.path.dirname(__path__[0]), 'data', 'phot')
SPECDIR = os.path.join(os.path.dirname(__path__[0]), 'data', 'spectra')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

LOCALBASEDIR = os.path.join(os.path.expanduser('~'), 'local')
if not os.path.exists(LOCALBASEDIR): os.mkdir(LOCALBASEDIR)

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'cdips_followup')
if not os.path.exists(LOCALDIR): os.mkdir(LOCALDIR)

FFICACHEDIR = os.path.join(os.path.expanduser('~'), 'local', 'unpopular_ffi_cache')
if not os.path.exists(FFICACHEDIR): os.mkdir(FFICACHEDIR)
