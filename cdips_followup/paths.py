import os
from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = os.path.join(os.path.dirname(__path__[0]), 'data', 'phot')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')

LOCALBASEDIR = os.path.join(os.path.expanduser('~'), 'local')
if not os.path.exists(LOCALBASEDIR):
    os.mkdir(LOCALBASEDIR)

LOCALDIR = os.path.join(os.path.expanduser('~'), 'local', 'cdips_followup')
if not os.path.exists(LOCALDIR):
    os.mkdir(LOCALDIR)
