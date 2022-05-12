import os
from cdips_followup import __path__
DATADIR = os.path.join(os.path.dirname(__path__[0]), 'data')
PHOTDIR = os.path.join(os.path.dirname(__path__[0]), 'data', 'phot')
RESULTSDIR = os.path.join(os.path.dirname(__path__[0]), 'results')
