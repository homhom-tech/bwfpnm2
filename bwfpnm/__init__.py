import scipy

if scipy.__version__ < '0.14.0':
	raise Exception('OpenPNM requires SciPy version 0.14.0 or greater')

__version__ = '1.0'

__requires__ = ['scipy', 'OpenPNM']

import openpnm

from . import Base
from . import Geometry
from . import Network
#from . import tests
from . import Phases
from . import Physics
from . import Algorithms
from . import Postprocessing
from . import Utilities
from . import routine
from .Utilities import isclose
