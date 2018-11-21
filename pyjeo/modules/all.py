"""Module containing all functions and methods from the others."""

from .properties import *
from .pjio import *
from .pixops import *
from .ngbops import *
from .geometry import *
from .ccops import *
from .clssfy import *
from .demops import *

from .properties import _Properties
from .pjio import _IO
from .pixops import _PixOps
from .ngbops import _NgbOps
from .geometry import _Geometry
from .ccops import _CCOps
from .clssfy import _Classify
from .demops import _DEMOps


class _All(_Properties, _IO, _PixOps, _NgbOps, _Geometry, _CCOps, _Classify,
           _DEMOps):

    pass
