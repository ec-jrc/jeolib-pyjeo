"""Module containing all functions and methods from the others."""

from .properties import *
from .pjio import *
from .pixops import *
from .ngbops import *
from .geometry import *
from .ccops import *
from .clssfy import *
from .demops import *

from .properties import _Properties, _PropertiesList
from .pjio import _IO, _IOList
from .pixops import _PixOps, _PixOpsList
from .ngbops import _NgbOps, _NgbOpsList
from .geometry import _Geometry, _GeometryList
from .ccops import _CCOps, _CCOpsList
from .clssfy import _Classify, _ClassifyList
from .demops import _DEMOps, _DEMOpsList


class _All(_Properties, _IO, _PixOps, _NgbOps, _Geometry, _CCOps, _Classify,
           _DEMOps):
    """Inherit all methods."""

    pass


class _AllList(_PropertiesList, _IOList, _PixOpsList, _NgbOpsList,
               _GeometryList, _CCOpsList, _ClassifyList, _DEMOpsList):
    """Inherit all methods for JimLists."""

    pass
