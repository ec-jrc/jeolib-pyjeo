"""Module containing all functions and methods from the others."""

from .properties import *
from .pjio import *
from .pixops import *
from .ngbops import *
from .geometry import *
from .ccops import *
from .clssfy import *
from .demops import *

from .properties import _Properties, _PropertiesList, _PropertiesVect
from .pjio import _IO, _IOList, _IOVect
from .pixops import _PixOps, _PixOpsList, _PixOpsVect
from .ngbops import _NgbOps, _NgbOpsList, _NgbOpsVect
from .geometry import _Geometry, _GeometryList, _GeometryVect
from .ccops import _CCOps, _CCOpsList, _CCOpsVect
from .clssfy import _Classify, _ClassifyList, _ClassifyVect
from .demops import _DEMOps, _DEMOpsList, _DEMOpsVect


class _All(_Properties, _IO, _PixOps, _NgbOps, _Geometry, _CCOps, _Classify,
           _DEMOps):
    """Inherit all methods."""

    pass


class _AllList(_PropertiesList, _IOList, _PixOpsList, _NgbOpsList,
               _GeometryList, _CCOpsList, _ClassifyList, _DEMOpsList):
    """Inherit all methods for JimLists."""

    pass


class _AllVect(_GeometryVect, _PropertiesVect, _ClassifyVect, _DEMOpsVect,
               _IOVect, _CCOpsVect, _PixOpsVect, _NgbOpsVect):
    """Inherit all methods for JimVects."""

    pass
