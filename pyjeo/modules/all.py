"""
Module containing all functions and methods from the others.

Author(s): Pieter.Kempeneers@ec.europa.eu,
           Ondrej Pesek,
           Pierre.Soille@ec.europa.eu

Copyright (C) 2018-2020 European Union (Joint Research Centre)

This file is part of pyjeo.

pyjeo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyjeo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.
"""

from .properties import *
from .pjio import *
from .pixops import *
from .ngbops import *
from .geometry import *
from .ccops import *
from .classify import *
from .demops import *

from .properties import _Properties, _PropertiesList, _PropertiesVect
from .pjio import _IO, _IOList, _IOVect
from .pixops import _PixOps, _PixOpsList, _PixOpsVect
from .ngbops import _NgbOps, _NgbOpsList, _NgbOpsVect
from .geometry import _Geometry, _GeometryList, _GeometryVect
from .ccops import _CCOps, _CCOpsList, _CCOpsVect
from .classify import _Classify, _ClassifyList, _ClassifyVect
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
