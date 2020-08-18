"""
Basic file containing Jim* objects and functions for numpy conversions.
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


class JimError(Exception):
    """Base class for Jim exceptions."""
    pass


class JimEmptyError(JimError):
    """Exception class for when Jim is empty and should not be."""
    pass


class JimInnerParametersError(JimError):
    """Exception class for when Jim has undefined parameters and should not."""
    pass


class JimBandsError(JimError):
    """Exception class for when Jim band out of bands."""
    pass


class JimIllegalArgumentError(JimError):
    """Exception class for when used arguments do not make sense."""
    pass


class JimListError(Exception):
    """Base class for JimList exceptions."""
    pass


class JimListIllegalArgumentError(JimListError):
    """Exception class for when used arguments do not make sense."""
    pass


class JimVectError(Exception):
    """Base class for JimVect exceptions."""
    pass


class JimVectEmptyError(JimVectError):
    """Exception class for when JimVect is empty and should not be."""
    pass


class JimVectIllegalArgumentError(JimVectError):
    """Exception class for when used arguments do not make sense."""
    pass
