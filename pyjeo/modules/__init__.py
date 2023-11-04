"""Define a Base class for all modules."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2023 European Union (Joint Research Centre)
#
# This file is part of pyjeo.
#
# pyjeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyjeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.


class JimModuleBase:
    """Base class for Jim modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original Jim object."""
        self._jim_object = caller


class JimListModuleBase:
    """Base class for JimList modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original JimList object."""
        self._jim_list = caller


class JimVectModuleBase:
    """Base class for JimVect modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original JimVect object."""
        self._jim_vect = caller
