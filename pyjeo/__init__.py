"""Initialization file containing version."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
#
# Copyright (C) 2018-2022 European Union (Joint Research Centre)
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
from .pyjeo import *

import os as _os
import tempfile as _tempfile
import random as _random
import string as _string


__version__ = '1.0.8'


def _check_graph(graph, allowed_values):
    """Check whether values used as graph are allowed.

    :param graph: an integer holding for the graph connectivity
    :param allowed_values: values allowed for the graph parameter
    """
    if graph not in allowed_values:
        raise exceptions.JimIllegalArgumentError(
            'Value {} not allow as a graph parameter. Only values {} are '
            'allowed.'.format(graph, allowed_values))


def _get_random_path():
    """Return path of non-existing file in the temp directory.

    Needed for intermediate products for JimVect destructive methods.
    """
    random_string = ''.join(_random.sample(_string.ascii_letters, 5))
    temp_dir = _tempfile.gettempdir()
    non_existing_path = _os.path.join(temp_dir, random_string)
    while _os.path.isfile(non_existing_path):
        random_string = ''.join(_random.sample(_string.ascii_letters, 5))
        non_existing_path = _os.path.join(temp_dir, random_string)

    return non_existing_path
