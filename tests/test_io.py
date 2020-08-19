"""Test suite for module pyjeo.io."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2020 European Union (Joint Research Centre)
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

import pyjeo as pj
import unittest

import os
import random
import string


tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/modis_ndvi_training.sqlite'


class BadIO(unittest.TestCase):
    """Test functions and methods from io module."""

    @staticmethod
    def test_close():
        """Test closing Jim objects."""
        jim0 = pj.Jim(tiles[0])
        nr_of_cols = jim0.properties.nrOfCol()

        assert nr_of_cols != 0, \
            'Error in creating a Jim object (wrong nr of col)'

        jim0.io.close()
        nr_of_cols = jim0.properties.nrOfCol()

        assert nr_of_cols == 0, \
            'Error in closing a Jim object (nr of col not equal 0 afterwards)'


class BadIOVects(unittest.TestCase):
    """Test JimVect funcs and methods from io module."""

    @staticmethod
    def test_write():
        """Test writing JimVect objects."""
        jimv = pj.JimVect(vector)

        output = os.path.join(
            '/tmp', ''.join(random.sample(string.ascii_letters, 5)))
        jimv.io.write(output)

        assert os.path.isfile(output), \
            'Error in io.write(filename) (file does not exist after writing)'

        output2 = os.path.join(
            '/tmp', ''.join(random.sample(string.ascii_letters, 5)))

        jimv_2 = pj.geometry.convexHull(jimv, output=output2)

        jimv.io.close()
        jimv_2.io.close()

        jimv_3 = pj.JimVect(output2)
        jimv_3.io.close()

        assert jimv_3.properties.getFeatureCount() == 1, \
            'Error in io.write() (changes not written to the original file)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadIO),
                  loader.loadTestsFromTestCase(BadIOVects)]
    return unittest.TestSuite(suite_list)
