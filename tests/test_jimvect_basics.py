"""Test suite for overriden basic methods for JimVect objects."""
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
import numpy as np


vector = 'tests/data/nuts_italy.sqlite'
v1     = 'tests/data/v1.json'


class BadBasicMethods(unittest.TestCase):
    """Test funcs and methods on the root level and operations for JimVects."""

    @staticmethod
    def test_numpy_conversions():
        """Test JimVect.np() method."""
        vect = pj.JimVect(vector)
        anp0 = vect.np(ln=0)  # layer 0
        anp1 = vect.np(ln=1)  # layer 1

        assert isinstance(anp0, np.ndarray), \
            'Error in numpy conversion layer 0: not an instance of ndarray'
        assert isinstance(anp1, np.ndarray), \
            'Error in numpy conversion layer 1: not an instance of ndarray'
        assert anp0.shape[0] == 1, \
            'Error in numpy conversion: first dimension (number of features)'
        assert anp0.shape[1] == 12, \
            'Error in numpy conversion: second dimension (number of fields)'
        assert np.append(anp0, anp1, axis=0).shape[0] == 2, \
            'Error in geometry.merge() ' \
            '(append first dimension)'
        assert np.append(anp0, anp1, axis=0).shape[1] == 12, \
            'Error in geometry.merge() ' \
            '(append second dimension)'

    @staticmethod
    def test_jimvect_creations():
        """Test creating of JimVect objects."""
        out_path = pj._get_random_path()

        # Test empty vector

        vect = pj.JimVect()

        assert vect.properties.isEmpty(), \
            'Error in creating an empty JimVect object'
        assert vect.properties.getLayerCount() == 0, \
            'Error in creating an empty JimVect object'

        # Test vector with filepath defined

        vect = pj.JimVect(vector)

        assert not vect.properties.isEmpty(), \
            'Error in creating a JimVect object with specified path'
        assert vect.properties.getLayerCount() == 2, \
            'Error in creating a JimVect object with specified path'

        # Test vector with filepath and kwargs

        vect = pj.JimVect(vector, ln='milano')

        assert not vect.properties.isEmpty(), \
            'Error in creating a JimVect object with specified path and kwargs'
        assert vect.properties.getLayerCount() == 1, \
            'Error in creating a JimVect object with specified path and kwargs'

        # Test with parent vector

        vect = pj.JimVect(vect, output=out_path)

        assert os.path.isfile(out_path), \
            'Error in creating JimVect object based on another JimVect ' \
            '(output does not exist)'
        assert not vect.properties.isEmpty(), \
            'Error in creating JimVect object based on another JimVect ' \
            '(JimVect empty)'
        assert vect.properties.getLayerCount() == 1, \
            'Error in creating JimVect object based on another JimVect ' \
            '(wrong number of layers)'

        # Copy constructor with new field

        # add field with unique key
        vect = pj.JimVect(v1)
        vect = pj.JimVect(vect, output=out_path, newfield='fid')
        assert 'fid' in vect.properties.getFieldNames(), \
            'Error in creating JimVect object based on another JimVect' \
            'with new fid (fid not found)'

        assert np.amax(vect.np(field=['fid'])) + 1 == \
            vect.properties.getFeatureCount(), \
            'Error in creating JimVect object based on another JimVect' \
            'with new fid (value)'

        vect.io.close()
        os.remove(out_path)

        # add field with fixed numeric label of default integer type
        vect = pj.JimVect(v1)
        vect = pj.JimVect(vect, output=out_path, newfield='label', newvalue=10)
        assert 'label' in vect.properties.getFieldNames(), \
            'Error in creating JimVect object based on another JimVect' \
            'with new field (label not found)'

        assert np.amax(vect.np(field=['label'])) == 10, \
            'Error in creating JimVect object based on another JimVect' \
            'with new field (value=10)'

        vect.io.close()
        os.remove(out_path)

        # add field with fixed numeric label of Real type
        vect = pj.JimVect(v1)
        vect = pj.JimVect(vect, output=out_path, newfield='label',
                          newtype='Real', newvalue=9.9, verbose=2)
        assert np.amax(vect.np(field=['label'])) == 9.9, \
            'Error in creating JimVect object based on another JimVect' \
            'with new field (value=9.9)'

        vect.io.close()
        os.remove(out_path)

        # Copy constructor with new field with two layers

        # add field with unique key
        vect = pj.JimVect(vector)
        vect = pj.JimVect(vect, output=out_path, newfield='fid')

        assert 'fid' in vect.properties.getFieldNames(), \
            'Error in creating JimVect object based on another JimVect' \
            'with new fid (fid not found)'
        max_value = 0
        for layer in range(0, vect.properties.getLayerCount()):
            max_value += np.amax(vect.np(field=['fid'])) + 1
        assert max_value == vect.properties.getFeatureCount(), \
            'Error in creating JimVect object based on another JimVect' \
            'with new fid (value)'

        vect.io.close()
        os.remove(out_path)

        # Test with JSON string

        jsonstring = \
            '{"polygons": ' \
                '{"type": "FeatureCollection", ' \
                '"crs": ' \
                    '{"type": "name", ' \
                    '"properties": ' \
                    '{"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}, ' \
                '"features": ' \
                    '[{"type": "Feature", ' \
                    '"properties": {"label": 1}, ' \
                    '"geometry": ' \
                        '{"type": "Polygon", ' \
                        '"coordinates": ' \
                            '[[[ 16.296883885037882, 48.07125730037879 ], ' \
                            '[ 16.29418254261364, 47.787616345833342 ], ' \
                            '[ 16.518393963825762, 47.814629770075761 ], ' \
                            '[ 16.413041609280306, 48.04424387613637 ], ' \
                            '[ 16.296883885037882, 48.07125730037879 ]]' \
            ']}}]}}'

        vect = pj.JimVect(jsonstring)

        assert not vect.properties.isEmpty(), \
            'Error in creating JimVect object based on json string' \
            '(JimVect empty)'
        assert vect.properties.getLayerCount() == 1, \
            'Error in creating JimVect object based on json string' \
            '(wrong number of layers)'
        assert vect.properties.getFeatureCount() == 1, \
            'Error in creating JimVect object based on json string' \
            '(wrong number of features)'

        vect = pj.JimVect(output=vector)

        assert not vect.properties.isEmpty(), \
            'Error in creating JimVect object based on kwarg output=... ' \
            '(JimVect empty)'
        assert vect.properties.getLayerCount() == 2, \
            'Error in creating JimVect object based on kwarg output=... ' \
            '(wrong number of layers)'

        # Test JimVect creation based on WKT

        wkt_string = 'POLYGON ((23.314208 37.768469, 24.039306 37.768469, ' \
                     '24.039306 38.214372, 23.314208 38.214372, ' \
                     '23.314208 37.768469))'

        vect = pj.JimVect(wkt=wkt_string)

        out_path = pj._get_random_path()
        wkt_jimvect = pj.JimVect(vect, output=out_path, oformat='GeoJSON')

        assert wkt_jimvect.properties.getFeatureCount() == 1, \
            'Error in creating JimVect object based on kwarg wkt=... '

        os.remove(out_path)

        # Test wrong JimVects catches

        try:
            _ = pj.JimVect(vect)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching wrong parameters for JimVect creation ' \
            '(no output, but initial vector)'

        try:
            _ = pj.JimVect(vect, ln='milano')
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching wrong parameters for JimVect creation ' \
            '(no output, but initial vector and kwargs)'

        try:
            _ = pj.JimVect(ulx=0, uly=1, lrx=1, lry=0, output=out_path)
            raised = False
        except pj.exceptions.JimVectIllegalArgumentError:
            raised = True

        assert raised, \
            'Error in catching non-existing output without defining ' \
            'a JimVect template'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadBasicMethods)]
    return unittest.TestSuite(suite_list)
