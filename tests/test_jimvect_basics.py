"""Test suite for overriden basic methods for JimVect objects."""

import pyjeo as pj
import unittest

import os
import random
import string


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


class BadBasicMethods(unittest.TestCase):
    """Test funcs and methods on the root level and operations for JimVects."""

    @staticmethod
    def test_jimvect_creations():
        """Test creating of JimVect objects."""
        out_path = os.sep.join(
            [os.sep + 'tmp',
             ''.join(random.choice(string.ascii_lowercase) for _ in range(11))
                + '.sqlite'])

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

        vect.io.close()
        os.remove(out_path)

        # Test with parent vector

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

        # Test wrong JimVects catches

        try:
            _ = pj.JimVect(vect)
            raised = False
        except AttributeError:
            raised = True

        assert raised, \
            'Error in catching wrong parameters for JimVect creation ' \
            '(no output, but initial vector)'

        try:
            _ = pj.JimVect(vect, ln='milano')
            raised = False
        except AttributeError:
            raised = True

        assert raised, \
            'Error in catching wrong parameters for JimVect creation ' \
            '(no output, but initial vector and kwargs)'

        try:
            _ = pj.JimVect(ulx=0, uly=1, lrx=1, lry=0, output=out_path)
            raised = False
        except AttributeError:
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
