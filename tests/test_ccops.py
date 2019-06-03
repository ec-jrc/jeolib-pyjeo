"""Test suite for module pyjeo.ccops."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadCCOps(unittest.TestCase):
    """Test functions and methods from ccops modules."""

    def test_distances(self):
        """Test the distance2d functions and methods."""
        jim = pj.Jim(tiles[0])

        jim.pixops.convert('Byte')
        distances = pj.ccops.distance2dEuclideanSquared(jim)
        jim.ccops.distance2dEuclideanSquared()

        assert jim.pixops.isEqual(distances), \
            'Error in ccops.distance2dEuclideanSquared()'

        stats = jim.stats.getStats()

        assert stats['min'] == 0, 'Error in ccops.distance2dEuclideanSquared()'
        assert stats['max'] <= \
               jim.properties.nrOfCol()*jim.properties.nrOfRow(), \
            'Error in ccops.distance2dEuclideanSquared()'

        jim = pj.Jim(tiles[0])

        distances = pj.ccops.distance2d4(jim)
        jim.ccops.distance2d4()
        stats = jim.stats.getStats(['max', 'min'])

        assert jim.pixops.isEqual(distances), \
            'Inconsistency in ccops.distance2d4() ' \
            '(method returns different result than function)'

        assert stats['max'] == jim.properties.nrOfRow() / 2 - 1, \
            'Error in Jim.ccops.distance2d4() (wrong maximum value)'
        assert stats['min'] == 0, \
            'Error in Jim.ccops.distance2d4() (wrong minimum value)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadCCOps)]
    return unittest.TestSuite(suite_list)
