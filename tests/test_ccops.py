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

        # Test distance2d4

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

        # Test distance2dChamfer

        jim = pj.Jim(tiles[0])[0:10, 0:10]

        chamfer_type = 11
        distances = pj.ccops.distance2dChamfer(jim, chamfer_type)
        jim.ccops.distance2dChamfer(chamfer_type)
        stats = jim.stats.getStats(['max', 'min'])

        assert jim.pixops.isEqual(distances), \
            'Inconsistency in ccops.distance2dChamfer() ' \
            '(method returns different result than function)'

        assert 0 < stats['min'] < 10, \
            'Error in Jim.ccops.distance2dChamfer() (suspicious minimum value)'
        assert stats['max'] <= chamfer_type, \
            'Error in Jim.ccops.distance2dChamfer() (wrong maximum value)'

        chamfer_type = 5711
        jim.ccops.distance2dChamfer(chamfer_type)

        assert not jim.pixops.isEqual(distances), \
            'Suspicious values for distance2dChamfer() ' \
            '(same results for different types of distance)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadCCOps)]
    return unittest.TestSuite(suite_list)
