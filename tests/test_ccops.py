"""Test suite for module pyjeo.ccops."""

import pyjeo as pj
import unittest

import numpy as np


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadCCOps(unittest.TestCase):
    """Test functions and methods from ccops modules."""

    def test_distances(self):
        """Test the distance functions and methods."""
        jim = pj.Jim(tiles[0])

        jim.pixops.convert('Byte')
        distances = pj.ccops.distance2dEuclideanSquared(jim)
        jim.ccops.distance2dEuclideanSquared()

        assert jim.pixops.isEqual(distances), \
            'Error in ccops.distance2dEuclideanSquared()'

        stats = jim.stats.getStats(band=0)

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

        # Test distance2dEuclideanConstrained

        jim1 = pj.Jim(tiles[0])
        jim1.pixops.convert('Byte')
        jim2 = pj.Jim(jim1)

        jim1.pixops.simpleThreshold(127, 250, 0, 1)
        jim2.pixops.simpleThreshold(127, 200, 0, 1)

        jim_byte = pj.Jim(jim2)

        distances = pj.ccops.distance2dEuclideanConstrained(jim2, jim1)
        jim_byte.ccops.distance2dEuclideanConstrained(jim1)
        stats = jim_byte.stats.getStats(band=0)
        max = stats['max']

        assert jim_byte.pixops.isEqual(distances), \
            'Inconsistency in ccops.distance2dEuclideanConstrained() ' \
            '(method returns different result than function)'

        assert not jim2.pixops.isEqual(distances), \
            'Error in ccops.distance2dEuclideanConstrained() ' \
            '(did not have any effect)'

        assert 0 < max < jim_byte.properties.nrOfCol() * \
               jim_byte.properties.nrOfRow(), \
            'Error in Jim.ccops.distance2dEuclideanConstrained() ' \
            '(maximum value not smaller than nrOfCol * nrOfRow or equal to 0)'
        assert stats['min'] == 0, \
            'Error in Jim.ccops.distance2dEuclideanConstrained() ' \
            '(minimum value not 0 for Byte)'

        # Test distanceInfluenceZones2dEuclidean
        nrow = ncol = 500
        jim = pj.Jim(nrow=nrow, ncol=ncol, otype='Byte')
        for i in range(15):
            jim[np.random.randint(0, 500), np.random.randint(0, 500)] = i

        jim_byte = pj.Jim(jim)

        distances = pj.ccops.distanceInfluenceZones2dEuclidean(jim_byte)
        jim_byte.ccops.distanceInfluenceZones2dEuclidean()

        stats = jim_byte.stats.getStats(['max', 'min'])

        assert jim_byte.pixops.isEqual(distances), \
            'Inconsistency in ccops.distanceInfluenceZones2dEuclidean() ' \
            '(method returns different result than function)'

        assert not jim.pixops.isEqual(distances), \
            'Error in ccops.distanceInfluenceZones2dEuclidean() ' \
            '(did not have any effect)'

        assert stats['max'] == 14, \
             'Error in Jim.ccops.distanceInfluenceZones2dEuclidean() ' \
             '(maximum value not 255 for Byte)'
        assert stats['min'] >= 1, \
            'Error in Jim.ccops.distanceInfluenceZones2dEuclidean() ' \
            '(minimum value not 1 for Byte)'

        # Test distanceGeodesic

        jim1 = pj.Jim(tiles[0])
        jim1.pixops.convert('Byte')
        jim2 = pj.Jim(jim1)
        # jim2.pixops.convert('Byte')

        jim1.pixops.simpleThreshold(127, 250, 0, 1)
        jim2.pixops.simpleThreshold(127, 200, 0, 1)

        jim_byte1 = pj.Jim(jim1)
        jim_byte2 = pj.Jim(jim1)

        distances = pj.ccops.distanceGeodesic(jim_byte1, jim2, 4)
        jim_byte1.ccops.distanceGeodesic(jim2, 4)
        stats = jim_byte1.stats.getStats(band=0)
        mean = stats['mean']

        assert jim_byte1.pixops.isEqual(distances), \
            'Inconsistency in ccops.distanceGeodesic() ' \
            '(method returns different result than function)'

        assert not jim1.pixops.isEqual(distances), \
            'Error in ccops.distanceGeodesic() ' \
            '(did not have any effect)'

        assert stats['max'] == 255, \
            'Error in Jim.ccops.distanceGeodesic() ' \
            '(maximum value not 255 for Byte)'
        assert stats['min'] == 0, \
            'Error in Jim.ccops.distanceGeodesic() ' \
            '(minimum value not 0 for Byte)'

        distances = pj.ccops.distanceGeodesic(jim_byte2, jim2, 8)
        jim_byte2.ccops.distanceGeodesic(jim2, 8)
        stats = jim_byte2.stats.getStats(band=0)

        assert jim_byte2.pixops.isEqual(distances), \
            'Inconsistency in ccops.distanceGeodesic() ' \
            '(method returns different result than function)'

        assert not jim_byte2.pixops.isEqual(jim_byte1), \
            'Error in ccops.distanceGeodesic() ' \
            '(the same result for graph=8 and graph=4)'

        assert stats['max'] == 255, \
            'Error in Jim.ccops.distanceGeodesic() ' \
            '(maximum value not 255 for Byte)'
        assert stats['min'] == 0, \
            'Error in Jim.ccops.distanceGeodesic() ' \
            '(minimum value not 0 for Byte)'
        assert stats['mean'] < mean, \
            'Error in Jim.ccops.distanceGeodesic() ' \
            '(mean value for graph=8 not smaller than for graph=4)'

    def test_labelling(self):
        """Test the distance functions and methods."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])

        jim1.pixops.convert('Byte')
        jim2.pixops.convert('Byte')

        # Test labelConstrainedCCsVariance()

        labelled = pj.ccops.labelConstrainedCCsVariance(
            jim1, 0, 0, 0, 0, 0, 0, 4)
        labelled_different = pj.ccops.labelConstrainedCCsVariance(
            jim1, 0, 0, 0, 0, 0, 0, 8)
        jim1_copy = pj.Jim(jim1)
        jim1_copy.ccops.labelConstrainedCCsVariance(0, 0, 0, 0, 0, 0, 4)

        stats = labelled.stats.getStats(band=0)

        assert jim1_copy.pixops.isEqual(labelled), \
            'Inconsistency in ccops.labelConstrainedCCsVariance() ' \
            '(method returns different result than function)'

        assert not labelled.pixops.isEqual(labelled_different), \
            'Error in ccops.labelConstrainedCCsVariance() ' \
            '(created the same object for different alpha value)'

        assert stats['min'] == 0, \
            'Error in Jim.ccops.labelConstrainedCCsVariance() ' \
            '(minimum value not equal to 0)'
        assert 0 < stats['max'] < jim1.properties.nrOfCol() * \
               jim1.properties.nrOfRow(), \
            'Error in Jim.ccops.labelConstrainedCCsVariance() ' \
            '(maximum value not smaller than nrOfCol * nrOfRow or ' \
            'equal to 0)'

        labelled_different = pj.ccops.labelConstrainedCCsVariance(
            jim1, 0, 0, 0, 1, 1, 5, 4)
        stats2 = labelled_different.stats.getStats(band=0)

        assert stats2['max'] < stats['max'], \
            'Error in Jim.ccops.labelConstrainedCCsVariance() ' \
            '(maximum value for rg2, rl2 and varmax2 > rg1, rl1 and varmax1 ' \
            'not smaller)'

        labelled_different = pj.ccops.labelConstrainedCCsVariance(
            jim1, 5, 5, 0, 0, 0, 0, 4)

        assert labelled_different[1, 1].np()[0, 0] == 0, \
            'Error in Jim.ccops.labelConstrainedCCsVariance() ' \
            '(some of parameters ox, oy, oz not applied)'

        non_zero1 = np.count_nonzero(jim1_copy.np())
        non_zero2 = np.count_nonzero(labelled_different.np())
        assert non_zero1 > non_zero2, \
            'Error in Jim.ccops.labelConstrainedCCsVariance() ' \
            '(some of parameters ox, oy, oz not applied)'

        # Test labelFlatZonesSeeded()
        nr_of_col = nr_of_row = 20
        jim = pj.Jim(jim1[0:nr_of_col, 0:nr_of_row])

        ngb = pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1

        seeds = pj.Jim(ncol=nr_of_col,
                       nrow=nr_of_row,
                       uniform=[0, 2],
                       otype='Byte')

        labelled = pj.ccops.labelFlatZonesSeeded(jim, ngb, seeds,
                                                 1, 1, 0)
        jim.ccops.labelFlatZonesSeeded(ngb, seeds, 1, 1, 0)

        labelled_different = pj.ccops.labelFlatZonesSeeded(
            jim, ngb, seeds, 0, 0, 0)

        stats = jim.stats.getStats(['min', 'max'])

        assert jim.pixops.isEqual(labelled), \
            'Inconsistency in ccops.labelFlatZonesSeeded() ' \
            '(method returns different result than function)'
        assert not labelled.pixops.isEqual(labelled_different), \
            'Error in ccops.labelFlatZonesSeeded() ' \
            '(created the same object for different ox and oy)'
        assert stats['min'] == 0, \
            'Error in Jim.ccops.labelFlatZonesSeeded() ' \
            '(minimum value not equal to 0)'
        assert labelled.np()[0, 0] == 0, \
            'Error in Jim.ccops.labelFlatZonesSeeded() ' \
            '(value at position [0, 0] with 3x3 jim_ngb not equal to 0)'
        assert 0 < stats['max'] < nr_of_col * nr_of_row, \
            'Error in Jim.ccops.labelFlatZonesSeeded() ' \
            '(maximum value not smaller than nrOfCol * nrOfRow or ' \
            'equal to 0)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadCCOps)]
    return unittest.TestSuite(suite_list)
