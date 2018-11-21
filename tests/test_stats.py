"""Test suite for module pyjeo.stats."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadStats(unittest.TestCase):
    """Test functions and methods for getting statistics."""

    def test_histograms(self):
        """Test that values of histograms are not suspicious."""
        jim1 = pj.Jim(tiles[0])

        histo1d = jim1.stats.getHisto1d()

        assert histo1d.stats.getStats()['min'] >= 0, \
            'Error in computing 1D histogram'

        jim2 = pj.Jim(tiles[0][:-8] + 'nir' + tiles[0][-5] + '.tif')

        histo2d = jim1.stats.getHisto2d(jim2)

        assert histo2d.stats.getStats()['min'] >= 0, \
            'Error in computing 2D histogram'

    def test_stretch(self):
        """Test stretching a Jim object."""
        jim = pj.Jim(testFile)

        jim_min = jim.stats.getStats()['min']
        jim_min_dict = jim.stats.getStats(function='min')

        for key, val in jim_min_dict.iteritems():
            assert jim_min == val, 'Error in getting statistics with getStats'

        jim_stretched = pj.stats.stretch(jim, dst_min=20)

        if jim_min != 20:
            assert jim_stretched.stats.getStats()['min'] != jim_min, \
                'Error in stretching Jim'
            assert not jim.pixops.isEqual(jim_stretched), \
                'Error in stretching Jim'

        jim.stats.stretch(dst_min=20)

        assert jim.pixops.isEqual(jim_stretched), 'Error in stretching Jim'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadStats)]
    return unittest.TestSuite(suite_list)
