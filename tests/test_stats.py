"""Test suite for module pyjeo.stats."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


class BadStats(unittest.TestCase):
    """Test functions and methods for getting statistics."""

    def test_getStats(self):
        """Test if values from getStats are not suspicious."""
        jim = pj.Jim(tiles[0])
        stats = jim.stats.getStats()

        max = stats['max']
        min = stats['min']
        mean = stats['mean']

        jim_min_dict = jim.stats.getStats(function='min')

        for key, val in jim_min_dict.iteritems():
            assert val == min, \
                'Error in getting statistics with stats.getStats()'

        assert min < mean < max, \
            'Error in getting statistics with stats.getStats()'

        assert jim_min_dict == pj.stats.getStats(jim, function='min'), \
            'Error in getting statistics with stats.getStats()'

    def test_histograms(self):
        """Test that values of histograms are not suspicious."""
        jim1 = pj.Jim(tiles[0])

        jim1_rows = jim1.properties.nrOfRow()
        jim1_cols = jim1.properties.nrOfCol()

        # Test getHisto1d()
        histo1d = jim1.stats.getHisto1d()

        assert histo1d.stats.getStats()['min'] >= 0, \
            'Error in computing 1D histogram'

        assert histo1d.stats.getStats()['max'] <= jim1_rows * jim1_cols, \
            'Error in computing 1D histogram'

        assert histo1d.pixops.isEqual(pj.stats.getHisto1d(jim1)), \
            'Function and method getHisto1d() return different results'

        jim2 = pj.Jim(tiles[0][:-8] + 'nir' + tiles[0][-5] + '.tif')

        # Test getHisto2d()
        histo2d = jim1.stats.getHisto2d(jim2)

        assert histo2d.stats.getStats()['min'] >= 0, \
            'Error in computing 2D histogram'

        assert histo2d.stats.getStats()['max'] <= jim1_rows * jim1_cols, \
            'Error in computing 2D histogram'

        assert histo2d.pixops.isEqual(pj.stats.getHisto2d(jim1, jim2)), \
            'Function and method getHisto2d() return different results'

        # TODO: Cover histo3d()

        # Test getHistoCumulative
        try:
            _ = jim1.stats.getHistoCumulative()
            failed = True
        except TypeError:
            failed = False
        assert not failed, 'Error in catching wrong data type in ' \
                           'stats.getHistoCumulative()'

        # histo_cumul = pj.pixops.convert(jim1, 4).stats.getHistoCumulative()
        jim2=pj.pixops.convert(jim1, 'Int32')
        histo_cumul=jim2.stats.getHistoCumulative()

        assert histo_cumul.properties.nrOfCol() == jim1_rows * jim1_cols + 1, \
            'Error in stats.getHistoCumulative() ' \
            '(nrOfCol != nrOfRow*Col of original)'

    def test_stretch(self):
        """Test stretching a Jim object."""
        jim = pj.Jim(testFile)

        jim_min = jim.stats.getStats()['min']

        jim_stretched = pj.stats.stretch(jim, dst_min=20)

        if jim_min != 20:
            assert jim_stretched.stats.getStats()['min'] != jim_min, \
                'Error in stretching Jim'
            assert not jim.pixops.isEqual(jim_stretched), \
                'Error in stretching Jim'

        jim.stats.stretch(dst_min=20)

        assert jim.pixops.isEqual(jim_stretched), 'Error in stretching Jim'

    def test_getStatProfile(self):
        """Test if values from getStatProfile are not wrong."""
        jim = pj.Jim(tiles[0])

        stats = jim.stats.getStats()

        min = stats['min']
        max = stats['max']

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        min_profile = pj.stats.getStatProfile(jim, 'min')
        jim.stats.statProfile('max')

        assert min_profile.stats.getStats()['min'] == min, \
            'Error in stats.getStatProfile()'

        assert jim.stats.getStats()['max'] == max, \
            'Error in stats.getStatProfile()'


class BadStatsLists(unittest.TestCase):
    """Test JimList functions and methods for getting statistics."""

    def test_getStatProfile(self):
        """Test if values from getStatProfile are not wrong."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])
        jiml = pj.JimList([jim1, jim2])

        stats1 = jiml.stats.getStats(['min', 'max'])
        stats2 = jim2.stats.getStats()

        min = stats1['min']
        if stats2['min'] < min:
            min = stats2['min']

        max = stats1['max']
        if stats2['max'] > max:
            max = stats2['max']

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        min_profile = jiml.stats.getStatProfile('min')
        max_profile = jiml.stats.getStatProfile('max')

        assert min_profile.stats.getStats()['min'] == min, \
            'Error in stats.getStatProfile()'

        assert max_profile.stats.getStats()['max'] == max, \
            'Error in stats.getStatProfile()'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadStats),
                  loader.loadTestsFromTestCase(BadStatsLists)]
    return unittest.TestSuite(suite_list)
