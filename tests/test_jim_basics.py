"""Test suite for overriden basic methods for Jim objects."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


class BadBasicMethods(unittest.TestCase):
    """Test functions and methods on the root level and operations for Jims."""

    def test_jim_creations(self):
        """Test creating of Jim objects."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(jim1, copyData=True)
        jim3 = pj.Jim(jim1, copyData=False)

        assert jim1.pixops.isEqual(jim2), 'Error in creating Jim object'
        assert not jim1.pixops.isEqual(jim3), 'Error in creating Jim object'

        jim4 = pj.Jim(jim1, nrow=5)

        assert jim1.pixops.isEqual(jim4), 'Error in ignoring kwargs when ' \
                                          'creating Jim object with ' \
                                          'Jim(jim, kwargs)'

        jim5 = pj.Jim(tiles[0], nrow=5)

        assert jim5.properties.nrOfRow() == 5, \
            'Error in creating Jim with Jim(filepath, kwargs)'
        assert jim5.properties.nrOfCol() == jim1.properties.nrOfCol(), \
            'Error in creating Jim with Jim(filepath, kwargs)'

        jim6 = pj.Jim(nrow=5, ncol=5)

        assert jim6.properties.nrOfRow() == 5, \
            'Error in creating Jim with Jim(kwargs)'
        assert jim6.all.nrOfCol() == 5, \
            'Error in creating Jim with Jim(kwargs)'

        jim7 = pj.Jim(nrow=500, ncol=500, otype='Float32', uniform=[0, 2],
                      seed=0)
        stats7 = jim7.stats.getStats(['min', 'max'])

        assert stats7['max'] < 2 and stats7['min'] > 0, \
            'Error in creating Jim with uniform distribution'

        jim8 = pj.Jim(nrow=500, ncol=500, otype='Float32', uniform=[0, 2],
                      seed=0)

        assert jim7.pixops.isEqual(jim8), \
            'Error in usage of seed keyword argument when creating Jim()'

        try:
            _ = pj.Jim(seed=5)
            failed = True
        except AttributeError:
            failed = False

        assert not failed, \
            'Error in catching a call of Jim creation with nonsense (kw)args'

    def test_numpy_conversions(self):
        """Test conversions to numpy and back."""
        jim = pj.Jim(tiles[0])

        jim_np = pj.jim2np(jim)
        jim_np2 = pj.np(jim)

        assert (jim_np == jim_np2).all(), 'Error in jim2np() or np()'

        new_jim = pj.np2jim(jim_np)

        assert jim.pixops.isEqual(new_jim), 'Error in jim2np() or np2jim()'

        try:
            jim_empty = pj.Jim()
            _ = jim_empty.np()
            failed = True
        except ValueError:
            failed = False

        assert not failed, \
            'Error in catching a call of Jim.np() for non-dimensional Jims'

    def test_getters_setters(self):
        """Test getters and setters."""
        jim1 = pj.Jim(tiles[0])
        stats1 = jim1.stats.getStats()

        jim1[0, 0] = stats1['mean']
        first = jim1[0, 0]
        # test
        stats = first.stats.getStats()
        assert stats['max'] == stats['min'] == int(stats1['mean']), \
            'Error in jim[int, int] (either get or set item)'

        jim1[-1, -1] = stats1['max'] + 1
        stats = jim1.stats.getStats()
        # test
        assert stats['max'] == stats1['max'] + 1, \
            'Error in jim[int, int] (either get or set item)'

        last = jim1[-1, -1]
        # test
        stats = last.stats.getStats()
        assert stats['max'] == stats['min'] == stats1['max'] + 1, \
            'Error in jim[-int, -int] (either get or set item)'

        last = jim1[-5::2, -5::2]
        # test
        stats = last.stats.getStats()
        assert stats['max'] == stats1['max'] + 1, \
            'Error in jim[-int:-int:stride, -int:-int:stride] or jim[slice] ' \
            '(either get or set item)'

        try:
            _ = jim1['a', 'a']
            failed = True
        except IndexError:
            failed = False
        assert not failed, \
            'Error in catching wrong indices in jim[index, index]'

        try:
            _ = jim1[1, 'a']
            failed = True
        except IndexError:
            failed = False
        assert not failed, \
            'Error in catching wrong indices in jim[index, index]'

        jim1 = pj.Jim(ncol=256, nrow=256, nband=2, nplane=2)
        jim1.properties.setProjection('epsg:5514')
        jim1.geometry.cropBand(0)
        jim1.pixops.setData(5)
        stats1 = jim1.stats.getStats()

        # first = jim1[0, 0, 0, 0]
        first = jim1[0, 0, 0]
        # test
        stats = first.stats.getStats()
        assert stats['max'] == stats['min'] == stats1['mean'] == 5, \
            'Error in jim[int, int, int, int] (either get or set item)'
        assert first.properties.nrOfBand() == 1, \
            'Error in jim[int, int, int, int] (either get or set item, ' \
            'wrong nrOfBand)'
        assert first.properties.nrOfPlane() == 1, \
            'Error in jim[int, int, int, int] (either get or set item, ' \
            'wrong nrOfPlane)'
        assert first.properties.getProjection() == \
               jim1.properties.getProjection(), \
            'Error in jim[int, int, int] (either get or set item, ' \
            'projection not transmitted)'

        # last = jim1[-1, -1, -2, -2]
        last = jim1[-2, -1, -1]
        # test
        stats = last.stats.getStats()
        assert stats['max'] == stats['min'] == stats['mean'] == 5, \
            'Error in jim[-int, -int, -int, -int] (either get or set item)'
        assert last.properties.nrOfBand() == 1, \
            'Error in jim[-int, -int, -int, -int] (either get or set item, ' \
            'wrong nrOfBand)'
        assert last.properties.nrOfPlane() == 1, \
            'Error in jim[-int, -int, -int, -int] (either get or set item, ' \
            'wrong nrOfPlane)'
        assert last.properties.getProjection() == \
               jim1.properties.getProjection(), \
            'Error in jim[-int, -int, -int] (either get or set item, ' \
            'projection not transmitted)'

        # last = jim1[-1, -1, -2:-1:1, -2]
        # last = jim1[-2,-1, -1, -2:-1:1, -2]
        # stats = last.stats.getStats()
        # assert stats['max'] == stats['min'] == stats['mean'] == 5, \
        #     'Error in jim[-int, -int] (either get or set item)'
        # assert last.properties.nrOfBand() == 1, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfBand)'
        # assert last.properties.nrOfPlane() == 2, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfPlane)'

        try:
            _ = jim1[0, 0, 'a']
            failed = True
        except IndexError:
            failed = False
        assert not failed, \
            'Error in catching wrong indices in jim[index, index]'

        # last = jim1[-1, -1, -2:-1:1, -2:-1:1]
        # stats = last.stats.getStats()
        # assert stats['max'] == stats['min'] == stats['mean'] == 5, \
        #     'Error in jim[-int, -int] (either get or set item)'
        # assert last.properties.nrOfBand() == 2, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfBand)'
        # assert last.properties.nrOfPlane() == 2, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfPlane)'

    def test_operators(self):
        """Test basic operators (+, -, *, /, =)."""
        jim1 = pj.Jim(tiles[0])
        # test
        stats1 = jim1.stats.getStats()
        jim2 = pj.Jim(tiles[1])
        # test
        stats2 = jim2.stats.getStats()

        jim3 = jim1 + jim2
        # test
        stats3 = jim3.stats.getStats()
        max = stats3['max']
        min = stats3['min']

        assert max <= stats1['max'] + stats2['max'], \
            'Error in operation type Jim + Jim'

        jim3 += 1
        stats3 = jim3.stats.getStats()

        assert stats3['max'] == max + 1, 'Error in operation type Jim += int'

        jim3 += jim3
        stats3 = jim3.stats.getStats()

        assert stats3['max'] == (max + 1) * 2, \
            'Error in operation type Jim += Jim'
        assert stats3['min'] == (min + 1) * 2, \
            'Error in operation type Jim += Jim'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadBasicMethods)]
    return unittest.TestSuite(suite_list)
