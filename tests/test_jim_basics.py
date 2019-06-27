"""Test suite for overriden basic methods for Jim objects."""

import pyjeo as pj
import unittest
import numpy as np


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/modis_ndvi_training.sqlite'


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

        anp = np.arange(2 * 100 * 256).reshape((2, 100, 256)).astype(np.float64)
        ajim = pj.np2jim(anp)

        assert ajim.properties.nrOfBand() == 1
        assert ajim.properties.nrOfPlane() == anp.shape[0]
        assert ajim.properties.nrOfRow() == anp.shape[1]
        assert ajim.properties.nrOfCol() == anp.shape[2]

        try:
            jim_empty = pj.Jim()
            _ = jim_empty.np()
            failed = True
        except ValueError:
            failed = False

        assert not failed, \
            'Error in catching a call of Jim.np() for non-dimensional Jims'

        try:
            _ = ajim.np(ajim.properties.nrOfBand() - 1)
            failed = False
        except:
            failed = True

        assert not failed, 'Error in Jim.np(band)'

        try:
            _ = ajim.np(ajim.properties.nrOfBand())
            failed = True
        except ValueError:
            failed = False

        assert not failed, 'Error in catching a call of Jim.np(band) with ' \
                           'band value greater than number of bands'

    def test_getters_setters(self):
        """Test getters and setters."""
        jim1 = pj.Jim(tiles[0])
        vect = pj.JimVect(vector)

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

        jim_same = jim1[:]

        assert jim1.pixops.isEqual(jim_same), \
            'Error in Jim[:] (get all items)'

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

        jim1[0, 0] += 1
        stats = jim1.stats.getStats()

        assert stats['max'] == stats['min'] + 1, \
            'Error in jim[int, int] when jim[int, int] returns a 1-D array' \
            '(either get or set item)'

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

        # Test JimVect usage in getters and setters as an argument

        try:
            _ = jim1[vect]
            failed = True
        except ValueError:
            failed = False
        assert not failed, 'Error in catching a JimVect used as an index ' \
                           'for a multiplanar Jim (get item)'

        try:
            jim1[vect] = 5
            failed = True
        except ValueError:
            failed = False
        assert not failed, 'Error in catching a JimVect used as an index ' \
                           'for a multiplanar Jim (set item)'

        modis = pj.Jim(testFile)
        modis.properties.clearNoData()
        modis_clipped = modis[vect]
        modis.properties.setNoDataVals(0)
        modis_clipped2 = modis[vect]

        bbox_vect = vect.properties.getBBox()
        bbox_clipped = modis_clipped.properties.getBBox()

        failed = False
        for i in range(len(bbox_vect)):
            if i % 2 == 0:
                delta = modis_clipped.properties.getDeltaX()
            else:
                delta = modis_clipped.properties.getDeltaY()

            if abs(bbox_vect[i] - bbox_clipped[i]) > delta:
                failed = True
                break

        assert not failed, 'Error in clipping a Jim by JimVect (Jim[JimVect])'

        assert modis_clipped.properties.getNoDataVals() == \
               modis_clipped2.properties.getNoDataVals() == [0], \
            'Error in clipping a Jim by JimVect (Jim[JimVect]) (noData not ' \
            'correctly transferred)'

        modis[vect] = 5
        fives = modis[vect]
        fives_stats = fives.stats.getStats()

        assert fives_stats['max'] == 5 and fives_stats['min'] in [0, 5], \
            'Error in using JimVect as an argument in Jim[JimVect] = number'

        modis[vect] = modis + 5

        assert modis[vect].pixops.isEqual(fives * 2) and \
               not modis[0, 0].pixops.isEqual(fives[0, 0]), \
            'Error in using JimVect as an argument in Jim[JimVect] = Jim'

        # Test Jim usage in getters and setters as an argument

        rand_jim = pj.Jim(nrow=50, ncol=50, uniform=[0, 2], otype='int8')
        rand_jim[0, 0] = 0
        rand_jim[0, 1] = 1
        stats = rand_jim.stats.getStats()
        twos = pj.Jim(nrow=50, ncol=50, uniform=[2, 2], otype='int8')

        twos_masked = twos[rand_jim]
        stats_masked = twos_masked.stats.getStats()

        assert stats_masked['max'] == stats['max'] * 2 and \
               stats_masked['mean'] == stats['mean'] * 2 and \
               stats_masked['min'] == stats['min'] == 0, \
            'Error in masking a Jim by Jim (Jim1[Jim2])'

        # Test a nonsense argument in [gs]etters
        try:
            rand_jim['a'] = 5
            failed = True
        except ValueError:
            failed = False

        assert not failed, 'Error in catching wrong indices like Jim["string"]'

    def test_operators(self):
        """Test basic operators (+, -, *, /, =, abs(), ~)."""
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

        # Test +=

        jim3 += 1
        stats3 = jim3.stats.getStats()

        assert stats3['max'] == max + 1, 'Error in operation type Jim += int'

        jim3 += jim3
        stats3 = jim3.stats.getStats()

        assert stats3['max'] == (max + 1) * 2, \
            'Error in operation type Jim += Jim'
        assert stats3['min'] == (min + 1) * 2, \
            'Error in operation type Jim += Jim'

        # Test specialities like __neg__, abs(), ~, ...

        zeros = jim3 + -jim3
        empty = pj.Jim(nrow=jim3.properties.nrOfRow(),
                       ncol=jim3.properties.nrOfCol(),
                       otype='int32')

        assert zeros.pixops.isEqual(empty), \
            'Error in -Jim (not returning negative values)'

        minus_ones = empty - 1

        jim3.pixops.convert('int32')

        jim3_plus_one = jim3 + abs(minus_ones)

        assert jim3_plus_one.pixops.isEqual(jim3 + 1), \
            'Error in abs(Jim)'

        assert (~jim3).pixops.isEqual(-1 - jim3), \
            'Error in a bit-wise inversion (~Jim)'

        # Test __radd__

        assert jim3_plus_one.pixops.isEqual(1 + jim3), \
            'Error in Jim.__radd__() (number + Jim)'

        # Test __rmul__ and __rsub__

        assert (-jim3 - jim3).pixops.isEqual(-2 * jim3), \
            'Error in operation of type Jim - Jim or number * Jim '

        # Test -=

        jim3_plus_one -= 1

        assert jim3.pixops.isEqual(jim3_plus_one), \
            'Error in operation of type Jim -= number'

        jim3_plus_one -= jim3

        assert jim3_plus_one.pixops.isEqual(empty), \
            'Error in operation of type Jim -= Jim'

        # Test divisions

        try:
            _ = jim1 / 2
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching the error from dividing an int Jim'

        minus_ones.pixops.convert('float32')

        halves = minus_ones / 2
        stats = halves.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
            'Error in Jim / number'

        halves /= minus_ones
        stats = halves.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
            'Error in Jim /= Jim'

        halves /= -1
        stats = halves.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
            'Error in Jim /= number'

        halves = halves / minus_ones
        stats = halves.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
            'Error in Jim / Jim'

        # Test modulo

        fifteens = pj.Jim(nrow=jim3.properties.nrOfRow(),
                          ncol=jim3.properties.nrOfCol(),
                          otype='float32')
        fifteens.pixops.setData(15)

        sevens = fifteens % 8
        stats = sevens.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 7, \
            'Error in Jim % number'

        test = pj.Jim(sevens)
        test %= 4
        stats = test.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 3, \
            'Error in Jim %= number'

        ones = fifteens % sevens
        stats = ones.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim % Jim'

        test2 = sevens
        test2 %= test
        stats = test2.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim %= Jim'

        # Test powering

        nines = test ** 2
        stats = nines.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 9, \
            'Error in Jim ** number'

        test **= 2

        assert test.pixops.isEqual(nines), 'Error in Jim **= number'

        # Test shifts
        test.pixops.convert('int32')

        seventy_twos = test << 3
        stats = seventy_twos.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 72, \
            'Error in Jim << number'

        test = seventy_twos >> 3

        assert test.pixops.isEqual(nines), \
            'Error in Jim >> number'

        test = pj.Jim(seventy_twos)
        test >>= 3

        assert test.pixops.isEqual(nines), \
            'Error in Jim >>= number'

        test <<= 3

        assert test.pixops.isEqual(seventy_twos), \
            'Error in Jim <<= number'

        try:
            _ = test << 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching wrong right side of << operation'

        try:
            _ = test >> 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching wrong right side of >> operation'

        try:
            test <<= 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching wrong right side of <<= operation'

        try:
            test >>= 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching wrong right side of >>= operation'

        # Test |
        ones.pixops.convert('int32')

        assert (seventy_twos | 1).pixops.isEqual(73 * ones), \
            'Error in operation of type Jim | number'

        assert (1 | seventy_twos).pixops.isEqual(73 * ones), \
            'Error in operation of type number | Jim'

        test |= 1

        assert test.pixops.isEqual(73 * ones), \
            'Error in operation of type Jim |= number'

        assert (seventy_twos | ones).pixops.isEqual(73 * ones), \
            'Error in operation of type Jim | Jim'

        test |= 9 * ones

        assert test.pixops.isEqual(73 * ones), \
            'Error in operation of type Jim |= Jim'

        try:
            _ = ones | 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of | operation'

        try:
            ones |= 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of |= operation'

        try:
            _ = 'a' | ones
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong left side of | operation'

        # Test ^

        assert (seventy_twos ^ 9).pixops.isEqual(65 * ones), \
            'Error in operation of type Jim ^ number'

        assert (9 ^ seventy_twos).pixops.isEqual(65 * ones), \
            'Error in operation of type number ^ Jim'

        test ^= 10

        assert test.pixops.isEqual(67 * ones), \
            'Error in operation of type Jim ^= number'

        assert (seventy_twos ^ (9 * ones)).pixops.isEqual(65 * ones), \
            'Error in operation of type Jim ^ Jim'

        test ^= (68 * ones)

        assert test.pixops.isEqual(7 * ones), \
            'Error in operation of type Jim ^= Jim'

        try:
            _ = ones ^ 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of ^ operation'

        try:
            ones ^= 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of ^= operation'

        try:
            _ = 'a' ^ ones
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong left side of ^ operation'

        # Test &

        assert (seventy_twos & 9).pixops.isEqual(8 * ones), \
            'Error in operation of type Jim & number'

        assert (9 & seventy_twos).pixops.isEqual(8 * ones), \
            'Error in operation of type number & Jim'

        test &= 9

        assert test.pixops.isEqual(ones), \
            'Error in operation of type Jim &= number'

        assert (seventy_twos & (9 * ones)).pixops.isEqual(8 * ones), \
            'Error in operation of type Jim & Jim'

        test += 10
        test &= (5 * ones)

        assert test.pixops.isEqual(ones), \
            'Error in operation of type Jim &= Jim'

        try:
            _ = ones & 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of & operation'

        try:
            ones &= 'a'
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong right side of &= operation'

        try:
            _ = 'a' & ones
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in catching wrong left side of & operation'

    def test_pixel_wise_conditions(self):
        """Test conditions like ==, !=, >, >=, <, <= for Jims."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(jim1)

        jim_equality = jim1 == jim2
        stats = jim_equality.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim == Jim'

        jim_not_equality = jim1 != jim2
        stats = jim_not_equality.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 0, \
            'Error in Jim != Jim'

        jim_greater = jim1 > jim2
        stats = jim_greater.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 0, \
            'Error in Jim > Jim'

        jim_greatere = jim1 >= jim2
        stats = jim_greatere.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim >= Jim'

        jim_lesser = jim1 < jim2
        stats = jim_lesser.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 0, \
            'Error in Jim < Jim'

        jim_lessere = jim1 <= jim2
        stats = jim_lessere.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim <= Jim'

        jim_equality_int = jim_equality == 1
        stats = jim_equality_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim == number'

        jim_not_equality_int = jim_not_equality != 1
        stats = jim_not_equality_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim != number'

        jim_greater_int = jim_greater > -1
        stats = jim_greater_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim > number'

        jim_greatere_int = jim_greatere >= 1
        stats = jim_greatere_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim >= number'

        jim_lesser_int = jim_lesser < 1
        stats = jim_lesser_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim < number'

        jim_lessere_int = jim_lessere <= 1
        stats = jim_lessere_int.stats.getStats()

        assert stats['max'] == stats['min'] == stats['mean'] == 1, \
            'Error in Jim <= number'

        # Now check it also when just some of the values are equal
        jim2[0:10, 0:10] += 1

        jim_equality = jim1 == jim2
        stats = jim_equality.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim == Jim'

        jim_not_equality = jim1 != jim2
        stats = jim_not_equality.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim != Jim'

        jim_greater = jim2 > jim1
        stats = jim_greater.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim > Jim'

        jim_greatere = jim1 >= jim2
        stats = jim_greatere.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim >= Jim'

        jim_lesser = jim1 < jim2
        stats = jim_lesser.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim < Jim'

        jim_lessere = jim2 <= jim1
        stats = jim_lessere.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim <= Jim'

        jim_equality_int = jim_equality == 1
        stats = jim_equality_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim == number'

        jim_not_equality_int = jim_not_equality != 1
        stats = jim_not_equality_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim != number'

        jim_greater_int = jim_greater > 0
        stats = jim_greater_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim > number'

        jim_greatere_int = jim_greatere >= 1
        stats = jim_greatere_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim >= number'

        jim_lesser_int = jim_lesser < 1
        stats = jim_lesser_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim < number'

        jim_lessere_int = jim_lessere <= 0
        stats = jim_lessere_int.stats.getStats()

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim <= number'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadBasicMethods)]
    return unittest.TestSuite(suite_list)
