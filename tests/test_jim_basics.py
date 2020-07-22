"""
test suite for overriden basic methods for Jim objects.
Author(s): Pieter.Kempeneers@ec.europa.eu,
           Ondrej Pesek,
           Pierre.Soille@ec.europa.eu
Copyright (C) 2018-2020 European Union (Joint Research Centre)

This file is part of pyjeo.

pyjeo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyjeo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.
"""

import pyjeo as pj
import unittest

import numpy as np
import warnings


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/modis_ndvi_training.sqlite'


class BadBasicMethods(unittest.TestCase):
    """Test functions and methods on the root level and operations for Jims."""

    @staticmethod
    def test_jim_creations():
        """Test creating of Jim objects."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(jim1, copy_data=True)
        jim3 = pj.Jim(jim1, copy_data=False)

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
        stats7 = jim7.stats.getStats(['min', 'max'], band=0)

        assert stats7['max'] < 2 and stats7['min'] > 0, \
            'Error in creating Jim with uniform distribution'

        jim8 = pj.Jim(nrow=500, ncol=500, otype='Float32', uniform=[0, 2],
                      seed=0)

        assert jim7.pixops.isEqual(jim8), \
            'Error in usage of seed keyword argument when creating Jim()'

        jim8_2 = pj.Jim(nrow=500, ncol=500, otype='Float32', uniform=2,
                        seed=0)

        assert jim7.pixops.isEqual(jim8_2), \
            'Error in using only one value for uniform creation of Jim()' \
            '(created Jim is not the same as the one created with[0, value])'

        try:
            _ = pj.Jim(nrow=5, ncol=5, otype='Float32', uniform=[0, 0, 0])
            raised = False
        except AttributeError:
            raised = True

        assert raised, \
            'Error in catching a call of Jim creation with wrong value ' \
            'parsed as the uniform argument (three values parsed)'

        # Test creation with mean and stdev

        jim9 = pj.Jim(nrow=500, ncol=500, otype='Float32', mean=500, stdev=50,
                      seed=0)

        stats = jim9.stats.getStats(['mean', 'stdev'])

        assert 500 - 50 < stats['mean'] < 500 + 50, \
            'Suspicious values in Jim when creating it with Jim(mean, stdev)' \
            '(mean seems too distant from the specified one)'
        assert 50 - 5 < stats['stdev'] < 50 + 5, \
            'Suspicious values in Jim when creating it with Jim(mean, stdev)' \
            '(stdev seems too different from the specified one)'

        # Test creation with only seed (mean=0, stdev=1)

        jim10 = pj.Jim(nrow=500, ncol=500, otype='Byte', seed=0)

        jim10_2 = pj.Jim(nrow=500, ncol=500, otype='Byte', mean=0, stdev=1,
                         seed=0)

        assert jim10.pixops.isEqual(jim10_2), \
            'Error when creating Jim with Jim(nrow, nncol, otype, seed) ' \
            '(not equal to Jim(nrow, nncol, otype, seed, mean=0, stdev=1))'

        try:
            _ = pj.Jim(seed=5)
            raised = False
        except AttributeError:
            raised = True

        assert raised, \
            'Error in catching a call of Jim creation with nonsense (kw)args'

        non_existing_path = pj._get_random_path()

        # TODO: Commented until issue #34 will be solved
        try:
            _ = pj.Jim(non_existing_path)
            raised = False
        except ValueError:
            raised = True

        assert raised, \
            'Error in catching a call of Jim creation with non-existing path'

    @staticmethod
    def test_numpy_conversions():
        """Test conversions to numpy and back."""
        jim = pj.Jim(tiles[0])

        jim_np = pj.jim2np(jim)
        jim_np2 = pj.np(jim)

        assert (jim_np == jim_np2).all(), 'Error in jim2np() or np()'

        new_jim = pj.np2jim(jim_np)

        assert jim.pixops.isEqual(new_jim), 'Error in jim2np() or np2jim()'

        anp = np.arange(2 * 100 * 256).reshape((2, 100, 256)).astype(
            np.float64)
        ajim = pj.np2jim(anp)
        nr_of_band = ajim.properties.nrOfBand()

        assert nr_of_band == 1
        assert ajim.properties.nrOfPlane() == anp.shape[0]
        assert ajim.properties.nrOfRow() == anp.shape[1]
        assert ajim.properties.nrOfCol() == anp.shape[2]

        try:
            jim_empty = pj.Jim()
            _ = jim_empty.np()
            raised = False
        except ValueError:
            raised = True

        assert raised, \
            'Error in catching a call of Jim.np() for non-dimensional Jims'

        nrow = ncol = nband = 5
        multib_jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, otype='Byte',
                            uniform=[0, 2], seed=0)

        band_last = multib_jim.np(nband - 1)
        cropped = pj.geometry.cropBand(multib_jim, nband - 1)

        assert (band_last == cropped.np()).all(), \
            'Error in Jim.np(band) ' \
            '(not returning the same object as cropBand(band).np())'

        band_last_minus = multib_jim.np(-1)

        assert (band_last == band_last_minus).all(), \
            'Error in Jim.np(-band) ' \
            '(Jim.np(-1) not returning the same object as ' \
            'Jim.np(Jim.properties.nrOfBand() - 1))'

        try:
            _ = ajim.np(ajim.properties.nrOfBand())
            raised = False
        except ValueError:
            raised = True

        assert raised, 'Error in catching a call of Jim.np(band) with ' \
                       'band value greater than number of bands'

    @staticmethod
    def test_getters_setters():
        """Test getters and setters."""
        jim1 = pj.Jim(tiles[0])
        vect = pj.JimVect(vector)

        stats1 = jim1.stats.getStats(band=0)

        jim1[0, 0] = stats1['mean']
        first = jim1[0, 0]
        # test
        stats = first.stats.getStats(band=0)
        assert stats['max'] == stats['min'] == int(stats1['mean']), \
            'Error in jim[int, int] (either get or set item)'

        jim1[-1, -1] = stats1['max'] + 1
        stats = jim1.stats.getStats(band=0)
        # test
        assert stats['max'] == stats1['max'] + 1, \
            'Error in jim[int, int] (either get or set item)'

        last = jim1[-1, -1]
        # test
        stats = last.stats.getStats(band=0)
        assert stats['max'] == stats['min'] == stats1['max'] + 1, \
            'Error in jim[-int, -int] (either get or set item)'

        last = jim1[-5::2, -5::2]
        # test
        stats = last.stats.getStats(band=0)
        assert stats['max'] == stats1['max'] + 1, \
            'Error in jim[-int:-int:stride, -int:-int:stride] or jim[slice] ' \
            '(either get or set item)'

        try:
            _ = jim1['a', 'a']
            raised = False
        except IndexError:
            raised = True
        assert raised, \
            'Error in catching wrong indices in jim[index, index]'

        try:
            _ = jim1[1, 'a']
            raised = False
        except IndexError:
            raised = True
        assert raised, \
            'Error in catching wrong indices in jim[index, index]'

        jim1 = pj.Jim(ncol=256, nrow=256, nband=2, nplane=2)
        jim1.properties.setProjection('epsg:5514')
        jim1.geometry.cropBand(0)
        jim1.pixops.setData(5)
        stats1 = jim1.stats.getStats(band=0)

        jim_same = jim1[:]

        assert jim1.pixops.isEqual(jim_same), \
            'Error in Jim[:] (get all items)'

        # first = jim1[0, 0, 0, 0]
        first = jim1[0, 0, 0]
        # test
        stats = first.stats.getStats(band=0)
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
        stats = last.stats.getStats(band=0)
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
        # stats = last.stats.getStats(band=0)
        # assert stats['max'] == stats['min'] == stats['mean'] == 5, \
        #     'Error in jim[-int, -int] (either get or set item)'
        # assert last.properties.nrOfBand() == 1, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfBand)'
        # assert last.properties.nrOfPlane() == 2, \
        #     'Error in jim[int, int, int, int] (either get or set item, ' \
        #     'wrong nrOfPlane)'

        jim1[0, 0] += 1
        stats = jim1.stats.getStats(band=0)

        assert stats['max'] == stats['min'] + 1, \
            'Error in jim[int, int] when jim[int, int] returns a 1-D array' \
            '(either get or set item)'

        try:
            _ = jim1[0, 0, 'a']
            raised = False
        except IndexError:
            raised = True
        assert raised, \
            'Error in catching wrong indices in jim[index, index]'

        # last = jim1[-1, -1, -2:-1:1, -2:-1:1]
        # stats = last.stats.getStats(band=0)
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
            raised = False
        except ValueError:
            raised = True
        assert raised, 'Error in catching a JimVect used as an index ' \
                       'for a multiplanar Jim (get item)'

        try:
            jim1[vect] = 5
            raised = False
        except ValueError:
            raised = True
        assert raised, 'Error in catching a JimVect used as an index ' \
                       'for a multiplanar Jim (set item)'

        modis = pj.Jim(testFile)
        modis.properties.clearNoData()
        modis_clipped = modis[vect]
        modis.properties.setNoDataVals(0)
        modis_clipped2 = modis[vect]

        bbox_vect = vect.properties.getBBox()
        bbox_clipped = modis_clipped.properties.getBBox()

        raised = True
        for i in range(len(bbox_vect)):
            if i % 2 == 0:
                delta = modis_clipped.properties.getDeltaX()
            else:
                delta = modis_clipped.properties.getDeltaY()

            if abs(bbox_vect[i] - bbox_clipped[i]) > delta:
                raised = False
                break

        assert raised, 'Error in clipping a Jim by JimVect (Jim[JimVect])'

        assert modis_clipped.properties.getNoDataVals() == \
               modis_clipped2.properties.getNoDataVals() == [0], \
            'Error in clipping a Jim by JimVect (Jim[JimVect]) (noData not ' \
            'correctly transferred)'

        modis[vect] = 5
        fives = modis[vect]
        fives_stats = fives.stats.getStats(band=0)

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
        stats = rand_jim.stats.getStats(band=0)
        twos = pj.Jim(nrow=50, ncol=50, uniform=[2, 2], otype='int8')

        twos_masked = twos[rand_jim]
        stats_masked = twos_masked.stats.getStats(band=0)

        assert stats_masked['max'] == stats['max'] * 2 and \
               stats_masked['mean'] == stats['mean'] * 2 and \
               stats_masked['min'] == stats['min'] == 0, \
            'Error in masking a Jim by Jim (Jim1[Jim2])'

        # Test a nonsense argument in [gs]etters
        try:
            rand_jim['a'] = 5
            raised = False
        except ValueError:
            raised = True

        assert raised, 'Error in catching wrong indices like Jim["string"]'

    @staticmethod
    def test_operators():
        """Test basic operators (+, -, *, /, =, abs(), ~)."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])

        for jim_one, jim_two in [(jim1, jim2), (pj.geometry.stackBand(jim1,
                                                                      jim1),
                                                pj.geometry.stackBand(jim2,
                                                                      jim2))]:
            nr_of_bands = jim_one.properties.nrOfBand()
            empty = pj.Jim(nrow=jim1.properties.nrOfRow(),
                           ncol=jim1.properties.nrOfCol(),
                           nband=nr_of_bands,
                           otype='int32')
            fifteens = pj.Jim(nrow=jim1.properties.nrOfRow(),
                              ncol=jim1.properties.nrOfCol(),
                              nband=nr_of_bands,
                              otype='float32')
            fifteens.pixops.setData(15)

            stats1 = jim_one.stats.getStats(band=0)
            stats2 = jim_two.stats.getStats(band=0)

            jim3 = jim_one + jim_two
            stats3 = jim3.stats.getStats(band=0)
            max = stats3['max']
            min = stats3['min']

            assert max <= stats1['max'] + stats2['max'], \
                'Error in operation type Jim + Jim'

            if nr_of_bands == 2:
                stats3 = jim3.stats.getStats(band=1)
                max = stats3['max']

                assert max <= stats1['max'] + stats2['max'], \
                    'Error in Jim + Jim (not performed for all bands)'

            # Test +=

            jim3 += 1
            stats3 = jim3.stats.getStats(band=0)

            assert stats3['max'] == max + 1, \
                'Error in operation type Jim += int'

            if nr_of_bands == 2:
                stats3 = jim3.stats.getStats(band=1)

                assert stats3['max'] == max + 1, \
                    'Error in Jim += int (not performed for all bands)'

            jim3 += jim3
            stats3 = jim3.stats.getStats(band=0)

            assert stats3['max'] == (max + 1) * 2, \
                'Error in operation type Jim += Jim'
            assert stats3['min'] == (min + 1) * 2, \
                'Error in operation type Jim += Jim'

            if nr_of_bands == 2:
                stats3 = jim3.stats.getStats(band=1)

                assert stats3['max'] == (max + 1) * 2, \
                    'Error in Jim += Jim (not performed for all bands)'
                assert stats3['min'] == (min + 1) * 2, \
                    'Error in Jim += Jim (not performed for all bands)'

            # Test * and *=

            fifteen_jim3 = pj.pixops.convert(jim3, 'float32') * fifteens
            stats3 = fifteen_jim3.stats.getStats(band=0)

            assert stats3['max'] == (max + 1) * 30, \
                'Error in operation type Jim * Jim'
            assert stats3['min'] == (min + 1) * 30, \
                'Error in operation type Jim * Jim'

            if nr_of_bands == 2:
                stats3 = fifteen_jim3.stats.getStats(band=1)

                assert stats3['max'] == (max + 1) * 30, \
                    'Error in Jim += Jim (not performed for all bands)'
                assert stats3['min'] == (min + 1) * 30, \
                    'Error in Jim += Jim (not performed for all bands)'

            fifteen_jim3 *= empty
            stats3 = fifteen_jim3.stats.getStats(band=0)

            assert stats3['min'] == stats3['max'] == stats3['mean'] == 0, \
                'Error in operation type Jim += Jim'

            if nr_of_bands == 2:
                stats3 = fifteen_jim3.stats.getStats(band=1)

                assert stats3['min'] == stats3['max'] == stats3['mean'] == 0, \
                    'Error in Jim += Jim (not performed for all bands)'

            # Test specialities like __neg__, abs(), ~, ...

            zeros = jim3 + -jim3

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
                _ = jim_one / 2
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching the error from dividing an int Jim'

            minus_ones.pixops.convert('float32')

            halves = minus_ones / 2
            stats = halves.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
                'Error in Jim / number'

            if nr_of_bands == 2:
                stats = halves.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
                    'Error in Jim / number (not performed for all bands)'

            halves /= minus_ones
            stats = halves.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
                'Error in Jim /= Jim'

            if nr_of_bands == 2:
                stats = halves.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
                    'Error in Jim /= Jim (not performed for all bands)'

            halves /= -1
            stats = halves.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
                'Error in Jim /= number'

            if nr_of_bands == 2:
                stats = halves.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == -0.5, \
                    'Error in Jim /= Jim (not performed for all bands)'

            halves = halves / minus_ones
            stats = halves.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
                'Error in Jim / Jim'

            if nr_of_bands == 2:
                stats = halves.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 0.5, \
                    'Error in Jim / Jim (not performed for all bands)'

            # Test modulo

            sevens = fifteens % 8
            stats = sevens.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 7, \
                'Error in Jim % number'

            if nr_of_bands == 2:
                stats = sevens.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 7, \
                    'Error in Jim % number (not performed for all bands)'

            test = pj.Jim(sevens)
            test %= 4
            stats = test.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 3, \
                'Error in Jim %= number'

            if nr_of_bands == 2:
                stats = test.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 3, \
                    'Error in Jim %= number (not performed for all bands)'

            ones = fifteens % sevens
            stats = ones.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim % Jim'

            if nr_of_bands == 2:
                stats = ones.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim % Jim (not performed for all bands)'

            test2 = sevens
            test2 %= test
            stats = test2.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim %= Jim'

            if nr_of_bands == 2:
                stats = test2.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim %= Jim (not performed for all bands)'

            # Test powering

            nines = test ** 2
            stats = nines.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 9, \
                'Error in Jim ** number'

            if nr_of_bands == 2:
                stats = nines.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 9, \
                    'Error in Jim ** number (not performed for all bands)'

            test **= 2

            assert test.pixops.isEqual(nines), 'Error in Jim **= number'

            # Test shifts
            test.pixops.convert('int32')

            seventy_twos = test << 3
            stats = seventy_twos.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 72, \
                'Error in Jim << number'

            if nr_of_bands == 2:
                stats = seventy_twos.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 72, \
                    'Error in Jim << number (not performed for all bands)'

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
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of << operation'

            try:
                _ = test >> 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of >> operation'

            try:
                test <<= 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of <<= operation'

            try:
                test >>= 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
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
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of | operation'

            try:
                ones |= 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of |= operation'

            try:
                _ = 'a' | ones
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong left side of | operation'

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
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of ^ operation'

            try:
                ones ^= 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of ^= operation'

            try:
                _ = 'a' ^ ones
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong left side of ^ operation'

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
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of & operation'

            try:
                ones &= 'a'
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong right side of &= operation'

            try:
                _ = 'a' & ones
                raised = False
            except TypeError:
                raised = True

            assert raised, \
                'Error in catching wrong left side of & operation'

    @staticmethod
    def test_pixel_wise_conditions():
        """Test conditions like ==, !=, >, >=, <, <= for Jims."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(jim1)

        for jim_one, jim_two in [(jim1, jim2), (pj.geometry.stackBand(jim1,
                                                                      jim1),
                                                pj.geometry.stackBand(jim2,
                                                                      jim2))]:
            nr_of_bands = jim_one.properties.nrOfBand()

            jim_equality = jim_one == jim_two
            stats = jim_equality.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim == Jim'

            if nr_of_bands == 2:
                stats = jim_equality.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim == Jim (not performed for all bands)'

            jim_not_equality = jim_one != jim_two
            stats = jim_not_equality.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                'Error in Jim != Jim'

            if nr_of_bands == 2:
                stats = jim_not_equality.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                    'Error in Jim != Jim (not performed for all bands)'

            jim_greater = jim_one > jim_two
            stats = jim_greater.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                'Error in Jim > Jim'

            if nr_of_bands == 2:
                stats = jim_greater.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                    'Error in Jim > Jim (not performed for all bands)'

            jim_greatere = jim_one >= jim_two
            stats = jim_greatere.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim >= Jim'

            if nr_of_bands == 2:
                stats = jim_greatere.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim >= Jim (not performed for all bands)'

            jim_lesser = jim_one < jim_two
            stats = jim_lesser.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                'Error in Jim < Jim'

            if nr_of_bands == 2:
                stats = jim_lesser.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 0, \
                    'Error in Jim < Jim (not performed for all bands)'

            jim_lessere = jim_one <= jim_two
            stats = jim_lessere.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim <= Jim'

            if nr_of_bands == 2:
                stats = jim_lessere.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim <= Jim (not performed for all bands)'

            jim_equality_int = jim_equality == 1
            stats = jim_equality_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim == number'

            if nr_of_bands == 2:
                stats = jim_equality_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim == Jim (not performed for all bands)'

            jim_not_equality_int = jim_not_equality != 1
            stats = jim_not_equality_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim != number'

            if nr_of_bands == 2:
                stats = jim_not_equality_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim != number (not performed for all bands)'

            jim_greater_int = jim_greater > -1
            stats = jim_greater_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim > number'

            if nr_of_bands == 2:
                stats = jim_greater_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim > number (not performed for all bands)'

            jim_greatere_int = jim_greatere >= 1
            stats = jim_greatere_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim >= number'

            if nr_of_bands == 2:
                stats = jim_greatere_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim >= number (not performed for all bands)'

            jim_lesser_int = jim_lesser < 1
            stats = jim_lesser_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim < number'

            if nr_of_bands == 2:
                stats = jim_lesser_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim < number (not performed for all bands)'

            jim_lessere_int = jim_lessere <= 1
            stats = jim_lessere_int.stats.getStats(band=0)

            assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                'Error in Jim <= number'

            if nr_of_bands == 2:
                stats = jim_lessere_int.stats.getStats(band=1)

                assert stats['max'] == stats['min'] == stats['mean'] == 1, \
                    'Error in Jim <= number (not performed for all bands)'

        # Now check it also when just some of the values are equal
        jim2[0:10, 0:10] += 1

        jim_equality = jim1 == jim2
        stats = jim_equality.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim == Jim'

        jim_not_equality = jim1 != jim2
        stats = jim_not_equality.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim != Jim'

        jim_greater = jim2 > jim1
        stats = jim_greater.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim > Jim'

        jim_greatere = jim1 >= jim2
        stats = jim_greatere.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim >= Jim'

        jim_lesser = jim1 < jim2
        stats = jim_lesser.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim < Jim'

        jim_lessere = jim2 <= jim1
        stats = jim_lessere.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim <= Jim'

        jim_equality_int = jim_equality == 1
        stats = jim_equality_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim == number'

        jim_not_equality_int = jim_not_equality != 1
        stats = jim_not_equality_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim != number'

        jim_greater_int = jim_greater > 0
        stats = jim_greater_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim > number'

        jim_greatere_int = jim_greatere >= 1
        stats = jim_greatere_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim >= number'

        jim_lesser_int = jim_lesser < 1
        stats = jim_lesser_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim < number'

        jim_lessere_int = jim_lessere <= 0
        stats = jim_lessere_int.stats.getStats(band=0)

        assert stats['max'] == 1 and stats['min'] == 0, \
            'Error in Jim <= number'

        # Test inner checks in operators

        jim_twobands = pj.Jim(ncol=5, nrow=5, nband=2, otype='GDT_Byte')
        jim_threebands = pj.Jim(ncol=5, nrow=5, nband=3, otype='GDT_Byte')

        try:
            jim_threebands == jim_twobands
            raised = False
        except IndexError:
            raised = True

        assert raised, 'Error in raising an error when left Jim has ' \
                       'more bands than the right one for basic conditions'

        warnings.filterwarnings('error', category=Warning)
        try:
            jim_twobands == jim_threebands
            raised = False
        except Warning:
            raised = True

        assert raised, 'Error in raising a warning when left Jim has ' \
                       'less bands than the right one for basic conditions'

        warnings.resetwarnings()
        warnings.filterwarnings('ignore', category=Warning)

        eq = jim_twobands == jim_threebands

        assert eq.properties.nrOfBand() == 2, \
            'Error in equations with the left Jim with less bands ' \
            '(the output Jim does not have the same number of bands)'

        warnings.resetwarnings()

    @staticmethod
    def test_checks():
        """Test checks of arguments appearing behind the scene."""
        jim1 = pj.Jim(tiles[0])
        jim1.pixops.convert('Byte')

        try:
            _ = pj.ccops.labelConstrainedCCsVariance(jim1, 0, 0, 0, 0, 0, 0,
                                                     pj.Jim(graph=0))
            raised = False
        except ValueError:
            raised = True

        assert raised, 'Error in raising an error when an invalid value ' \
                       'is passed as a graph in a function (for example ' \
                       'ccops.labelConstrainedCCsVariance())'

    @staticmethod
    def test_args_different_formats():
        """Test the parsing of arguments as bytes, unicode, etc."""
        jim0 = pj.Jim(tiles[0])
        jim1 = pj.Jim(u'tests/data/red1.tif')
        jim2 = pj.Jim(b'tests/data/red1.tif')

        assert jim1.pixops.isEqual(jim0), \
            'Error in creating Jim object with unicode path'
        assert jim2.pixops.isEqual(jim0), \
            'Error in creating Jim object with byte path'

        jim0 = pj.Jim(nrow=500, ncol=500, otype='Float32', uniform=[0, 2],
                      seed=0)
        jim1 = pj.Jim(nrow=500, ncol=500, otype=u'Float32', uniform=[0, 2],
                      seed=0)
        jim2 = pj.Jim(nrow=500, ncol=500, otype=b'Float32', uniform=[0, 2],
                      seed=0)

        assert jim1.pixops.isEqual(jim0), \
            'Error in the usage of unicode strings inside arguments'
        assert jim2.pixops.isEqual(jim0), \
            'Error in the usage of byte strings inside arguments'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadBasicMethods)]
    return unittest.TestSuite(suite_list)
