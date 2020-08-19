"""Test suite for module pyjeo.stats."""
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


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


class BadStats(unittest.TestCase):
    """Test functions and methods for getting statistics."""

    @staticmethod
    def test_getStats():
        """Test if values from getStats are not suspicious."""
        jim = pj.Jim(tiles[0])
        stats = jim.stats.getStats(band=0)
        stats2 = pj.stats.getStats(jim, band=0)

        assert stats == stats2, 'Inconsistency in getStats() (method returns' \
                                ' different result than function)'

        max = stats['max']
        min = stats['min']
        mean = stats['mean']

        jim_min_dict = jim.stats.getStats(function='min', band=0)

        for key, val in jim_min_dict.items():
            assert val == min, \
                'Error in getting statistics with stats.getStats()'

        assert min < mean < max, \
            'Error in getting statistics with stats.getStats()'

        jim.properties.clearNoData()

        stats = jim.stats.getStats('median,invalid', band=0)
        stats2 = pj.stats.getStats(jim, 'median,invalid', band=0)

        assert stats == stats2, 'Inconsistency in getStats() (method returns' \
                                ' different result than function)'

        assert stats['ninvalid'] == 0, \
            'Error in properties.clearNoData() or ' \
            'getStats(function="invalid") (nodata detected after clearing)'

        stats = jim.stats.getStats(['max', 'median', 'invalid'], nodata=max,
                                   band=0)
        stats2 = pj.stats.getStats(jim, ['max', 'median', 'invalid'],
                                   nodata=max, band=0)

        assert stats == stats2, 'Inconsistency in getStats() (method returns' \
                                ' different result than function)'

        assert stats['max'] < max, \
            'Error in using nodata kwarg for stats.getStats()'

        assert stats['max'] > stats['median'] > min, \
            'Suspicious value of median with getStats() ' \
            '(not smaller than max or not bigger than min)'

        assert stats['ninvalid'] > 0, 'Error in getStats(function="invalid")' \
                                      ' (no nodata detected, but should be)'

    @staticmethod
    def test_histograms():
        """Test that values of histograms are not suspicious."""
        jim1 = pj.Jim(tiles[0])

        jim1_rows = jim1.properties.nrOfRow()
        jim1_cols = jim1.properties.nrOfCol()

        # Test getHisto1d()
        histo1d = jim1.stats.getHisto1d()

        assert histo1d.stats.getStats(band=0)['min'] >= 0, \
            'Error in computing 1D histogram'

        assert histo1d.stats.getStats(band=0)['max'] <= jim1_rows * jim1_cols,\
            'Error in computing 1D histogram'

        assert histo1d.properties.isEqual(pj.stats.getHisto1d(jim1)), \
            'Function and method getHisto1d() return different results'

        jim2 = pj.Jim(tiles[0][:-8] + 'nir' + tiles[0][-5] + '.tif')

        # Test getHisto2d()
        histo2d = jim1.stats.getHisto2d(jim2)

        assert histo2d.stats.getStats(band=0)['min'] >= 0, \
            'Error in computing 2D histogram'

        assert histo2d.stats.getStats(band=0)['max'] <= jim1_rows * jim1_cols,\
            'Error in computing 2D histogram'

        assert histo2d.properties.isEqual(pj.stats.getHisto2d(jim1, jim2)), \
            'Function and method getHisto2d() return different results'

        # TODO: Cover histo3d()

        # Test getHistoCumulative
        try:
            _ = jim1.stats.getHistoCumulative()
            raised = False
        except pj.exceptions.JimInnerParametersError:
            raised = True
        assert raised, 'Error in catching wrong data type in ' \
                       'stats.getHistoCumulative()'

        # histo_cumul = pj.pixops.convert(jim1, 4).stats.getHistoCumulative()
        jim2 = pj.pixops.convert(jim1, 'Int32')
        histo_cumul = jim2.stats.getHistoCumulative()

        assert histo_cumul.properties.nrOfCol() == jim1_rows * jim1_cols + 1, \
            'Error in stats.getHistoCumulative() ' \
            '(nrOfCol != nrOfRow*Col of original)'

    @staticmethod
    def test_stretch():
        """Test stretching a Jim object."""
        jim = pj.Jim(testFile, band=[0, 1])

        jim_min = jim.stats.getStats(band=0)['min']

        jim_stretched = pj.stats.stretch(jim, dst_min=20)

        if jim_min != 20:
            assert jim_stretched.stats.getStats(band=0)['min'] != jim_min, \
                'Error in stretching Jim'
            assert not jim.properties.isEqual(jim_stretched), \
                'Error in stretching Jim'

        jim.stats.stretch(dst_min=20)

        assert jim.properties.isEqual(jim_stretched), 'Error in stretching Jim'

    @staticmethod
    def test_getStatProfile():
        """Test if values from getStatProfile are not wrong."""
        jim = pj.Jim(tiles[0])

        stats = jim.stats.getStats(band=0)

        min = stats['min']
        max = stats['max']

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        min_profile = pj.stats.getStatProfile(jim, 'min')
        jim.stats.statProfile('max')

        assert min_profile.stats.getStats(band=0)['min'] == min, \
            'Error in stats.getStatProfile()'

        assert jim.stats.getStats(band=0)['max'] == max, \
            'Error in stats.getStatProfile()'


class BadStatsLists(unittest.TestCase):
    """Test JimList functions and methods for getting statistics."""

    @staticmethod
    def test_getStats():
        """Test if values from getStats are not suspicious."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])
        jiml = pj.JimList([jim1, jim2])

        stats = jiml.stats.getStats(band=0)
        stats2 = pj.stats.getStats(jiml, band=0)

        assert stats == stats2, \
            'Inconsistency in JimList.stats.getStats() ' \
            '(method returns different result than function)'

        max = stats['max']
        min = stats['min']
        mean = stats['mean']

        jim_minmax_dict = jiml.stats.getStats(function=['min', 'max'], band=0)

        assert jim_minmax_dict['min'] == min, \
            'Error in getting statistics with JimList.stats.getStats()' \
            '(min value from getStats() and ' \
            'getStats(function=["min", "max"]) do not correspond)'

        assert jim_minmax_dict['max'] == max, \
            'Error in getting statistics with JimList.stats.getStats()' \
            '(max value from getStats() and ' \
            'getStats(function=["min", "max"]) do not correspond)'

        assert min < mean < max, \
            'Error in getting statistics with JimList.stats.getStats()'

        jiml.properties.clearNoData()

        stats = jiml.stats.getStats('median,invalid', band=0)
        stats2 = pj.stats.getStats(jiml, 'median,invalid', band=0)

        assert stats == stats2, \
            'Inconsistency in JimList.stats.getStats() ' \
            '(method returns different result than function)'

        assert stats['ninvalid'] == 0, \
            'Error in JimList.properties.clearNoData() or ' \
            'JimList.stats.getStats(function="invalid") ' \
            '(nodata detected after clearing)'

        stats = jiml.stats.getStats(['max', 'median', 'invalid'], nodata=max,
                                    band=0)
        stats2 = pj.stats.getStats(jiml, ['max', 'median', 'invalid'],
                                   nodata=max, band=0)

        assert stats == stats2, \
            'Inconsistency in JimList.stats.getStats() ' \
            '(method returns different result than function)'

        assert stats['max'] < max, \
            'Error in using nodata kwarg for JimList.stats.getStats()'

        assert stats['max'] > stats['median'] > min, \
            'Suspicious value of median with JimList.stats.getStats() ' \
            '(not smaller than max or not bigger than min)'

        assert stats['ninvalid'] > 0, \
            'Error in JimList.stats.getStats(function="invalid") ' \
            '(no nodata detected, but should be)'

    @staticmethod
    def test_getStatProfile():
        """Test if values from getStatProfile are not wrong."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])
        jiml = pj.JimList([jim1, jim2])

        stats1 = jiml.stats.getStats(['min', 'max'], band=0)
        stats2 = jim2.stats.getStats(band=0)

        min = stats1['min']
        if stats2['min'] < min:
            min = stats2['min']

        max = stats1['max']
        if stats2['max'] > max:
            max = stats2['max']

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        min_profile = jiml.stats.getStatProfile('min')
        max_profile = jiml.stats.getStatProfile('max')

        assert min_profile.stats.getStats(band=0)['min'] == min, \
            'Error in stats.getStatProfile()'

        assert max_profile.stats.getStats(band=0)['max'] == max, \
            'Error in stats.getStatProfile()'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadStats),
                  loader.loadTestsFromTestCase(BadStatsLists)]
    return unittest.TestSuite(suite_list)
