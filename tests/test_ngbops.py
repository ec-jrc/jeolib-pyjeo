"""Test suite for module pyjeo.ngbops."""

import pyjeo as pj
import unittest

import numpy as np


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadNgbOps(unittest.TestCase):
    """Test functions and methods from ngbops modules."""

    def test_filters(self):
        """Test filter1d() and filter2d() functions and methods."""
        jim = pj.Jim(testFile)
        jim.geometry.cropBand(band=[0, 1, 2])

        stats_max_jim = max([stats['max'] for stats in [
            jim.stats.getStats('max', band=i) for i in [0, 1, 2]]])
        stats_min_jim = min([stats['min'] for stats in [
            jim.stats.getStats('min', band=i) for i in [0, 1, 2]]])

        min_jim = pj.ngbops.filter2d(jim, 'min', dx=3)
        jim.ngbops.filter2d('min', dx=3)

        stats_max_min_jim = max([stats['max'] for stats in [
            jim.stats.getStats('max', band=i) for i in [0, 1, 2]]])
        stats_min_min_jim = min([stats['min'] for stats in [
            jim.stats.getStats('min', band=i) for i in [0, 1, 2]]])

        assert min_jim.pixops.isEqual(jim), \
            'Inconsistency in ngbops.filter2d() ' \
            '(method returns different result than function)'
        assert stats_min_min_jim == stats_min_jim, \
            'Error in ngbops.filter2d() (wrong values)'
        assert stats_max_min_jim < stats_max_jim, \
            'Error in ngbops.filter2d() (wrong values)'

        max_jim = pj.ngbops.filter1d(jim, 'max', dz=3, pad='zeros',
                                     otype='Int16')
        jim.ngbops.filter1d('max', dz=3, pad='zeros', otype='Int16')

        stats_max_max_jim = max([stats['max'] for stats in [
            jim.stats.getStats('max', band=i) for i in [0, 1, 2]]])
        stats_min_max_jim = min([stats['min'] for stats in [
            jim.stats.getStats('min', band=i) for i in [0, 1, 2]]])

        assert max_jim.pixops.isEqual(jim), \
            'Inconsistency in ngbops.filter1d() ' \
            '(method returns different result than function)'
        assert stats_max_max_jim == stats_max_min_jim, \
            'Error in ngbops.filter1d() (wrong values)'
        assert stats_min_max_jim >= stats_min_min_jim, \
            'Error in ngbops.filter1d() (wrong values)'

        # TODO: Test numpy array as filter for filter2d()
        #       (must wait until the bug #16 in filter2d() in jiplib will be
        #       fixed)

    def test_erode_dilate(self):
        """Test morphoDilate... and morphoErode...() functions and methods."""
        jim = pj.Jim(nrow=500, ncol=500, otype='Byte', uniform=[0, 2], seed=0)

        arr = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
        jim[0:3, 0:3] = arr

        stats = jim.stats.getStats()

        # Test morphoDilateDiamond()
        dilated = pj.ngbops.morphoDilateDiamond(jim)
        jim.ngbops.morphoDilateDiamond()

        stats_dilated = jim.stats.getStats()

        assert jim.pixops.isEqual(dilated), \
            'Inconsistency in ngbops.morphoDilateDiamond() ' \
            '(method returns different result than function)'
        assert stats_dilated['max'] == 1, \
            'Error in ngbops.morphoDilateDiamond() (max value is not equal ' \
            'to 1, but is {})'.format(stats_dilated['max'])
        assert stats_dilated['min'] == 0, \
            'Error in ngbops.morphoDilateDiamond() ' \
            '(min value is not equal to 0)'
        assert stats_dilated['mean'] > stats['mean'], \
            'Error in ngbops.morphoDilateDiamond() ' \
            '(mean value is not greater than the one of the original Jim)'
        assert jim[0, 0] == 0, \
            'Error in ngbops.morphoDilateDiamond() ' \
            '(value whose 4-ngbs are of value 0 not equal to 0)'
        assert jim[0, 1] == jim[1, 0] == jim[1, 1] == jim[1, 2] == \
               jim[2, 1] == 1, \
            'Error in ngbops.morphoDilateDiamond() ' \
            '(values in 4-ngb of value 1 not equal to 1)'

        # Test morphoErodeDiamond()
        jim[0:3, 0:3] = arr
        jim[0:3, 1] = np.array([1, 1, 1])
        stats = stats_dilated
        eroded = pj.ngbops.morphoErodeDiamond(jim)
        jim.ngbops.morphoErodeDiamond()

        stats_eroded = jim.stats.getStats()

        assert jim.pixops.isEqual(eroded), \
            'Inconsistency in ngbops.morphoErodeDiamond() ' \
            '(method returns different result than function)'
        assert stats_eroded['max'] == 1, \
            'Error in ngbops.morphoErodeDiamond() (max value is not equal ' \
            'to 1, but is {})'.format(stats_eroded['max'])
        assert stats_eroded['min'] == 0, \
            'Error in ngbops.morphoErodeDiamond() ' \
            '(min value is not equal to 0)'
        assert stats_eroded['mean'] < stats['mean'], \
            'Error in ngbops.morphoErodeDiamond() ' \
            '(mean value is not smaller than the one of the original Jim)'
        assert jim[1, 0] == 1, \
            'Error in ngbops.morphoErodeDiamond() ' \
            '(value whose 4-ngbs are of value 1 not equal to 1)'
        assert jim[0, 0] == jim[0, 1] == jim[0, 2] == jim[1, 1] == \
               jim[1, 2] == jim[2, 0] == jim[2, 1] == jim[2, 2] == 0, \
            'Error in ngbops.morphoErodeDiamond() ' \
            '(values in 4-ngb of value 0 not equal to 0)'

        # Test morphoDilateLine()
        jim[0:3, 0:3] = arr
        stats = stats_eroded

        dilated = pj.ngbops.morphoDilateLine(jim, 1, 0, 3, 0)
        jim.ngbops.morphoDilateLine(1, 0, 3, 0)

        stats_dilated = jim.stats.getStats()

        assert jim.pixops.isEqual(dilated), \
            'Inconsistency in ngbops.morphoDilateLine() ' \
            '(method returns different result than function)'
        assert stats_dilated['max'] == 1, \
            'Error in ngbops.morphoDilateLine() ' \
            '(max value is not the same as of the original Jim)'
        assert stats_dilated['min'] == 0, \
            'Error in ngbops.morphoDilateLine() ' \
            '(min value is not equal to 0)'
        assert stats_dilated['mean'] > stats['mean'], \
            'Error in ngbops.morphoDilateLine() ' \
            '(mean value is not greater than the one of the original Jim)'
        assert jim[0, 0] == jim[0, 1] == jim[2, 0] == jim[2, 1] == 0, \
            'Error in ngbops.morphoDilateLine() ' \
            '(values surrounded in left-right ngb by value 0 not equal to 0)'
        assert jim[1, 0] == jim[1, 1] == jim[1, 2] == 1, \
            'Error in ngbops.morphoDilateLine() ' \
            '(values in left-right ngb of value 1 not equal to 1)'

        # Test morphoErodeLine()
        jim[0:3, 0:3] = arr
        jim[1, 0] = 1
        stats = stats_dilated

        dilated = pj.ngbops.morphoErodeLine(jim, 1, 0, 3, 0)
        jim.ngbops.morphoErodeLine(1, 0, 3, 0)

        stats_dilated = jim.stats.getStats()

        assert jim.pixops.isEqual(dilated), \
            'Inconsistency in ngbops.morphoErodeLine() ' \
            '(method returns different result than function)'
        assert stats_dilated['max'] == 1, \
            'Error in ngbops.morphoErodeLine() ' \
            '(max value is not the same as of the original Jim)'
        assert stats_dilated['min'] == 0, \
            'Error in ngbops.morphoErodeLine() ' \
            '(min value is not equal to 0)'
        assert stats_dilated['mean'] < stats['mean'], \
            'Error in ngbops.morphoErodeLine() ' \
            '(mean value is not lower than the one of the original Jim)'
        assert jim[0, 0] == jim[0, 1] == jim[0, 2] == jim[1, 1] == \
               jim[1, 2] == jim[2, 0] == jim[2, 1] == jim[2, 2] == 0, \
            'Error in ngbops.morphoErodeLine() ' \
            '(values in left-right ngb of value 0 not equal to 0)'
        assert jim[1, 0] == 1, \
            'Error in ngbops.morphoErodeLine() ' \
            '(value 1 not in left-right ngb of value 0 not equal to 1)'

        # TODO: Test morphoErode() and morphoDilate()
        #       (must wait until the issue #18 in jiplib will be
        #       fixed)

        # Test morphoGradientByDilationDiamond()
        jim_copy = pj.Jim(jim)

        dilated = pj.ngbops.morphoGradientByDilationDiamond(jim)
        jim.ngbops.morphoGradientByDilationDiamond()

        assert jim.pixops.isEqual(dilated), \
            'Inconsistency in ngbops.morphoGradientByDilationDiamond() ' \
            '(method returns different result than function)'
        assert jim.pixops.isEqual(
            pj.ngbops.morphoDilateDiamond(jim_copy)) - jim_copy, \
            'Error in ngbops.morphoGradientByDilationDiamond() ' \
            '(jim.ngbops.morphoGradientByDilationDiamond() not equal to ' \
            'pj.ngbops.morphoDilateDiamond(jim)) - jim'

        # Test morphoGradientByErosionDiamond()
        jim_copy = pj.Jim(jim)

        eroded = pj.ngbops.morphoGradientByErosionDiamond(jim)
        jim.ngbops.morphoGradientByErosionDiamond()

        assert jim.pixops.isEqual(eroded), \
            'Inconsistency in ngbops.morphoGradientByErosionDiamond() ' \
            '(method returns different result than function)'
        assert jim.pixops.isEqual(
            jim_copy - pj.ngbops.morphoErodeDiamond(jim_copy)), \
            'Error in ngbops.morphoGradientByErosionDiamond() ' \
            '(jim.ngbops.morphoGradientByErosionDiamond() not equal to ' \
            'jim - pj.ngbops.morphoErodeDiamond(jim))'

    def test_edgeWeight(self):
        """Test edgeWeight() function and method."""
        jim = pj.Jim(nrow=500, ncol=500, otype='Byte', uniform=[0, 2], seed=0)

        max_value = 2
        arr = np.array([[0, 1, 0], [0, max_value, 1], [0, 0, 0]])
        jim[0:3, 0:3] = arr

        # Test edgeWeight() with dir=0 and type=0
        hori_diff = pj.ngbops.edgeWeight(jim)
        jim.ngbops.edgeWeight()
        stats = jim.stats.getStats(['max', 'min'])

        assert jim.pixops.isEqual(hori_diff), \
            'Inconsistency in ngbops.edgeWeight() ' \
            '(method returns different result than function)'
        assert stats['max'] <= max_value, \
            'Error in ngbops.edgeWeight() ' \
            '(max value not <= max value before edgeWight())'
        assert stats['min'] >= 0, \
            'Error in ngbops.edgeWeight() ' \
            '(min value not >= min value before edgeWight())'
        assert jim[2, 0] == jim[2, 1] == 0, \
            'Error in ngbops.edgeWeight() ' \
            '(processing of [0, 0, 0] did not return [0, 0, x] for type=0)'
        assert jim[0, 0] == jim[0, 1] == 1, \
            'Error in ngbops.edgeWeight() ' \
            '(processing of [0, 1, 0] did not return [1, 1, x] for type=0)'
        assert jim[1, 0] == jim[1, 1] == 2, \
            'Error in ngbops.edgeWeight() ' \
            '(processing of [0, 2, 1] did not return [2, 1, x] for type=0)'

        # Test edgeWeight() with dir=1 and type=1
        jim[0:3, 0:3] = arr
        verti_max = pj.ngbops.edgeWeight(jim, 1, 1)
        jim.ngbops.edgeWeight(1, 1)
        stats = jim.stats.getStats(['max', 'min'])

        assert jim.pixops.isEqual(verti_max), \
            'Inconsistency in ngbops.edgeWeight(dir=1, type=1) ' \
            '(method returns different result than function)'
        assert stats['max'] <= max_value, \
            'Error in ngbops.edgeWeight(dir=1, type=1) ' \
            '(max value not <= max value before edgeWight())'
        assert stats['min'] >= 0, \
            'Error in ngbops.edgeWeight(dir=1, type=1) ' \
            '(min value not >= min value before edgeWight())'
        assert jim[0, 0] == jim[1, 0] == 0, \
            'Error in ngbops.edgeWeight(dir=1, type=1) ' \
            '(processing of [[0], [0], [0]] did not return [[0], [0], [x]])'
        assert jim[0, 2] == jim[1, 2] == 1, \
            'Error in ngbops.edgeWeight(dir=1, type=1) ' \
            '(processing of [[0], [1], [0]] did not return [[1], [1], [x]])'
        assert jim[0, 1] == jim[1, 1] == 2, \
            'Error in ngbops.edgeWeight(dir=1, type=1) ' \
            '(processing of [[1], [2], [0]] did not return [[2], [2], [x]])'

        # Test edgeWeight() with dir=0 and type=2
        jim[0:3, 0:3] = arr
        hori_min = pj.ngbops.edgeWeight(jim, type=2)
        jim.ngbops.edgeWeight(type=2)
        stats = jim.stats.getStats(['max', 'min'])

        assert jim.pixops.isEqual(hori_min), \
            'Inconsistency in ngbops.edgeWeight(type=2) ' \
            '(method returns different result than function)'
        assert stats['max'] <= max_value, \
            'Error in ngbops.edgeWeight(type=2) ' \
            '(max value not <= max value before edgeWight())'
        assert stats['min'] >= 0, \
            'Error in ngbops.edgeWeight(type=2) ' \
            '(min value not >= min value before edgeWight())'
        assert jim[0, 0] == jim[0, 1] == jim[1, 0] == jim[2, 0] == \
               jim[2, 1] == 0, \
            'Error in ngbops.edgeWeight(type=2) ' \
            '(processing of [0, 0, 1, 0, 2] did not return [0, 0, 0, 0, x]'
        assert jim[1, 1] == 1, \
            'Error in ngbops.edgeWeight(type=2) ' \
            '(processing of [2, 1] did not return [1, x]'

    def test_getDissim(self):
        """Test getDissim() function and method."""
        jim1 = pj.Jim(nrow=500, ncol=500, otype='Byte', uniform=[0, 3])
        jim2 = pj.Jim(nrow=500, ncol=500, otype='Byte', uniform=[0, 2])

        max_value = 3
        arr1 = np.array([[0, 1, 0], [0, 2, 1], [0, 0, 0]])
        jim1[0:3, 0:3] = arr1
        arr2 = np.array([[0, 0, 0], [0, max_value, 1], [0, 0, 0]])
        jim2[0:3, 0:3] = arr2

        # Test getDissim() with dissimType=0
        hori_diss1, vert_diss1 = pj.ngbops.getDissim(jim1)
        hori_diss2, vert_diss2 = jim1.ngbops.getDissim()
        stats = jim1.stats.getStats(['max', 'min'])

        assert hori_diss1 == hori_diss2 and vert_diss1 == vert_diss2, \
            'Inconsistency in ngbops.getDissim() ' \
            '(method returns different result than function)'
        assert stats['max'] <= 2, \
            'Error in ngbops.getDissim() ' \
            '(max value not <= max value before getDissim())'
        assert stats['min'] >= 0, \
            'Error in ngbops.getDissim() ' \
            '(min value not >= min value before getDissim())'
        assert hori_diss1[2, 0] == hori_diss1[2, 1] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=0)'
        assert hori_diss1[0, 0] == hori_diss1[0, 1] == hori_diss1[1, 1] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=0)'
        assert hori_diss1[1, 0] == 2, \
            'Error in ngbops.getDissim() ' \
            '(processing of [3, 1] and [2, 1] did not return ' \
            '[2, x] for dissimType=0)'
        assert vert_diss1[0, 0] == vert_diss1[1, 0] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=0)'
        assert vert_diss1[0, 1] == vert_diss1[0, 2] == vert_diss1[1, 2] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=0)'
        assert vert_diss1[1, 1] == 2, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 3, 0] and [1, 2, 0] did not return ' \
            '[3, 3, x] for dissimType=0)'

        # Test getDissim() with dissimType=0 and multiple Jims
        hori_diss1, vert_diss1 = pj.ngbops.getDissim([jim1, jim2])
        hori_diss2, vert_diss2 = jim1.ngbops.getDissim(jim2)
        hori_diss3, vert_diss3 = jim1.ngbops.getDissim([jim2])
        stats = jim1.stats.getStats(['max', 'min'])

        assert hori_diss1 == hori_diss2 and vert_diss1 == vert_diss2, \
            'Inconsistency in ngbops.getDissim() ' \
            '(method returns different result than function)'
        assert hori_diss3 == hori_diss2 and vert_diss3 == vert_diss2, \
            'Error in ngbops.getDissim(jimo=x) ' \
            '(when one Jim parsed into jimo as Jim returns different result ' \
            'than parsed as a list)'
        assert stats['max'] <= max_value, \
            'Error in ngbops.getDissim() ' \
            '(max value not <= max value before getDissim())'
        assert stats['min'] >= 0, \
            'Error in ngbops.getDissim() ' \
            '(min value not >= min value before getDissim())'
        assert hori_diss1[2, 0] == hori_diss1[2, 1] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=0)'
        assert hori_diss1[0, 0] == hori_diss1[0, 1] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=0)'
        assert hori_diss1[1, 1] == 2, \
            'Error in ngbops.getDissim() ' \
            '(processing of [3, 1] and [2, 1] did not return ' \
            '[2, x] for dissimType=0)'
        assert hori_diss1[1, 0] == 3, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 3] and [0, 2] did not return ' \
            '[3, x] for dissimType=0)'
        assert vert_diss1[0, 0] == vert_diss1[1, 0] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=0)'
        assert vert_diss1[0, 2] == vert_diss1[1, 2] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=0)'
        assert vert_diss1[0, 1] == vert_diss1[1, 1] == 3, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 3, 0] and [1, 2, 0] did not return ' \
            '[3, 3, x] for dissimType=0)'

        # Test getDissim() with dissimType=1
        hori_diss1, vert_diss1 = pj.ngbops.getDissim([jim1, jim2], 1)
        hori_diss2, vert_diss2 = jim1.ngbops.getDissim(jim2, 1)
        stats = jim1.stats.getStats(['max', 'min'])

        assert hori_diss1 == hori_diss2 and vert_diss1 == vert_diss2, \
            'Inconsistency in ngbops.getDissim() ' \
            '(method returns different result than function)'
        assert stats['max'] <= max_value, \
            'Error in ngbops.getDissim() ' \
            '(max value not <= max value before getDissim())'
        assert stats['min'] >= 0, \
            'Error in ngbops.getDissim() ' \
            '(min value not >= min value before getDissim())'
        assert hori_diss1[2, 0] == hori_diss1[2, 1] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=1)'
        assert hori_diss1[0, 0] == hori_diss1[0, 1] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=1)'
        assert hori_diss1[1, 1] == 2, \
            'Error in ngbops.getDissim() ' \
            '(processing of [3, 1] and [2, 1] did not return ' \
            '[2, x] for dissimType=1)'
        assert hori_diss1[1, 0] == 3, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 3] and [0, 2] did not return ' \
            '[3, x] for dissimType=1)'
        assert vert_diss1[0, 0] == vert_diss1[1, 0] == 0, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 0, 0] did not return ' \
            '[0, 0, x] for dissimType=1)'
        assert vert_diss1[0, 2] == vert_diss1[1, 2] == 1, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 0, 0] and [0, 1, 0] did not return ' \
            '[1, 1, x] for dissimType=1)'
        assert vert_diss1[0, 1] == vert_diss1[1, 1] == 3, \
            'Error in ngbops.getDissim() ' \
            '(processing of [0, 3, 0] and [1, 2, 0] did not return ' \
            '[3, 3, x] for dissimType=1)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadNgbOps)]
    return unittest.TestSuite(suite_list)
