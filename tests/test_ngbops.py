"""Test suite for module pyjeo.ngbops."""

import pyjeo as pj
import unittest


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
        #       (must wait until the bug in filter2d() in jiplib will be fixed)


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadNgbOps)]
    return unittest.TestSuite(suite_list)
