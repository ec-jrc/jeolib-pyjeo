"""Test suite for module pyjeo.pixops."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadPixOps(unittest.TestCase):
    """Test functions and methods from pixops modules."""

    def test_NDVI(self):
        """Test computing NDVI in different ways."""
        jim = pj.Jim(testFile)

        jim_red = pj.geometry.cropBand(jim, 0)
        jim_nir = pj.geometry.cropBand(jim, 1)

        ndvi = pj.pixops.NDVI(jim, 0, 1)
        jim.pixops.NDVI(0, 1)

        assert jim.pixops.isEqual(ndvi), 'Error in computing NDVI'

        ndvi = pj.pixops.NDVISeparateBands(jim_red, jim_nir)
        jim_red.pixops.NDVISeparateBands(jim_nir)

        assert jim_red.pixops.isEqual(ndvi), 'Error in computing ' \
                                             'NDVISeparateBands'

        assert jim_red.pixops.isEqual(jim), 'Error in computing NDVI or ' \
                                            'NDVISeparateBands'

        ndvi = pj.pixops.NDVISeparateBands(jim_red, jim_nir.convertToFloat32())

        assert not jim.pixops.isEqual(ndvi), 'Error in computing NDVI'

    def test_supremum(self):
        """Test picking up supremum from different computed NDVIs."""
        for tile in tiles:
            jim4 = pj.Jim(tile)
            jim8 = pj.Jim(tile[:-8] + 'nir' + tile[-5] + '.tif')

            if 'max_ndvi' in locals():
                b = pj.pixops.NDVISeparateBands(jim4, jim8)

                max_ndvi[b > max_ndvi] = b
                max_ndvi_func = pj.pixops.supremum(max_ndvi, b)
                assert max_ndvi.pixops.isEqual(max_ndvi_func), \
                    'Error in computing supremum'
                max_ndvi2.pixops.supremum(b)
                assert max_ndvi.pixops.isEqual(max_ndvi2), \
                    'Error in computing supremum'
            else:
                max_ndvi = pj.pixops.NDVISeparateBands(jim4, jim8)
                max_ndvi2 = pj.Jim(max_ndvi)

    def test_setFunctions(self):
        """Test setData and setThreshold functions and methods."""
        jim = pj.Jim(tiles[0])

        stats = jim.stats.getStats()
        min = stats['min']
        max = stats['max']

        thresholded = pj.pixops.setThreshold(jim, min=min+1, max=max-1,
                                             nodata=(max+min)/2)

        thresholded_stats = thresholded.stats.getStats()
        nodata = thresholded.properties.getNoDataVals()

        assert thresholded_stats != stats, \
            'Error in pixops.setThreshold() or stats.getStats()'
        assert thresholded_stats['min'] >= min + 1, \
            'Error in pixops.setThreshold() or stats.getStats()'
        assert thresholded_stats['max'] <= max - 1, \
            'Error in pixops.setThreshold() or stats.getStats()'
        assert thresholded.properties.getNoDataVals() == [((max+min)/2)], \
            'Error in pixops.setThreshold() or properties.getNoDataVals()'

        jim.pixops.setThreshold(min=min+1, max=max-1, nodata=(max+min)/2)

        stats = jim.stats.getStats()

        assert stats == thresholded_stats, \
            'Error in pixops.setThreshold() or stats.getStats()'
        assert jim.properties.getNoDataVals() == nodata, \
            'Error in pixops.setThreshold() or properties.getNoDataVals()'

        jim.setData(5)

        stats = jim.stats.getStats()

        assert all(v==5 for v in [stats['min'], stats['max'], stats['mean']]),\
            'Error in pixops.setData() or stats.getStats()'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadPixOps)]
    return unittest.TestSuite(suite_list)
