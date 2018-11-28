"""Test suite for module pyjeo.properties."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']


class BadProps(unittest.TestCase):
    """Test functions and methods for getting and setting properties."""

    jim = pj.Jim(tiles[0])
    ulx = jim.properties.getUlx()
    uly = jim.properties.getUly()
    lrx = jim.properties.getLrx()
    lry = jim.properties.getLry()

    def test_no_data_vals(self):
        """Test functions connected to no data values."""
        self.jim = pj.Jim(tiles[0])

        self.jim.properties.clearNoData()
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [], \
            'Error in properties.clearNoData() or getNoDataVals()'

        self.jim.properties.pushNoDataVal(5)
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [5], 'Error in properties.pushNoData()'

        self.jim.properties.pushNoDataVal(10)
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [5, 10], 'Error in properties.pushNoData()'

        self.jim.properties.setNoDataVals(22)
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [22], \
            'Error in properties.setNoDataVals() for one value'

        self.jim.properties.setNoDataVals([1, 2, 3])
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [1, 2, 3], \
            'Error in properties.setNoDataVals() for multiple values'

        self.jim.properties.clearNoData()
        no_data = self.jim.properties.getNoDataVals()

        assert no_data == [], 'Error in properties.setNoDataVals()'

    def test_covers(self):
        """Test covers() method."""
        assert self.jim.properties.covers(self.ulx + 0.1, self.uly - 0.1), \
            'Error in properties.covers(), getUlx() or getUly()'
        assert self.jim.properties.covers(self.lrx - 0.1, self.lry + 0.1), \
            'Error in properties.covers(), getUlx() or getUly()'
        assert not self.jim.properties.covers(self.lrx + 0.1, self.lry), \
            'Error in properties.covers(), getUlx() or getUly()'

    def test_bbox(self):
        """Test getBoundingBox() method."""
        assert self.jim.properties.getBoundingBox() == [self.ulx, self.uly,
                                                        self.lrx, self.lry], \
            'Error in properties.getBoundingBox()'

    def test_geo_transform(self):
        """Test GeoTransform() methods."""
        gt = self.jim.properties.getGeoTransform()

        assert gt[0] == self.ulx, 'Error in properties.getGeoTransform()'
        assert gt[3] == self.uly, 'Error in properties.getGeoTransform()'

        assert gt[1] == self.jim.properties.getDeltaX(), \
            'Error in properties.getDeltaX() or getGeoTransform()'
        assert abs(gt[5]) == self.jim.properties.getDeltaY(), \
            'Error in properties.getDeltaY() or getGeoTransform()'

        new_gt = [0, 1, 0, 0, 0, 1]

        self.jim.properties.setGeoTransform(new_gt)

        assert self.jim.properties.getGeoTransform() == new_gt, \
            'Error in properties.setGeoTransform()'

        jim = pj.Jim(tiles[0])
        if gt == new_gt:
            self.jim.properties.setGeoTransform([1, 2, 0, 1, 0, 2])

        self.jim.properties.copyGeoTransform(jim)

        assert self.jim.properties.getGeoTransform() == gt, \
            'Error in properties.copyGeoTransform()'

    def test_projection(self):
        """Test Projection() methods."""
        proj = self.jim.properties.getProjection()
        self.jim.properties.setProjection('EPSG:5514')

        assert self.jim.properties.getProjection() != proj, \
            'Error in properties.getProjection() or setProjection()'

    def test_data_types(self):
        """Test getDataType() and conversion methods."""
        assert pj.Jim(
            self.jim.convertToUchar8()).properties.getDataType() == 1, \
            'Error in properties.getDataType() or converting methods'

        assert pj.Jim(
            self.jim.convertToUint32()).properties.getDataType() == 5, \
            'Error in properties.getDataType() or converting methods'

        assert pj.Jim(
            self.jim.convertToFloat32()).properties.getDataType() == 6, \
            'Error in properties.getDataType() or converting methods'

        assert pj.Jim(
            self.jim.convertToDouble64()).properties.getDataType() == 7, \
            'Error in properties.getDataType() or converting methods'

        # TODO: Types 0, 2, 3, 5, 8, 9, 10, 11
        # TODO: Change the second assertion to 4 after fixing the bug in jiplib


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadProps)]
    return unittest.TestSuite(suite_list)
