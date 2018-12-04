"""Test suite for module pyjeo.properties."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


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

        try:
            self.jim.properties.setNoDataVals('string')
            failed = True
        except TypeError:
            failed = False

        assert not failed, 'Error in type checks in properties.setNoDataVals()'

    def test_covers(self):
        """Test covers() method."""
        assert self.jim.properties.covers(self.ulx + 0.1, self.uly - 0.1), \
            'Error in properties.covers(), getUlx() or getUly()'
        assert self.jim.properties.covers(self.lrx - 0.1, self.lry + 0.1), \
            'Error in properties.covers(), getUlx() or getUly()'
        assert not self.jim.properties.covers(self.lrx + 0.1, self.lry), \
            'Error in properties.covers(), getUlx() or getUly()'

    def test_bbox(self):
        """Test getBBox() method."""
        assert self.jim.properties.getBBox() == [self.ulx, self.uly,
                                                        self.lrx, self.lry], \
            'Error in properties.getBBox()'

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

        assert self.jim.properties.nrOfPlane() == \
               self.jim.properties.nrOfBand(), \
            'Error in properties.nrOfPlains() or properties.nrOfBands()'

        center_pos = self.jim.properties.getCenterPos()

        assert center_pos[0] > self.ulx and center_pos[0] < self.lrx, \
            'Error in properties.getCenterPos()'
        assert center_pos[1] < self.uly and center_pos[1] > self.lry, \
            'Error in properties.getCenterPos()'

    def test_projection(self):
        """Test Projection() methods."""
        proj = self.jim.properties.getProjection()
        self.jim.properties.setProjection('EPSG:5514')

        assert self.jim.properties.getProjection() != proj, \
            'Error in properties.getProjection() or setProjection()'


class BadPropsLists(unittest.TestCase):
    """Test JimList funcs and methods for getting and setting properties."""

    jim1 = pj.Jim(testFile)
    jim2 = pj.Jim(tiles[1])

    jiml = pj.JimList([jim1, jim2])
    ulx = jiml.properties.getUlx()
    uly = jiml.properties.getUly()
    lrx = jiml.properties.getLrx()
    lry = jiml.properties.getLry()

    def test_no_data_vals(self):
        """Test JimList methods connected to no data values."""
        self.jiml.properties.clearNoData()

        self.jiml.properties.pushNoDataVal(0)

        assert self.jiml.properties.getNoDataVals() == [0], \
            'Error in properties.pushNoDataVal() or getNoDataVals for JimList'

        self.jiml.properties.pushNoDataVal(1)

        assert self.jiml.properties.getNoDataVals() == [0, 1], \
            'Error in properties.pushNoDataVal() for JimList'

        self.jiml.properties.clearNoData()

        assert self.jiml.properties.getNoDataVals() == [], \
            'Error in properties.clearNoData() for JimList'

    def test_covers(self):
        """Test JimList properties.covers() method."""
        assert self.jiml.properties.covers(self.ulx + 0.1, self.uly - 0.1), \
            'Error in properties.covers(), getUlx() or getUly()'
        assert self.jiml.properties.covers(self.lrx - 0.1, self.lry + 0.1), \
            'Error in properties.covers(), getLrx() or getLry()'
        assert not self.jiml.properties.covers(self.lrx + 0.1, self.lry), \
            'Error in properties.covers(), getLrx() or getLry()'

    def test_bbox(self):
        """Test JimList properties.getBBox() method."""
        assert self.jiml.properties.getBBox() == [self.ulx, self.uly,
                                                  self.lrx, self.lry], \
            'Error in properties.getBBox()'

    def test_selectGeo(self):
        """Test JimList properties.selectGeo() function and method."""
        self.jiml.properties.selectGeo(self.ulx + 0.1, self.uly - 0.1)
        assert len(self.jiml) == 1, 'Error in properties.selectGeo()'

        self.jiml.properties.selectGeo(0, 0)
        assert len(self.jiml) == 0, 'Error in properties.selectGeo()'


class BadPropsVects(unittest.TestCase):
    """Test JimVect funcs and methods for getting and setting properties."""

    jimv = pj.JimVect(vector)

    def test_geospatial_infos(self):
        """Test JimVect methods connected to geospatial informations."""
        bbox = self.jimv.properties.getBBox()

        assert self.jimv.properties.getUlx() == bbox[0], \
            'Error in properties.getBBox() or properties.getUlx() for JimVects'
        assert self.jimv.properties.getUly() == bbox[1], \
            'Error in properties.getBBox() or properties.getUly() for JimVects'
        assert self.jimv.properties.getLrx() == bbox[2], \
            'Error in properties.getBBox() or properties.getLrx() for JimVects'
        assert self.jimv.properties.getLry() == bbox[3], \
            'Error in properties.getBBox() or properties.getLry() for JimVects'

    def test_feature_layer_counts(self):
        """Test JimVect methods getLayerCount() and getFeatureType()"""
        assert self.jimv.properties.getLayerCount() == \
               self.jimv.properties.getFeatureCount() == 2, \
            'Error in properties.getLayerCount or getFeatureCount for JimVects'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadProps),
                  loader.loadTestsFromTestCase(BadPropsLists),
                  loader.loadTestsFromTestCase(BadPropsVects)]
    return unittest.TestSuite(suite_list)
