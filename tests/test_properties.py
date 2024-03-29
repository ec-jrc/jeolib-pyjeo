"""Test suite for module pyjeo.properties."""
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
import numpy as np
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

    def test_isEqual(self):
        """Test isEqual() method."""
        jim2 = pj.Jim(tiles[0])
        jim3 = pj.Jim(tiles[1])

        assert self.jim.properties.isEqual(jim2), \
            'Error in properties.isEqual() method (for two Jims)'
        assert not self.jim.properties.isEqual(jim3), \
            'Error in properties.isEqual() method (non-equality)'
        assert not self.jim.properties.isEqual(1), \
            'Error in properties.isEqual() method (for not Jim object)'

        assert pj.properties.isEqual(self.jim, jim2), \
            'Error in properties.isEqual() function (for two Jims)'
        assert not pj.properties.isEqual(self.jim, jim3), \
            'Error in properties.isEqual() function (non-equality)'
        assert not pj.properties.isEqual(self.jim, 1), \
            'Error in properties.isEqual() function (for not Jim object)'

        # Test multi-plane objects

        jim_oneplane = pj.Jim(ncol=5, nrow=5, nband=2, otype='GDT_Byte')
        jim_twoplanes = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                               otype='GDT_Byte')

        jim_oneplane.pixops.setData(5)
        jim_twoplanes.pixops.setData(5)

        assert not jim_oneplane.properties.isEqual(jim_twoplanes), \
            'Error in Jim.properties.isEqual() ' \
            '(not False for different number of planes)'

        assert not pj.properties.isEqual(jim_oneplane, jim_twoplanes), \
            'Error in pj.properties.isEqual() ' \
            '(not False for different number of planes)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(5)

        assert jim_twoplanes.properties.isEqual(jim_twoplanes2), \
            'Error in pj.properties.isEqual() ' \
            '(wrong result for multiplane Jim object)'

        assert pj.properties.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.properties.isEqual() ' \
            '(wrong result for multiplane Jim object)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=3, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(5)

        assert not jim_twoplanes.properties.isEqual(jim_twoplanes2), \
            'Error in Jim.properties.isEqual() ' \
            '(not False for different number of bands for multiplane Jim)'

        assert not pj.properties.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.properties.isEqual() ' \
            '(not False for different number of bands for multiplane Jim)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(6)

        assert not jim_twoplanes.properties.isEqual(jim_twoplanes2), \
            'Error in pj.properties.isEqual() ' \
            '(wrong result for multiplane multiband Jim object)'

        assert not pj.properties.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.properties.isEqual() ' \
            '(wrong result for multiplane multiband Jim object)'

        # Test different band number recognition

        jim_oneplane_oneband = pj.geometry.cropBand(jim_oneplane, 0)

        assert not jim_oneplane_oneband.properties.isEqual(jim_oneplane), \
            'Error in properties.isEqual() ' \
            '(did not recognize different number of bands)'

    def test_no_data_vals(self):
        """Test functions connected to no data values."""
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
            raised = False
        except pj.exceptions.JimIllegalArgumentError:
            raised = True

        assert raised, 'Error in type checks in properties.setNoDataVals()'

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
        bbox = self.jim.properties.getBBox()
        assert bbox == [self.ulx, self.uly, self.lrx, self.lry], \
            'Error in properties.getBBox()'

        """Test Jim properties.getBBox() method at target projection."""
        assert bbox == self.jim.properties.getBBox(
            t_srs=self.jim.properties.getProjection()), \
            'Error in properties getBBox() with projection '

        """Test Jim properties.getBBox() method at target projection."""
        bbox = self.jim.properties.getBBox(t_srs = 4326)
        assert bbox[0] > 7.70 and bbox[0] < 7.71, \
            "Error in properties.getBBox() with projection [0]: {}".format(bbox[0])
        assert bbox[1] > 46.04 and bbox[1] < 46.05, \
            "Error in properties.getBBox() with projection [1]: {}".format(bbox[1])
        assert bbox[2] > 7.77 and bbox[2] < 7.78, \
            "Error in properties.getBBox() with projection [2]: {}".format(bbox[2])
        assert bbox[3] > 46.00 and bbox[3] < 46.01, \
            "Error in properties.getBBox() with projection [3]: {}".format(bbox[3])

    def test_refpix(self):
        """Test getRefPix() method."""
        assert self.jim.properties.getRefPix() == [402520.0,
                                                   5097440.0],\
            'Error in properties.getRefPix()'

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

        assert self.ulx < center_pos[0] < self.lrx, \
            'Error in properties.getCenterPos()'
        assert self.uly > center_pos[1] > self.lry, \
            'Error in properties.getCenterPos()'

    def test_projection(self):
        """Test Projection() methods."""
        proj = self.jim.properties.getProjection()
        self.jim.properties.setProjection('EPSG:5514')

        assert self.jim.properties.getProjection() != proj, \
            'Error in properties.getProjection() or setProjection()'

    @staticmethod
    def test_georeference():
        """Test georeference (projection + geotransform) methods."""
        jim1 = pj.Jim(testFile, band=[0, 1])
        np1 = jim1.np()[:]
        jim_1 = pj.np2jim(np1)
        jim_1.properties.copyGeoReference(jim1)

        assert jim_1.properties.getProjection() == \
               jim1.properties.getProjection(),\
            'Error in properties.copyGeoReference getProjection()'

        assert jim_1.properties.getGeoTransform() == \
               jim1.properties.getGeoTransform(),\
            'Error in properties.copyGeoReference getGeoTransform()'


class BadPropsLists(unittest.TestCase):
    """Test JimList funcs and methods for getting and setting properties."""

    jim1 = pj.Jim(testFile, band=[0, 1])
    jim2 = pj.Jim(tiles[1])
    jim3 = pj.Jim(tiles[0])

    jiml = pj.JimList([jim1, jim2])
    ulx = jiml.properties.getUlx()
    uly = jiml.properties.getUly()
    lrx = jiml.properties.getLrx()
    lry = jiml.properties.getLry()

    jiml2 = pj.JimList([jim2, jim3])
    ulx2 = jiml2.properties.getUlx()
    uly2 = jiml2.properties.getUly()
    lrx2 = jiml2.properties.getLrx()
    lry2 = jiml2.properties.getLry()

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
        bbox = self.jiml2.properties.getBBox()
        assert bbox == [self.ulx2, self.uly2, self.lrx2, self.lry2], \
            'Error in properties.getBBox()'

        """Test JimList properties.getBBox() method at target projection."""
        assert bbox == self.jiml2.properties.getBBox(
            t_srs=self.jiml2[0].properties.getProjection()),\
            'Error in properties getBBox() with projection '

        """Test JimList properties.getBBox() method at target projection."""
        bbox = self.jiml2.properties.getBBox(t_srs = 4326)
        assert bbox[0] > 7.70 and bbox[0] < 7.71, \
            "Error in properties.getBBox() with projection [0]: {}".format(bbox[0])
        assert bbox[1] > 46.04 and bbox[1] < 46.05, \
            "Error in properties.getBBox() with projection [1]: {}".format(bbox[1])
        assert bbox[2] > 7.77 and bbox[2] < 7.78, \
            "Error in properties.getBBox() with projection [2]: {}".format(bbox[2])
        assert bbox[3] > 46.00 and bbox[3] < 46.01, \
            "Error in properties.getBBox() with projection [3]: {}".format(bbox[3])

    def test_selectGeo(self):
        """Test JimList properties.selectGeo() function and method."""
        self.jiml.properties.selectGeo(self.ulx + 0.1, self.uly - 0.1)
        assert len(self.jiml) == 1, 'Error in properties.selectGeo()'

        self.jiml.properties.selectGeo(0, 0)
        assert len(self.jiml) == 0, 'Error in properties.selectGeo()'


class BadPropsVects(unittest.TestCase):
    """Test JimVect funcs and methods for getting and setting properties."""

    jimv = pj.JimVect(vector)

    def test_isEqual(self):
        """Test JimVect isEqual method."""
        assert self.jimv.properties.isEqual(self.jimv) is True, \
            'Error in properties.isEqual(): not equal to itself (2 layers)'
        assert pj.properties.isEqual(self.jimv, self.jimv) is True, \
            'Error in function properties.isEqual(): not equal to itself'

        jimv1 = pj.JimVect(vector, ln='milano')
        jimv2 = pj.JimVect(vector, ln='lodi')

        assert self.jimv.properties.isEqual(jimv1) is False, \
            'Error in properties.isEqual(): \
            2 layer vector should be different from 1 layer vector'
        assert pj.properties.isEqual(self.jimv, jimv1) is False, \
            'Error in function properties.isEqual(): \
            2 layer vector should be different from 1 layer vector'
        assert jimv1.properties.isEqual(jimv1) is True, \
            'Error in properties.isEqual(): not equal to itself (1 layer)'
        assert pj.properties.isEqual(jimv1, jimv1) is True, \
            'Error in function properties.isEqual(): \
            not equal to itself (1 layer)'
        assert jimv1.properties.isEqual(jimv2) is False, \
            'Error in properties.isEqual(): \
            layer milano should be different from lodi'
        assert pj.properties.isEqual(jimv1, jimv2) is False, \
            'Error in function properties.isEqual(): \
            layer milano should be different from lodi'

        jimv1_double = pj.geometry.append(jimv1, jimv1,
                                          output=pj._get_random_path())

        assert not jimv1_double.properties.isEqual(jimv1), \
            'Error in properties.isEqual() method ' \
            '(JimVects with different number of features considered equal)'
        assert not pj.properties.isEqual(jimv1_double, jimv1), \
            'Error in properties.isEqual() function ' \
            '(JimVects with different number of features considered equal)'

        # Test with the argument Jim empty

        jimv_empty = pj.JimVect()

        equal = jimv1.properties.isEqual(jimv_empty)

        assert equal is False, \
            'Error in properties.isEqual(JimVect) method ' \
            '(did not return False when the argument JimVect was empty)'

        equal = pj.properties.isEqual(jimv1, jimv_empty)

        assert equal is False, \
            'Error in properties.isEqual(JimVect, JimVect) function ' \
            '(did not return False when the argument JimVect was empty)'

        equal = jimv_empty.properties.isEqual(jimv1)
        assert equal is False, \
            'Error in properties.isEqual(JimVect, JimVect) method ' \
            '(did not return False when the argument JimVect was empty)'

        equal = pj.properties.isEqual(jimv_empty, jimv1)
        assert equal is False, \
            'Error in properties.isEqual(JimVect, JimVect) function ' \
            '(did not return False when the argument JimVect was empty)'

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

        """Test JimVect properties.getBBox() method at target projection."""
        assert bbox == self.jimv.properties.getBBox(
            t_srs=self.jimv.properties.getProjection()), \
            'Error in properties getBBox() with projection '

        """Test JimVect properties.getBBox() method at target projection."""
        bbox = self.jimv.properties.getBBox(t_srs = 4326)
        assert bbox[0] > 8.70 and bbox[0] < 8.71, \
            "Error in properties.getBBox() with projection [0]: {}".format(bbox[0])
        assert bbox[1] > 45.73 and bbox[1] < 45.74, \
            "Error in properties.getBBox() with projection [1]: {}".format(bbox[1])
        assert bbox[2] > 9.89 and bbox[2] < 9.90, \
            "Error in properties.getBBox() with projection [2]: {}".format(bbox[2])
        assert bbox[3] > 45.05 and bbox[3] < 45.06, \
            "Error in properties.getBBox() with projection [3]: {}".format(bbox[3])

        # assert bbox[0] > 7.70 and bbox[0] < 7.71, \
        # assert bbox[1] > 46.04 and bbox[1] < 46.05, \
        # assert bbox[2] > 7.77 and bbox[2] < 7.78, \
        # assert bbox[3] > 46.00 and bbox[3] < 46.01, \
        projection = self.jimv.properties.getProjection()

        assert (projection.split('GEOGCS["')[1][:6] == 'ETRS89' or
                projection.split('GEOGCS["')[1][:13] == 'GCS_ETRS_1989'), \
                'Error in properties.getProjection() for JimVects (wrong GEOGCS)'

    def test_feature_layer_counts(self):
        """Test JimVect methods getLayerCount() and getFeatureType()."""
        assert self.jimv.properties.getLayerCount() == 2, \
            'Error in properties.getLayerCount() for JimVects'

        assert self.jimv.properties.getFeatureCount() == \
               self.jimv.properties.getFeatureCount(0) + \
               self.jimv.properties.getFeatureCount(1), \
            'Error in properties.getFeatureCount() for JimVects'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadProps),
                  loader.loadTestsFromTestCase(BadPropsLists),
                  loader.loadTestsFromTestCase(BadPropsVects)]
    return unittest.TestSuite(suite_list)
