"""Test suite for module pyjeo.geometry."""

import pyjeo as pj
import unittest

import random
import string
import os

tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']

rasterfn = 'tests/data/modis_ndvi_2010.tif'
vectorfn = 'tests/data/modis_ndvi_training.sqlite'
outputfn = os.path.join('/tmp',
                        ''.join(random.sample(string.ascii_letters, 5)))
warpedfn = os.path.join('/tmp',
                        ''.join(random.sample(string.ascii_letters, 5)))


class BadGeometry(unittest.TestCase):
    """Test functions and methods from geometry module."""

    def test_stack(self):
        """Test the stack band method."""
        jim0 = pj.Jim(rasterfn, band=0)
        jim1 = pj.Jim(rasterfn, band=1)
        jim2 = pj.Jim(rasterfn, band=2)
        jimlist = pj.JimList([jim0, jim1, jim2])
        jimliststack = pj.geometry.stackBand(jimlist)
        jimstack = pj.geometry.stackBand(pj.geometry.stackBand(jim0, jim1),
                                         jim2)
        assert jimliststack.pixops.isEqual(jimstack), \
            'Error in geometry.stackBand() ' \
            '(jimliststack not equal to jimstack)'
        jimliststack = pj.geometry.stackBand(jimlist, jim0)
        jimliststack.geometry.cropBand([0, 1, 2])
        assert jimliststack.pixops.isEqual(jimstack), \
            'Error in geometry.stackBand() ' \
            '(jimliststack not equal to jimstack after crop)'
        jim3 = pj.JimList([jim0, jim1, jim2]).geometry.stackBand()
        jim3.geometry.cropBand([0, 1])
        assert pj.geometry.cropBand(jim3, 0).pixops.isEqual(jim0), \
            'Error jim3 not equal to jim0'
        jim4 = pj.JimList([jim0, jim1, jim2]).geometry.stackBand()
        jim4.geometry.cropBand([1, 2])
        assert pj.geometry.cropBand(jim4, 0).pixops.isEqual(jim1), \
            'Error jim3 not equal to jim1'
        jim = pj.JimList([jim3, jim4]).geometry.stackBand([0, 1])
        jim.geometry.cropBand([1, 2])
        jim.geometry.cropBand(0)
        assert jim.pixops.isEqual(jim1), \
            'Error jim not equal to jim1'

        # Test the stack plane method
        jim0 = pj.Jim(rasterfn, band=0)
        jim1 = pj.Jim(rasterfn, band=1)
        jim2 = pj.Jim(rasterfn, band=2)
        jimlist = pj.JimList([jim0, jim1, jim2])
        jimliststack = pj.geometry.stackPlane(jimlist)
        jimstack = pj.geometry.stackPlane(pj.geometry.stackPlane(jim0, jim1),
                                          jim2)
        assert jimliststack.pixops.isEqual(jimstack), \
            'Error in geometry.stackPlane() ' \
            '(jimliststack not equal to jimstack)'
        jim3 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim3.geometry.cropPlane([0, 1])
        assert pj.geometry.cropPlane(jim3, 0).pixops.isEqual(jim0), \
            'Error jim3 not equal to jim0'
        jim4 = pj.JimList([jim0, jim1, jim2]).geometry.stackPlane()
        jim4.geometry.cropPlane([1, 2])
        assert pj.geometry.cropPlane(jim4, 0).pixops.isEqual(jim1), \
            'Error jim3 not equal to jim1'
        jim = pj.JimList([jim3, jim4]).geometry.stackPlane()
        jim.geometry.cropPlane([1, 2])
        jim.geometry.cropPlane(0)
        assert jim.pixops.isEqual(jim1), \
            'Error jim not equal to jim1'

        # Test the reduce plane method
        jimstack = pj.Jim(rasterfn, band=[0, 1, 2], band2plane=True)
        jimreduce = pj.geometry.reducePlane(jimstack, rule='max')
        assert jimreduce.properties.nrOfPlane() == 1
        jimstack.geometry.reducePlane(rule='max')
        assert jimstack.properties.nrOfPlane() == 1
        assert jimreduce.pixops.isEqual(jimstack), \
            'Error in geometry.reducePlane()'

    def test_warp(self):
        """Test the warp method."""
        jim0 = pj.Jim(rasterfn)
        jim_warped = pj.geometry.warp(jim0, 'epsg:4326')

        assert jim_warped.properties.getProjection()[-7:-3] == '4326', \
            'Error in geometry.warp(): EPSG not changed'

        jim_warped.geometry.warp('epsg:3035')

        jim_warped.geometry.crop(ulx=jim0.properties.getUlx(),
                                 uly=jim0.properties.getUly(),
                                 lrx=jim0.properties.getLrx(),
                                 lry=jim0.properties.getLry(),
                                 dx=jim0.properties.getDeltaX(),
                                 dy=jim0.properties.getDeltaY())
        jim0.geometry.crop(ulx=jim0.properties.getUlx(),
                           uly=jim0.properties.getUly(),
                           lrx=jim0.properties.getLrx(),
                           lry=jim0.properties.getLry(),
                           dx=jim0.properties.getDeltaX(),
                           dy=jim0.properties.getDeltaY())
        jim_warped.io.write(warpedfn)
        assert jim0.properties.getBBox() == jim_warped.properties.getBBox(), \
            'Error in geometry.warp(): BBox'
        assert jim0.properties.nrOfCol() == jim_warped.properties.nrOfCol(), \
            'Error in geometry.warp(): nrOfCol'
        assert jim0.properties.nrOfRow() == jim_warped.properties.nrOfRow(), \
            'Error in geometry.warp(): nrOfRow'
        assert jim0.properties.nrOfBand() == jim_warped.properties.nrOfBand(),\
            'Error in geometry.warp(): nrOfBand'

    def test_extractOgr(self):
        """Test the extractOgr method."""
        jim0 = pj.Jim(rasterfn)
        sample = pj.JimVect(vectorfn)
        for band in range(0, 12):
            jl0 = pj.JimList([pj.geometry.cropBand(jim0, band)])
            bandname = 'B' + str(band)
            if not band:
                v = jl0.geometry.extractOgr(sample, rule='mean',
                                            output=outputfn, oformat='SQLite',
                                            co=['OVERWRITE=YES'],
                                            bandname=bandname, fid='fid')
                v.io.write()
                assert v.properties.getFeatureCount() == 11, \
                    'Error in geometry.extractOgr() feature count (1)'
                assert 'fid' in v.properties.getFieldNames(), \
                    'Error in geometry.extractOgr() field names (1)'
                v.io.close()
            else:
                v1 = pj.JimVect(outputfn)
                v2 = jl0.geometry.extractOgr(
                    sample, rule='mean', output='/vsimem/v2.sqlite',
                    oformat='SQLite', co=['OVERWRITE=YES'],
                    bandname=bandname, fid='fid')
                v = pj.geometry.join(v1, v2, output=outputfn, oformat='SQLite',
                                     co=['OVERWRITE=YES'], key=['fid'])
                assert v.properties.getFeatureCount() == 11, \
                    'Error in geometry.extractOgr() feature count (2)'
                assert 'fid' in v.properties.getFieldNames(), \
                    'Error in geometry.extractOgr() field names (1)'
                if band == 11:
                    assert len(v.properties.getFieldNames()) == 14, \
                        'Error in geometry.extractOgr() number of field names'
                v1.io.close()
                v2.io.close()
                v.io.write()
                v.io.close()
            jl0.io.close()
        sample.io.close()

    def test_crop(self):
        """Test crop...() functions and methods."""
        raster = pj.geometry.stackPlane(pj.Jim(rasterfn), pj.Jim(rasterfn))
        vector = pj.JimVect(vectorfn)
        raster_bbox = raster.properties.getBBox()
        raster_dx = raster.properties.getDeltaX()
        raster_dy = raster.properties.getDeltaY()
        vector_bbox = vector.properties.getBBox()

        # Test cropOgr()
        cropped = pj.geometry.cropOgr(raster, vector)
        raster.geometry.cropOgr(vector)
        raster_bbox_cropped = raster.properties.getBBox()

        mod_x = (raster_bbox_cropped[2] - raster_bbox_cropped[0]) % raster_dx
        mod_y = (raster_bbox_cropped[3] - raster_bbox_cropped[1]) % raster_dy

        assert raster.pixops.isEqual(cropped), \
            'Inconsistency in geometry.cropOgr() ' \
            '(method returns different result than function)'
        assert raster_bbox != raster_bbox_cropped, \
            'Error in geometry.cropOgr() ' \
            '(BBox not changed after crop)'
        assert raster_bbox_cropped[:2] == vector_bbox[:2] and mod_x == 0 and \
               mod_y == 0, \
            'Error in geometry.cropOgr() ' \
            '(new BBox values not equal to the vector one)'

    def test_coords_transformations(self):
        """Test geo2image() and image2geo() functions and methods."""
        jim = pj.Jim(tiles[0])
        geo_ulx, geo_uly, geo_lrx, geo_lry = jim.properties.getBBox()
        delta_x = jim.properties.getDeltaX()
        delta_y = jim.properties.getDeltaY()
        ulx_pix_center = geo_ulx + delta_x / 2
        uly_pix_center = geo_uly - delta_y / 2

        im_x, im_y = pj.geometry.geo2image(jim, geo_ulx, geo_uly)
        im_x_m, im_y_m = jim.geometry.geo2image(geo_ulx, geo_uly)

        assert im_x == im_x_m and im_y == im_y_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert im_x == 0 and im_y == 0, \
            'Error in geometry.geo2image()' \
            '(geo2image(BBox[0], BBox[1]) did not return zeros)'

        im_x, im_y = pj.geometry.geo2image(jim, 0, 0)
        im_x_m, im_y_m = jim.geometry.geo2image(0, 0)

        assert im_x == im_x_m and im_y == im_y_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert im_x == -geo_ulx / delta_x and im_y == geo_uly / delta_y, \
            'Error in geometry.geo2image()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'

        geo_ulx_2, geo_uly_2 = pj.geometry.image2geo(jim, 0, 0)
        geo_ulx_2_m, geo_uly_2_m = jim.geometry.image2geo(0, 0)

        assert geo_ulx_2 == geo_ulx_2_m and geo_uly_2 == geo_uly_2_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert geo_ulx_2 == ulx_pix_center and geo_uly_2 == uly_pix_center, \
            'Error in geometry.image2geo()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'

        geo_ulx_2, geo_uly_2 = pj.geometry.image2geo(jim, im_x, im_y)
        geo_ulx_2_m, geo_uly_2_m = jim.geometry.image2geo(im_x, im_y)

        assert geo_ulx_2 == geo_ulx_2_m and geo_uly_2 == geo_uly_2_m, \
            'Inconsistency in geometry.geo2image() ' \
            '(method returns different result than function)'
        assert geo_ulx_2 == 5 and geo_uly_2 == -5, \
            'Error in geometry.image2geo()' \
            '(geo2image(0, 0) did not return original values divided by dX/dY)'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadGeometry)]
    return unittest.TestSuite(suite_list)
