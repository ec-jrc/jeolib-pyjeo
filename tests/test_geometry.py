"""Test suite for module pyjeo.geometry."""

import pyjeo as pj
import unittest

tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']

rasterfn = 'tests/data/modis_ndvi_2010.tif'
vectorfn = 'tests/data/modis_ndvi_training.sqlite'
outputfn = pj._get_random_path()
warpedfn = pj._get_random_path()


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

        # Test the band2plane method
        jimband = pj.Jim(rasterfn, band=[0,1,2])
        jimplane = pj.Jim(rasterfn, band=[0,1,2], band2plane=True)
        jimband.geometry.band2plane()
        assert jimband.pixops.isEqual(jimplane), \
            'Error in geometry.band2plane() ' \
            '(jimband not equal to jimsplane)'

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

    def test_image_frames(self):
        """Test imageFrame...() functions and methods."""
        nrow = ncol = 500
        nband = nplane = 2
        jim = pj.Jim(nrow=nrow, ncol=ncol, nband=nband, nplane=nplane,
                     otype='Byte', uniform=[0, 2], seed=0)

        # Test imageFrameAdd()
        #      (for 1-band Jim, see test below)

        added = pj.geometry.imageFrameAdd(jim, 1, 2, 1, 2, 1, 2, 5)
        jim.geometry.imageFrameAdd(1, 2, 1, 2, 1, 2, 5)

        assert jim.pixops.isEqual(added), \
            'Inconsistency in geometry.imageFrameAdd() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of bands changed)'
        assert jim[0, 100, 100].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for new planes)'
        assert jim[1, 0, 0].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for the frame)'
        assert 0 <= jim[1, 1, 1].np()[0] <= 1, \
            'Error in geometry.imageFrameAdd() ' \
            '(values in the original image changed)'

        # Test imageFrameSubtract()
        #      (for 1-band Jim, see test below)

        subtracted = pj.geometry.imageFrameSubtract(jim, 1, 2, 1, 2, 1, 2)
        jim.geometry.imageFrameSubtract(1, 2, 1, 2, 1, 2)

        assert jim.pixops.isEqual(subtracted), \
            'Inconsistency in geometry.imageFrameSubtract() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of bands changed)'
        assert (jim.np() == added.np()[1:3, 1:-2, 1:-2]).all(), \
            'Error in geometry.imageFrameSubtract() ' \
            '(changed values in the original Jim)'

        # Test imageFrameSet()

        setted = pj.geometry.imageFrameSet(jim, 1, 2, 1, 2, 1, 0, 5)
        jim.geometry.imageFrameSet(1, 2, 1, 2, 1, 0, 5)

        assert jim.pixops.isEqual(setted), \
            'Inconsistency in geometry.imageFrameSet() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSet() ' \
            '(number of cols changed)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSet() ' \
            '(number of rows changed)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSet() ' \
            '(number of planes changed)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSet() ' \
            '(number of bands changed)'
        assert jim[1, 0, 0].np() == 5, \
            'Error in geometry.imageFrameSet() ' \
            '(value not used for the frame)'
        assert 0 <= jim[1, 1, 1].np()[0] <= 1,\
            'Error in geometry.imageFrameSet() ' \
            '(values outside the frame changed)'

        # Test imageFrameSet() with a specified band

        setted = pj.geometry.imageFrameSet(jim, 1, 2, 1, 2, 0, 1, 10, band=0)
        jim.geometry.imageFrameSet(1, 2, 1, 2, 0, 1, 10, band=0)

        assert jim.pixops.isEqual(setted), \
            'Inconsistency in geometry.imageFrameSet() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSet() ' \
            '(number of cols changed)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSet() ' \
            '(number of rows changed)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSet() ' \
            '(number of planes changed)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSet() ' \
            '(number of bands changed)'
        assert jim[1, 0, 0].np() == 10, \
            'Error in geometry.imageFrameSet() ' \
            '(value not used for the frame)'
        assert jim[0, 1, 1].np() == 5,\
            'Error in geometry.imageFrameSet() ' \
            '(values outside the frame changed)'

        # Test imageFrameAdd() for 1-band Jim

        jim.geometry.cropBand(0)
        nband = 1

        added = pj.geometry.imageFrameAdd(jim, 1, 2, 1, 2, 1, 2, 10)
        jim.geometry.imageFrameAdd(1, 2, 1, 2, 1, 2, 10)

        assert jim.pixops.isEqual(added), \
            'Inconsistency in geometry.imageFrameAdd() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane + 3, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameAdd() ' \
            '(number of bands changed)'
        assert jim[0, 100, 100].np() == 10, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for new planes)'
        assert jim[1, 0, 0].np() == 10, \
            'Error in geometry.imageFrameAdd() ' \
            '(value not used for the frame)'
        assert jim[1, 2, 2].np() == 5, \
            'Error in geometry.imageFrameAdd() ' \
            '(values in the original image changed)'

        # Test imageFrameSubtract() for 1-band Jim

        subtracted = pj.geometry.imageFrameSubtract(jim, 1, 2, 1, 2, 1, 2)
        jim.geometry.imageFrameSubtract(1, 2, 1, 2, 1, 2)

        assert jim.pixops.isEqual(subtracted), \
            'Inconsistency in geometry.imageFrameSubtract() ' \
            '(method returns different result than function)'
        assert jim.properties.nrOfCol() == ncol, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of cols not raised or not raised to the right number)'
        assert jim.properties.nrOfRow() == nrow, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of rows not raised or not raised to the right number)'
        assert jim.properties.nrOfPlane() == nplane, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of planes not raised or not raised to the right number)'
        assert jim.properties.nrOfBand() == nband, \
            'Error in geometry.imageFrameSubtract() ' \
            '(number of bands changed)'
        assert (jim.np() == added.np()[1:3, 1:-2, 1:-2]).all(), \
            'Error in geometry.imageFrameSubtract() ' \
            '(changed values in the original Jim)'

class BadGeometryVects(unittest.TestCase):
    """Test functions and methods from geometry module."""

    def test_intersect(self):
        """Test the stack band method."""
        jim = pj.Jim(rasterfn, band=0)
        jimv = pj.JimVect(vectorfn)

        bbox = jim.properties.getBBox()
        new_ulx = (bbox[0] + bbox[2]) / 2.0
        jim_cropped = pj.geometry.crop(jim, ulx=new_ulx, uly=bbox[1],
                                       lrx=bbox[2], lry=bbox[3])

        nr_of_features = jimv.properties.getFeatureCount()

        non_existing_path = pj._get_random_path()

        intersected = pj.geometry.intersect(jimv, jim_cropped,
                                            non_existing_path)
        jimv.geometry.intersect(jim_cropped)

        intersected.io.write()
        jimv.io.write()

        feature_count_func = intersected.properties.getFeatureCount()
        feature_count_meth = jimv.properties.getFeatureCount()

        assert feature_count_func == feature_count_meth, \
            'Inconsistency in geometry.intersect() ' \
            '(method returns different result than function)'
        assert feature_count_meth < nr_of_features, \
            'Error in geometry.intersect() ' \
            '(not less features on intersected area than on the whole)'

        try:
            _ = pj.geometry.intersect(jimv, jimv, pj._get_random_path())
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching a call of geometry.intersect(Jim) ' \
            'function where the argument is not an instance of a Jim object'

        try:
            _ = jimv.geometry.intersect(jimv)
            failed = True
        except TypeError:
            failed = False

        assert not failed, \
            'Error in catching a call of geometry.intersect(Jim) method ' \
            'where the argument is not an instance of a Jim object'


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadGeometry),
                  loader.loadTestsFromTestCase(BadGeometryVects)]
    return unittest.TestSuite(suite_list)
