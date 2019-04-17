import pyjeo as pj
import unittest

rasterfn = 'tests/data/modis_ndvi_2010.tif'
vectorfn = 'tests/data/modis_ndvi_training.sqlite'
outputfn = 'tests/data/modis_ndvi_features.sqlite'
warpedfn = 'tests/data/modis_ndvi_warped.tif'


class BadGeometry(unittest.TestCase):

    def test_stack(self):
        """Test the stack method."""
        jim0 = pj.Jim(rasterfn, band=0)
        jim1 = pj.Jim(rasterfn, band=1)
        jim2 = pj.Jim(rasterfn, band=2)
        jimlist = pj.JimList([jim0,jim1,jim2])
        jimliststack = pj.geometry.stackBand(jimlist)
        jimstack = pj.geometry.stackBand(pj.geometry.stackBand(jim0,jim1),jim2)
        assert jimliststack.pixops.isEqual(jimstack), \
            'Error in geometry.stackBand(): jimliststack not equal to jimstack'
        jimliststack = pj.geometry.stackBand(jimlist, jim0)
        jimliststack.geometry.cropBand([0,1,2])
        assert jimliststack.pixops.isEqual(jimstack), \
            'Error in geometry.stackBand(): jimliststack not equal to jimstack after crop'
        jim3=pj.JimList([jim0,jim1,jim2]).geometry.stackBand()
        jim3.geometry.cropBand([0,1])
        assert pj.geometry.cropBand(jim3,0).pixops.isEqual(jim0), \
            'Error jim3 not equal to jim0'
        jim4=pj.JimList([jim0,jim1,jim2]).geometry.stackBand()
        jim4.geometry.cropBand([1,2])
        assert pj.geometry.cropBand(jim4,0).pixops.isEqual(jim1), \
            'Error jim3 not equal to jim1'
        jim=pj.JimList([jim3,jim4]).geometry.stackBand([0,1])
        jim.geometry.cropBand([1,2])
        jim.geometry.cropBand(0)
        assert jim.pixops.isEqual(jim1), \
            'Error jim not equal to jim1'


    def test_warp(self):
        """Test the warp method."""
        jim0 = pj.Jim(rasterfn)
        jim_warped = pj.geometry.warp(jim0,'epsg:4326')

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
                v2=jl0.geometry.extractOgr(sample,rule='mean',output='/vsimem/v2.sqlite',oformat='SQLite',co=['OVERWRITE=YES'],bandname=bandname,fid='fid')
                v=pj.geometry.join(v1, v2,output=outputfn,oformat='SQLite',co=['OVERWRITE=YES'],key=['fid']);
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

def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadGeometry)]
    return unittest.TestSuite(suite_list)
