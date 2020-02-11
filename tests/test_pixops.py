"""Test suite for module pyjeo.pixops."""

import pyjeo as pj
import unittest


testFile = 'tests/data/modis_ndvi_2010.tif'
tiles = ['tests/data/red1.tif', 'tests/data/red2.tif']
vector = 'tests/data/nuts_italy.sqlite'


class BadPixOps(unittest.TestCase):
    """Test functions and methods from pixops modules."""

    def test_isEqual(self):
        """Test isEqual() function."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[0])
        jim3 = pj.Jim(tiles[1])

        assert jim1.pixops.isEqual(jim2), \
            'Error in pixops.isEqual() method (for two Jims)'
        assert not jim1.pixops.isEqual(jim3), \
            'Error in pixops.isEqual() method (non-equality)'
        assert not jim1.pixops.isEqual(1), \
            'Error in pixops.isEqual() method (for not Jim object)'

        assert pj.pixops.isEqual(jim1, jim2), \
            'Error in pixops.isEqual() function (for two Jims)'
        assert not pj.pixops.isEqual(jim1, jim3), \
            'Error in pixops.isEqual() function (non-equality)'
        assert not pj.pixops.isEqual(jim1, 1), \
            'Error in pixops.isEqual() function (for not Jim object)'

        jim_oneplane = pj.Jim(ncol=5, nrow=5, nband=2, otype='GDT_Byte')
        jim_twoplanes = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                               otype='GDT_Byte')

        jim_oneplane.pixops.setData(5)
        jim_twoplanes.pixops.setData(5)

        assert not jim_oneplane.pixops.isEqual(jim_twoplanes), \
            'Error in Jim.pixops.isEqual() ' \
            '(not False for different number of planes)'

        assert not pj.pixops.isEqual(jim_oneplane, jim_twoplanes), \
            'Error in pj.pixops.isEqual() ' \
            '(not False for different number of planes)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(5)

        assert jim_twoplanes.pixops.isEqual(jim_twoplanes2), \
            'Error in pj.pixops.isEqual() ' \
            '(wrong result for multiplane Jim object)'

        assert pj.pixops.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.pixops.isEqual() ' \
            '(wrong result for multiplane Jim object)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=3, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(5)

        assert not jim_twoplanes.pixops.isEqual(jim_twoplanes2), \
            'Error in Jim.pixops.isEqual() ' \
            '(not False for different number of bands for multiplane Jim)'

        assert not pj.pixops.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.pixops.isEqual() ' \
            '(not False for different number of bands for multiplane Jim)'

        jim_twoplanes2 = pj.Jim(ncol=5, nrow=5, nband=2, nplane=2,
                                otype='GDT_Byte')
        jim_twoplanes2.pixops.setData(6)

        assert not jim_twoplanes.pixops.isEqual(jim_twoplanes2), \
            'Error in pj.pixops.isEqual() ' \
            '(wrong result for multiplane multiband Jim object)'

        assert not pj.pixops.isEqual(jim_twoplanes, jim_twoplanes2), \
            'Error in pj.pixops.isEqual() ' \
            '(wrong result for multiplane multiband Jim object)'

    def test_NDVI(self):
        """Test computing NDVI in different ways."""
        jim = pj.Jim(testFile)

        jim_red = pj.geometry.cropBand(jim, 0)
        jim_nir = pj.geometry.cropBand(jim, 1)

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        ndvi = pj.pixops.NDVI(jim, 0, 1)
        jim.pixops.NDVI(0, 1)

        assert jim.pixops.isEqual(ndvi), 'Error in computing NDVI'

        ndvi = pj.pixops.NDVISeparateBands(jim_red, jim_nir)
        jim_red.pixops.NDVISeparateBands(jim_nir)

        assert jim_red.pixops.isEqual(ndvi), 'Error in computing ' \
                                             'NDVISeparateBands'

        assert jim_red.pixops.isEqual(jim), 'Error in computing NDVI or ' \
                                            'NDVISeparateBands'

        jim_nir.pixops.convert('Float32')
        ndvi = pj.pixops.NDVISeparateBands(jim_red, jim_nir)

        assert not jim.pixops.isEqual(ndvi), 'Error in computing NDVI'

    def test_supremum_infimum(self):
        """Test picking up supremum and infimum from computed NDVIs."""
        for tile in tiles:
            jim4 = pj.Jim(tile)
            jim8 = pj.Jim(tile[:-8] + 'nir' + tile[-5] + '.tif')

            if 'max_ndvi' and 'min_ndvi' in locals():
                b = pj.pixops.NDVISeparateBands(jim4, jim8)

                max_ndvi[b > max_ndvi] = b
                max_ndvi_func = pj.pixops.supremum(max_ndvi, b)
                max_ndvi2.pixops.supremum(b)
                max_ndvi3 = pj.pixops.supremum(pj.JimList([max_ndvi2, b]))
                max_ndvi4 = pj.Jim(max_ndvi2)
                max_ndvi4.pixops.supremum(pj.JimList([b, max_ndvi2]))
                max_ndvi5 = pj.JimList([max_ndvi2, b]).pixops.supremum()

                assert max_ndvi.pixops.isEqual(max_ndvi2), \
                    'Inconsistency in pixops.supremum() ' \
                    '(method returns different result than function)'
                assert max_ndvi.pixops.isEqual(max_ndvi_func), \
                    'Error in pixops.supremum or jim[a>jim] = a'
                assert max_ndvi.pixops.isEqual(max_ndvi3), \
                    'Error in pixops.supremum JimList 3'
                assert max_ndvi.pixops.isEqual(max_ndvi4), \
                    'Error in pixops.supremum JimList 4'
                assert max_ndvi.pixops.isEqual(max_ndvi5), \
                    'Error in pixops.supremum JimList 5'

                min_ndvi[b < min_ndvi] = b
                min_ndvi_func = pj.pixops.infimum(min_ndvi, b)
                min_ndvi2.pixops.infimum(b)
                min_ndvi3 = pj.pixops.infimum(pj.JimList([min_ndvi2, b]))
                min_ndvi4 = pj.Jim(min_ndvi2)
                min_ndvi4.pixops.infimum(pj.JimList([b, min_ndvi2]))
                min_ndvi5 = pj.JimList([min_ndvi2, b]).pixops.infimum()

                assert min_ndvi.pixops.isEqual(min_ndvi2), \
                    'Inconsistency in pixops.infimum() ' \
                    '(method returns different result than function)'
                assert min_ndvi.pixops.isEqual(min_ndvi_func), \
                    'Error in pixops.infimum or jim[a<jim] = a'
                assert min_ndvi.pixops.isEqual(min_ndvi3), \
                    'Error in pixops.infimum JimList 3'
                assert min_ndvi.pixops.isEqual(min_ndvi4), \
                    'Error in pixops.infimum JimList 4'
                assert min_ndvi.pixops.isEqual(min_ndvi5), \
                    'Error in pixops.infimum JimList 5'
            else:
                max_ndvi = pj.pixops.NDVISeparateBands(jim4, jim8)
                max_ndvi2 = pj.Jim(max_ndvi)

                min_ndvi = pj.Jim(max_ndvi)
                min_ndvi2 = pj.Jim(min_ndvi)

    def test_stretch(self):
        """Test stretch function."""
        for tile in tiles:
            jim = pj.Jim(tile)

            stretched = pj.pixops.stretch(jim, otype='GDT_Byte', dst_min=0,
                                          dst_max=255, cc_min=2, cc_max=98)
            stretched_eq = pj.pixops.stretch(jim, otype='GDT_Byte', dst_min=0,
                                             dst_max=255, cc_min=2, cc_max=98,
                                             eq=True)
            jim.pixops.stretch(otype='GDT_Byte', dst_min=0, dst_max=255,
                               cc_min=2, cc_max=98)
            assert jim.pixops.isEqual(stretched), \
                'Inconsistency in pixops.stretch() ' \
                '(method returns different result than function)'

            theStats = stretched_eq.stats.getStats(['min', 'max', 'histogram'],
                                                   src_min=1, src_max=254)
            assert theStats['min'] == 1, \
                'Error in pixops.stretch(): min is not 1'
            assert theStats['max'] == 254, \
                'Error in pixops.stretch(): max is not 254'
            assert max(theStats['histogram']) <= 1150, \
                'Error in pixops.stretch(): max is not < 1150'
            assert min(theStats['histogram']) >= 850, \
                'Error in pixops.stretch(): max is not >= 850'

    def test_setFunctions(self):
        """Test setData, setLevel and setThreshold functions and methods."""
        jim = pj.Jim(testFile)

        stats = jim.stats.getStats(band=0)
        min = stats['min']
        max = stats['max']

        thresholded = pj.pixops.setThreshold(jim, min=min+1, max=max-1,
                                             nodata=(max+min)/2)

        thresholded_stats = thresholded.stats.getStats(band=0)
        nodata = thresholded.properties.getNoDataVals()

        assert thresholded_stats != stats, \
            'Error in pixops.setThreshold() or stats.getStats(band=0)'
        assert thresholded_stats['min'] >= min + 1, \
            'Error in pixops.setThreshold() or stats.getStats(band=0)'
        assert thresholded_stats['max'] <= max - 1, \
            'Error in pixops.setThreshold() or stats.getStats(band=0)'
        assert thresholded.properties.getNoDataVals() == [((max+min)/2)], \
            'Error in pixops.setThreshold() or properties.getNoDataVals()'

        jim.pixops.setThreshold(min=min+1, max=max-1, nodata=(max+min)/2)

        stats_jim_thresholded = jim.stats.getStats(band=0)

        assert stats_jim_thresholded == thresholded_stats, \
            'Error in pixops.setThreshold() or stats.getStats(band=0)'
        assert jim.properties.getNoDataVals() == nodata, \
            'Error in pixops.setThreshold() or properties.getNoDataVals()'

        jim = pj.Jim(testFile)

        levelled = pj.pixops.setLevel(jim, min=min+1, max=max-1, val=max+1)
        jim.pixops.setLevel(min=min+1, max=max-1, val=max+1)

        levelled_stats = levelled.stats.getStats(band=0)

        assert jim.pixops.isEqual(levelled), \
            'Inconsistency in pixops.setLevel() ' \
            '(method returns different result than function)'
        assert levelled_stats != stats, \
            'Error in pixops.setLevel() (object not changed)'
        assert levelled_stats['min'] == min, \
            'Error in pixops.setLevel() (min value changed)'
        assert levelled_stats['max'] == max + 1, \
            'Error in pixops.setLevel() (value not setted)'

        fives = pj.pixops.setData(jim, 5, bands=[0, 1])
        jim.pixops.setData(5, bands=[0, 1])

        assert jim.pixops.isEqual(fives), \
            'Inconsistency in pixops.setData() ' \
            '(method returns different result than function)'

        fives = pj.pixops.setData(jim, 5)
        jim.pixops.setData(5)

        assert jim.pixops.isEqual(fives), \
            'Inconsistency in pixops.setData() ' \
            '(method returns different result than function)'

        stats = jim.stats.getStats(band=0)

        assert all(v == 5 for v in [stats['min'],
                                    stats['max'],
                                    stats['mean']]),\
            'Error in pixops.setData() or stats.getStats(band=0)'

        jim.pixops.setData(10, dx=jim.properties.getDeltaX()+0.1, bands=[0, 1])

        stats = jim.stats.getStats(band=0)

        assert all(v == 10 for v in [stats['min'],
                                     stats['max'],
                                     stats['mean']]), \
            'Error in pixops.setData() or stats.getStats(band=0)'

    def test_convert(self):
        """Test data type conversions."""
        jim = pj.Jim(tiles[0])

        jim.pixops.convert('Byte')
        assert jim.properties.getDataType() == 'Byte', \
            'Error in pixops.convert()'

        jim.pixops.convert('UInt16')
        assert jim.properties.getDataType() == 'UInt16', \
            'Error in pixops.convert()'

        jim.pixops.convert('Int16')
        assert jim.properties.getDataType() == 'Int16', \
            'Error in pixops.convert()'

        jim.pixops.convert('UInt32', scale=1)
        assert jim.properties.getDataType() == 'UInt32', \
            'Error in pixops.convert()'
        # TODO: Change the assertion to 4 after fixing the bug in jiplib

        jim.pixops.convert('Int32')
        assert jim.properties.getDataType() == 'Int32', \
            'Error in pixops.convert()'

        jim.pixops.convert('Float32')
        assert jim.properties.getDataType() == 'Float32', \
            'Error in pixops.convert()'

        jim.pixops.convert('Float64')
        assert jim.properties.getDataType() == 'Float64', \
            'Error in pixops.convert()'

        # jim.pixops.convert('CInt16')
        # assert jim.properties.getDataType() == 8, 'Error in pixops.convert()'
        #
        # jim.pixops.convert('CInt32')
        # assert jim.properties.getDataType() == 9, \
        #   'Error in pixops.convert()'
        #
        # jim.pixops.convert('CFloat32')
        # assert jim.properties.getDataType() == 10, \
        # 'Error in ixops.convert()'
        #
        # jim.pixops.convert('CFloat64')
        # assert jim.properties.getDataType() == 11, \
        # 'Error in pixops.convert()'
        # # TODO: Uncomment when bug fixed in jiplib

        jim.pixops.convert('Byte', a_srs='EPSG:5514')
        assert jim.properties.getDataType() == 'Byte', \
            'Error in pixops.convert()'

        jim.pixops.convert('UInt16', a_srs='EPSG:5514')
        assert jim.properties.getDataType() == 'UInt16', \
            'Error in pixops.convert()'

        jim.pixops.convert('Float32', a_srs='EPSG:5514')
        assert jim.properties.getDataType() == 'Float32', \
            'Error in pixops.convert()'

        jim.pixops.convert('Float64', a_srs='EPSG:5514')
        assert jim.properties.getDataType() == 'Float64', \
            'Error in pixops.convert()'

        try:
            jim.pixops.convert('string')
            raised = False
        except TypeError:
            raised = True

        assert raised, \
            'Error in checks for non-supported data types in pixops.convert()'

        a = pj.pixops.convert(jim, 'Byte')
        assert a.properties.getDataType() == 'Byte', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'UInt16')
        assert a.properties.getDataType() == 'UInt16', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'Int16')
        assert a.properties.getDataType() == 'Int16', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'UInt32')
        assert a.properties.getDataType() == 'UInt32', \
            'Error in pixops.convert()'
        # TODO: Change the assertion to 4 after fixing the bug in jiplib

        a = pj.pixops.convert(a, 'Int32')
        assert a.properties.getDataType() == 'Int32', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'Float32')
        assert a.properties.getDataType() == 'Float32', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'Float64')
        assert a.properties.getDataType() == 'Float64', \
            'Error in pixops.convert()'

        # a = pj.pixops.convert(jim, 'CInt16')
        # assert a.properties.getDataType() == 8, 'Error in pixops.convert()'
        #
        # a = pj.pixops.convert(jim, 'CInt32')
        # assert a.properties.getDataType() == 9, 'Error in pixops.convert()'
        #
        # a = pj.pixops.convert(jim, 'CFloat32')
        # assert a.properties.getDataType() == 10, 'Error in pixops.convert()'
        #
        # a = pj.pixops.convert(jim, 'CFloat64')
        # assert a.properties.getDataType() == 11, 'Error in pixops.convert()'
        # # TODO: Uncomment when bug fixed in jiplib

        a = pj.pixops.convert(a, 'Byte', a_srs='EPSG:5514')
        assert a.properties.getDataType() == 'Byte', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'UInt16', a_srs='EPSG:5514')
        assert a.properties.getDataType() == 'UInt16', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'Float32', a_srs='EPSG:5514')
        assert a.properties.getDataType() == 'Float32', \
            'Error in pixops.convert()'

        a = pj.pixops.convert(a, 'Float64', a_srs='EPSG:5514')
        assert a.properties.getDataType() == 'Float64', \
            'Error in pixops.convert()'

        try:
            pj.pixops.convert(a, 'string')
            raised = False
        except TypeError:
            raised = True

        assert raised, \
            'Error in checks for non-supported data types in pixops.convert()'

    def test_histoCompress(self):
        """Test histoCompress() function and method."""
        jim = pj.Jim(tiles[0])

        destructive_object = pj.Jim(jim)
        compressed = pj.pixops.histoCompress(destructive_object)
        destructive_object.pixops.histoCompress()

        assert destructive_object.pixops.isEqual(compressed), \
            'Error in pixops.histoCompress()'
        assert destructive_object.stats.getStats(band=0)['min'] == 0, \
            'Error in pixops.histoCompress() (minimum value not 0)'

        compressed = pj.pixops.histoCompress(jim, 0)
        jim.pixops.histoCompress(0)

        assert jim.pixops.isEqual(compressed), \
            'Error in pixops.histoCompress()'
        assert jim.stats.getStats(band=0)['min'] == 0, \
            'Error in pixops.histoCompress() (minimum value not 0)'

    def test_simple_op(self):
        """Test simpleArithOp() and simpleBitwiseOp() functions and methods."""
        jim = pj.Jim(tiles[0])

        # Test simpleArithOp()
        ones = pj.pixops.setData(jim, 1)
        jim_plus_one = jim + 1
        jim_plus_one_arith = pj.pixops.simpleArithOp(jim, ones, 0)
        jim.pixops.simpleArithOp(ones, 0)

        assert jim.pixops.isEqual(jim_plus_one_arith), \
            'Inconsistency in pixops.simpleArithOp() ' \
            '(method returns different result than function)'

        assert jim.pixops.isEqual(jim_plus_one), \
            'Error in pixops.simpleArithOp(op=0) or Jim + int ' \
            '(Results not equal)'

        jim_minus_one = jim - 1
        jim_minus_one_arith = pj.pixops.simpleArithOp(jim, ones, 1)
        jim.pixops.simpleArithOp(ones, 1)

        assert jim.pixops.isEqual(jim_minus_one_arith), \
            'Inconsistency in pixops.simpleArithOp() ' \
            '(method returns different result than function)'

        assert jim.pixops.isEqual(jim_minus_one), \
            'Error in pixops.simpleArithOp(op=1) or Jim - int ' \
            '(Results not equal)'

        # Test simpleBitwiseOp()
        jim_and = jim & ones
        jim_and_bitwise = pj.pixops.simpleBitwiseOp(jim, ones, 10)
        jim.pixops.simpleBitwiseOp(ones, 10)

        assert jim.pixops.isEqual(jim_and_bitwise), \
            'Inconsistency in pixops.simpleBitwiseOp() ' \
            '(method returns different result than function)'

        assert jim.pixops.isEqual(jim_and), \
            'Error in pixops.simpleBitwiseOp(op=10) or Jim & Jim ' \
            '(Results not equal)'

        jim_or = jim | jim_minus_one
        jim_or_bitwise = pj.pixops.simpleBitwiseOp(jim, jim_minus_one, 11)
        jim.pixops.simpleBitwiseOp(jim_minus_one, 11)

        assert jim.pixops.isEqual(jim_or_bitwise), \
            'Inconsistency in pixops.simpleBitwiseOp() ' \
            '(method returns different result than function)'

        assert jim.pixops.isEqual(jim_or), \
            'Error in pixops.simpleBitwiseOp(op=11) or Jim & Jim ' \
            '(Results not equal)'


class BadPixOpsLists(unittest.TestCase):
    """Test JimList functions and methods from pixops modules."""

    def test_composite(self):
        """Test composite() function."""
        jim1 = pj.Jim(tiles[0])
        jim2 = pj.Jim(tiles[1])
        jim3 = pj.Jim(tiles[0][:-8] + 'nir' + tiles[0][-5] + '.tif')
        jiml1 = pj.JimList([jim1, jim2])
        jiml2 = pj.JimList([jim1, jim3])

        # TODO: Suppress output originating in jiplib (flag `quiet`, please?)
        comp1 = pj.pixops.composite(jiml1)
        comp2 = jiml1.pixops.composite()
        comp3 = jiml2.pixops.composite()
        comp4 = jiml1.pixops.composite(otype='Byte')

        assert comp1.pixops.isEqual(comp2)
        assert not comp1.pixops.isEqual(comp3)
        assert not comp2.pixops.isEqual(comp4)


def load_tests(loader=None, tests=None, pattern=None):
    """Load tests."""
    if not loader:
        loader = unittest.TestLoader()
    suite_list = [loader.loadTestsFromTestCase(BadPixOps),
                  loader.loadTestsFromTestCase(BadPixOpsLists)]
    return unittest.TestSuite(suite_list)
