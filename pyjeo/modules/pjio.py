try:
    import jiplib as _jl
except:
    from jeodpp import jiplib as _jl
import pyjeo as _pj


# def GDALRead(fn, band=0, nXOff=0, nYOff=0, nXSize=None, nYSize=None,
#              nBufXSize=None, nBufYSize=None):
#     """Read a GDAL compatible image stored in the filename.

#     :param fn: a string for the name of an image file (possibly its path)
#     :param band: an integer for the band number, 0 for first band
#     :param nXOff: integer for the pixel offset to the top left corner
#     :param nYOff: integer for the line offset to the top left corner
#     :param nXSize: integer for the width of the region
#     :param nYSize: integer for the height of the region
#     :param nBufXSize: integer for the number of columns of output image
#     :param nBufYSize: integer for the number of lines of output image
#     :return: a Jim object
#     """
#     if nXSize is None:
#         nXSize = nXOff
#     if nYSize is None:
#         nYSize = nYOff
#     if nBufXSize is None:
#         nBufXSize = nXSize
#     if nBufYSize is None:
#         nBufYSize = nYSize
#     return _pj.Jim(_jl.GDALRead(fn, band, nXOff, nYOff, nXSize, nYSize,
#                                nBufXSize, nBufYSize))


def createJim(filename=None, **kwargs):
    """
    Create a new Jim object, either :ref:`from file <create_Jim_from_file>` or :ref:`create new <create_Jim_new>`


    :param filename: Path to a raster dataset or another Jim object
    :param: see supported keys in table below
    :return: a Jim object

    .. _create_Jim_from_file:

    :Create Jim object from file:

    Supported keys as arguments:

    ======== ===================================================
    band     Bands to open, index starts from 0
    ulx      Upper left x value bounding box
    uly      Upper left y value bounding box
    lrx      Lower right x value bounding box
    lry      Lower right y value bounding box
    dx       Resolution in x
    dy       Resolution in y
    resample Resample algorithm used for reading pixel data in case of interpolation
    extent   get boundary from extent from polygons in vector dataset
    nodata   Nodata value to put in image
    noread   Set this flag to True to not read data when opening
    ======== ===================================================

    .. note:: You can specify a different spatial reference system to define the region of interest to read set with keys ulx, uly, lrx, and lry with the extra key 't_srs'. Notice this will not re-project the resulting image. You can use the function :py:func:`geometry:warp` for this.
    ..

    .. note:: resample values: please check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a

    Example:

    Create Jim image object by opening an existing file (file content will automatically be read in memory)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=jl.createJim(ifn)
        #do stuff with jim ...
        jim.close()

    Create Jim image object by opening an existing file for specific region of interest and spatial resolution using cubic convolution resampling::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim0=jl.createJim(ifn,'noread'=True)
        ULX=jim0.getUlx()
        ULY=jim0.getUly()
        LRX=jim0.getUlx()+100*jim0.getDeltaX()
        LRY=jim0.getUly()-100*jim0.getDeltaY()
        jim=jl.Jim.createImg(ifn,ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,dx=5,dy=5,resample='GRIORA_Cubic')
        #do stuff with jim ...
        jim.close()

    .. _create_Jim_new:

    :Create a new Jim image object by defining image attributes (not read from file):

    Supported keys as arguments:

    ===== =================
    ncol  Number of columns
    nrow  Number of rows
    nband (default: 1) Number of bands
    otype (default: Byte) Data type ({Byte/Int16/UInt16/UInt32/Int32/Int64/UInt64/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64})
    a_srs Assign the spatial reference for the output file, e.g., psg:3035 to use European projection and force to European grid
    ===== =================

    Supported keys used to initialize random pixel values in new Jim image object:

    ======= ============================================
    seed    (default: 0) seed value for random generator
    mean    (default: 0) Mean value for random generator
    stdev   (default: 0) Standard deviation for Gaussian random generator
    uniform (default: 0) Start and end values for random value with uniform distribution
    ======= ============================================

    Create a new georeferenced Jim image object by defining the projection epsg code, bounding box, and pixel size::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=jl.Jim.createImg(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',dx=100,dy=100)
        #do stuff with jim ...
        jim.close()

    Create a new georeferenced Jim image object for writing by defining the projection epsg code, bounding box and number of rows and columns::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=jl.Jim.createImg(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',ncol=1098,nrow=1098)
        #do stuff with jim ...
        jim.close()
    """
    return _pj.Jim(_jl.createJim(filename, **kwargs))


def createVector(filename, **kwargs):
    """Create an empty VectorOgr object.

    Created object is an instance of the basis vector class of the pyJEO
    library

    :param filename: Path to a vector dataset or another VectorOgr object
    :param: see supported keys in table below
    :return: a VectorOgr object
    """
    return _jl.createVector(filename, **kwargs)


# def delJimObject(jim_object):
#     if isinstance(jim_object, _pj.Jim):
#         del jim_object.properties, jim_object.io, jim_object.pixops, \
#             jim_object.ngbops, jim_object.geometry, jim_object.ccops, \
#             jim_object.clssfy, jim_object.demops, jim_object.all


class _IO():
    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def write(self, filename, **kwargs):
        """ Write the raster dataset to file in a GDAL supported format

        : param filename: output filename to write to

        Supported keys as arguments:

        ======== ===================================================
        oformat  (default: GTiff) Output image (GDAL supported) format
        co       Creation option for output file. Multiple options can be specified as a list
        nodata   Nodata value to put in image
        ======== ===================================================

        .. note::

           Supported GDAL output formats are restricted to those that support creation (see http://www.gdal.org/formats_list.html#footnote1)

        Example:

        Create Jim image object by opening an existing file in jp2 format. Then write to a compressed and tiled file in the default GeoTIFF format::

           ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
           jim=pj.io.createJim({'filename':ifn})
           jim.io.write('/tmp/test.tif','co':['COMPRESS=LZW','TILED=YES']})
        """

        kwargs.update({'filename': filename})
        self._jim_object.write(kwargs)

    def dumpImg(self, **kwargs):
        """ Dump the raster dataset to output (standard output or ASCII file).

        Supported keys as arguments:

        =========  =============================================================
        output     Output ascii file (Default is empty: dump to standard output)
        oformat    Output format: matrix or list (x,y,z) form. Default is matrix format
        geo        (bool) Set to True to dump x and y in spatial reference system of raster dataset (for list form only). Default is to dump column and row index (starting from 0)
        band       Band index to dump
        srcnodata  Do not dump these no data values (for list form only)
        force      (bool) Set to True to force full dump even for large images
        =========  =============================================================
        """
        self._jim_object.dumpImg(kwargs)
