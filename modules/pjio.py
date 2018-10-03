import jiplib as _jl
import pyjeo as _pj


def GDALRead(fn, band=0, nXOff=0, nYOff=0, nXSize=None, nYSize=None,
             nBufXSize=None, nBufYSize=None):
    """Read a GDAL compatible image stored in the filename.

    :param fn: a string for the name of an image file (possibly its path)
    :param band: an integer for the band number, 0 for first band
    :param nXOff: integer for the pixel offset to the top left corner
    :param nYOff: integer for the line offset to the top left corner
    :param nXSize: integer for the width of the region
    :param nYSize: integer for the height of the region
    :param nBufXSize: integer for the number of columns of output image
    :param nBufYSize: integer for the number of lines of output image
    :return: a Jim object
    """
    if nXSize is None:
        nXSize = nXOff
    if nYSize is None:
        nYSize = nYOff
    if nBufXSize is None:
        nBufXSize = nXSize
    if nBufYSize is None:
        nBufYSize = nYSize
    return _pj.Jim(_jl.GDALRead(fn, band, nXOff, nYOff, nXSize, nYSize,
                               nBufXSize, nBufYSize))


def createJim(filename, readData=True, **kwargs):
    """

    Supported keys as arguments:

    ======== ===================================================
    filename input filename to read from (GDAL supported format)
    nodata   Nodata value to put in image
    band     Bands to open, index starts from 0
    ulx      Upper left x value bounding box
    uly      Upper left y value bounding box
    lrx      Lower right x value bounding box
    lry      Lower right y value bounding box
    dx       Resolution in x
    dy       Resolution in y
    resample Resample algorithm used for reading pixel data in case of interpolation (default: GRIORA_NearestNeighbour). Check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a for available options.
    extent   get boundary from extent from polygons in vector dataset
    noread   Set this flag to True to not read data when opening
    ======== ===================================================

    .. note::
    You can specify a different spatial reference system to define the region of interest to read set with keys ulx, uly, lrx, and lry with the extra key 't_srs'. Notice this will not re-project the resulting image. You can use the function :py:func:Jim:`warp` for this.
    ..

    resample: (default: GRIORA_NearestNeighbour) Resample algorithm used for reading pixel data in case of interpolation GRIORA_NearestNeighbour | GRIORA_Bilinear | GRIORA_Cubic | GRIORA_CubicSpline | GRIORA_Lanczos | GRIORA_Average | GRIORA_Average | GRIORA_Gauss (check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a)

    Supported keys when creating new Jim image object not read from file:

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

    Returns:
    This instance of Jim object (self)

    Example:

    Create Jim image object by opening an existing file (file content will automatically be read in memory)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=jl.createJim('filename'=ifn)
        #do stuff with jim ...
        jim.close()

    The key 'filename' is optional::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=jl.createJim(ifn)
        #do stuff with jim ...
        jim.close()

    Create Jim image object by opening an existing file for specific region of interest and spatial resolution using cubic convolution resampling::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim0=jl.createJim(filename=ifn,'noread'=True)
        ULX=jim0.getUlx()
        ULY=jim0.getUly()
        LRX=jim0.getUlx()+100*jim0.getDeltaX()
        LRY=jim0.getUly()-100*jim0.getDeltaY()
        jim=jl.Jim.createImg(filename=ifn,ulx:ULX,'uly':ULY,'lrx':LRX,'lry':LRY,'dx':5,'dy':5,'resample':'GRIORA_Cubic'})
        #do stuff with jim ...
        jim.close()

    Create a new georeferenced Jim image object by defining the projection epsg code, bounding box, and pixel size::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        dict={'ulx':ULX,'uly':ULY,'lrx':LRX,'lry':LRY,'a_srs':projection}
        dict.update({'otype':'GDT_UInt16'})
        dict.update({'dy':100,'dx':100})
        jim=jl.Jim.createImg(dict)
        #do stuff with jim ...
        jim.close()

    Create a new georeferenced Jim image object for writing by defining the projection epsg code, bounding box and number of rows and columns::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        dict={'ulx':ULX,'uly':ULY,'lrx':LRX,'lry':LRY,'a_srs':projection}
        dict.update({'otype':'GDT_UInt16'})
        nrow=1098
        ncol=1098
        dict.update({'nrow':nrow,'ncol':ncol})
        jim=jl.Jim.createImg(dict)
        #do stuff with jim ...
        jim.close()
    """
    return _pj.Jim(_jl.createJim(filename, readData, **kwargs))


def createVector(filename, readData=True, **kwargs):
    """Create an empty VectorOgr object.

    Created object is an instance of the basis vector class of the pyJEO
    library

    :param filename: Path to a vector file or another VectorOgr object
    :param readData: boolean parameter whether read data or create just shuck
    :return: a VectorOgr object
    """
    return _jl.createVector(filename, readData, **kwargs)


class _IO():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object
