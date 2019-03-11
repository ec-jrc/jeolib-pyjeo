"""Module for input-output operations."""

import pyjeo as _pj
try:
    import jiplib as _jl
except ImportError:
    from jeodpp import jiplib as _jl


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


# def createJim(filename=None, **kwargs):
#     """
#     :param filename: Path to a raster dataset or another Jim object
#     :param: see supported keys in table below
#     :return: a Jim object
#     """

#     return _pj.Jim(_jl.createJim(filename, **kwargs))


# def createJimList(jims_list, **kwargs):
#     return _pj.JimList(_jl.createJimList(jims_list, **kwargs))


# def createVector(filename, **kwargs):
#     """Create an empty VectorOgr object.

#     Created object is an instance of the basis vector class of the pyJEO
#     library

#     :param filename: Path to a vector dataset or another VectorOgr object
#     :param: see supported keys in table below
#     :return: a VectorOgr object
#     """
#     return _pj.JimVect(_jl.createVector(filename, **kwargs))


class _IO():
    """Define all IO methods."""

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

           Supported GDAL output formats are restricted to those that support `creation <http://www.gdal.org/formats_list.html#footnote1>`_.

        Example:

        Create Jim image object by opening an existing file in jp2 format. Then write to a compressed and tiled file in the default GeoTIFF format::

           ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
           jim=pj.Jim({'filename':ifn})
           jim.io.write('/tmp/test.tif','co':['COMPRESS=LZW','TILED=YES']})
        """
        kwargs.update({'filename': filename})
        self._jim_object._jipjim.write(kwargs)

    def dumpImg(self, **kwargs):
        """ Dump the raster dataset to output (standard output or ASCII file).

        Supported keys as arguments:

        =========  ============================================================
        output     Output ascii file (Default is empty: dump to standard
                   output)
        oformat    Output format: matrix or list (x,y,z) form. Default is
                   matrix format
        geo        (bool) Set to True to dump x and y in spatial reference
                   system of raster dataset (for list form only). Default is
                   to dump column and row index (starting from 0)
        band       Band index to dump
        srcnodata  Do not dump these no data values (for list form only)
        force      (bool) Set to True to force full dump even for large images
        =========  ============================================================
        """
        self._jim_object._jipjim.dumpImg(kwargs)

    def dumpImg3D(self, x, y, z, nx, ny):
        """
        Dump on screen a dx*dy window with the image values around
        coordinates (x,y) and within the plane z.

        :param x: x coordinate
        :param y: y coordinate
        :param z: z coordinate
        :param nx: integer for size of window along x-axis
        :param ny: integer for size of window along y-axis
        """
        self._jim_object._jipjim.imageDump(x, y, z, nx, ny)

    def close(self):
        """
        Close Jim object

        """
        self._jim_object._jipjim.close()

class _IOList():
    """Define all IO methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def close(self):
        """
        Close all Jim object in the JimList object
        """
        self._jim_list._jipjimlist.close()

class _IOVect():
    """Define all IO methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    def write(self, filename=None):
        """
        Write JimVect object to file. If no filename is provided, the original filename with which the JimVect object was created will be used.

        :param filename: path to a raster dataset or another Jim object
        """
        if filename:
            self._jim_vect._jipjimvect.write(filename)
        else:
            self._jim_vect._jipjimvect.write()

    def close(self):
        """
        Close JimVect object.

        """
        self._jim_vect._jipjimvect.close()

    def dumpVect(self, name=None, **kwargs):
        """
        Dump vector content to screen or file (if output argument is provided)

        :param name: the list of field name(s) to dump (default is empty: dump all fields)
        :param output: output ascii file (default is empty: dump to standard output)
        """
        if name:
            kwargs.update({'name': name})
        self._jim_vect._jipjimvect.dumpOgr(**kwargs)
