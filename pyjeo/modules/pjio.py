"""Module for input-output operations."""

import os as _os
from osgeo import ogr as _ogr
import pyjeo as _pj
import jiplib as _jl

class _IO(_pj.modules.JimModuleBase):
    """Define all IO methods."""

    def close(self):
        """Close Jim object."""
        self._jim_object._jipjim.close()

    def dumpImg(self,
                **kwargs):
        """Dump the raster dataset to output (standard output or ASCII file).

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

    def dumpImg3D(self,
                  x: float,
                  y: float,
                  z: float,
                  nx: int,
                  ny: int):
        """
        Dump on screen a dx*dy window with the image values around coordinates.

        Dumped within the plane z.

        :param x: x coordinate
        :param y: y coordinate
        :param z: z coordinate
        :param nx: integer for size of window along x-axis
        :param ny: integer for size of window along y-axis
        """
        self._jim_object._jipjim.imageDump(x, y, z, nx, ny)

    def open(self, image=None, **kwargs):
        """Open Jim object.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        ajim=_pj.Jim(image, **kwargs)
        self._jim_object._set(ajim._jipjim)

    def write(self,
              filename: str,
              **kwargs):
        """Write the raster dataset to file in a GDAL supported format.

        :param filename: output filename to write to

        Supported keys as arguments:

        ======== ===================================================
        oformat  (default: GTiff) Output image (GDAL supported) format
        co       Creation option for output file. Multiple options can be
                 specified as a list
        nodata   Nodata value to put in image
        ======== ===================================================

        .. note::

           Supported GDAL output formats are restricted to those that support
           `creation <http://www.gdal.org/formats_list.html#footnote1>`_.

        Example:

        Create Jim image object by opening an existing file in jp2 format.
        Then write to a compressed and tiled file in the default GeoTIFF
        format::

           ifn = '/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
           jim = pj.Jim({'filename': ifn})
           jim.io.write('/tmp/test.tif', 'co': ['COMPRESS=LZW', 'TILED=YES']})
        """
        kwargs.update({'filename': filename})
        self._jim_object._jipjim.write(kwargs)


class _IOList(_pj.modules.JimListModuleBase):
    """Define all IO methods for JimLists."""

    def close(self):
        """Close all Jim object in the JimList object."""
        self._jim_list._jipjimlist.close()


class _IOVect(_pj.modules.JimVectModuleBase):
    """Define all IO methods for JimVects."""

    def close(self):
        """Close JimVect object."""
        self._jim_vect._jipjimvect.close()

    def dumpVect(self,
                 name=None,
                 output: str = None,
                 **kwargs):
        """
        Dump vector content to screen or file (if output argument is provided).

        :param name: the list of field name(s) to dump
            (default is empty: dump all fields)
        :param output: output ascii file
            (default is empty: dump to standard output)
        """
        if not kwargs:
            kwargs = {}

        if name:
            kwargs.update({'name': name})
        if output is not None:
            kwargs.update({'output': output})

        self._jim_vect._jipjimvect.dumpOgr(kwargs)

    def open(self, vector=None, **kwargs):
        """Open JimVect object.

        :param vector: path to a vector dataset or another JimVect object
        """

        keys = kwargs.keys()

        if isinstance(vector, _pj.JimVect):
            if 'output' not in keys:
                raise AttributeError(
                    "Parameter output required for copy constructor")
        elif not vector:
            if 'output' in keys and not _os.path.isfile(kwargs['output']):
                raise AttributeError('Output path does not exist and the '
                                     'template vector is not specified')
        if kwargs:
            if vector:
                if isinstance(vector, _pj.JimVect):
                    kwargs.update({'filename': kwargs.pop('output', None)})
                    kwargs.update({'access':1})
                    print(kwargs)
                    self._jim_vect._jipjimvect.open(kwargs)
                    self._jim_vect._jipjimvect.copy(vector._jipjimvect, kwargs)
                    self._jim_vect.io.write()

                else:
                    kwargs.update({'filename': vector})
                    self._jim_vect._jipjimvect.open(kwargs)
            elif 'wkt' in kwargs:
                geom = _ogr.CreateGeometryFromWkt(kwargs.pop('wkt'))
                self._jim_vect._jipjimvect.open(geom.ExportToJson())
            else:
                kwargs.update({'filename': kwargs.pop('output', None)})
                self._jim_vect._jipjimvect.open(kwargs)
        else:
            if vector:
                self._jim_vect._jipjimvect.open(vector, **kwargs)
            else:
                self._jim_vect._jipjimvect.open()

        # if isinstance(vector,str):
        #     kwargs.update({'filename':kwargs.pop('output')})
        #     self._jim_vect._jipjimvect.open(kwargs)

        # elif isinstance(vector,_pj.JimVect):
        #     # self._jim_vect._jipjimvect.open(kwargs.pop('output'), kwargs.pop('oformat',None))
        #     kwargs.update({'filename':kwargs.pop('output')})
        #     kwargs.update({'access':1})
        #     print(kwargs)
        #     self._jim_vect._jipjimvect.open(kwargs)
        #     self._jim_vect._jipjimvect.copy(vector._jipjimvect, kwargs)
        #     self._jim_vect.io.write()

    def write(self,
              filename: str = None):
        """Write JimVect object to file.

        If no filename is provided, the original filename with which
        the JimVect object was created will be used.

        :param filename: path to a raster dataset or another Jim object
        """
        if filename:
            self._jim_vect._jipjimvect.write(filename)
        else:
            self._jim_vect._jipjimvect.write()
