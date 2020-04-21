"""Module for pixel-wise operations."""

import numpy as _np

import pyjeo as _pj
import jiplib as _jl


def composite(jim_list,
              crule: str = 'overwrite',
              **kwargs):
    """Composite Jims in a JimList using a composite rule.

    :param jim_list: List of Jims to composite
    :param crule: Rule for the composition
    :return: composited Jim object
    """
    kwargs.update({'crule': crule})
    return _pj.Jim(jim_list._jipjimlist.composite(kwargs))


def convert(jim_object,
            otype):
    """Convert Jim image with respect to data type.

    :param jim_object: Jim object to be used for the conversion
    :param otype: Data type for output image
    :return: a Jim object


    Example:

    Convert data type of input image to byte::

        jim0 = pj.Jim('/path/to/raster.tif')
        jim0.pixops.convert(otype='Byte')

    Clip raster dataset between 0 and 255 (set all other values to 0),
    then convert data type to byte::

        jim1 = pj.Jim('/path/to/raster.tif')
        jim1.pixops.setThreshold(min=0, max=255, nodata=0)
        jim1.pixops.convert('Byte')
    """
    if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
        otype = 'GDT_Byte'
    elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
        otype = 'GDT_UInt16'
    elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
        otype = 'GDT_Int16'
    elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
        otype = 'GDT_UInt32'
    elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
        otype = 'GDT_Int32'
    elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
        otype = 'GDT_Float32'
    elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
        otype = 'GDT_Float64'
    elif otype in ['int64', 'Int64', 'JDT_Int64']:
        otype = 'JDT_Int64'
    elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
        otype = 'JDT_UInt64'
    elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
        otype = 'JDT_UInt64'
    else:
        raise TypeError("Output type {} not supported".format(otype))
    # TODO: Support CTypes

    return _pj.Jim(jim_object._jipjim.convertDataType(otype))


def histoCompress(jim_object,
                  band: int = None):
    """Redistribute the intensity of histogram to fit full range of values.

    Redistributes the intensity values of im in such a way that the
    minimum value is zero and that all subsequent intensity values are
    used until the maximum value (i.e. no gaps).

    :param jim_object: Jim object to be used for the histoCompress
    :param band: from which to get the histogram
    :return: a Jim object
    """
    if band is not None:
        return _pj.Jim(jim_object._jipjim.histoCompress(band))
    else:
        return _pj.Jim(jim_object._jipjim.histoCompress())


def infimum(jim,
            *args):
    """Create Jim composed using minimum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of smalles values from provided Jim objects
    """
    if isinstance(jim, _pj.JimList):
        inf = None
        for newJim in jim:
            if inf is None:
                inf = _pj.Jim(jim[0])
            else:
                inf._jipjim.d_pointOpArith(newJim._jipjim, 4)
    else:
        inf = _pj.Jim(jim)
    for newJim in args:
        inf._jipjim.d_pointOpArith(newJim._jipjim, 4)

    return inf


# def modulo(jim_object, val):
#     """Set all pixels to their value modulo val.

#     :param val:  modulo value (integer)

#     Modifies the instance on which the method was called.

#     """
#     return _pj.Jim(jim_object._jipjim.pointOpModulo(val))


def NDVI(jim_object,
         band_red: int,
         band_nir: int):
    """Compute NDVI on the Jim object.

    :param jim_object: Jim object from which the red and NIR bands are to be
        derived
    :param band_red: index of band with values of red
    :param band_nir: index of band with values of NIR
    :return: a Jim object with values of NDVI
    """
    red = _pj.geometry.cropBand(jim_object, band_red)
    nir = _pj.geometry.cropBand(jim_object, band_nir)

    return _pj.Jim(nir._jipjim.pointOpNDI(red._jipjim))


def NDVISeparateBands(jim_red,
                      jim_nir):
    """Compute NDVI from two Jim objects.

    Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

    :param jim_red: Jim object with values of red
    :param jim_nir: Jim object with values of NIR
    :return: a Jim object with values of NDVI
    """
    return _pj.Jim(jim_nir._jipjim.pointOpNDI(jim_red._jipjim))


def setData(jim,
            value: float,
            ulx: float = None,
            uly: float = None,
            lrx: float = None,
            lry: float = None,
            bands=None,
            dx: int = 0,
            dy: int = 0,
            nogeo: bool = False):
    """Set range of pixels to value.

    :param jim: a Jim object
    :param value: new value for pixels of Jim object
    :param ulx: upper left corner x coordinate (in projected coordinates
        if geo is True, else in image coordinates)
    :param uly: upper left corner y coordinate (in projected coordinates
        if geo is True, else in image coordinates)
    :param lrx: lower right corner x coordinate (in projected coordinates
        if geo is True, else in image coordinates)
    :param lry: lower right corner y coordinate (in projected coordinates
        if geo is True, else in image coordinates)
    :param bands: List of band indices to crop (index is 0 based)
    :param dx: spatial resolution in x to crop (stride if geo is False)
    :param dy: spatial resolution in y to crop (stride if geo is False)
    :param nogeo: use geospatial coordinates if False, image coordinates
        if True
    :return: a Jim object
    """
    jout = _pj.Jim(jim)

    if bands is None:
        bands = range(jim.properties.nrOfBand())

    jout.pixops.setData(value, ulx, uly, lrx, lry, bands, dx, dy, nogeo)
    return jout


def setLevel(jim_object,
             min: float,
             max: float,
             val: float):
    """Set all pixels within min and max values to val.

    :param jim_object: a Jim object used fas a basis for the level setting
    :param min:  Minimum threshold value
    :param max:  Maximum threshold value
    :param val:  All pixels within [min,max] are set to val
    """
    return _pj.Jim(jim_object._jipjim.pointOpSetLevel(min, max, val))


def setThreshold(jim_object,
                 **kwargs):
    """Apply min and max threshold to pixel values in raster dataset.

    :param jim_object: the Jim object on which to set threshold
    :param kwargs: See table :py:meth:`~pixops._PixOps.setThreshold`.

    for help, please refer to the corresponding
    method :py:meth:`~pixops._PixOps.setThreshold`.
    """
    return _pj.Jim(jim_object._jipjim.setThreshold(kwargs))


def simpleArithOp(jim1,
                  jim2,
                  op: int,
                  *args):
    """Create Jim composed using a simple arithmetic operation (coded with op).

    :param jim1: Jim object
    :param jim2: Jim object (to be sure that at least one is provided)
    :param op: integer for operation type
    :param args: Jim objects
    :return: Jim holding specified arithmetic operation with from provided
        Jim objects
    """
    jout = _pj.Jim(jim1)
    jims = [jim2]
    jims.extend(args)
    for jim in jims:
        jout._jipjim.d_pointOpArith(jim._jipjim, op)

    return _pj.Jim(jout)


def simpleBitwiseOp(jim,
                    another_jim,
                    op: int,
                    *args):
    """Create Jim composed using a simple bitwise operation (coded with op).

    :param jim: Jim object
    :param another_jim: Jim object (to be sure that at least one is provided)
    :param op: integer for operation type
    :param args: Jim objects
    :return: Jim holding specified bitwise operation with from provided
        Jim objects
    """
    jout = _pj.Jim(jim)
    jims = [another_jim]
    jims.extend(args)
    for newJim in jims:
        jout._jipjim.d_pointOpBitwise(newJim._jipjim, op)

    return _pj.Jim(jout)


def simpleThreshold(jim_object,
                    min: float,
                    max: float,
                    bg_val: float,
                    fg_val: float):
    """Set all pixels within min and max values to val.

    :param jim_object: Jim object to be used for the threshold
    :param min: Minimum threshold value
    :param max: Maximum threshold value
    :param bg_val: All pixels outside [min,max] are set to bg_val
    :param fg_val: All pixels within [min,max] are set to fg_val
    """
    return _pj.Jim(jim_object._jipjim.pointOpThresh(min, max, fg_val, bg_val))


def stretch(jim_object,
            **kwargs):
    """Stretch pixel values

    :param jim_object: Jim object to be used for the stretch
    :param kwargs: See table below

    Modifies the instance on which the method was called.

    +---------------+--------------------------------------------------+
    | key           | value                                            |
    +==================+===============================================+
    | nodata        | do not consider these values                     |
    +---------------+--------------------------------------------------+
    | src_min       | clip source below this minimum value             |
    +---------------+--------------------------------------------------+
    | src_max       | clip source above this maximum value             |
    +---------------+--------------------------------------------------+
    | dst_min       | mininum value in output image                    |
    +---------------+--------------------------------------------------+
    | dst_max       | maximum value in output image                    |
    +---------------+--------------------------------------------------+
    | cc_min        | cumulative count cut from                        |
    +---------------+--------------------------------------------------+
    | cc_max        | cumulative count cut to                          |
    +---------------+--------------------------------------------------+
    | band          | band to stretch                                  |
    +---------------+--------------------------------------------------+
    | eq            | Histogram equalization                           |
    +---------------+--------------------------------------------------+
    | otype         | Output data type                                 |
    +---------------+--------------------------------------------------+

    Example:

    Convert data type of input image to byte with min 0 and max 255
    while stretching between cumulative counts 2 and 98 pct::

        jim = pj.Jim('/path/to/raster.tif')
        jim_stretched = pj.pixops.stretch(jim, otype='GDT_Byte', dst_min=0,
                                          dst_max=255, cc_min=2, cc_max=98)
    """

    return _pj.Jim(jim_object._jipjim.stretch(kwargs))


def supremum(jim,
             *args):
    """Create Jim composed using maximum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    if isinstance(jim, _pj.JimList):
        sup = None
        for newJim in jim:
            if sup is None:
                sup = _pj.Jim(jim[0])
            else:
                sup._jipjim.d_pointOpArith(newJim._jipjim, 5)
    else:
        sup = _pj.Jim(jim)
    for newJim in args:
        sup._jipjim.d_pointOpArith(newJim._jipjim, 5)

    return sup


class _PixOps(_pj.modules.JimModuleBase):
    """Define all PixOps methods."""

    def convert(self,
                otype):
        """Convert Jim image with respect to data type.

        :param otype: Data type for output image

        Modifies the instance on which the method was called.

        Example:

        Convert data type of input image to byte::

            jim0 = pj.Jim('/path/to/raster.tif')
            jim0.pixops.convert(otype='Byte')

        Clip raster dataset between 0 and 255 (set all other values to 0),
        then convert data type to byte::

            jim1 = pj.Jim('/path/to/raster.tif')
            jim1.setThreshold(min=0, max=255, nodata=0)
            jim1.pixops.convert('Byte')
        """
        if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
            otype = 'GDT_Byte'
        elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
            otype = 'GDT_UInt16'
        elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
            otype = 'GDT_Int16'
        elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
            otype = 'GDT_UInt32'
        elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
            otype = 'GDT_Int32'
        elif otype in [6, 'float32', 'Float32', 'GDT_Float32',
                       _jl.GDT_Float32]:
            otype = 'GDT_Float32'
        elif otype in [7, 'float64', 'Float64', 'GDT_Float64',
                       _jl.GDT_Float64]:
            otype = 'GDT_Float64'
        elif otype in ['int64', 'Int64', 'JDT_Int64']:
            otype = 'JDT_Int64'
        elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
            otype = 'JDT_UInt64'
        elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
            otype = 'JDT_UInt64'
        else:
            raise TypeError("Output type {} not supported".format(otype))
        # TODO: Support CTypes

        self._jim_object._set(
            self._jim_object._jipjim.convertDataType(otype))

    def histoCompress(self,
                      band: int = None):
        """Redistribute the intensity of histogram to fit full range of values.

        Redistributes the intensity values of im in such a way that the
        minimum value is zero and that all subsequent intensity values are
        used until the maximum value (i.e. no gaps).

        Modifies the instance on which the method was called.

        :param band: from which to get the histogram
        """
        if band is not None:
            self._jim_object._jipjim.d_histoCompress(band)
        else:
            self._jim_object._jipjim.d_histoCompress()

    def infimum(self,
                *args):
        """Change values of Jim using mimimum composition rule.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            if isinstance(jim, _pj.JimList):
                for newJim in jim:
                    self._jim_object._jipjim.d_pointOpArith(newJim._jipjim, 4)
            else:
                self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 4)

    # def modulo(self, val):
    #     """Set all pixels to their value modulo val

    #     :param val:  modulo value (integer)

    #      Modifies the instance on which the method was called.

    #     """
    #     self._jim_object._jipjim.d_pointOpModulo(val)

    def NDVI(self,
             band_red: int,
             band_nir: int):
        """Compute NDVI on the Jim object.

        Modifies the instance on which the method was called.

        :param band_red: index of band with values of red
        :param band_nir: index of band with values of NIR
        """
        red = _pj.geometry.cropBand(self._jim_object, band_red)
        nir = _pj.geometry.cropBand(self._jim_object, band_nir)

        self._jim_object._set(nir._jipjim.pointOpNDI(red._jipjim))

    def NDVISeparateBands(self,
                          jim_nir):
        """Compute NDVI from two Jims (call on red band, use NIR as param).

        Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

        Modifies the instance on which the method was called.

        :param jim_nir: Jim object with values of NIR
        """
        self._jim_object._set(
            jim_nir._jipjim.pointOpNDI(self._jim_object._jipjim))

    def setData(self,
                value: float,
                ulx: float = None,
                uly: float = None,
                lrx: float = None,
                lry: float = None,
                bands=None,
                dx: int = 0,
                dy: int = 0,
                nogeo: bool = False):
        """Set range of pixels to value.

        :param value: new value for pixels of Jim object
        :param ulx: upper left corner x coordinate (in projected coordinates
            if geo is True, else in image coordinates)
        :param uly: upper left corner y coordinate (in projected coordinates
            if geo is True, else in image coordinates)
        :param lrx: lower right corner x coordinate (in projected coordinates
            if geo is True, else in image coordinates)
        :param lry: lower right corner y coordinate (in projected coordinates
            if geo is True, else in image coordinates)
        :param bands: List of band indices to crop (index is 0 based)
        :param dx: spatial resolution in x to crop (stride if geo is False)
        :param dy: spatial resolution in y to crop (stride if geo is False)
        :param nogeo: use geospatial coordinates if False, image coordinates
            if True
        :return: a Jim object
        """
        if bands is None:
            bands = range(self._jim_object.properties.nrOfBand())

        if all(v is None for v in [ulx, uly, lrx, lry]) and dx == 0 and \
                dy == 0 and not nogeo:
            for band in bands:
                self._jim_object._jipjim.setData(value, band)
        else:
            if ulx is None:
                ulx = self._jim_object.properties.getUlx()
            if uly is None:
                uly = self._jim_object.properties.getUly()
            if lrx is None:
                lrx = self._jim_object.properties.getLrx()
            if lry is None:
                lry = self._jim_object.properties.getLry()
            for band in bands:
                self._jim_object._jipjim.setData(value, ulx, uly, lrx, lry,
                                                 band, dx, dy, nogeo)

    def setLevel(self,
                 min: float,
                 max: float,
                 val: float):
        """Set all pixels within min and max values to val.

        :param min:  Minimum threshold value
        :param max:  Maximum threshold value
        :param val:  All pixels within [min,max] are set to val

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_pointOpSetLevel(min, max, val)

    def setThreshold(self,
                     **kwargs):
        """Apply min and max threshold to pixel values in raster dataset.

        :param kwargs: See table below

        Modifies the instance on which the method was called.


        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | min              | Minimum threshold value (if pixel value < min    |
        |                  | set pixel value to no data)                      |
        +------------------+--------------------------------------------------+
        | max              | Maximum threshold value                          |
        |                  | (if pixel value < max set pixel value to no data)|
        +------------------+--------------------------------------------------+
        | value            | value to be set if within min and max            |
        |                  | (if not set, valid pixels will remain their input|
        |                  | value)                                           |
        +------------------+--------------------------------------------------+
        | abs              | Set to True to perform threshold test to absolute|
        |                  | pixel values                                     |
        +------------------+--------------------------------------------------+
        | nodata           | Set pixel value to this no data if pixel value   |
        |                  | < min or > max                                   |
        +------------------+--------------------------------------------------+

        .. note::

            A simplified interface to set a threshold is provided
            via :ref:`indexing <indexing>` (see also example below).

        .. _setThreshold_example:

        Example:

        Mask all values not within [0,250] and set to 255::

            jim0[(jim0<0) | (jim0>250)]=255

            Mask all values not within [0,250] and set to 255 (no data)::

            jim_threshold=jim.setThreshold(min=0,max=250,nodata=255)
        """
        self._jim_object._set(self._jim_object._jipjim.setThreshold(kwargs))

    def simpleArithOp(self,
                      jim,
                      op: int,
                      *args):
        """Change values of Jim using an arithmetic operation (coded by op).

        Modifies the instance on which the method was called.

        :param jim: Jim object (to be sure that at least one is provided)
        :param op: integer coding operation type (see table below)
        :param args: Jim objects
        """
        jims = [jim]
        jims.extend(args)
        for jim in jims:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, op)

    def simpleBitwiseOp(self,
                        jim,
                        op: int,
                        *args):
        """Change values of Jim using an bitwise operation (coded by op).

        Modifies the instance on which the method was called.

        :param jim: Jim object (to be sure that at least one is provided)
        :param op: integer coding operation type (see table below)
        :param args: Jim objects
        """
        jims = [jim]
        jims.extend(args)
        for jim in jims:
            self._jim_object._jipjim.d_pointOpBitwise(jim._jipjim, op)

    def simpleThreshold(self,
                        min: float,
                        max: float,
                        bg_val: float,
                        fg_val: float):
        """Set all pixels within min and max values to val.

        Modifies the instance on which the method was called.

        :param min:  Minimum threshold value
        :param max:  Maximum threshold value
        :param bg_val:  All pixels outside [min,max] are set to bg_val
        :param fg_val:  All pixels within [min,max] are set to fg_val
        """
        self._jim_object._jipjim.d_pointOpThresh(min, max, bg_val, fg_val)

    def stretch(self,
                **kwargs):
        """Stretch pixel values

        Modifies the instance on which the method was called.

        :param kwargs: See table below

        Modifies the instance on which the method was called.

        +---------------+--------------------------------------------------+
        | key           | value                                            |
        +==================+===============================================+
        | nodata        | do not consider these values                     |
        +---------------+--------------------------------------------------+
        | src_min       | clip source below this minimum value             |
        +---------------+--------------------------------------------------+
        | src_max       | clip source above this maximum value             |
        +---------------+--------------------------------------------------+
        | dst_min       | mininum value in output image                    |
        +---------------+--------------------------------------------------+
        | dst_max       | maximum value in output image                    |
        +---------------+--------------------------------------------------+
        | cc_min        | cumulative count cut from                        |
        +---------------+--------------------------------------------------+
        | cc_max        | cumulative count cut to                          |
        +---------------+--------------------------------------------------+
        | band          | band to stretch                                  |
        +---------------+--------------------------------------------------+
        | eq            | Histogram equalization                           |
        +---------------+--------------------------------------------------+
        | otype         | Output data type                                 |
        +---------------+--------------------------------------------------+

        Example:

        Convert data type of input image to byte with min 0 and max 255
        while stretching between cumulative counts 2 and 98 pct::

            jim = pj.Jim('/path/to/raster.tif')
            jim_stretched = pj.pixops.stretch(jim, otype='GDT_Byte', dst_min=0,
                                              dst_max=255, cc_min=2, cc_max=98)
        """

        self._jim_object._set(self._jim_object._jipjim.stretch(kwargs))

    def supremum(self,
                 *args):
        """Change values of Jim using maximum composition rule.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            if isinstance(jim, _pj.JimList):
                for newJim in jim:
                    self._jim_object._jipjim.d_pointOpArith(newJim._jipjim, 5)
            else:
                self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 5)


class _PixOpsList(_pj.modules.JimListModuleBase):
    """Define all PixOps methods for JimLists."""

    def composite(self,
                  crule: str = 'overwrite',
                  **kwargs):
        """Composite Jims in a JimList using a composite rule.

        :param crule: Rule for the composition
        :return: a composited Jim object
        """
        kwargs.update({'crule': crule})
        return _pj.Jim(self._jim_list._jipjimlist.composite(kwargs))

    def infimum(self):
        """Create Jim composed using minimum rule from Jim objects in JimList.

        :return: Jim composed of smallest values from provided Jim objects
        """
        inf = None
        for newJim in self._jim_list:
            if inf is None:
                inf = _pj.Jim(self._jim_list[0])
            else:
                inf._jipjim.d_pointOpArith(newJim._jipjim, 4)
        return inf

    def supremum(self):
        """Create Jim composed using maximum rule from Jim objects in JimList.

        :return: Jim composed of biggest values from provided Jim objects
        """
        sup = None
        for newJim in self._jim_list:
            if sup is None:
                sup = _pj.Jim(self._jim_list[0])
            else:
                sup._jipjim.d_pointOpArith(newJim._jipjim, 5)
        return sup


class _PixOpsVect(_pj.modules.JimVectModuleBase):
    """Define all PixOps methods for JimVects."""

    pass
