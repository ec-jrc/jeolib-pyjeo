"""Module for pixel-wise operations."""

import pyjeo as _pj
try:
    import jiplib as _jl
except ImportError:
    from jeodpp import jiplib as _jl

import numpy

def convert(jim_object, otype, **kwargs):
    """Convert Jim image with respect to data type.

    :param jim_object: a Jim object
    :param otype: Data type for output image
    :param kwargs: See table below
    :return: a Jim object

    +------------------+------------------------------------------------------+
    | key              | value                                                |
    +==================+======================================================+
    | scale            | Scale output: output=scale*input+offset              |
    +------------------+------------------------------------------------------+
    | offset           | Apply offset: output=scale*input+offset              |
    +------------------+------------------------------------------------------+
    | autoscale        | Scale output to min and max, e.g., [0,255]           |
    +------------------+------------------------------------------------------+
    | a_srs            | Override the projection for the output file          |
    +------------------+------------------------------------------------------+

    .. note::
        To ignore some pixels from the extraction process, see list
        of :ref:`mask <extract_mask>` key values:

    Example:

    Convert data type of input image to byte using autoscale::

        jim0=jl.io.createJim('/path/to/raster.tif')
        jim0.convert(Byte,autoscale=[0,255])

        Clip raster dataset between 0 and 255 (set all other values to 0), then
        convert data type to byte::

        jim1=jl.io.createJim('/path/to/raster.tif')
        jim1.setThreshold(min=0,max=255,nodata=0)
        jim1.convert('Byte')
    """

    #commented out for now
    # if len(kwargs) == 0:
    #     if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
    #         return _pj.Jim(jim_object._jipjim.convertToUchar8())
    #     elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
    #         return _pj.Jim(jim_object._jipjim.convertToUint16())
    #     elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
    #         return _pj.Jim(jim_object._jipjim.convertToUint32())
    #     elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
    #         return _pj.Jim(jim_object._jipjim.convertToFloat32())
    #     elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
    #         return _pj.Jim(jim_object._jipjim.convertToDouble64())

    if otype in [1, 'int8', 'uint8', 'UInt8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
        kwargs.update({'otype': 'GDT_Byte'})
    elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
        kwargs.update({'otype': 'GDT_UInt16'})
    elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
        kwargs.update({'otype': 'GDT_Int16'})
    elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
        kwargs.update({'otype': 'GDT_UInt32'})
    elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
        kwargs.update({'otype': 'GDT_Int32'})
    elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
        kwargs.update({'otype': 'GDT_Float32'})
    elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
        kwargs.update({'otype': 'GDT_Float64'})
    else:
        raise TypeError("Output type {} not supported".format(otype))
        # TODO: Support CTypes when bug fixed in jiplib

    return _pj.Jim(jim_object._jipjim.convert(kwargs))


def isEqual(first_jim, second_jim):
    if isinstance(second_jim, _pj.Jim) and isinstance(first_jim, _pj.Jim):
        return numpy.array_equal(first_jim.np(),second_jim.np())
        # return first_jim._jipjim.isEqual(second_jim._jipjim)
    else:
        return False


def histoCompress(jim_object, band=None):
    """Redistribute the intensity of histogram to fit full range of values.

    Redistributes the intensity values of im in such a way that the
    minimum value is zero and that all subsequent intensity values are
    used until the maximum value (i.e. no gaps).

    :return: a  Jim object
    """
    if band is not None:
        return _pj.Jim(jim_object._jipjim.histoCompress(band))
    else:
        return _pj.Jim(jim_object._jipjim.histoCompress())


def NDVI(jim_object, redBand, nirBand):
    """Compute NDVI on the Jim object.

    :param redBand: index of band with values of red
    :param nirBand: index of band with values of NIR
    :return: a Jim object with values of NDVI
    """
    red = _pj.geometry.cropBand(jim_object, redBand)
    nir = _pj.geometry.cropBand(jim_object, nirBand)

    return _pj.Jim(nir._jipjim.pointOpNDI(red._jipjim))


def NDVISeparateBands(redJim, nirJim):
    """Compute NDVI from two Jim objects.

    :param redJim: Jim object with values of red
    :param nirJim: Jim object with values of NIR
    :return: a Jim object with values of NDVI
    """
    return _pj.Jim(nirJim._jipjim.pointOpNDI(redJim._jipjim))


def modulo(jim_object, val):
    """Set all pixels to their value modulo val

    :param val:  modulo value (integer)

    Modifies the instance on which the method was called.

    """
    return _pj.Jim(jim_object._jipjim.pointOpModulo(val))


def blank(jim_object, val):
    """Set all pixels to val

    :param val:  All pixels within [min,max] are set to val

    Modifies the instance on which the method was called.

    """
    return _pj.Jim(jim_object._jipjim.pointOpBlank(val))


def setLevel(jim_object, min, max, val):
    """Set all pixels within min and max values to val

    :param min:  Minimum threshold value       
    :param max:  Maximum threshold value
    :param val:  All pixels within [min,max] are set to val

    Modifies the instance on which the method was called.

    """
    return _pj.Jim(jim_object._jipjim.pointOpSetLevel(min, max, val))


def simpleThreshold(jim_object, min, max, bg_val, fg_val):
    """Set all pixels within min and max values to val

    :param min:  Minimum threshold value       
    :param max:  Maximum threshold value
    :param bg_val:  All pixels outside [min,max] are set to bg_val
    :param fg_val:  All pixels within [min,max] are set to fg_val

    Modifies the instance on which the method was called.

    """
    return _pj.Jim(jim_object._jipjim.pointOpThresh(min, max, fg_val, bg_val))


def setThreshold(jim_object, **kwargs):
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

        A simplified interface to set a threshold is provided via the index
        operator [] (see :ref:`example <setitem_example>` below).

    .. _setitem_example:

    Example:

    Mask all values not within [0,250] and set to 255::

        jim0[(jim0<0) | (jim0>250)]=255

        Mask all values not within [0,250] and set to 255 (no data)::

        jim_threshold=jim.setThreshold(min=0,max=250,nodata=255)
    """
    return _pj.Jim(jim_object._jipjim.setThreshold(kwargs))


def setData(jim, value, ulx=None, uly=None, lrx=None, lry=None, bands=[0],
            dx=0, dy=0, nogeo=False):
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
    jout.pixops.setData(value)
    return _pj.Jim(jout)


def simpleArithOp(jim, op, *args):
    """Create Jim composed using a simple arithmetic operation (coded with op) from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param op: integer for operation type
    :param args: Jim objects
    :return: Jim holding specified arithmetic operation with from provided Jim objects
    """
    jout = _pj.Jim(jim)
    for newJim in args:
        jout._jipjim.d_pointOpArith(newJim._jipjim, op)

    return _pj.Jim(jout)


def simpleBitwiseOp(jim, op, *args):
    """Create Jim composed using a simple bitwise operation (coded with op) from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param op: integer for operation type
    :param args: Jim objects
    :return: Jim holding specified bitwise operation with from provided Jim objects
    """
    jout = _pj.Jim(jim)
    for newJim in args:
        jout._jipjim.d_pointOpBitwise(newJim._jipjim, op)

    return _pj.Jim(jout)


def supremum(jim, *args):
    """Create Jim composed using maximum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    supremum = _pj.Jim(jim)
    for newJim in args:
        supremum._jipjim.d_pointOpArith(newJim._jipjim, 5)

    return _pj.Jim(supremum)


def infimum(jim, *args):
    """Create Jim composed using minimum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    infimum = _pj.Jim(jim)
    for newJim in args:
        infimum._jipjim.d_pointOpArith(newJim._jipjim, 4)

    return _pj.Jim(infimum)


def composite(jim_list, crule='overwrite', **kwargs):
    kwargs.update({'crule': crule})
    return _pj.Jim(jim_list._jipjimlist.composite(kwargs))


class _PixOps():
    """Define all PixOps methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def convert(self, otype, **kwargs):
        """Convert Jim image with respect to data type.

        :param otype: Data type for output image
        :param kwargs: See table below

        Modifies the instance on which the method was called.

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | scale            | Scale output: output=scale*input+offset          |
        +------------------+--------------------------------------------------+
        | offset           | Apply offset: output=scale*input+offset          |
        +------------------+--------------------------------------------------+
        | autoscale        | Scale output to min and max, e.g., [0,255]       |
        +------------------+--------------------------------------------------+
        | a_srs            | Override the projection for the output file      |
        +------------------+--------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list
            of :ref:`mask <extract_mask>` key values:

        Example:

        Convert data type of input image to byte using autoscale::

            jim0=jl.io.createJim('/path/to/raster.tif')
            jim0.convert(otype=Byte,autoscale=[0,255])

            Clip raster dataset between 0 and 255 (set all other values to 0),
            then convert data type to byte::

            jim1=jl.io.createJim('/path/to/raster.tif')
            jim1.setThreshold(min=0,max=255,nodata=0)
            jim1.convert(Byte)
        """

        #commented out for now
        # if len(kwargs) == 0:
        #     if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
        #         self._jim_object._jipjim.d_convertToUchar8()
        #         return None
        #     elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
        #         self._jim_object._set(
        #             self._jim_object._jipjim.convertToUint16())
        #         return None
        #     elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
        #         self._jim_object._set(
        #             self._jim_object._jipjim.convertToUint32())
        #         return None
        #     elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
        #         self._jim_object._set(
        #             self._jim_object._jipjim.convertToFloat32())
        #         return None
        #     elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
        #         self._jim_object._set(
        #             self._jim_object._jipjim.convertToDouble64())
        #         return None

        if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
            kwargs.update({'otype': 'GDT_Byte'})
        elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
            kwargs.update({'otype': 'GDT_UInt16'})
        elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
            kwargs.update({'otype': 'GDT_Int16'})
        elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
            kwargs.update({'otype': 'GDT_UInt32'})
        elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
            kwargs.update({'otype': 'GDT_Int32'})
        elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
            kwargs.update({'otype': 'GDT_Float32'})
        elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
            kwargs.update({'otype': 'GDT_Float64'})
        else:
            raise TypeError("Output type {} not supported".format(otype))
        # TODO: Support CTypes when bug fixed in jiplib

        self._jim_object._set(self._jim_object._jipjim.convert(kwargs))

    def histoCompress(self, band=None):
        """Redistribute the intensity of histogram to fit full range of values.

        Redistributes the intensity values of im in such a way that the
        minimum value is zero and that all subsequent intensity values are
        used until the maximum value (i.e. no gaps).

        Modifies the instance on which the method was called.
        """
        if band is not None:
            self._jim_object._jipjim.d_histoCompress(band)
        else:
            self._jim_object._jipjim.d_histoCompress()

    def isEqual(self, other):
        if isinstance(other, _pj.Jim):
            # return self._jim_object._jipjim.isEqual(other._jipjim)
            return numpy.array_equal(self._jim_object.np(),other.np())
        else:
            return False

    def NDVI(self, redBand, nirBand):
        """Compute NDVI on the Jim object.

        Modifies the instance on which the method was called.

        :param redBand: index of band with values of red
        :param nirBand: index of band with values of NIR
        """
        red = _pj.geometry.cropBand(self._jim_object, redBand)
        nir = _pj.geometry.cropBand(self._jim_object, nirBand)

        self._jim_object._set(nir._jipjim.pointOpNDI(red._jipjim))

    def NDVISeparateBands(self, nirJim):
        """Compute NDVI from two Jims (call on red band, use NIR as param).

        Modifies the instance on which the method was called.

        :param nirJim: Jim object with values of NIR
        """
        self._jim_object._set(
            nirJim._jipjim.pointOpNDI(self._jim_object._jipjim))


    def modulo(self, val):
        """Set all pixels to their value modulo val

        :param val:  modulo value (integer)

         Modifies the instance on which the method was called.

        """
        self._jim_object._jipjim.d_pointOpModulo(val)


    def blank(self, val):
        """Set all pixels to val

        :param val:  All pixels are set to val

         Modifies the instance on which the method was called.

        """
        self._jim_object._jipjim.d_pointOpBlank(val)

    def setData(self, value, ulx=None, uly=None, lrx=None, lry=None, bands=[0],
                dx=0, dy=0, nogeo=False):
        """Set range of pixels to value.

        :param jim_object: a Jim object
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
        if all(v is None for v in [ulx, uly, lrx, lry]) and dx==0 and dy==0 \
                and not nogeo:
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



    def simpleThreshold(self, min, max, bg_val, fg_val):
        """Set all pixels within min and max values to val

        :param min:  Minimum threshold value       
        :param max:  Maximum threshold value
        :param bg_val:  All pixels outside [min,max] are set to bg_val
        :param fg_val:  All pixels within [min,max] are set to fg_val

        Modifies the instance on which the method was called.

        """
        self._jim_object._jipjim.d_pointOpThresh(min, max, bg_val, fg_val)



    def setThreshold(self, **kwargs):
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

            A simplified interface to set a threshold is provided via the index
            operator [] (see :ref:`example <setitem_example>` below).

        .. _setitem_example:

        Example:

        Mask all values not within [0,250] and set to 255::

            jim0[(jim0<0) | (jim0>250)]=255

            Mask all values not within [0,250] and set to 255 (no data)::

            jim_threshold=jim.setThreshold(min=0,max=250,nodata=255)
        """
        self._jim_object._set(self._jim_object._jipjim.setThreshold(kwargs))


    def simpleArithOp(self, op, *args):
        """Change values of Jim using an arithmetic operation (coded by op) from provided Jim objects.

        Modifies the instance on which the method was called.

        :param op: integer coding operation type (see table below)
        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, op)
            
    def simpleBitwiseOp(self, op, *args):
        """Change values of Jim using an bitwise operation (coded by op) from provided Jim objects.

        Modifies the instance on which the method was called.

        :param op: integer coding operation type (see table below)
        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpBitwise(jim._jipjim, op)
    
    def supremum(self, *args):
        """Change values of Jim using maximum composition rule from provided Jim objects.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 5)

    def infimum(self, *args):
        """Change values of Jim using mimimum composition rule from provided Jim objects.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 4)


class _PixOpsList():
    """Define all PixOps methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def composite(self, crule='overwrite', **kwargs):
        kwargs.update({'crule': crule})
        return _pj.Jim(self._jim_list._jipjimlist.composite(kwargs))


class _PixOpsVect():
    """Define all PixOps methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
