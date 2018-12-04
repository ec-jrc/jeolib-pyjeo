"""Module for pixel-wise operations."""

try:
    import pyjeo as _pj
except ImportError:
    try:
        from jeodpp import pyjeo as _pj
    except ImportError:
        import jeodpp.pyjeo as _pj


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

    if len(kwargs) == 0:
        if otype in ['Byte', 'GDT_Byte', _jl.GDT_Byte]:
            return _pj.Jim(jim_object.convertToUchar8)
        elif otype in ['UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
            return _pj.Jim(jim_object.convertToUint16)
        elif otype in ['UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
            return _pj.Jim(jim_object.convertToUint32)
        elif otype in ['Float32', 'GDT_Float32', _jl.GDT_Float32]:
            return _pj.Jim(jim_object.convertToFloat32())
        elif otype in ['Float64', 'GDT_Float64', _jl.GDT_Float64]:
            return _pj.Jim(jim_object.convertToFloat64())
    if otype in ['Byte', 'GDT_Byte', _jl.GDT_Byte]:
        kwargs.update({'otype': 'GDT_Byte'})
    elif otype in ['UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
        kwargs.update({'otype': 'GDT_UInt16'})
    elif otype in ['Int16', 'GDT_Int16', _jl.GDT_Int16]:
        kwargs.update({'otype': 'GDT_Int16'})
    elif otype in ['UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
        kwargs.update({'otype': 'GDT_UInt32'})
    elif otype in ['Int32', 'GDT_Int32', _jl.GDT_Int32]:
        kwargs.update({'otype': 'GDT_Int32'})
    elif otype in ['Float32', 'GDT_Float32', _jl.GDT_Float32]:
        kwargs.update({'otype': 'GDT_Float32'})
    elif otype in ['Float64', 'GDT_Float64', _jl.GDT_Float64]:
        kwargs.update({'otype': 'GDT_Float32'})
    else:
        print("Warning: output type {} not supported".format(otype))
    kwargs.update({'otype': otype})
    # retJim = _pj.io.createJim(jim_object)
    # return _pj.Jim(retJim.convert(kwargs))
    return _pj.Jim(jim_object.convert(kwargs))


def isEqual(first_jim, second_jim):
    if isinstance(other, _pj.Jim):
        return first_jim.isEqual(second_jim)
    else:
        return False


def NDVI(jim_object, redBand, nirBand):
    """Compute NDVI on the Jim object.

    :param redBand: index of band with values of red
    :param nirBand: index of band with values of NIR
    :return: a Jim object with values of NDVI
    """
    red = _pj.geometry.cropBand(jim_object, redBand)
    nir = _pj.geometry.cropBand(jim_object, nirBand)

    return _pj.Jim(nir.pointOpNDI(red))


def NDVISeparateBands(redJim, nirJim):
    """Compute NDVI from two Jim objects.

    :param redJim: Jim object with values of red
    :param nirJim: Jim object with values of NIR
    :return: a Jim object with values of NDVI
    """
    return _pj.Jim(nirJim.pointOpNDI(redJim))


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
    return _pj.Jim(jim_object.setThreshold(kwargs))


def supremum(jim, *args):
    """Create Jim composed of biggest values from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    supremum = _pj.Jim(jim)
    for newJim in args:
        supremum = supremum.pointOpArith(newJim, 5)

    return _pj.Jim(supremum)


def composite(jim_list, **kwargs):
    return _pj.Jim(jim_list.composite(kwargs))


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
        kwargs.update({'otype': otype})
        self._jim_object._set(self._jim_object.convert(kwargs))

    def isEqual(self, other):
        if isinstance(other, _pj.Jim):
            return self._jim_object.isEqual(other)
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

        self._jim_object._set(nir.pointOpNDI(red))

    def NDVISeparateBands(self, nirJim):
        """Compute NDVI from two Jims (call on red band, use NIR as param).

        Modifies the instance on which the method was called.

        :param nirJim: Jim object with values of NIR
        """
        self._jim_object._set(nirJim.pointOpNDI(self._jim_object))

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
        if all(v is None for v in [ulx, uly, lrx, lry]) and bands==[0] and \
                dx==0 and dy==0 and not nogeo:
            for band in bands:
                self._jim_object.setData(value)
        else:
            for band in bands:
                self._jim_object.setData(value, ulx, uly, lrx, lry, band, dx,
                                         dy, nogeo)

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
        self._jim_object._set(self._jim_object.setThreshold(kwargs))

    def supremum(self, *args):
        """Change values of Jim for the biggest ones from provided Jim objects.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            self._jim_object.d_pointOpArith(jim, 5)


class _PixOpsList():
    """Define all PixOps methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def composite(self, **kwargs):
        return _pj.Jim(self._jim_list.composite(kwargs))


class _PixOpsVect():
    """Define all PixOps methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
