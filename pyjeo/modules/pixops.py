"""Module for pixel-wise operations."""

import numpy

import pyjeo as _pj
import jiplib as _jl
from . import JimModuleBase as _JimModuleBase
from . import JimListModuleBase as _JimListModuleBase
from . import JimVectModuleBase as _JimVectModuleBase


# def blank(jim_object, val):
#     """Set all pixels to val.

#     :param val:  All pixels within [min,max] are set to val

#     Modifies the instance on which the method was called.

#     """
#     return _pj.Jim(jim_object._jipjim.pointOpBlank(val))


def composite(jim_list, crule='overwrite', **kwargs):
    """Composite Jims in a JimList using a composite rule."""
    kwargs.update({'crule': crule})
    return _pj.Jim(jim_list._jipjimlist.composite(kwargs))


def convert(jim_object, otype, **kwargs):
    """Convert Jim image with respect to data type.

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
    if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
        otype = 'GDT_Byte'
        nptype = numpy.int8
    elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
        otype = 'GDT_UInt16'
        nptype = numpy.uint16
    elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
        otype = 'GDT_Int16'
        nptype = numpy.int16
    elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
        otype = 'GDT_UInt32'
        nptype = numpy.uint32
    elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
        otype = 'GDT_Int32'
        nptype = numpy.int32
    elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
        otype = 'GDT_Float32'
        nptype = numpy.float32
    elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
        otype = 'GDT_Float64'
        nptype = numpy.float64
    elif otype in ['int64', 'Int64', 'JDT_Int64']:
        otype = 'JDT_Int64'
        nptype = numpy.int64
    elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
        otype = 'JDT_UInt64'
        nptype = numpy.uint64
    elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
        otype = 'JDT_UInt64'
        nptype = numpy.uint64
    else:
        raise TypeError("Output type {} not supported".format(otype))
    # TODO: Support CTypes

    if len(kwargs) and jim_object.properties.nrOfPlane() == 1:
        kwargs.update({'otype': otype})
        return _pj.Jim(jim_object._jipjim.convert(kwargs))
    else:
        jimnew = _pj.Jim(ncol=jim_object.properties.nrOfCol(),
                         nrow=jim_object.properties.nrOfRow(),
                         nband=jim_object.properties.nrOfBand(),
                         nplane=jim_object.properties.nrOfPlane(),
                         otype=otype)
        jimnew.properties.setProjection(jim_object.properties.getProjection())
        jimnew.properties.setGeoTransform(
            jim_object.properties.getGeoTransform())
        for iband in range(0, jim_object.properties.nrOfBand()):
            jimnew.np(iband)[:] = jim_object.np(iband).astype(nptype)
        return jimnew

# def convert(jim_object, otype, **kwargs):
#     """Convert Jim image with respect to data type.

#     :param jim_object: a Jim object
#     :param otype: Data type for output image
#     :param kwargs: See table below
#     :return: a Jim object

#     +------------------+--------------------------------------------------+
#     | key              | value                                            |
#     +==================+==================================================+
#     | scale            | Scale output: output=scale*input+offset          |
#     +------------------+--------------------------------------------------+
#     | offset           | Apply offset: output=scale*input+offset          |
#     +------------------+--------------------------------------------------+
#     | autoscale        | Scale output to min and max, e.g., [0,255]       |
#     +------------------+--------------------------------------------------+
#     | a_srs            | Override the projection for the output file      |
#     +------------------+--------------------------------------------------+

#     Example:

#     Convert data type of input image to byte using autoscale::

#         jim0=jl.io.createJim('/path/to/raster.tif')
#         jim0.convert(Byte,autoscale=[0,255])

#         Clip raster dataset between 0 and 255 (set all other values to 0),
#         then convert data type to byte::

#         jim1=jl.io.createJim('/path/to/raster.tif')
#         jim1.setThreshold(min=0,max=255,nodata=0)
#         jim1.convert('Byte')
#     """
#     if otype in [1, 'int8', 'uint8', 'UInt8', 'Byte', 'GDT_Byte',
#                  _jl.GDT_Byte]:
#         kwargs.update({'otype': 'GDT_Byte'})
#     elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
#         kwargs.update({'otype': 'GDT_UInt16'})
#     elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
#         kwargs.update({'otype': 'GDT_Int16'})
#     elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
#         kwargs.update({'otype': 'GDT_UInt32'})
#     elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
#         kwargs.update({'otype': 'GDT_Int32'})
#     elif otype in [6, 'float32', 'Float32', 'GDT_Float32', _jl.GDT_Float32]:
#         kwargs.update({'otype': 'GDT_Float32'})
#     elif otype in [7, 'float64', 'Float64', 'GDT_Float64', _jl.GDT_Float64]:
#         kwargs.update({'otype': 'GDT_Float64'})
#     else:
#         raise TypeError("Output type {} not supported".format(otype))
#         # TODO: Support CTypes when bug fixed in jiplib

#     return _pj.Jim(jim_object._jipjim.convert(kwargs))


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


def infimum(jim, *args):
    """Create Jim composed using minimum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    infimum = _pj.Jim(jim)
    for newJim in args:
        infimum._jipjim.d_pointOpArith(newJim._jipjim, 4)

    return infimum


def isEqual(first_jim, second_jim):
    """Check if the values of one Jim object are the same as in another one.

    :param first_jim: a Jim object
    :param second_jim: a Jim object
    :return: True if the values are equal, zero otherwise
    """
    if isinstance(second_jim, _pj.Jim) and isinstance(first_jim, _pj.Jim):
        if first_jim.properties.nrOfPlane() != \
                second_jim.properties.nrOfPlane() or \
                first_jim.properties.nrOfBand() != \
                second_jim.properties.nrOfBand():
            return False
        if first_jim.properties.nrOfPlane() == 1:
            for iband in range(0, first_jim.properties.nrOfBand()):
                if not numpy.array_equal(first_jim.np(iband),
                                         second_jim.np(iband)):
                    return False
            return True
        else:
            for iplane in range(0, first_jim.properties.nrOfPlane()):
                first_plane = _pj.geometry.cropPlane(first_jim, iplane)
                second_plane = _pj.geometry.cropPlane(second_jim, iplane)
                if first_plane.properties.nrOfBand() != \
                        second_plane.properties.nrOfBand():
                    return False
                for iband in range(0, first_plane.properties.nrOfBand()):
                    if not numpy.array_equal(first_plane.np(iband),
                                             second_plane.np(iband)):
                        return False
            return True
    else:
        return False


# def modulo(jim_object, val):
#     """Set all pixels to their value modulo val.

#     :param val:  modulo value (integer)

#     Modifies the instance on which the method was called.

#     """
#     return _pj.Jim(jim_object._jipjim.pointOpModulo(val))


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

    Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

    :param redJim: Jim object with values of red
    :param nirJim: Jim object with values of NIR
    :return: a Jim object with values of NDVI
    """
    return _pj.Jim(nirJim._jipjim.pointOpNDI(redJim._jipjim))


def setData(jim, value, ulx=None, uly=None, lrx=None, lry=None, bands=None,
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

    if bands is None:
        bands = range(jim.properties.nrOfBand())

    jout.pixops.setData(value, ulx, uly, lrx, lry, bands, dx, dy, nogeo)
    return jout


def setLevel(jim_object, min, max, val):
    """Set all pixels within min and max values to val.

    :param min:  Minimum threshold value
    :param max:  Maximum threshold value
    :param val:  All pixels within [min,max] are set to val
    """
    return _pj.Jim(jim_object._jipjim.pointOpSetLevel(min, max, val))


def setThreshold(jim_object, **kwargs):
    """Apply min and max threshold to pixel values in raster dataset.

    :jim_object: the Jim object on which to set threshold
    :param kwargs: See table :py:meth:`~pixops._PixOps.setThreshold`.

    for help, please refer to the corresponding
    method :py:meth:`~pixops._PixOps.setThreshold`.
    """
    return _pj.Jim(jim_object._jipjim.setThreshold(kwargs))


def simpleArithOp(jim1, jim2, op, *args):
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


def simpleBitwiseOp(jim, another_jim, op, *args):
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


def simpleThreshold(jim_object, min, max, bg_val, fg_val):
    """Set all pixels within min and max values to val.

    :param min:  Minimum threshold value
    :param max:  Maximum threshold value
    :param bg_val:  All pixels outside [min,max] are set to bg_val
    :param fg_val:  All pixels within [min,max] are set to fg_val
    """
    return _pj.Jim(jim_object._jipjim.pointOpThresh(min, max, fg_val, bg_val))


def supremum(jim, *args):
    """Create Jim composed using maximum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of biggest values from provided Jim objects
    """
    supremum = _pj.Jim(jim)
    for newJim in args:
        supremum._jipjim.d_pointOpArith(newJim._jipjim, 5)

    return supremum


class _PixOps(_JimModuleBase):
    """Define all PixOps methods."""

    # def blank(self, val):
    #     """Set all pixels to val

    #     :param val:  All pixels are set to val

    #      Modifies the instance on which the method was called.

    #     """
    #     self._jim_object._jipjim.d_pointOpBlank(val)

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
        if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
            kwargs.update({'otype': 'GDT_Byte'})
            otype = 'GDT_Byte'
            nptype = numpy.int8
        elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
            otype = 'GDT_UInt16'
            nptype = numpy.uint16
        elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
            otype = 'GDT_Int16'
            nptype = numpy.int16
        elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
            otype = 'GDT_UInt32'
            nptype = numpy.uint32
        elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
            otype = 'GDT_Int32'
            nptype = numpy.int32
        elif otype in [6, 'float32', 'Float32', 'GDT_Float32',
                       _jl.GDT_Float32]:
            otype = 'GDT_Float32'
            nptype = numpy.float32
        elif otype in [7, 'float64', 'Float64', 'GDT_Float64',
                       _jl.GDT_Float64]:
            otype = 'GDT_Float64'
            nptype = numpy.float64
        elif otype in ['int64', 'Int64', 'JDT_Int64']:
            otype = 'JDT_Int64'
            nptype = numpy.int64
        elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
            otype = 'JDT_UInt64'
            nptype = numpy.uint64
        elif otype in ['uint64', 'UInt64', 'JDT_UInt64']:
            otype = 'JDT_UInt64'
            nptype = numpy.uint64
        else:
            raise TypeError("Output type {} not supported".format(otype))
        # TODO: Support CTypes
        if len(kwargs) and self._jim_object.properties.nrOfPlane() == 1:
            kwargs.update({'otype': otype})
            self._jim_object._set(self._jim_object._jipjim.convert(kwargs))
        else:
            jimnew = _pj.Jim(ncol=self._jim_object.properties.nrOfCol(),
                             nrow=self._jim_object.properties.nrOfRow(),
                             nband=self._jim_object.properties.nrOfBand(),
                             nplane=self._jim_object.properties.nrOfPlane(),
                             otype=otype)
            jimnew.properties.setProjection(
                self._jim_object.properties.getProjection())
            jimnew.properties.setGeoTransform(
                self._jim_object.properties.getGeoTransform())
            for iband in range(0, self._jim_object.properties.nrOfBand()):
                jimnew.np(iband)[:] = self._jim_object.np(iband).astype(nptype)
            self._jim_object._set(jimnew._jipjim)

    # def convert(self, otype, **kwargs):
    #     """Convert Jim image with respect to data type.

    #     :param otype: Data type for output image
    #     :param kwargs: See table below

    #     Modifies the instance on which the method was called.

    #     +------------------+----------------------------------------------+
    #     | key              | value                                        |
    #     +==================+==============================================+
    #     | scale            | Scale output: output=scale*input+offset      |
    #     +------------------+----------------------------------------------+
    #     | offset           | Apply offset: output=scale*input+offset      |
    #     +------------------+----------------------------------------------+
    #     | autoscale        | Scale output to min and max, e.g., [0,255]   |
    #     +------------------+----------------------------------------------+
    #     | a_srs            | Override the projection for the output file  |
    #     +------------------+----------------------------------------------+

    #     Example:

    #     Convert data type of input image to byte using autoscale::

    #         jim0=jl.io.createJim('/path/to/raster.tif')
    #         jim0.convert(otype=Byte,autoscale=[0,255])

    #         Clip raster dataset between 0 and 255 (set all other values to
    #         0), then convert data type to byte::

    #         jim1=jl.io.createJim('/path/to/raster.tif')
    #         jim1.setThreshold(min=0,max=255,nodata=0)
    #         jim1.convert(Byte)
    #     """
    #     if otype in [1, 'int8', 'uint8', 'Byte', 'GDT_Byte', _jl.GDT_Byte]:
    #         kwargs.update({'otype': 'GDT_Byte'})
    #     elif otype in [2, 'uint16', 'UInt16', 'GDT_UInt16', _jl.GDT_UInt16]:
    #         kwargs.update({'otype': 'GDT_UInt16'})
    #     elif otype in [3, 'int16', 'Int16', 'GDT_Int16', _jl.GDT_Int16]:
    #         kwargs.update({'otype': 'GDT_Int16'})
    #     elif otype in [4, 'uint32', 'UInt32', 'GDT_UInt32', _jl.GDT_UInt32]:
    #         kwargs.update({'otype': 'GDT_UInt32'})
    #     elif otype in [5, 'int32', 'Int32', 'GDT_Int32', _jl.GDT_Int32]:
    #         kwargs.update({'otype': 'GDT_Int32'})
    #     elif otype in [6, 'float32', 'Float32', 'GDT_Float32',
    #                    _jl.GDT_Float32]:
    #         kwargs.update({'otype': 'GDT_Float32'})
    #     elif otype in [7, 'float64', 'Float64', 'GDT_Float64',
    #                    _jl.GDT_Float64]:
    #         kwargs.update({'otype': 'GDT_Float64'})
    #     else:
    #         raise TypeError("Output type {} not supported".format(otype))
    #     # TODO: Support CTypes when bug fixed in jiplib

    #     self._jim_object._set(self._jim_object._jipjim.convert(kwargs))

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

    def infimum(self, *args):
        """Change values of Jim using mimimum composition rule.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 4)

    def isEqual(self, other):
        """Check if the values of one Jim object are the same as in another.

        :param other: a Jim object
        :return: True if the values are equal, zero otherwise
        """
        if isinstance(other, _pj.Jim):
            if self._jim_object.properties.nrOfPlane() != \
                    other.properties.nrOfPlane() or \
                    self._jim_object.properties.nrOfBand() != \
                    other.properties.nrOfBand():
                return False
            if self._jim_object.properties.nrOfPlane() == 1:
                for iband in range(0, self._jim_object.properties.nrOfBand()):
                    if not numpy.array_equal(self._jim_object.np(iband),
                                             other.np(iband)):
                        return False
                return True
            else:
                for iplane in range(0,
                                    self._jim_object.properties.nrOfPlane()):
                    first_plane = _pj.geometry.cropPlane(self._jim_object,
                                                         iplane)
                    second_plane = _pj.geometry.cropPlane(other, iplane)
                    if first_plane.properties.nrOfBand() != \
                            second_plane.properties.nrOfBand():
                        return False
                    for iband in range(0, first_plane.properties.nrOfBand()):
                        if not numpy.array_equal(first_plane.np(iband),
                                                 second_plane.np(iband)):
                            return False
                return True
        else:
            return False

    # def modulo(self, val):
    #     """Set all pixels to their value modulo val

    #     :param val:  modulo value (integer)

    #      Modifies the instance on which the method was called.

    #     """
    #     self._jim_object._jipjim.d_pointOpModulo(val)

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

        Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

        Modifies the instance on which the method was called.

        :param nirJim: Jim object with values of NIR
        """
        self._jim_object._set(
            nirJim._jipjim.pointOpNDI(self._jim_object._jipjim))

    def setData(self, value, ulx=None, uly=None, lrx=None, lry=None,
                bands=None, dx=0, dy=0, nogeo=False):
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

    def setLevel(self, min, max, val):
        """Set all pixels within min and max values to val.

        :param min:  Minimum threshold value
        :param max:  Maximum threshold value
        :param val:  All pixels within [min,max] are set to val

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_pointOpSetLevel(min, max, val)

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

    def simpleArithOp(self, jim, op, *args):
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

    def simpleBitwiseOp(self, jim, op, *args):
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

    def simpleThreshold(self, min, max, bg_val, fg_val):
        """Set all pixels within min and max values to val.

        Modifies the instance on which the method was called.

        :param min:  Minimum threshold value
        :param max:  Maximum threshold value
        :param bg_val:  All pixels outside [min,max] are set to bg_val
        :param fg_val:  All pixels within [min,max] are set to fg_val
        """
        self._jim_object._jipjim.d_pointOpThresh(min, max, bg_val, fg_val)

    def supremum(self, *args):
        """Change values of Jim using maximum composition rule.

        Modifies the instance on which the method was called.

        :param args: Jim objects
        """
        for jim in args:
            self._jim_object._jipjim.d_pointOpArith(jim._jipjim, 5)


class _PixOpsList(_JimListModuleBase):
    """Define all PixOps methods for JimLists."""

    def composite(self, crule='overwrite', **kwargs):
        """Composite Jims in a JimList using a composite rule."""
        kwargs.update({'crule': crule})
        return _pj.Jim(self._jim_list._jipjimlist.composite(kwargs))


class _PixOpsVect(_JimVectModuleBase):
    """Define all PixOps methods for JimVects."""

    pass
