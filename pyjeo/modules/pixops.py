"""Module for pixel-wise operations."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2023 European Union (Joint Research Centre)
#
# This file is part of pyjeo.
#
# pyjeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyjeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.

import numpy as _np

import pyjeo as _pj
import jiplib as _jl


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
        raise _pj.exceptions.JimIllegalArgumentError(
            "Output type {} not supported".format(otype))
    # TODO: Support CTypes

    ret_jim = _pj.Jim(jim_object._jipjim.convertDataType(otype))
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


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
        ret_jim = _pj.Jim(jim_object._jipjim.histoCompress(band))
    else:
        ret_jim = _pj.Jim(jim_object._jipjim.histoCompress())

    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim

def infimum(jim_object,
            *args):
    """Create Jim composed using minimum rule from provided Jim objects.

    :param jim: Jim object (to be sure that at least one is provided)
    :param args: Jim objects
    :return: Jim composed of smalles values from provided Jim objects
    """
    if isinstance(jim_object, _pj.JimList):
        ret_jim = None
        for newJim in jim_object:
            if ret_jim is None:
                ret_jim = _pj.Jim(jim_object[0])
            else:
                ret_jim._jipjim.d_pointOpArith(newJim._jipjim, 4)
    else:
        ret_jim = _pj.Jim(jim_object)
    for newJim in args:
        ret_jim._jipjim.d_pointOpArith(newJim._jipjim, 4)

    return ret_jim


# def modulo(jim_object, val):
#     """Set all pixels to their value modulo val.

#     :param val:  modulo value (integer)

#     Modifies the instance on which the method was called.

#     """
#     return _pj.Jim(jim_object._jipjim.pointOpModulo(val))


def NDVI(jim_object,
         red,
         nir,
         name: str = None,
         nodata=-2,
         scale: float = 1.0,
         offset: float = 0,
         addBand: bool = False):
        """Compute NDVI on the Jim object.

        :param jim_object: Jim object from which the red and NIR bands are
            to be derived
        :param red: index of band with values of red
        :param nir: index of band with values of NIR
        :param name: name of band for NDVI
        :param scale: scale to multiply NDVI
        :param offset: offset to add NDVI
        :param addBand: add NDVI as a new band
        """
        if jim_object.properties.getDataType() not in (
                'Float32','Float64'):
            if scale <= 1.0:
                raise UserWarning("Warning data type is not floating point\
                , use scale and offset or change data type")
        nptype = jim_object.np().dtype.type
        if isinstance(red, int):
            redindex = red
        else:
            redindex = jim_object.properties.getDimension(
                'band').index(red)
        if isinstance(nir, int):
            nirindex = nir
        else:
            nirindex = jim_object.properties.getDimension(
                'band').index(nir)
        if addBand:
            ret_jim = _pj.Jim(jim_object)
            rednp = ret_jim.np(redindex).astype(float)
            nirnp = ret_jim.np(nirindex).astype(float)
            ndvi = (nirnp - rednp) / (rednp + nirnp)
            if offset != 0:
                ndvi += offset
            if scale != 1:
                ndvi *= scale
            ndvi[_np.isnan(ndvi)] = nodata
            ndvi = ndvi.astype(nptype)
            jimndvi = _pj.np2jim(ndvi)
            jimndvi.properties.copyGeoReference(ret_jim)
            if name is not None:
                jimndvi.properties.setDimension(name, 'band')
            ret_jim.geometry.stackBand(jimndvi)
        else:
            ret_jim = _pj.geometry.cropBand(jim_object, [red, nir])
            ret_jim.pixops.convert('GDT_Float32')
            rednp = ret_jim.np(0).astype(float)
            nirnp = ret_jim.np(1).astype(float)
            ndvi = (nirnp - rednp) / (rednp + nirnp)
            if offset != 0:
                ndvi += offset
            if scale != 1:
                ndvi *= scale
            ndvi[_np.isnan(ndvi)] = nodata
            ndvi = ndvi.astype(nptype)
            ret_jim.np(0)[:] = ndvi
            ret_jim.geometry.cropBand(0)
            if name is not None:
                ret_jim.properties.setDimension(name, 'band')
        return ret_jim


def NDVISeparateBands(jim_red,
                      jim_nir,
                      name: str = None):
    """Compute NDVI from two Jim objects.

    Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

    :param jim_red: Jim object with values of red
    :param jim_nir: Jim object with values of NIR
    :param name: name of band for NDVI
    :return: a Jim object with values of NDVI
    """
    if name is not None:
        ndvi = _pj.Jim(jim_nir._jipjim.pointOpNDI(jim_red._jipjim))
        ndvi.properties.setDimension(name, 'band')
        return ndvi
    else:
        return _pj.Jim(jim_nir._jipjim.pointOpNDI(jim_red._jipjim))


def setData(jim,
            value: float,
            bbox: list = None,
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
    :param bbox: bounding box (instead of ulx, uly, lrx, lry)
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

    jout.pixops.setData(value, bbox, ulx, uly, lrx, lry, bands, dx, dy, nogeo)
    jout.properties.setDimension(jim.properties.getDimension())
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
    ret_jim = _pj.Jim(jim_object._jipjim.pointOpSetLevel(min, max, val))
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


def setThreshold(jim_object,
                 **kwargs):
    """Apply min and max threshold to pixel values in raster dataset.

    :param jim_object: the Jim object on which to set threshold
    :param kwargs: See table :py:meth:`~pixops._PixOps.setThreshold`.

    for help, please refer to the corresponding
    method :py:meth:`~pixops._PixOps.setThreshold`.
    """
    ret_jim = _pj.Jim(jim_object._jipjim.setThreshold(kwargs))
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


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

    ret_jim = _pj.Jim(jout)
    ret_jim.properties.setDimension(jim1.properties.getDimension())
    return ret_jim


def simpleBitwiseOp(jim_object,
                    another_jim,
                    op: int,
                    *args):
    """Create Jim composed using a simple bitwise operation (coded with op).

    :param jim_object: Jim object
    :param another_jim: Jim object (to be sure that at least one is provided)
    :param op: integer for operation type
    :param args: Jim objects
    :return: Jim holding specified bitwise operation with from provided
        Jim objects
    """
    jout = _pj.Jim(jim_object)
    jims = [another_jim]
    jims.extend(args)
    for newJim in jims:
        jout._jipjim.d_pointOpBitwise(newJim._jipjim, op)

    ret_jim = _pj.Jim(jout)
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


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
    ret_jim = _pj.Jim(jim_object._jipjim.pointOpThresh(min, max, fg_val, bg_val))
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


def stretch(jim_object,
            **kwargs):
    """Stretch pixel values.

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

    ret_jim = _pj.Jim(jim_object._jipjim.stretch(kwargs))
    ret_jim.properties.setDimension(jim_object.properties.getDimension())
    return ret_jim


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
            raise _pj.exceptions.JimIllegalArgumentError(
                "Output type {} not supported".format(otype))
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
             red,
             nir,
             name: str = None,
             nodata=-2,
             scale: float = 1.0,
             offset: float = 0,
             addBand: bool = False):
        """Compute NDVI on the Jim object.

        Modifies the instance on which the method was called.

        :param red: index of band with values of red
        :param nir: index of band with values of NIR
        :param name: name of band for NDVI
        :param scale: scale to multiply NDVI
        :param offset: offset to add NDVI
        :param addBand: add NDVI as a new band
        """
        if self._jim_object.properties.getDataType() not in (
                'Float32','Float64'):
            if scale <= 1.0:
                raise UserWarning("Warning data type is not floating point\
                , use scale and offset or change data type")
        nptype = self._jim_object.np().dtype.type
        if isinstance(red, int):
            redindex = red
        else:
            redindex = self._jim_object.properties.getDimension(
                'band').index(red)
        if isinstance(nir, int):
            nirindex = nir
        else:
            nirindex = self._jim_object.properties.getDimension(
                'band').index(nir)
        if addBand:
            rednp = self._jim_object.np(redindex).astype(float)
            nirnp = self._jim_object.np(nirindex).astype(float)
            ndvi = (nirnp - rednp) / (rednp + nirnp)
            if offset != 0:
                ndvi += offset
            if scale != 1:
                ndvi *= scale
            ndvi[_np.isnan(ndvi)] = nodata
            ndvi = ndvi.astype(nptype)
            jimndvi = _pj.np2jim(ndvi)
            jimndvi.properties.copyGeoReference(self._jim_object)
            if name is not None:
                jimndvi.properties.setDimension(name, 'band')
            self._jim_object.geometry.stackBand(jimndvi)
        else:
            self._jim_object.geometry.cropBand([red, nir])
            self._jim_object.pixops.convert('GDT_Float32')
            rednp = self._jim_object.np(0).astype(float)
            nirnp = self._jim_object.np(1).astype(float)
            ndvi = (nirnp - rednp) / (rednp + nirnp)
            if offset != 0:
                ndvi += offset
            if scale != 1:
                ndvi *= scale
            ndvi[_np.isnan(ndvi)] = nodata
            ndvi = ndvi.astype(nptype)
            self._jim_object.np(0)[:] = ndvi
            self._jim_object.geometry.cropBand(0)
            if name is not None:
                self._jim_object.properties.setDimension(name, 'band')

    def NDVISeparateBands(self,
                          jim_nir,
                          name: str = None):
        """Compute NDVI from two Jims (call on red band, use NIR as param).

        Values in both red and NIR equal to 0 will obtain an NDVI value of -2)

        Modifies the instance on which the method was called.

        :param jim_nir: Jim object with values of NIR
        :param name: name of band for NDVI
        """

        if name is not None:
            ndvi = _pj.Jim(jim_nir._jipjim.pointOpNDI(self._jim_object._jipjim))
            ndvi.properties.setDimension(name, 'band')
            self._jim_object._set(ndvi._jim_object)
        else:
            self._jim_object._set(
                jim_nir._jipjim.pointOpNDI(self._jim_object._jipjim))

    def setData(self,
                value: float,
                bbox: list = None,
                ulx: float = None,
                uly: float = None,
                lrx: float = None,
                lry: float = None,
                bands = None,
                dx: int = 0,
                dy: int = 0,
                nogeo: bool = False):
        """Set range of pixels to value.

        :param value: new value for pixels of Jim object
        :param bbox: bounding box (instead of ulx, uly, lrx, lry)
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

        if bbox is not None:
            ulx = bbox[0]
            uly = bbox[1]
            lrx = bbox[2]
            lry = bbox[3]
        if bands is None:
            bands = list(range(self._jim_object.properties.nrOfBand()))
        elif not isinstance(bands,list):
            bands=[bands]

        if all(v is None for v in [ulx, uly, lrx, lry]) and dx == 0 and \
                dy == 0 and not nogeo:
            for band in bands:
                self._jim_object._jipjim.setData(value, band)
        else:
            if nogeo:
                if ulx is None:
                    ulx = 0
                if uly is None:
                    uly = 0
                if lrx is None:
                    lrx = self._jim_object.properties.nrOfCol()
                if lry is None:
                    lry = self._jim_object.properties.nrOfRow()
                uli = ulx
                ulj = uly
                lri = lrx
                lrj = lry
            else:
                if ulx is None:
                    ulx = self._jim_object.properties.getUlx()
                if uly is None:
                    uly = self._jim_object.properties.getUly()
                if lrx is None:
                    lrx = self._jim_object.properties.getLrx()
                if lry is None:
                    lry = self._jim_object.properties.getLry()
                upperLeftImage = self._jim_object.geometry.geo2image(ulx, uly)
                uli = upperLeftImage[0]
                ulj = upperLeftImage[1]
                lowerRightImage = self._jim_object.geometry.geo2image(lrx, lry)
                lri = lowerRightImage[0]
                lrj = lowerRightImage[1]
            if lri < 0 or lrj < 0:
                return None
            overlapColFrom = max(0, uli)
            overlapColTo = min(lri, self._jim_object.properties.nrOfCol())
            overlapRowFrom = max(0, ulj)
            overlapRowTo = min(lrj, self._jim_object.properties.nrOfRow())
            strideX = 1
            if dx > self._jim_object.properties.getDeltaX():
                strideX = int(dx // self._jim_object.properties.getDeltaX())
                assert strideX == dx / self._jim_object.properties.getDeltaX()
            strideY = 1
            if dy > self._jim_object.properties.getDeltaY():
                strideY = int(dy // self._jim_object.properties.getDeltaY())
                assert strideY == dy / self._jim_object.properties.getDeltaY()
            if self._jim_object.properties.nrOfPlane() > 1:
                self._jim_object[:, overlapRowFrom:overlapRowTo:strideY, overlapColFrom:overlapColTo:strideX] = value
            else:
                self._jim_object[overlapRowFrom:overlapRowTo:strideY, overlapColFrom:overlapColTo:strideX] = value

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
        |                  | (if pixel value > max set pixel value to no data)|
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

        Mask all values not within [0, 250] and set to 255::

            jim0[(jim0<0) | (jim0>250)] = 255

        Mask all values not within [0, 250] and set to 255 (no data)::

            jim_threshold = jim.setThreshold(min=0, max=250, nodata=255)
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
        """Stretch pixel values.

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
