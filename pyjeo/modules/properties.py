"""Module for accessing Jim attributes and geospatial informations."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2020 European Union (Joint Research Centre)
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
import osgeo
import pyjeo as _pj

# def imageInfo(jim_object):
#     """Return image information (number of lines, columns, etc.)

#     """
#     return _pj.Jim(jim_object._jipjim.imageInfo())


def isEqual(first_jim,
            second_jim):
    """Check if the values of one Jim object are the same as in another one.

    :param first_jim: a Jim or JimVect object
    :param second_jim: a Jim or JimVect object
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
                if not _np.array_equal(first_jim.np(iband),
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
                    if not _np.array_equal(first_plane.np(iband),
                                           second_plane.np(iband)):
                        return False
            return True
    elif isinstance(first_jim, _pj.JimVect):
        if first_jim.properties.isEmpty():
            raise _pj.exceptions.JimVectEmptyError('first_jim is empty')
        if second_jim.properties.isEmpty():
            return False

        if not isinstance(second_jim, _pj.JimVect):
            return False

        if first_jim.properties.getLayerCount() != \
                second_jim.properties.getLayerCount():
            return False
        if first_jim.properties.getFeatureCount() != \
                second_jim.properties.getFeatureCount():
            return False
        if first_jim.properties.getProjection() != \
                second_jim.properties.getProjection():
            return False
        if first_jim.properties.getBBox() != second_jim.properties.getBBox():
            return False

        for ilayer in range(0, first_jim.properties.getLayerCount()):
            #todo: check geometry equality for each feature
            jim1_fieldnames = first_jim.properties.getFieldNames()
            jim2_fieldnames = second_jim.properties.getFieldNames()
            if jim1_fieldnames != jim2_fieldnames:
                return False
            if not _np.array_equal(first_jim.np(ln=ilayer),
                                   second_jim.np(ln=ilayer)):
                return False

        return True
    else:
        return False


class _Properties(_pj.modules.JimModuleBase):
    """Define all properties methods."""

    def clearNoData(self):
        """Clear the list of no data values for this raster dataset."""
        self._jim_object._jipjim.clearNoData()

    def copyGeoTransform(self,
                         sec_jim_object):
        """Copy geotransform information from another georeferenced image.

        :param sec_jim_object: Jim from which the geotransform information
            should be copied.
        """
        self._jim_object._jipjim.copyGeoTransform(sec_jim_object._jipjim)

    def copyGeoReference(self,
                         sec_jim_object):
        """Copy projection and geotransform information from another image.

        :param sec_jim_object: Jim from which the projection and geotransform
            information should be copied.
        """
        self._jim_object._jipjim.copyGeoReference(sec_jim_object._jipjim)

    def covers(self,
               *args):
        """Check if a geolocation is covered by this dataset.

        You can check the coverage for a :ref:`point of interest
        <covers1>` or :ref:`region of interest <covers2>`.

        .. _covers1:

        :Using a point coordinates:

        * ``x`` (float): x coordinate in spatial reference system of \
                         the rasterdataset
        * ``y`` (float): y coordinate in spatial reference system of \
                         the raster dataset

        .. _covers2:

        :Using a region of interest:

        * ``ulx`` (float): upper left x coordinate in spatial reference \
                           system of the raster dataset
        * ``uly`` (float): upper left y coordinate in spatial reference \
                           system of the raster dataset
        * ``lrx`` (float): lower right x coordinate in spatial reference \
                           system of the raster dataset
        * ``lry`` (float): lower right x coordinate in spatial reference \
                           system of the raster dataset
        * ``all`` (bool): set to True if the entire bounding box must be \
                          covered by the raster dataset

        Returns:
        True if the raster dataset covers the point or region of interest.
        """
        return self._jim_object._jipjim.covers(*args)

    def getBBox(self, t_srs = None):
        """Get the bounding box (georeferenced) coordinates of this dataset.

        :t_srs: A target reference system (e.g., 'epsg:3035, 3035, or WKT)
        :return: A list with upper left x, upper left y, lower right x, and
            lower right y
        """
        bbox = self._jim_object._jipjim.getBoundingBox()
        if t_srs is not None:
            # create coordinate transformation
            inSpatialRef = osgeo.osr.SpatialReference()
            outSpatialRef = osgeo.osr.SpatialReference()

            inSpatialRef.ImportFromWkt(self._jim_object.properties.getProjection())
            if int(osgeo.__version__[0]) >= 3:
                #hanges axis order: https://github.com/OSGeo/gdal/issues/1546
                inSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                outSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            if isinstance(t_srs, int):
                outSpatialRef.ImportFromEPSG(t_srs)
            elif isinstance(t_srs, str):
                if 'EPSG:' in t_srs:
                    t_srs=int(t_srs.split("EPSG:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                elif 'epsg:' in t_srs:
                    t_srs=int(t_srs.split("epsg:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                else:
                    outSpatialRef.ImportFromWkt(t_srs)
            else:
                raise JimIllegalArgumentError("Error: coordinate reference "
                                              "system must be integer "
                                              "representing epsg code or string")
            coordTransform = osgeo.osr.CoordinateTransformation(inSpatialRef,
                                                                outSpatialRef)
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[0], bbox[1])
            point.Transform(coordTransform)
            ulx = point.GetX()
            uly = point.GetY()
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[2], bbox[3])
            point.Transform(coordTransform)
            lrx = point.GetX()
            lry = point.GetY()
            bbox = [ulx, uly, lrx, lry]
        return bbox

    def getCenterPos(self):
        """
        Get the center position of this dataset in georeferenced coordinates.

        :return: A list with the center position of this dataset in
            georeferenced coordinates.
        """
        return self._jim_object._jipjim.getCenterPos()

    def getDataType(self):
        """Get the internal datatype for this raster dataset.

        :return: The datatype id of this Jim object
        """
        otype = self._jim_object._jipjim.getDataType()
        if otype == 1:
            return 'Byte'
        elif otype == 2:
            return 'UInt16'
        elif otype == 3:
            return 'Int16'
        elif otype == 4:
            return 'UInt32'
        elif otype == 5:
            return 'Int32'
        elif otype == 6:
            return 'Float32'
        elif otype == 7:
            return 'Float64'
        elif otype == 8:
            return 'CInt16'
        elif otype == 9:
            return 'CInt32'
        elif otype == 10:
            return 'CFloat32'
        elif otype == 11:
            return 'CFloat64'
        elif otype == 12:
            return 'TypeCount'
        elif otype == 14:
            return 'Int64'
        elif otype == 15:
            return 'UInt64'
        elif otype == 16:
            return 'JDT_Word'
        else:
            raise _pj.exceptions.JimInnerParametersError(
                "Unknown data format".format(otype))

    def getDeltaX(self):
        """Get the pixel cell spacing in x.

        :return: The pixel cell spacing in x.
        """
        return self._jim_object._jipjim.getDeltaX()

    def getDeltaY(self):
        """Get the piyel cell spacing in y.

        :return: The piyel cell spacing in y.
        """
        return self._jim_object._jipjim.getDeltaY()

    def getGeoTransform(self):
        """Get the geotransform data for this dataset as a list of floats.

        :returns: List of floats with geotransform data:

        * [0] top left x
        * [1] w-e pixel resolution
        * [2] rotation, 0 if image is "north up"
        * [3] top left y
        * [4] rotation, 0 if image is "north up"
        * [5] n-s pixel resolution
        """
        return self._jim_object._jipjim.getGeoTransform()

    def getLrx(self):
        """Get the lower left corner x coordinate of this dataset.

        :return: The lower left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_object._jipjim.getLrx()

    def getLry(self):
        """Get the lower left corner y coordinate of this dataset.

        :return: The lower left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_object._jipjim.getLry()

    def getNoDataVals(self):
        """Get the list of no data values."""
        return self._jim_object._jipjim.getNoDataValues()

    def getProjection(self):
        """Get the projection for this dataset in well known text (wkt) format.

        :return: The projection string in well known text format.
        """
        return self._jim_object._jipjim.getProjection()

    def getRefPix(self,
                  *args):
        """Calculate the reference pixel as the center of gravity pixel.

        Calculated as the weighted average of all values not taking into
        account no data values for a specific band (start counting from 0).

        :return: The reference pixel as the centre of gravity pixel (weighted
            average of all values not taking into account no data values) for
            a specific band (start counting from 0).
        """
        return self._jim_object._jipjim.getRefPix(*args)

    def getUlx(self):
        """Get the upper left corner x coordinate of this dataset.

        :return: The upper left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_object._jipjim.getUlx()

    def getUly(self):
        """Get the upper left corner y coordinate of this dataset.

        :return: The upper left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_object._jipjim.getUly()

    def imageInfo(self):
        """Return image information (number of lines, columns, etc.)."""
        self._jim_object._jipjim.imageInfo()

    def isEqual(self,
                other):
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
                    if not _np.array_equal(self._jim_object.np(iband),
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
                        if not _np.array_equal(first_plane.np(iband),
                                               second_plane.np(iband)):
                            return False
                return True
        else:
            return False

    def nrOfBand(self):
        """Get number of bands in this raster dataset.

        :return: The number of bands in this raster dataset
        """
        return self._jim_object._jipjim.nrOfBand()

    def nrOfCol(self):
        """Get number of columns in this raster dataset.

        :return: The number of columns in this raster dataset
        """
        return self._jim_object._jipjim.nrOfCol()

    def nrOfPlane(self):
        """Get number of planes in this raster dataset.

        :return: The number of planes in this raster dataset
        """
        return self._jim_object._jipjim.nrOfPlane()

    def nrOfRow(self):
        """Get number of rows in this raster dataset.

        :return: The number of rows in this raster dataset
        """
        return self._jim_object._jipjim.nrOfRow()

    def printNoDataVals(self):
        """Print the list of no data values of this raster dataset."""
        self._jim_object._jipjim.printNoDataValues()

    def pushNoDataVal(self,
                      value: float):
        """Push a no data value for this raster dataset."""
        self._jim_object._jipjim.pushNoDataValue(value)

    def setGeoTransform(self,
                        *args):
        """Set the geotransform data for this dataset.

        :args: List of floats with geotransform data:

        * [0] top left x
        * [1] w-e pixel resolution
        * [2] rotation, 0 if image is "north up"
        * [3] top left y
        * [4] rotation, 0 if image is "north up"
        * [5] n-s pixel resolution
        """
        self._jim_object._jipjim.setGeoTransform(*args)

    def setNoDataVals(self,
                      value: float):
        """Set the list of no data value for this raster dataset."""
        if isinstance(value, list):
            self._jim_object._jipjim.setNoData(value)
        elif type(value) in (float, int):
            self._jim_object._jipjim.setNoDataValue(value)
        else:
            raise _pj.exceptions.JimIllegalArgumentError(
                'setNoDataVals not implemented for value type {}'.format(
                    type(value)))

    def setProjection(self,
                      *args):
        """Set the projection for this dataset.

        Set it in well known text (wkt) format or as a string epsg:<epsgcode>.

        :args: the projection for this dataset in well known text (wkt) format
            or as a string epsg:<epsgcode>.
        """
        self._jim_object._jipjim.setProjection(*args)


class _PropertiesList(_pj.modules.JimListModuleBase):
    """Define all properties methods for JimLists."""

    def clearNoData(self):
        """Clear the list of no data values for this JimList object."""
        self._jim_list._jipjimlist.clearNoData()

    def covers(self,
               *args):
        """Check if a geolocation is covered by this dataset.

        You can check the coverage for a :ref:`point of interest
        <coversl1>` or :ref:`region of interest <coversl2>`.

        .. _coversl1:

        :Using a point coordinates:

        * ``x`` (float): x coordinate in spatial reference system of \
                         the raster dataset
        * ``y`` (float): y coordinate in spatial reference system of \
                         the raster dataset

        .. _coversl2:

        :Using a region of interest:

        * ``ulx`` (float): upper left x coordinate in spatial reference \
                           system of the raster dataset
        * ``uly`` (float): upper left y coordinate in spatial reference \
                           system of the raster dataset
        * ``lrx`` (float): lower right x coordinate in spatial reference \
                           system of the raster dataset
        * ``lry`` (float): lower right x coordinate in spatial reference \
                           system of the raster dataset
        * ``all`` (bool): set to True if the entire bounding box must be \
                          covered by the raster dataset

        Returns:
        True if the raster dataset covers the point or region of interest.
        """
        return self._jim_list._jipjimlist.covers(*args)

    def getBBox(self, t_srs = None):
        """Get the bounding box (georeferenced) coordinates of this dataset.

        :t_srs: A target reference system (e.g., 'epsg:3035, 3035, or WKT)
        :return: A list with upper left x, upper left y, lower right x, and
            lower right y
        """
        bbox = self._jim_list._jipjimlist.getBoundingBox()
        if t_srs is not None:
            # create coordinate transformation
            inSpatialRef = osgeo.osr.SpatialReference()
            outSpatialRef = osgeo.osr.SpatialReference()

            inSpatialRef.ImportFromWkt(self._jim_list.properties.getProjection())
            if int(osgeo.__version__[0]) >= 3:
                #hanges axis order: https://github.com/OSGeo/gdal/issues/1546
                inSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                outSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            if isinstance(t_srs, int):
                outSpatialRef.ImportFromEPSG(t_srs)
            elif isinstance(t_srs, str):
                if 'EPSG:' in t_srs:
                    t_srs=int(t_srs.split("EPSG:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                elif 'epsg:' in t_srs:
                    t_srs=int(t_srs.split("epsg:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                else:
                    outSpatialRef.ImportFromWkt(t_srs)
            else:
                raise JimIllegalArgumentError("Error: coordinate reference "
                                              "system must be integer "
                                              "representing epsg code or string")
            coordTransform = osgeo.osr.CoordinateTransformation(inSpatialRef,
                                                                outSpatialRef)
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[0], bbox[1])
            point.Transform(coordTransform)
            ulx = point.GetX()
            uly = point.GetY()
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[2], bbox[3])
            point.Transform(coordTransform)
            lrx = point.GetX()
            lry = point.GetY()
            bbox = [ulx, uly, lrx, lry]
        return bbox

    def getLrx(self):
        """Get the lower left corner x coordinate of this dataset.

        :return: The lower left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_list._jipjimlist.getLrx()

    def getLry(self):
        """Get the lower left corner y coordinate of this dataset.

        :return: The lower left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_list._jipjimlist.getLry()

    def getNoDataVals(self):
        """Get the list of no data values."""
        return self._jim_list._jipjimlist.getNoDataValues()

    def getUlx(self):
        """Get the upper left corner x coordinate of this dataset.

        :return: The upper left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_list._jipjimlist.getUlx()

    def getUly(self):
        """Get the upper left corner y coordinate of this dataset.

        :return: The upper left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_list._jipjimlist.getUly()

    def pushNoDataVal(self,
                      value: float):
        """Push a no data value for this raster JimList object."""
        self._jim_list._jipjimlist.pushNoDataValue(value)

    def selectGeo(self, *args):
        """Select geographical properties (ulx, uly, ...)."""
        self._jim_list._jipjimlist.selectGeo(*args)
        self._jim_list._set(self._jim_list)


class _PropertiesVect(_pj.modules.JimVectModuleBase):
    """Define all properties methods for JimVects."""

    def getBBox(self, t_srs = None):
        """Get the bounding box (georeferenced) coordinates of this dataset.

        :return: A list with upper left x, upper left y, lower right x, and
            lower right y
        """
        bbox = self._jim_vect._jipjimvect.getBoundingBox()
        if t_srs is not None:
            # create coordinate transformation
            inSpatialRef = osgeo.osr.SpatialReference()
            outSpatialRef = osgeo.osr.SpatialReference()

            inSpatialRef.ImportFromWkt(self._jim_vect.properties.getProjection())
            if int(osgeo.__version__[0]) >= 3:
                #hanges axis order: https://github.com/OSGeo/gdal/issues/1546
                inSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
                outSpatialRef.SetAxisMappingStrategy(
                    osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            if isinstance(t_srs, int):
                outSpatialRef.ImportFromEPSG(t_srs)
            elif isinstance(t_srs, str):
                if 'EPSG:' in t_srs:
                    t_srs=int(t_srs.split("EPSG:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                elif 'epsg:' in t_srs:
                    t_srs=int(t_srs.split("epsg:",1)[1])
                    outSpatialRef.ImportFromEPSG(t_srs)
                else:
                    outSpatialRef.ImportFromWkt(t_srs)
            else:
                raise JimIllegalArgumentError("Error: coordinate reference "
                                              "system must be integer "
                                              "representing epsg code or string")
            coordTransform = osgeo.osr.CoordinateTransformation(inSpatialRef,
                                                                outSpatialRef)
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[0], bbox[1])
            point.Transform(coordTransform)
            ulx = point.GetX()
            uly = point.GetY()
            point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
            point.AddPoint(bbox[2], bbox[3])
            point.Transform(coordTransform)
            lrx = point.GetX()
            lry = point.GetY()
            bbox = [ulx, uly, lrx, lry]
        return bbox

    def getFeatureCount(self,
                        layer: int = None):
        """Get the number of features of this dataset.

        :layer: the layer number (default is all layers)
        :return: the number of features of this layer
        """
        if layer is not None:
            return self._jim_vect._jipjimvect.getFeatureCount(layer)
        else:
            return self._jim_vect._jipjimvect.getFeatureCount()

    def getFieldNames(self,
                      layer: int = 0):
        """Get the list of field names for this dataset.

        :param layer: The layer to get the field names, index starting from 0
            (default is 0: first layer)
        :return: The list of field names
            dataset
        """
        return self._jim_vect._jipjimvect.getFieldNames(layer)

    def getLayerCount(self):
        """Get the number of layers of this dataset.

        :return: the number of layers of this dataset
        """
        return self._jim_vect._jipjimvect.getLayerCount()

    def getLrx(self):
        """Get the lower left corner x coordinate of this dataset.

        :return: The lower left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_vect._jipjimvect.getLrx()

    def getLry(self):
        """Get the lower left corner y coordinate of this dataset.

        :return: The lower left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_vect._jipjimvect.getLry()

    def getProjection(self,
                      layer: int = 0):
        """Get the projection for this dataset in well known text (wkt) format.

        :param layer: The layer to get the projection from, index starting
            from 0 (default is 0: first layer)
        :return: The projection string in well known text format.
        """
        return self._jim_vect._jipjimvect.getProjection(layer)

    def getUlx(self):
        """Get the upper left corner x coordinate of this dataset.

        :return: The upper left corner x (georeferenced) coordinate of this
            dataset
        """
        return self._jim_vect._jipjimvect.getUlx()

    def getUly(self):
        """Get the upper left corner y coordinate of this dataset.

        :return: The upper left corner y (georeferenced) coordinate of this
            dataset
        """
        return self._jim_vect._jipjimvect.getUly()

    def isEmpty(self):
        """Check if object contains features (non-emtpy).

        :return: True if empty, False if not
        """
        return self._jim_vect._jipjimvect.isEmpty()

    def isEqual(self,
                other):
        """Check if the values of one Jim object are the same as in another.

        :param other: a JimVect object
        :return: True if the values are equal, zero otherwise
        """
        if self._jim_vect.properties.isEmpty():
            raise _pj.exceptions.JimVectEmptyError(
                'JimVect on which the method was called is empty')

        if not isinstance(other, _pj.JimVect):
            return False
        if other.properties.isEmpty():
            return False

        if self._jim_vect.properties.getLayerCount() != \
                other.properties.getLayerCount():
            return False
        if self._jim_vect.properties.getFeatureCount() != \
                other.properties.getFeatureCount():
            return False
        if self._jim_vect.properties.getProjection() != \
                other.properties.getProjection():
            return False
        if self._jim_vect.properties.getBBox() != other.properties.getBBox():
            return False

        for ilayer in range(0, self._jim_vect.properties.getLayerCount()):
            #todo: check geometry equality for each feature
            jim1_fieldnames = self._jim_vect.properties.getFieldNames()
            jim2_fieldnames = other.properties.getFieldNames()
            if jim1_fieldnames != jim2_fieldnames:
                return False
            if not _np.array_equal(self._jim_vect.np(ln=ilayer),
                                   other.np(ln=ilayer)):
                return False

        return True
