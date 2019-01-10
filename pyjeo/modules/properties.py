"""Module for accessing Jim attributes and geospatial informations."""

import pyjeo as _pj




def imageInfo(jim_object):
    """Return image information (number of lines, columns, etc.)

    """
    return _pj.Jim(jim_object._jipjim.imageInfo())



class _Properties():
    """Define all properties methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def imageInfo(self):
        """Return image information (number of lines, columns, etc.)
        """
        self._jim_object._jipjim.imageInfo()

    def nrOfCol(self):
        """Get number of columns in this raster dataset.

        :return: The number of columns in this raster dataset
        """
        return self._jim_object._jipjim.nrOfCol()

    def nrOfRow(self):
        """Get number of rows in this raster dataset.

        :return: The number of rows in this raster dataset
        """
        return self._jim_object._jipjim.nrOfRow()

    def nrOfPlane(self):
        """Get number of planes in this raster dataset.

        :return: The number of planes in this raster dataset
        """
        return self._jim_object._jipjim.nrOfPlane()

    def nrOfBand(self):
        """Get number of bands in this raster dataset.

        :return: The number of bands in this raster dataset
        """
        return self._jim_object._jipjim.nrOfBand()

    def printNoDataVals(self):
        """Print the list of no data values of this raster dataset."""
        self._jim_object._jipjim.printNoDataValues()

    def pushNoDataVal(self, value):
        """Push a no data value for this raster dataset."""
        self._jim_object._jipjim.pushNoDataValue(value)

    def setNoDataVals(self, value):
        """Set the list of no data value for this raster dataset."""
        if isinstance(value, list):
            self._jim_object._jipjim.setNoData(value)
        elif type(value) in (float, int):
            self._jim_object._jipjim.setNoDataValue(value)
        else:
            raise TypeError('Error: setNoDataVals not implemented for value '
                            'type {}'.format(type(value)))

    def getNoDataVals(self):
        """Get the list of no data values."""
        return self._jim_object._jipjim.getNoDataValues()

    def clearNoData(self):
        """Clear the list of no data values for this raster dataset."""
        self._jim_object._jipjim.clearNoData()

    def getDataType(self):
        """Get the internal datatype for this raster dataset.

        :return: The datatype id of this Jim object
        """
        return self._jim_object._jipjim.getDataType()

    def covers(self, *args):
        """Check if a geolocation is covered by this dataset.

        You can check the coverage for a :ref:`point of interest
        <covers1>` or :ref:`region of interest <covers2>`.

        .. _covers1:

        :Using a point coordinates:

        * ``x`` (float): x coordinate in spatial reference system of the raster dataset
        * ``y`` (float): y coordinate in spatial reference system of the raster dataset

        .. _covers2:

        :Using a region of interest:

        * ``ulx`` (float): upper left x coordinate in spatial reference system of the raster dataset
        * ``uly`` (float): upper left y coordinate in spatial reference system of the raster dataset
        * ``lrx`` (float): lower right x coordinate in spatial reference system of the raster dataset
        * ``lry`` (float): lower right x coordinate in spatial reference system of the raster dataset
        * ``all`` (bool): set to True if the entire bounding box must be covered by the raster dataset

        Returns:
        True if the raster dataset covers the point or region of interest.
        """
        return self._jim_object._jipjim.covers(*args)

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

    def setGeoTransform(self, *args):
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

    def copyGeoTransform(self, sec_jim_object):
        """Copy geotransform information from another georeferenced image."""
        self._jim_object._jipjim.copyGeoTransform(sec_jim_object._jipjim)

    def getProjection(self):
        """Get the projection for this dataget in well known text (wkt) format.

        :return: The projection string in well known text format.
        """
        return self._jim_object._jipjim.getProjection()

    def setProjection(self, *args):
        """Set the projection for this dataset.

        Set it in well known text (wkt) format or as a string epsg:<epsgcode>.

        :args: the projection for this dataset in well known text (wkt) format
            or as a string epsg:<epsgcode>.
        """
        self._jim_object._jipjim.setProjection(*args)

    def getCenterPos(self):
        """
        Get the center position of this dataset in georeferenced coordinates.

        :return: A list with the center position of this dataset in
            georeferenced coordinates.
        """
        return self._jim_object._jipjim.getCenterPos()

    def getBBox(self):
        """Get the bounding box (georeferenced) coordinates of this dataset.

        :return: A list with upper left x, upper left y, lower right x, and
            lower right y
        """
        return self._jim_object._jipjim.getBoundingBox()

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

    def getRefPix(self, *args):
        """Calculate the reference pixel as the center of gravity pixel.

        Calculated as the weighted average of all values not taking into
        account no data values for a specific band (start counting from 0).

        :return: The reference pixel as the centre of gravity pixel (weighted
            average of all values not taking into account no data values) for
            a specific band (start counting from 0).
        """
        return self._jim_object._jipjim.getRefPix(*args)


class _PropertiesList():
    """Define all properties methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def pushNoDataVal(self, value):
        """Push a no data value for this raster JimList object."""
        self._jim_list._jipjimlist.pushNoDataValue(value)

    def clearNoData(self):
        """Clear the list of no data values for this JimList object."""
        self._jim_list._jipjimlist.clearNoData()

    def getNoDataVals(self):
        """Get the list of no data values."""
        return self._jim_list._jipjimlist.getNoDataValues()

    def covers(self, *args):
        """Check if a geolocation is covered by this dataset.

        You can check the coverage for a :ref:`point of interest
        <covers1>` or :ref:`region of interest <covers2>`.

        .. _covers1:

        :Using a point coordinates:

        * ``x`` (float): x coordinate in spatial reference system of the raster dataset
        * ``y`` (float): y coordinate in spatial reference system of the raster dataset

        .. _covers2:

        :Using a region of interest:

        * ``ulx`` (float): upper left x coordinate in spatial reference system of the raster dataset
        * ``uly`` (float): upper left y coordinate in spatial reference system of the raster dataset
        * ``lrx`` (float): lower right x coordinate in spatial reference system of the raster dataset
        * ``lry`` (float): lower right x coordinate in spatial reference system of the raster dataset
        * ``all`` (bool): set to True if the entire bounding box must be covered by the raster dataset

        Returns:
        True if the raster dataset covers the point or region of interest.
        """
        return self._jim_list._jipjimlist.covers(*args)

    def getBBox(self):
        """Get the bounding box (georeferenced) coordinates of this dataset.

        :return: A list with upper left x, upper left y, lower right x, and
            lower right y
        """
        return self._jim_list._jipjimlist.getBoundingBox()

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

    def selectGeo(self, *args):
        self._jim_list._jipjimlist.selectGeo(*args)
        self._jim_list._set(self._jim_list)


class _PropertiesVect():
    """Define all properties methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    def getLayerCount(self):
        return self._jim_vect._jipjimvect.getLayerCount()

    def getFeatureCount(self):
        return self._jim_vect._jipjimvect.getFeatureCount()

    def getBBox(self):
        return self._jim_vect._jipjimvect.getBoundingBox()

    def getUlx(self):
        return self._jim_vect._jipjimvect.getUlx()

    def getUly(self):
        return self._jim_vect._jipjimvect.getUly()

    def getLrx(self):
        return self._jim_vect._jipjimvect.getLrx()

    def getLry(self):
        return self._jim_vect._jipjimvect.getLry()

    def getFields(self, layer=0):
        return self._jim_vect._jipjimvect.getFields(layer)
