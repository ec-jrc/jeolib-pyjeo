import jiplib as _jl

class _Properties():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

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
        _jl.dumpxyz(self._jim_object, x, y, z, nx, ny)

    def nrOfCol(self):
        """Get number of columns in this raster dataset.

        :return: The number of columns in this raster dataset
        """
        return self._jim_object.nrOfCol()

    def nrOfRow(self):
        """Get number of rows in this raster dataset.

        :return: The number of rows in this raster dataset
        """
        return self._jim_object.nrOfRow()

    def nrOfPlane(self):
        """Get number of planes in this raster dataset.

        :return: The number of planes in this raster dataset
        """
        return self._jim_object.nrOfPlane()

    def nrOfBand(self):
        """Get number of bands in this raster dataset.

        :return: The number of bands in this raster dataset
        """
        return self._jim_object.nrOfBand()

    def printNoDataVals(self):
        """
        Print the list of no data values of this raster dataset
        """
        self._jim_object.printNoDataValues()

    def pushNoDataVal(self, value):
        """
        Push a no data value for this raster dataset
        """
        self._jim_object.pushNoDataValue(value)

    def setNoDataVals(self, value):
        """
        Set the list of no data value for this raster dataset
        """
        if isinstance(value, list):
            self._jim_object.setNoData(value)
        else:
            self._jim_object.setNoDataValue(value)

    def getNoDataVals(self, value):
        """
        Get the list of no data values
        """
        return self._jim_object.getNoDataValues()

    def clearNoData(self):
        """
        Clear the list of no data values for this raster dataset
        """
        self._jim_object.clearNoData()

    def getDataType(self):
        """
        Get the internal datatype for this raster dataset
        :return: The datatype id of this Jim object
        """
        return self._jim_object.getDataType()

    def covers(self, *args):
        """
        Check if a geolocation is covered by this dataset. You can check the coverage for a :ref:`point of interest <covers1>` or :ref:`region of interest <covers2>`).

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
        return self._jim_object.covers(*args)

    def getGeoTransform(self):
        """
        Get the geotransform data for this dataset as a list of floats.

        :returns: List of floats with geotransform data:

        * [0] top left x
        * [1] w-e pixel resolution
        * [2] rotation, 0 if image is "north up"
        * [3] top left y
        * [4] rotation, 0 if image is "north up"
        * [5] n-s pixel resolution
        """
        return self._jim_object.getGeoTransform()

    def setGeoTransform(self, *args):
        """
        Set the geotransform data for this dataset.

        :args: List of floats with geotransform data:

        * [0] top left x
        * [1] w-e pixel resolution
        * [2] rotation, 0 if image is "north up"
        * [3] top left y
        * [4] rotation, 0 if image is "north up"
        * [5] n-s pixel resolution
        """
        self._jim_object.setGeoTransform(*args)

    def copyGeoTransform(self, sec_jim_object):
        """
        Copy geotransform information from another georeferenced image.
        """
        self._jim_object.copyGeoTransform(sec_jim_object)

    def getProjection(self):
        """
        Get the projection for this dataget in well known text (wkt) format.

        :return: The projection string in well known text format.
        """
        return self._jim_object.getProjection()

    def setProjection(self, *args):
        """
        Set the projection for this dataset in well known text (wkt) format or as a string epsg:<epsgcode>.

        :args: the projection for this dataset in well known text (wkt) format or as a string epsg:<epsgcode>.
        """
        self._jim_object.setProjection(*args)

    def getCenterPos(self):
        """
        Get the center position of this dataset in georeferenced coordinates

        :return: A list with the center position of this dataset in georeferenced coordinates.
        """
        return self._jim_object.getCenterPos()

    def getUlx(self):
        """
        Get the upper left corner x (georeferenced) coordinate of this dataset

        :return: The upper left corner x (georeferenced) coordinate of this dataset
        """
        return self._jim_object.getUlx()

    def getUly(self):
        """
        Get the upper left corner y (georeferenced) coordinate of this dataset

        :return: The upper left corner y (georeferenced) coordinate of this dataset
        """
        return self._jim_object.getUly()

    def getLrx(self):
        """
        Get the lower left corner x (georeferenced) coordinate of this dataset

        :return: The lower left corner x (georeferenced) coordinate of this dataset
        """
        return self._jim_object.getLrx()

    def getLry(self):
        """
        Get the lower left corner y (georeferenced) coordinate of this dataset

        :return: The lower left corner y (georeferenced) coordinate of this dataset
        """
        return self._jim_object.getLry()

    def getDeltaX(self):
        """
        Get the pixel cell spacing in x.

        :return: The pixel cell spacing in x.
        """
        return self._jim_object.getDeltaX()

    def getDeltaY(self):
        """
        Get the piyel cell spacing in y.

        :return: The piyel cell spacing in y.
        """
        return self._jim_object.getDeltaY()

    def getRefPix(self, *args):
        """
        Calculate the reference pixel as the center of gravity pixel (weighted average of all values not taking into account no data values) for a specific band (start counting from 0).

        :return: The reference pixel as the centre of gravity pixel (weighted average of all values not taking into account no data values) for a specific band (start counting from 0).
        """
        return self._jim_object.getRefPix(*args)
