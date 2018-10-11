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

        :return:The number of columns in this raster dataset
        """
        return self._jim_object.nrOfCol()

    def nrOfRow(self):
        """Get number of rows in this raster dataset.

        :return:The number of rows in this raster dataset
        """
        return self._jim_object.nrOfRow()

    def nrOfPlane(self):
        """Get number of planes in this raster dataset.

        :return:The number of planes in this raster dataset
        """
        return self._jim_object.nrOfPlane()

    def nrOfBand(self):
        """Get number of bands in this raster dataset.

        :return:The number of bands in this raster dataset
        """
        return self._jim_object.nrOfBand()

    def printNoDataVals(self):
        self._jim_object.printNoDataValues()

    def pushNoDataVal(self, value):
        self._jim_object.pushNoDataValue(value)

    def setNoDataVals(self, value):
        if isinstance(value, list):
            self._jim_object.setNoData(value)
        else:
            self._jim_object.setNoDataValue(value)

    def clearNoData(self):
        self._jim_object.clearNoData()
