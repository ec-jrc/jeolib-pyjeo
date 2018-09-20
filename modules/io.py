import jiplib as _jl


def GDALRead(fn, band=0, nXOff=0, nYOff=0, nXSize=None, nYSize=None,
             nBufXSize=None, nBufYSize=None):
    """Read a GDAL compatible image stored in the filename.

    :param fn: a string for the name of an image file (possibly its path)
    :param band: an integer for the band number, 0 for first band
    :param nXOff: integer for the pixel offset to the top left corner
    :param nYOff: integer for the line offset to the top left corner
    :param nXSize: integer for the width of the region
    :param nYSize: integer for the height of the region
    :param nBufXSize: integer for the number of columns of output image
    :param nBufYSize: integer for the number of lines of output image
    :return: a Jim object
    """
    if nXSize is None:
        nXSize = nXOff
    if nYSize is None:
        nYSize = nYOff
    if nBufXSize is None:
        nBufXSize = nXSize
    if nBufYSize is None:
        nBufYSize = nYSize
    return _jl.GDALRead(fn, band, nXOff, nYOff, nXSize, nYSize,
                        nBufXSize, nBufYSize)


def createPyJim(filename, **kwargs):
    """Create an empty Jim object.

    Created object is an instance of the basis image class of the pyJEO library

    :param filename: Path to a raster file
    :return: a Jim object
    """
    return _jl.createJim(filename, **kwargs)


def createPyVector(filename):
    """Create an empty VectorOgr object.

    Created object is an instance of the basis vector class of the pyJEO
    library

    :param filename: Path to a vector file
    :return: a VectorOgr object
    """
    return _jl.createVector(filename)


class _IO():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object
