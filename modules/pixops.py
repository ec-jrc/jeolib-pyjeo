import jiplib as _jl
import pyjeo as _pj

def convert(jim_object, **kwargs):
    """Convert Jim image with respect to data type.

    :param jim_object: a Jim object
    :param kwargs: See table below
    :return: a Jim object
    Modifies the instance on which the method was called.


    +------------------+---------------------------------------------------------------------------------+
    | key              | value                                                                           |
    +==================+=================================================================================+
    | otype            | Data type for output image                                                      |
    +------------------+---------------------------------------------------------------------------------+
    | scale            | Scale output: output=scale*input+offset                                         |
    +------------------+---------------------------------------------------------------------------------+
    | offset           | Apply offset: output=scale*input+offset                                         |
    +------------------+---------------------------------------------------------------------------------+
    | autoscale        | Scale output to min and max, e.g., [0,255]                                      |
    +------------------+---------------------------------------------------------------------------------+
    | a_srs            | Override the projection for the output file                                     |
    +------------------+---------------------------------------------------------------------------------+

    .. note::
        To ignore some pixels from the extraction process, see list of :ref:`mask <extract_mask>` key values:

    Example:

    Convert data type of input image to byte using autoscale::

    jim0=jl.io.createJim('/path/to/raster.tif')
    jim0.convert(otype=Byte,autoscale=[0,255])

    Clip raster dataset between 0 and 255 (set all other values to 0), then convert data type to byte::

    jim1=jl.io.createJim('/path/to/raster.tif')
    jim1.setThreshold(min=0,max=255,nodata=0)
    jim1.convert({'otype':'Byte'})
    """
    # return _pj.Jim(jim_object.convert(kwargs))
    retJim=_pj.io.createJim(jim_object)
    return _pj.Jim(retJim.convert(kwargs))


def NDVI(redJim, nirJim):
    """Compute NDVI from two Jim objects.

    :param redJim: Jim object with values of red
    :param nirJim: Jim object with values of NIR
    :return: a Jim object with values of NDVI
    """
    numerator = _pj.Jim(nirJim.convertToFloat32()) - \
                _pj.Jim(redJim.convertToFloat32())
    denominator = _pj.Jim(nirJim.convertToFloat32()) + \
                  _pj.Jim(redJim.convertToFloat32())

    denominator[denominator==0] = 1
    ndvi = numerator / denominator
    ndvi[denominator==0] = 1

    return _pj.Jim(ndvi)
    # TODO: Check NODATA values


class _PixOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def isEqual(self, other):
        if isinstance(other, _pj.Jim):
            return self._jim_object.isEqual(other)
        else:
            return False

    def setData(self, value, ulx, uly, lrx, lry, bands=0, dx=None, dy=None, geo=True):
        """Set range of pixels to value.

        :param jim_object: a Jim object
        :param value: new value for pixels of Jim object
        :param ulx: upper left corner x coordinate (in projected coordinates if geo is True, else in image coordinates)
        :param uly: upper left corner y coordinate (in projected coordinates if geo is True, else in image coordinates)
        :param lrx: lower right corner x coordinate (in projected coordinates if geo is True, else in image coordinates)
        :param lry: lower right corner y coordinate (in projected coordinates if geo is True, else in image coordinates)
        :param bands: List of band indices to crop (index is 0 based)
        :param dx: spatial resolution in x to crop (stride if geo is False)
        :param dy: spatial resolution in y to crop (stride if geo is False)
        :param geo: use geospatial coordinates if True, image coordinates if False
        :return: a Jim object
        """
        if not dx:
            dx=0
        if not dy:
            dy=0
        for band in bands:
            self._jim_object.setData(value, ulx, uly, lrx, lry, band, dx, dy, geo)

    # def setData(self, value, bands=0):
    #     """Set all pixels to value.

    #     :param jim_object: a Jim object
    #     :param value: new value for pixels of Jim object
    #     :param band: List of band indices to crop (index is 0 based)
    #     :return: a Jim object
    #     """
    #     for band in bands:
    #         self._jim_object.setData(value,band)

    def convert(self, **kwargs):
        """Convert Jim image with respect to data type.

        :param kwargs: See table below
        Modifies the instance on which the method was called.


        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | otype            | Data type for output image                                                      |
        +------------------+---------------------------------------------------------------------------------+
        | scale            | Scale output: output=scale*input+offset                                         |
        +------------------+---------------------------------------------------------------------------------+
        | offset           | Apply offset: output=scale*input+offset                                         |
        +------------------+---------------------------------------------------------------------------------+
        | autoscale        | Scale output to min and max, e.g., [0,255]                                      |
        +------------------+---------------------------------------------------------------------------------+
        | a_srs            | Override the projection for the output file                                     |
        +------------------+---------------------------------------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list of :ref:`mask <extract_mask>` key values:

        Example:

        Convert data type of input image to byte using autoscale::

        jim0=jl.io.createJim('/path/to/raster.tif')
        jim0.convert(otype=Byte,autoscale=[0,255])

        Clip raster dataset between 0 and 255 (set all other values to 0), then convert data type to byte::

        jim1=jl.io.createJim('/path/to/raster.tif')
        jim1.setThreshold(min=0,max=255,nodata=0)
        jim1.convert({'otype':'Byte'})
        """
        self._jim_object._set(self._jim_object.convert(kwargs))

    def setThreshold(self, **kwargs):
        """Apply minimum and maximum threshold to pixel values in raster dataset.

        :param kwargs: See table below
        Modifies the instance on which the method was called.


        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | min              | Minimum threshold value (if pixel value < min set pixel value to no data)       |
        +------------------+---------------------------------------------------------------------------------+
        | max              | Maximum threshold value (if pixel value < max set pixel value to no data)       |
        +------------------+---------------------------------------------------------------------------------+
        | value            | value to be set if within min and max                                           |
        |                  | (if not set, valid pixels will remain their input value)                        |
        +------------------+---------------------------------------------------------------------------------+
        | abs              | Set to True to perform threshold test to absolute pixel values                  |
        +------------------+---------------------------------------------------------------------------------+
        | nodata           | Set pixel value to this no data if pixel value < min or > max                   |
        +------------------+---------------------------------------------------------------------------------+

        Example:

        Mask all values not within [0,250] and set to 255 (no data)::

        jim_threshold=jim.setThreshold(min=0,max=250,nodata=255)
        """
        self._jim_object._set(self._jim_object.setThreshold(kwargs))

    def NDVI(self, redBand, nirBand):
        """Compute NDVI on the Jim object.

        Modifies the instance on which the method was called.

        :param redBand: index of band with values of red
        :param nirBand: index of band with values of NIR
        """
        red = _pj.geometry.cropBand(self._jim_object, redBand)
        nir = _pj.geometry.cropBand(self._jim_object, nirBand)

        numerator = _pj.Jim(nir.convertToFloat32()) - \
                    _pj.Jim(red.convertToFloat32())
        denominator = _pj.Jim(nir.convertToFloat32()) + \
                      _pj.Jim(red.convertToFloat32())

        denominator[denominator==0] = 1
        ndvi = numerator / denominator
        ndvi[denominator==0] = 1

        self._jim_object._set(ndvi)
        # TODO: Check NODATA values
