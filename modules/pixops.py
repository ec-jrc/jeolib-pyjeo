import jiplib as _jl
import pyjeo as _pj


def pointOpBitWise(jim_object, sec_jim_object, operation_code):
    """Bitwise operation between two images.

    :param jim_object: a Jim object
    :param sec_jim_object: a Jim object
    :param operation_code: 10 or AND op, 11 or OR op, and 12 or XOR op
    :return: a Jim object
    """
    return _pj.Jim(jim_object.pointOpBitwise(sec_jim_object, operation_code))


def pointOpBlank(jim_object, value):
    """Set all pixels of image to value.

    :param jim_object: a Jim object
    :param value: new value for pixels of Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object.pointOpBlank(value))


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
    _pj.Jim(jim_object.convert(kwargs))


class _PixOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def equal(self, aJim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            anp=_jl.jim2np(aJim)
            selfnp=np.uint8(1)*(selfnp==anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim

    def pointOpBitWise(self, sec_jim_object, operation_code):
        """Bitwise operation between two images.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param operation_code: 10 or AND op, 11 or OR op, and 12 or XOR op
        """
        self._jim_object._set(self._jim_object.pointOpBitwise(sec_jim_object,
                                                              operation_code))

    def pointOpArithCst(self, value, operation_code):
        """Bitwise operation between two images.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param operation_code: todo
        """
        self._jim_object.d_pointOpArithCst(value,operation_code)

    def pointOpArith(self, sec_jim_object, operation_code):
        """Bitwise operation between two images.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param operation_code: todo
        """
        self._jim_object.d_pointOpArith(sec_jim_object,operation_code)

    def pointOpBlank(self, value):
        """Set all pixels of image to value.

        Modifies the instance on which the method was called.

        :param value: new value for pixels of Jim object
        """
        self._jim_object.d_pointOpBlank(value)

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
