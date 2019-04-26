"""Module for neighbourhood operations."""

import pyjeo as _pj
import numpy
from scipy import signal
from scipy import ndimage

def filter1d(jim_object, filter, dz=None, pad=None, otype=None, **kwargs):
    """Subset raster dataset in spectral/temporal domain.

    Filter Jim object in spectral/temporal domain performed on multi-band
    raster dataset.

    :param jim_object: a Jim object
    :param filter: filter function (see values for different filter
        types :ref:`supported filters <filters1d>`)
    :param dz: filter kernel size in z (spectral/temporal dimension), must be
        odd (example: 3)
    :param pad: Padding method for filtering (how to handle edge effects).
        Possible values are: symmetric (default), replicate, circular,
        zero (pad with 0)
    :param otype: Data type for output image
    :return: filtered Jim object

    see the corresponding method :py:meth:`.filter1d` for more information
    """
    kwargs.update({'filter': filter})
    if dz:
        kwargs.update({'dz': dz})
    if pad:
        kwargs.update({'pad': pad})
    if otype:
        kwargs.update({'otype': otype})
    return _pj.Jim(jim_object._jipjim.filter1d(kwargs))


def filter2d(jim_object, filter, **kwargs):
    """Subset raster dataset in spectral/temporal domain.

    Filter Jim object in spatial domain performed on single or multi-band
    raster dataset.

    :param jim_object: a Jim object
    :param filter: filter function (see values for different filter
        types :ref:`supported filters <filters2d>`)
    :param dx: filter kernel size in x, use odd values only (default is 3)
    :param dy: filter kernel size in y, use odd values only (default is 3)
    :param pad: Padding method for filtering (how to handle edge effects).
        Possible values are: symmetric (default), replicate, circular,
        zero (pad with 0)
    :param otype: Data type for output image
    :return: filtered Jim object

    see the corresponding method :py:meth:`.filter2d` for more information
    """
    if isinstance(filter,numpy.ndarray ):
        kwargs.update({'dy': filter.shape[0]})
        kwargs.update({'dx': filter.shape[1]})
        kwargs.update({'tap':filter.flatten().tolist()})
    else:
        kwargs.update({'filter': filter})
    return _pj.Jim(jim_object._jipjim.filter2d(kwargs))


def morphoDilateDiamond(jim_object, ox=1, oy=1):
    """Output the dilation of im using the elementary diamond shaped SE

    :param ox: integer for origin along x-axis of SE (default is 1 for centred)
    :param oy: integer for origin along y-axis of SE (default is 1 for centred)
    """
    return _pj.Jim(jim_object._jipjim.morphoDilateNgb4(ox, oy))


def morphoErodeDiamond(jim_object, ox=1, oy=1):
    """Output the erosion of im using the elementary diamond shaped SE

    :param ox: integer for origin along x-axis of SE (default is 1 for centred)
    :param oy: integer for origin along y-axis of SE (default is 1 for centred)
    """
    return _pj.Jim(jim_object._jipjim.morphoErodeNgb4(ox, oy))


def morphoErodeLine(jim_object, dx, dy, k, o, type=0):
    """Output the erosion of im using the line SE with slope dy/dx, length k, origin o, and line type (see details at  :cite:`soille-breen-jones96`)

    :param jim_object: image on which to perform the erosion
    :param dx: integer for displacement along x-axis to set slope
    :param dy: integer for displacement along y-axis to set slope
    :param k: integer for number of pixels of line SE
    :param o: integer for origin of SE
    :param type: integer for line type (0 for plain and 1 for periodic).  0 is the default value
    """
    return _pj.Jim(jim_object._jipjim.morphoErodeLine(dx, dy, k, o, type))


def morphoDilateLine(jim_object, dx, dy, k, o, type=0):
    """Output the dilation of im using the line SE with slope dy/dx, length k, origin o, and line type (see details at  :cite:`soille-breen-jones96`)

    :param jim_object: image on which to perform the dilation
    :param dx: integer for displacement along x-axis to set slope
    :param dy: integer for displacement along y-axis to set slope
    :param k: integer for number of pixels of line SE
    :param o: integer for origin of SE
    :param type: integer for line type (0 for plain and 1 for periodic).  0 is the default value
    """
    return _pj.Jim(jim_object._jipjim.morphoDilateLine(dx, dy, k, o, type))


def morphoErode(jim_object, sec_jim_object, ox, oy, oz, trFlag=0):
    """Output the erosion of im using the SE defined by imse.

    Its origin is set at coordinates (x,y,z). The reflection of the SE
    is considered if trflag equals 1 (no reflection by default). Points of
    the SE are points with a non zero value in imse.

    :param jim_object: image on which to perform the erosion
    :param sec_jim_object: an image node for SE (UCHAR type)
    :param ox: x coordinate
    :param oy: y coordinate
    :param oz: z coordinate
    :param trFlag: optional parameter (0 or 1)
    """
    return _pj.Jim(jim_object._jipjim.morphoErode(sec_jim_object._jipjim,
                                           ox, oy, oz, trFlag))

def morphoDilate(jim_object, sec_jim_object, ox, oy, oz, trFlag=0):
    """Output the dilation of im using the SE defined by imse.

    Its origin is set at coordinates (x,y,z). The reflection of the SE
    is considered if trflag equals 1 (no reflection by default). Points of
    the SE are points with a non zero value in imse.

    :param jim_object: image on which to perform the dilation
    :param sec_jim_object: an image node for SE (UCHAR type)
    :param ox: x coordinate
    :param oy: y coordinate
    :param oz: z coordinate
    :param trFlag: optional parameter (0 or 1)
    """
    return _pj.Jim(jim_object._jipjim.morphoDilate(sec_jim_object._jipjim,
                                           ox, oy, oz, trFlag))

def morphoGradientByErosionDiamond(jim_object):
    """Output the gradient by erosion of im using the elementary diamond shaped SE
    """
    return jim_object - _pj.Jim(jim_object._jipjim.morphoErodeNgb4(1, 1))


def morphoGradientByErosionDiamondFrame(jim_object):
    """Output the gradient by erosion of im using the elementary diamond shaped SE
    """
    jim0 = _pj.Jim(jim_object)
    pixmax = jim0.stats.getStats('max')['max']
    jim0.geometry.imageFrameAdd(1, 1, 1, 1, 0, 0, pixmax)
    jim0.ngbops.morphoErodeDiamond()
    jim0.geometry.imageFrameSubtract(1, 1, 1, 1, 0, 0)
    jim0.pixops.simpleArithOp(jim_object, 16) # 16 for SUBSWAP_op
    return jim0

def morphoGradientByDilationDiamondFrame(jim_object):
    """Output the gradient by erosion of im using the elementary diamond shaped SE
    """
    jim0 = _pj.Jim(jim_object)
    pixmin = jim0.stats.getStats('min')['min']
    jim0.geometry.imageFrameAdd(1, 1, 1, 1, 0, 0, pixmin)
    jim0.ngbops.morphoDilateDiamond()
    jim0.geometry.imageFrameSubtract(1, 1, 1, 1, 0, 0)
    jim0.pixops.simpleArithOp(jim_object, 1) #  for SUB_op
    return jim0

def edgeWeight(jim_object, dir=0, type=0):
    """Computes the weights of the horizontal or vertical edges linking any pair of horizontally or vertically adjacent pixels.

    :param dir:  integer for coding edge direction (horizontal if 0, vertical otherwise).
    :param type: integer determining how the edge weights are computed: 0 for absolute difference (default), 1 for maximum value, 2 for minimum value.
    """
    return _pj.Jim(jim_object._jipjim.edgeWeight(dir, type))

def getDissim(jim_object_list, dissimType=0):
    """Compute the dissimilarities between horizontal and vertical pairs of adjacent pixels.

    :param jim_object_list: a list of grey level Jim objects with the same definition domain.  The dissimilarities are calculated for each image separately and composed using the point-wise maximum rule.
    :param dissimType: integer value indicating the type of dissimilarity measure
                       0 (default) for absolute difference
                       1 for dissimilarity measure countering the chaining effect as described in :cite:`soille2011ismm`
    :return: a list of 2 Jim objects holding the horizontal and vertical dissimilarities respectively
    """
    DIR_HORI    = 0
    DIR_VERT    = 1
    ABS_DIFF_op = 0
    MAX_op      = 1
    MIN_op      = 2

    if dissimType==0:
        h_dissim=_pj.ngbops.edgeWeight(jim_object_list[0], DIR_HORI, ABS_DIFF_op)
        v_dissim=_pj.ngbops.edgeWeight(jim_object_list[0], DIR_VERT, ABS_DIFF_op)

        for im in jim_object_list[1:]:
            h_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_HORI, ABS_DIFF_op))
            v_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_VERT, ABS_DIFF_op))
            
    elif dissimType==1:
        print("toto1")
        mingraderograddil = _pj.pixops.infimum(_pj.ngbops.morphoGradientByDilationDiamondFrame(jim_object_list[0]), _pj.ngbops.morphoGradientByErosionDiamondFrame(jim_object_list[0]))
        print("toto2")
        h_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_HORI, MAX_op)
        print("toto3")
        v_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_VERT, MAX_op)
        print("toto4")
        mingraderograddil = 0
        h_dissim.pixops.supremum(_pj.ngbops.edgeWeight(jim_object_list[0], DIR_HORI, ABS_DIFF_op))
        v_dissim.pixops.supremum(_pj.ngbops.edgeWeight(jim_object_list[0], DIR_VERT, ABS_DIFF_op))
        print("toto5")

        for im in jim_object_list[1:]:
            mingraderograddil = _pj.pixops.infimum(_pj.ngbops.morphoGradientByDilationDiamondFrame(im), _pj.ngbops.morphoGradientByErosionDiamondFrame(im))
            h_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil, DIR_HORI, MAX_op)
            v_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil, DIR_VERT, MAX_op)
            mingraderograddil = 0
            h_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_HORI, ABS_DIFF_op))
            v_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_VERT, ABS_DIFF_op))
           
            h_dissim.pixops.supremum(h_dissim_crt)
            v_dissim.pixops.supremum(v_dissim_crt)

    return [h_dissim, v_dissim]
        

class _NgbOps():
    """Define all NgbOps methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def filter1d(self, filter, dz=None, pad=None, otype=None, **kwargs):
        """Subset raster dataset in spectral/temporal domain.

        Filter Jim image in spectral/temporal domain performed on multi-band
        raster dataset.

        :param filter: filter function (see values for different filter
            types :ref:`supported filters <filters1d>`)
        :param dz: filter kernel size in z (spectral/temporal dimension), must
            be odd (example: 3)
        :param pad: Padding method for filtering (how to handle edge effects).
            Possible values are: symmetric (default), replicate, circular,
            zero (pad with 0)
        :param otype: Data type for output image

        .. _filters1d:

        :Morphological filters:

        +---------------------+-----------------------------------------------+
        | filter              | description                                   |
        +=====================+===============================================+
        | dilate              | morphological dilation                        |
        +---------------------+-----------------------------------------------+
        | erode               | morphological erosion                         |
        +---------------------+-----------------------------------------------+
        | close               | morpholigical closing (dilate+erode)          |
        +---------------------+-----------------------------------------------+
        | open                | morpholigical opening (erode+dilate)          |
        +---------------------+-----------------------------------------------+

        .. note::

            The morphological filter uses a linear structural element with
            a size defined by the key 'dz'

        Example:

        Perform a morphological dilation with a linear structural element
        of size 5::

            jim_filtered=jim.filter1d('dilate',dz=5)

        :Statistical filters:

        +--------------+------------------------------------------------------+
        | filter       | description                                          |
        +==============+======================================================+
        | smoothnodata | smooth nodata values (set nodata option!)            |
        +--------------+------------------------------------------------------+
        | nvalid       | report number of valid (not nodata) values in window |
        +--------------+------------------------------------------------------+
        | median       | perform a median filter                              |
        +--------------+------------------------------------------------------+
        | var          | calculate variance in window                         |
        +--------------+------------------------------------------------------+
        | min          | calculate minimum in window                          |
        +--------------+------------------------------------------------------+
        | max          | calculate maximum in window                          |
        +--------------+------------------------------------------------------+
        | sum          | calculate sum in window                              |
        +--------------+------------------------------------------------------+
        | mean         | calculate mean in window                             |
        +--------------+------------------------------------------------------+
        | stdev        | calculate standard deviation in window               |
        +--------------+------------------------------------------------------+
        | percentile   | calculate percentile value in window                 |
        +--------------+------------------------------------------------------+
        | proportion   | calculate proportion in window                       |
        +--------------+------------------------------------------------------+

        .. note::
            You can specify the no data value for the smoothnodata filter
            with the extra key 'nodata' and a list of no data values.
            The interpolation type can be set with the key 'interp'
            (check complete list of
            `values <http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html>`_,
            removing the leading "gsl_interp").

        Example:

        Smooth the 0 valued pixel values using a linear interpolation in
        a spectral/temporal neighborhood of 5 bands::

            jim.filter1d('smoothnodata',dz=5,nodata=0,interp='linear')

        :Wavelet filters:

        Perform a wavelet transform (or inverse) in spectral/temporal domain.

        .. note::
            The wavelet coefficients can be positive and negative. If the input
            raster dataset has an unsigned data type, make sure to set
            the output to a signed data type using the key 'otype'.

        You can use set the wavelet family with the key 'family' in
        the dictionary. The following wavelets are supported as values:

        * daubechies
        * daubechies_centered
        * haar
        * haar_centered
        * bspline
        * bspline_centered

        +----------+--------------------------------------+
        | filter   | description                          |
        +==========+======================================+
        | dwt      | discrete wavelet transform           |
        +----------+--------------------------------------+
        | dwti     | discrete inverse wavelet transform   |
        +----------+--------------------------------------+
        | dwt_cut  | DWT approximation in spectral domain |
        +----------+--------------------------------------+

        .. note::
            The filter 'dwt_cut' performs a forward and inverse transform,
            approximating the input signal. The approximation is performed by
            discarding a percentile of the wavelet coefficients that can be
            set with the key 'threshold'. A threshold of 0 (default) retains
            all and a threshold of 50 discards the lower half of the wavelet
            coefficients.

        Example:

        Approximate the multi-temporal raster dataset by discarding the lower
        20 percent of the coefficients after a discrete wavelet transform.
        The input dataset has a Byte data type. We wavelet transform is
        calculated using an Int16 data type. The approximated image is then
        converted to a Byte dataset, making sure all values below 0 and
        above 255 are set to 0::

            jim_multitemp.filter1d('dwt_cut',threshold=20, otype=Int16)
            jim_multitemp[(jim_multitemp<0) | (jim_multitemp>255)]=0
            jim_multitemp.convert(otype='Byte')

        :Hyperspectral filters:

        Hyperspectral filters assume the bands in the input raster dataset
        correspond to contiguous spectral bands. Full width half max (FWHM)
        and spectral response filters are supported. They converts an N band
        input raster dataset to an M (< N) band output raster dataset.

        The full width half max (FWHM) filter expects a list of M center
        wavelenghts and a corresponding list of M FWHM values. The M center
        wavelenghts define the output wavelenghts and must be provided with
        the key 'wavelengthOut'. For the FHWM, use the key 'fwhm' and a list
        of M values. The algorithm needs to know the N wavelenghts that
        correspond to the N bands of the input raster dataset. Use the key
        'wavelengthIn' and a list of N values. The units of input, output and
        FWHM are arbitrary, but should be identical (e.g., nm).

        Example:

        Covert the hyperspectral input raster dataset, with the wavelengths
        defined in the list wavelenghts_in to a multispectral raster dataset
        with three bands, corresponding to Red, Green, and Blue::

            wavelengths_in=[]
            #define the wavelenghts of the input raster dataset

            if len(wavelengths_in) == jim_hyperspectral.nrOfBand():
                jim_hyperspectral.filter1d(wavelengthIn=wavelenghts_in,wavelengthOut=[650,510,475],fwhm=[50,50,50])
            else:
                print("Error: number of input wavelengths must be equal to number of bands in input raster dataset")

        .. note::
                The input wavelenghts are automatically interpolated. You can
                specify the interpolation using the key 'interp' and values as
                listed interpolation
                http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html

        The spectral response filter (SRF)

        The input raster dataset is filtered with M of spectral response
        functions (SRF).  Each spectral response function must be provided by
        the user in an ASCII file that consists of two columns: wavelengths
        and response. Use the key 'srf' and a list of paths to the ASCII
        file(s). The algorithm automatically takes care of the normalization
        of the SRF.

        Example:

        Covert the hyperspectral input raster dataset, to a multispectral
        raster dataset with three bands, corresponding to Red, Green, and Blue
        as defined in the ASCII text files 'srf_red.txt', 'srf_green.txt',
        'srf_blue.txt'::

            wavelengths_in=[]
            #specify the wavelenghts of the input raster dataset

            if len(wavelengths_in) == jim_hyperspectral.nrOfBand():
                rgb=jim_hyperspectral.filter1d(wavelengthIn=wavelenghts_in,srf=['srf_red.txt','srf_green.txt','srf_blue.txt'])
            else:
                print("Error: number of input wavelengths must be equal to number of bands in input raster dataset")

        .. note::
            The input wavelenghts are automatically interpolated. You can
            specify the interpolation using the key 'interp' and values as
            listed interpolation
            http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html


        :Custom filters:

        For the custom filter, you can specify your own taps using the keyword
        'tapz' and a list of filter tap values. The tap values are
        automatically normalized by the algorithm.

        Example:

        Perform a simple smoothing filter by defining three identical tap
        values::

            jim.filter1d(tapz=[1,1,1])
        """
        if dz:
            kwargs.update({'dz': dz})
        if pad:
            kwargs.update({'pad': pad})
        if otype:
            kwargs.update({'otype': otype})
        self._jim_object._set(self._jim_object._jipjim.filter1d(kwargs))

    def filter2d(self, filter, **kwargs):
        """Subset raster dataset in spectral/temporal domain.

        Filter Jim object in spatial domain performed on single or multi-band
        raster dataset.

        :param jim_object: a Jim object
        :param filter: filter function (see values for different filter
            types :ref:`supported filters <filters2d>`)
        :param dx: filter kernel size in x, use odd values only (default is 3)
        :param dy: filter kernel size in y, use odd values only (default is 3)
        :param pad: Padding method for filtering (how to handle edge effects).
            Possible values are: symmetric (default), replicate, circular,
            zero (pad with 0)
        :param otype: Data type for output image

        .. _filters2d:

        **Edge detection**

        +---------------------+-----------------------------------------------+
        | filter              | description                                   |
        +=====================+===============================================+
        | sobelx              | Sobel operator in x direction                 |
        +---------------------+-----------------------------------------------+
        | sobely              | Sobel operator in y direction                 |
        +---------------------+-----------------------------------------------+
        | sobelxy             | Sobel operator in x and y direction           |
        +---------------------+-----------------------------------------------+
        | homog               | binary value indicating if pixel is identical |
        |                     | to all pixels in kernel                       |
        +---------------------+-----------------------------------------------+
        | heterog             | binary value indicating if pixel is different |
        |                     | than all pixels in kernel                     |
        +---------------------+-----------------------------------------------+

        Example:

        Perform Sobel edge detection in both x and direction::

            jim_filtered=jim.filter2d('sobelxy')

        **Morphological filters**

        +---------------------+-----------------------------------------------+
        | filter              | description                                   |
        +=====================+===============================================+
        | dilate              | morphological dilation                        |
        +---------------------+-----------------------------------------------+
        | erode               | morphological erosion                         |
        +---------------------+-----------------------------------------------+
        | close               | morpholigical closing (dilate+erode)          |
        +---------------------+-----------------------------------------------+
        | open                | morpholigical opening (erode+dilate)          |
        +---------------------+-----------------------------------------------+

        .. note::
            You can use the optional key 'class' with a list value to take only
            these pixel values into account. For instance, use 'class':[255] to
            dilate clouds in the raster dataset that have been flagged with
            value 255. In addition, you can use a circular disc kernel (set
            the key 'circular' to True).

        Example:

        Perform a morphological dilation using a circular kernel with size
        (diameter) of 5 pixels::

            jim.filter2d('dilate',dx=5,dy=5,circular=True)


        .. note::
            For a more comprehensive list of morphological operators, please refer
            to the corresponding methods, e.g., :py:meth:`~._NgbOps.morphoDilate`

        **Statistical filters**

        +--------------+------------------------------------------------------+
        | filter       | description                                          |
        +==============+======================================================+
        | smoothnodata | smooth nodata values (set nodata option!)            |
        +--------------+------------------------------------------------------+
        | nvalid       | report number of valid (not nodata) values in window |
        +--------------+------------------------------------------------------+
        | median       | perform a median filter                              |
        +--------------+------------------------------------------------------+
        | var          | calculate variance in window                         |
        +--------------+------------------------------------------------------+
        | min          | calculate minimum in window                          |
        +--------------+------------------------------------------------------+
        | max          | calculate maximum in window                          |
        +--------------+------------------------------------------------------+
        | ismin        | binary value indicating if pixel is minimum in kernel|
        +--------------+------------------------------------------------------+
        | ismax        | binary value indicating if pixel is maximum in kernel|
        +--------------+------------------------------------------------------+
        | sum          | calculate sum in window                              |
        +--------------+------------------------------------------------------+
        | mode         | calculate the mode (only for categorical values)     |
        +--------------+------------------------------------------------------+
        | mean         | calculate mean in window                             |
        +--------------+------------------------------------------------------+
        | stdev        | calculate standard deviation in window               |
        +--------------+------------------------------------------------------+
        | percentile   | calculate percentile value in window                 |
        +--------------+------------------------------------------------------+
        | proportion   | calculate proportion in window                       |
        +--------------+------------------------------------------------------+

        .. note::
            You can specify the no data value for the smoothnodata filter with
            the extra key 'nodata' and a list of no data values.
            The interpolation type can be set with the key 'interp' (check
            complete list of
            `values <http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html>`_,
            removing the leading "gsl_interp").

        Example:

        Perform a median filter with kernel size of 3x3 pixels::

            jim.filter2d('median',dx=5, dy=5)

        **Wavelet filters**

        Perform a wavelet transform (or inverse) in spatial domain.

        .. note::
            The wavelet coefficients can be positive and negative. If the input
            raster dataset has an unsigned data type, make sure to set
            the output to a signed data type using the key 'otype'.

        You can use set the wavelet family with the key 'family' in
        the dictionary. The following wavelets are supported as values:

        * daubechies
        * daubechies_centered
        * haar
        * haar_centered
        * bspline
        * bspline_centered

        +----------+--------------------------------------+
        | filter   | description                          |
        +==========+======================================+
        | dwt      | discrete wavelet transform           |
        +----------+--------------------------------------+
        | dwti     | discrete inverse wavelet transform   |
        +----------+--------------------------------------+
        | dwt_cut  | DWT approximation in spectral domain |
        +----------+--------------------------------------+

        .. note::
            The filter 'dwt_cut' performs a forward and inverse transform,
            approximating the input signal. The approximation is performed by
            discarding a percentile of the wavelet coefficients that can be
            set with the key 'threshold'. A threshold of 0 (default) retains
            all and a threshold of 50 discards the lower half of the wavelet
            coefficients.

        Example:

        Approximate the multi-temporal raster dataset by discarding the lower
        20 percent of the coefficients after a discrete wavelet transform.
        The input dataset has a Byte data type. We wavelet transform is
        calculated using an Int16 data type. The approximated image is then
        converted to a Byte dataset, making sure all values below 0 and above
        255 are set to 0::

            jim_multitemp.filter2d('dwt_cut',threshold=20, otype=Int16)
            jim_multitemp[(jim_multitemp<0) | (jim_multitemp>255)]=0
            jim_multitemp.convert(otype='Byte')
        """
        if isinstance(filter,numpy.ndarray ):
            kwargs.update({'dy': filter.shape[0]})
            kwargs.update({'dx': filter.shape[1]})
            kwargs.update({'tap':filter.flatten().tolist()})
        else:
            kwargs.update({'filter': filter})
        self._jim_object._set(self._jim_object._jipjim.filter2d(kwargs))

    def morphoErodeDiamond(self, ox=1, oy=1):
        """Output the erosion of im using the elementary diamond shaped SE

        Its origin is set at coordinates (x,y).

        :param ox: x coordinate
        :param oy: y coordinate
        """
        self._jim_object._jipjim.d_morphoErodeNgb4(ox, oy)

    def morphoDilateDiamond(self, ox=1, oy=1):
        """Output the dilation of im using the elementary diamond shaped SE

        Its origin is set at coordinates (x,y).

        :param ox: x coordinate
        :param oy: y coordinate
        """
        self._jim_object._jipjim.d_morphoDilateNgb4(ox, oy)

    def morphoErode(self, sec_jim_object, ox, oy, oz, trFlag=0):
        """Output the erosion of im using the SE defined by imse.

        Its origin is set at coordinates (x,y,z). The reflection of the SE
        is considered if trflag equals 1 (no reflection by default). Points of
        the SE are points with a non zero value in imse.

        :param sec_jim_object: an image node for SE (UCHAR type)
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate
        :param trFlag: optional parameter (0 or 1)
        """
        self._jim_object._jipjim.morphoErode(sec_jim_object._jipjim,
                                             ox, oy, oz, trFlag)

    def morphoDilate(self, sec_jim_object, ox, oy, oz, trFlag=0):
        """Output the dilation of im using the SE defined by imse.

        Its origin is set at coordinates (x,y,z). The reflection of the SE
        is considered if trflag equals 1 (no reflection by default). Points of
        the SE are points with a non zero value in imse.

        :param sec_jim_object: an image node for SE (UCHAR type)
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate
        :param trFlag: optional parameter (0 or 1)
        """
        self._jim_object._jipjim.morphoDilate(sec_jim_object._jipjim,
                                             ox, oy, oz, trFlag)


class _NgbOpsList():
    """Define all NgbOps methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class _NgbOpsVect():
    """Define all NgbOps methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
