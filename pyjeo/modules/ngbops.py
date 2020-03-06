"""Module for neighbourhood operations."""

import numpy as _np

import pyjeo as _pj


def dwt1d(jim_object, wavelet=None, family=None, **kwargs):
    """Compute discrete forward wavelet transform in time-spectral domain.

    :param jim_object: a Jim object of data type GDT_Float64
    :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
        haar_centered, bspline, bspline_centered
    :param family: wavelet family,
        see also https://www.gnu.org/software/gsl/doc/html/dwt.html
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        dwt = pj.ngbops.dwt1d(jim)
        dwt.ngbops.dwti1d()
    """
    if wavelet is not None:
        kwargs.update({'wavelet': wavelet})

    if family is not None:
        kwargs.update({'family': family})

    return _pj.Jim(jim_object._jipjim.dwt1d(kwargs))


def dwt2d(jim_object, wavelet=None, family=None, **kwargs):
    """Compute forward discrete wavelet transform in the spatial domain.

    :param jim_object: a Jim object of data type GDT_Float64
    :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
        haar_centered, bspline, bspline_centered
    :param family: wavelet family, see also
        https://www.gnu.org/software/gsl/doc/html/dwt.html
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        dwt = pj.ngbops.dwt2d(jim)
        dwt.ngbops.dwti2d()
    """
    if wavelet is not None:
        kwargs.update({'wavelet': wavelet})

    if family is not None:
        kwargs.update({'family': family})

    return _pj.Jim(jim_object._jipjim.dwt2d(kwargs))


def dwti1d(jim_object, wavelet=None, family=None, **kwargs):
    """Compute inverse discrete wavelet transform in time-spectral domain.

    :param jim_object: a Jim object of data type GDT_Float64
    :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
        haar_centered, bspline, bspline_centered
    :param family: wavelet family, see also
        https://www.gnu.org/software/gsl/doc/html/dwt.html
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        dwt = pj.ngbops.dwt1d(jim)
        dwt.ngbops.dwti1d()

    Approximate a 3D image by setting all wavelet coefficients below
    some percentile value (e.g., 10) to 0::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        jim.ngbops.dwt1d()
        jimabs = pj.Jim(jim)
        jimabs = abs(jimabs)
        thresholds = np.percentile(jimabs.np(), 90, axis=0)
        jim[jimabs < thresholds] = 0
        jim.ngbops.dwti1d()
    """
    if wavelet is not None:
        kwargs.update({'wavelet': wavelet})

    if family is not None:
        kwargs.update({'family': family})

    return _pj.Jim(jim_object._jipjim.dwt1d(kwargs))


def dwti2d(jim_object, wavelet=None, family=None, **kwargs):
    """Compute inverse discrete wavelet transform in the spatial domain.

    :param jim_object: a Jim object of data type GDT_Float64
    :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
        haar_centered, bspline, bspline_centered
    :param family: wavelet family, see also
        https://www.gnu.org/software/gsl/doc/html/dwt.html
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        dwt = pj.ngbops.dwt2d(jim)
        dwt.ngbops.dwti2d()
    """
    if wavelet is not None:
        kwargs.update({'wavelet': wavelet})

    if family is not None:
        kwargs.update({'family': family})

    return _pj.Jim(jim_object._jipjim.dwt2d(kwargs))


def edgeWeight(jim_object, dir=0, type=0):
    """Compute the weights of the horizontal or vertical edges.

    Linking any pair of horizontally or vertically adjacent pixels.

    :param jim_object: Jim object on which to perform the computation
    :param dir:  integer for coding edge direction
        (horizontal if 0, vertical otherwise).
    :param type: integer determining how the edge weights are computed:
        0 for absolute difference (default),
        1 for maximum value,
        2 for minimum value.
    """
    return _pj.Jim(jim_object._jipjim.edgeWeight(dir, type))


def filter1d(jim_object, filter, dz=None, pad='symmetric', otype=None,
             **kwargs):
    """Subset raster dataset in spectral/temporal domain.

    This function is deprecated

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
    kwargs.update({'pad': pad})

    if dz:
        kwargs.update({'dz': dz})
    if otype:
        kwargs.update({'otype': otype})

    return _pj.Jim(jim_object._jipjim.filter1d(kwargs))


def filter2d(jim_object, filter, dx=3, dy=3, pad='symmetric', otype=None,
             **kwargs):
    """Subset raster dataset in spectral/temporal domain.

    #This function is deprecated

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
    kwargs.update({'dx': dx, 'dy': dy, 'pad': pad})

    if otype is not None:
        kwargs.update({'otype': otype})

    if isinstance(filter, _np.ndarray):
        taps = kwargs.pop('filter')
        kwargs.update({'taps': taps})
        return firfilter2d(jim_object, kwargs)
    else:
        kwargs.update({'filter': filter})
    return _pj.Jim(jim_object._jipjim.filter2d(kwargs))


def firfilter1d(jim_object, taps, pad='symmetric', **kwargs):
    """Compute the finite impulse response filter in time-spectral domain.

    :param jim_object: a Jim object
        (the same data type will be used for output)
    :param taps: 1D array of filter taps
    :param pad: Padding method for filtering (how to handle edge effects).
        Choose between: symmetric, replicate, circular, zero (pad with 0)
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/image.tif')

        filtered = pj.ngbops.firfilter1d(jim, taps=[1, 2, 1], pad='symmetric')
    """
    if len(taps.shape) != 1:
        raise ValueError('Error: taps should be 1D array')

    taps = _np.array(taps).tolist()
    kwargs.update({'taps': taps})
    kwargs.update({'pad': pad})
    return _pj.Jim(jim_object._jipjim.firfilter1d(kwargs))


def firfilter2d(jim_object, taps, nodata=None, norm=None, **kwargs):
    """Compute the finite impulse response filter in spatial domain.

    :param jim_object: a Jim object
        (the same data type will be used for output)
    :param taps: 2D array of filter taps
    :param nodata: list of no data values not to take into account when
        calculating the filter response
    :param norm: normalize tap values
    :return: filtered Jim object

    Example::

        jim = pj.Jim('/path/to/image.tif')

        filtered = pj.ngbops.firfilter2d(jim, taps=[1, 2, 1], norm=True)
    """
    if len(taps.shape) != 2:
        raise ValueError('Error: taps should be 2D array')

    taps = _np.array(taps)
    kwargs.update({'taps': taps.flatten().tolist()})
    kwargs.update({'dimx': taps.shape[1]})
    kwargs.update({'dimy': taps.shape[0]})

    if nodata is not None:
        kwargs.update({'nodata': nodata})
    if norm is not None:
        kwargs.update({'norm': norm})

    return _pj.Jim(jim_object._jipjim.firfilter2d(kwargs))


def getDissim(jimo, dissim_type=0):
    """Compute the dissimilarities.

    Compute the dissimilarities between horizontal and vertical pairs of
    adjacent pixels.

    :param jimo: a list of grey level Jim objects with the same
        definition domain. The dissimilarities are calculated for each
        image separately and composed using the point-wise maximum rule.
    :param dissim_type: integer value indicating the type of dissimilarity \
                            measure.
                        0 (default) for absolute difference.
                        1 for dissimilarity measure countering \
                          the chaining effect as described in \
                          :cite:`soille2011ismm`
    :return: a list of 2 Jim objects holding the horizontal and vertical
        dissimilarities respectively
    """
    DIR_HORI = 0
    DIR_VERT = 1
    ABS_DIFF_op = 0
    MAX_op = 1
    MIN_op = 2

    if isinstance(jimo, _pj.Jim):
        jim_object_list = _pj.JimList([jimo])
    else:
        jim_object_list = jimo

    if dissim_type == 0:
        # TODO: Check if looping through everything with edgeWeight() and
        #       then one call of supremum is not faster
        h_dissim = _pj.ngbops.edgeWeight(jim_object_list[0], DIR_HORI,
                                         ABS_DIFF_op)
        v_dissim = _pj.ngbops.edgeWeight(jim_object_list[0], DIR_VERT,
                                         ABS_DIFF_op)

        for im in jim_object_list[1:]:
            h_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_HORI,
                                                           ABS_DIFF_op))
            v_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_VERT,
                                                           ABS_DIFF_op))
    elif dissim_type == 1:
        mingraderograddil = _pj.pixops.infimum(
            _pj.ngbops.morphoGradientByDilationDiamond(
                jim_object_list[0]),
            _pj.ngbops.morphoGradientByErosionDiamond(jim_object_list[0]))
        h_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_HORI, MAX_op)
        v_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_VERT, MAX_op)
        mingraderograddil = 0
        h_dissim.pixops.supremum(
            _pj.ngbops.edgeWeight(jim_object_list[0], DIR_HORI, ABS_DIFF_op))
        v_dissim.pixops.supremum(
            _pj.ngbops.edgeWeight(jim_object_list[0], DIR_VERT, ABS_DIFF_op))

        for im in jim_object_list[1:]:
            mingraderograddil = _pj.pixops.infimum(
                _pj.ngbops.morphoGradientByDilationDiamond(im),
                _pj.ngbops.morphoGradientByErosionDiamond(im))
            h_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil, DIR_HORI,
                                                 MAX_op)
            v_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil, DIR_VERT,
                                                 MAX_op)
            mingraderograddil = 0
            h_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_HORI,
                                                               ABS_DIFF_op))
            v_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_VERT,
                                                               ABS_DIFF_op))

            h_dissim.pixops.supremum(h_dissim_crt)
            v_dissim.pixops.supremum(v_dissim_crt)
    else:
        raise ValueError('dissim_type {} not supported. Supported only values '
                         '0 and 1'.format(dissim_type))

    return [h_dissim, v_dissim]


def labelPixNgb(jim_object, sec_jim_object, ox, oy, oz):
    """Label pix ngb.

    :param jim_object: Jim object on which to perform the labelling
    :param sec_jim_object: a Jim object
    :param ox: x coordinate
    :param oy: y coordinate
    :param oz: z coordinate

    :return: labelled Jim object
    """
    return _pj.Jim(jim_object._jipjim.labelPixNgb(sec_jim_object._jipjim,
                                                  ox, oy, oz))


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


def morphoDilateDiamond(jim_object):
    """Output the dilation of im using the elementary diamond shaped SE.

    :param jim_object: a Jim object
    """
    ox = 1
    oy = 1
    return _pj.Jim(jim_object._jipjim.morphoDilateNgb4(ox, oy))


def morphoDilateLine(jim_object, dx, dy, k, o, type=0):
    """Output the dilation of im.

    Uses the line SE with slope dy/dx, length k, origin o, and line type
    (see details at  :cite:`soille-breen-jones96`).

    :param jim_object: image on which to perform the dilation
    :param dx: integer for displacement along x-axis to set slope
    :param dy: integer for displacement along y-axis to set slope
    :param k: integer for number of pixels of line SE
        (must be odd; if not, it is extended by one pixel)
    :param o: integer for origin of SE
    :param type: integer for line type (0 for plain and 1 for periodic).
        0 is the default value
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


def morphoErodeDiamond(jim_object):
    """Output the erosion of im using the elementary diamond shaped SE.

    :param jim_object: a Jim object
    """
    ox = 1
    oy = 1
    return _pj.Jim(jim_object._jipjim.morphoErodeNgb4(ox, oy))


def morphoErodeLine(jim_object, dx, dy, k, o, type=0):
    """Output the erosion of im using a line segment.

    Uses the line SE with slope dy/dx, length k, origin o, and line type. See
    details at  :cite:`soille-breen-jones96`.

    :param jim_object: image on which to perform the erosion
    :param dx: integer for displacement along x-axis to set slope
    :param dy: integer for displacement along y-axis to set slope
    :param k: integer for number of pixels of line SE
        (must be odd; if not, it is extended by one pixel)
    :param o: integer for origin of SE
    :param type: integer for line type (0 for plain and 1 for periodic).
        0 is the default value
    """
    return _pj.Jim(jim_object._jipjim.morphoErodeLine(dx, dy, k, o, type))


def morphoGradientByDilationDiamond(jim_object):
    """Output the gradient by dilation of Jim.

    Uses the elementary diamond shaped SE.
    """
    return _pj.Jim(jim_object._jipjim.morphoDilateNgb4(1, 1)) - jim_object


def morphoGradientByErosionDiamond(jim_object):
    """Output the gradient by erosion of Jim.

    Uses the elementary diamond shaped SE.
    """
    return jim_object - _pj.Jim(jim_object._jipjim.morphoErodeNgb4(1, 1))


def savgolay(jim_object, **kwargs):
    """Compute the Savitzky-Golay filter in the time-spectral domain.

    :param jim_object: a Jim object of data type GDT_Float64
    :param nl: Number of leftward (past) data points used in Savitzky-
        Golay filter)
    :param nr: Number of rightward (future) data points used in Savitzky-
        Golay filter)
    :param ld: order of the derivative desired in Savitzky-Golay filter
        (e.g., ld=0 for smoothed function)
    :param m: order of the smoothing polynomial in Savitzky-Golay filter,
        also equal to the highest conserved moment; usual values are
        m=2 or m=4)
    :return: filtered Jim object

    Example:

    Perform a Savitzky-Golay filter to reconstruct a time series data set
    as in `J. Chen 2004 <https://doi.org/10.1016/j.rse.2004.03.014>`_::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        savgol = pj.ngbops.savgolay(jim, nl=7, nr=7, m=2, pad='replicate')
        for loop in range(0, 10):
            savgol.pixops.convert(otype=jim.properties.getDataType())
            savgol[savgol < jim] = jim
            savgol = pj.ngbops.savgolay(savgol, nl=4, nr=4, m=6,
                                        pad='replicate')
    """
    return _pj.Jim(jim_object._jipjim.savgolay(kwargs))


def smoothNoData1d(jim_object, nodata=0, **kwargs):
    """Smooth nodata in spectral/temporal domain.

    :param jim_object: input Jim object
    :param nodata: no data value to interpolate
    :param interpolationType: type of interpolation for spectral filtering
        (see https://www.gnu.org/software/gsl/doc/html/interp.html)

    Example::

        jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
        jim.pixops.convert('GDT_Float64')

        pj.ngbops.smoothNoData1d(jim, 0)
    """
    kwargs.update({'nodata': nodata})
    return _pj.Jim(jim_object._jipjim.smoothNoData1d(kwargs))


# works but slower than directly in Numpy:  jim.np()[0]=jim.np().min(axis=0)
# def stats1d(jim_object, method, nodata=0, **kwargs):
#     """Calculate statistics for each pixel in the spectral/temporal domain

#     :param jim_object: input Jim object
#     :param method: (list of) methods to calculate
#     :param nodata: no data value to interpolate
#     :return: Jim object with statistics

#     Example:

#     jim=pj.Jim('/path/to/multi-band/image.tif',band2plane=True)

#     jimstats=pj.ngbops.stats1d(jim,method=['min','max])
#     """

#     kwargs.update({'method': method})
#     kwargs.update({'nodata':nodata})
#     print("method is: {}".format(method))
#     print("kwargs is: {}".format(kwargs))
#     return _pj.Jim(jim_object._jipjim.stats1d(kwargs))

class _NgbOps(_pj.modules.JimModuleBase):
    """Define all NgbOps methods."""

    def dwt1d(self, wavelet=None, family=None, **kwargs):
        """Compute discrete forward wavelet transform in time-spectral domain.

        :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
            haar_centered, bspline, bspline_centered
        :param family: wavelet family, see also
            https://www.gnu.org/software/gsl/doc/html/dwt.html

        Example::

            jim=pj.Jim('/path/to/multi-band/image.tif',band2plane=True)
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.dwt1d()
            jim.ngbops.dwti1d()
        """
        if wavelet is not None:
            kwargs.update({'wavelet': wavelet})

        if family is not None:
            kwargs.update({'family': family})

        self._jim_object._jipjim.d_dwt1d(kwargs)

    def dwt2d(self, wavelet=None, family=None, **kwargs):
        """Compute forward discrete wavelet transform in the spatial domain.

        :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
            haar_centered, bspline, bspline_centered
        :param family: wavelet family, see also
            https://www.gnu.org/software/gsl/doc/html/dwt.html

        Example::

            jim=pj.Jim('/path/to/image.tif')
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.dwt2d()
            jim.ngbops.dwti2d()
        """
        if wavelet is not None:
            kwargs.update({'wavelet': wavelet})

        if family is not None:
            kwargs.update({'family': family})

        self._jim_object._jipjim.d_dwt2d(kwargs)

    def dwti1d(self, wavelet=None, family=None, **kwargs):
        """Compute inverse discrete wavelet transform in time-spectral domain.

        :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
            haar_centered, bspline, bspline_centered
        :param family: wavelet family, see also
            https://www.gnu.org/software/gsl/doc/html/dwt.html

        Example::

            jim=pj.Jim('/path/to/multi-band/image.tif',band2plane=True)
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.dwt1d()
            jim.ngbops.dwti1d()

        Approximate a 3D image by setting all wavelet coefficients below
        some percentile value (e.g., 10) to 0::

            jim=pj.Jim('/path/to/multi-band/image.tif',band2plane=True)
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.dwt1d()
            jimabs=pj.Jim(jim)
            jimabs=abs(jimabs)
            thresholds=np.percentile(jimabs.np(),90,axis=0)
            jim[jimabs<thresholds]=0
            jim.ngbops.dwti1d()
        """
        if wavelet is not None:
            kwargs.update({'wavelet': wavelet})

        if family is not None:
            kwargs.update({'family': family})

        self._jim_object._jipjim.d_dwti1d(kwargs)

    def dwti2d(self, wavelet=None, family=None, **kwargs):
        """Compute inverse discrete wavelet transform in the spatial domain.

        :param wavelet: wavelet type: daubechies,daubechies_centered, haar,
            haar_centered, bspline, bspline_centered
        :param family: wavelet family, see also
            https://www.gnu.org/software/gsl/doc/html/dwt.html

        Example::

            jim=pj.Jim('/path/to/image.tif')
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.dwt2d()
            jim.ngbops.dwti2d()
        """
        if wavelet is not None:
            kwargs.update({'wavelet': wavelet})

        if family is not None:
            kwargs.update({'family': family})

        self._jim_object._jipjim.d_dwti2d(kwargs)

    def edgeWeight(self, dir=0, type=0):
        """Compute the weights of the horizontal or vertical edges.

        Linking any pair of horizontally or vertically adjacent pixels.

        Modifies the instance on which the method was called.

        :param dir:  integer for coding edge direction
            (horizontal if 0, vertical otherwise).
        :param type: integer determining how the edge weights are computed:
            0 for absolute difference (default),
            1 for maximum value,
            2 for minimum value.
        """
        self._jim_object._set(self._jim_object._jipjim.edgeWeight(dir, type))

    def filter1d(self, filter, dz=None, pad='symmetric', otype=None, **kwargs):
        """Subset raster dataset in spectral/temporal domain.

        This function is deprecated

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

            jim_filtered = jim.ngbops.filter1d('dilate', dz=5)

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
            `values <https://www.gnu.org/software/gsl/doc/html/interp.html>`_,
            removing the leading "gsl_interp").

        Example:

        Smooth the 0 valued pixel values using a linear interpolation in
        a spectral/temporal neighborhood of 5 bands::

            jim.ngbops.filter1d('smoothnodata', dz=5, nodata=0,
                                interp='linear')

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

            jim_multitemp.ngbops.filter1d('dwt_cut', threshold=20, otype=Int16)
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

            wavelengths_in = []
            # define the wavelenghts of the input raster dataset

            if len(wavelengths_in) == jim_hyperspectral.nrOfBand():
                jim_hyperspectral.ngbops.filter1d(
                    wavelengthIn=wavelenghts_in,
                    wavelengthOut=[650, 510, 475],
                    fwhm=[50, 50, 50])
            else:
                print("Error: number of input wavelengths must be equal to "
                      "number of bands in input raster dataset")

        .. note::
                The input wavelenghts are automatically interpolated. You can
                specify the interpolation using the key 'interp' and values as
                listed interpolation
                https://www.gnu.org/software/gsl/doc/html/interp.html

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

            wavelengths_in = []
            # specify the wavelenghts of the input raster dataset

            if len(wavelengths_in) == jim_hyperspectral.nrOfBand():
                rgb=jim_hyperspectral.ngbops.filter1d(
                wavelengthIn=wavelenghts_in,
                srf=['srf_red.txt', 'srf_green.txt', 'srf_blue.txt'])
            else:
                print("Error: number of input wavelengths must be equal to "
                      "number of bands in input raster dataset")

        .. note::
            The input wavelenghts are automatically interpolated. You can
            specify the interpolation using the key 'interp' and values as
            listed interpolation
            https://www.gnu.org/software/gsl/doc/html/interp.html


        :Custom filters:

        For the custom filter, you can specify your own taps using the keyword
        'tapz' and a list of filter tap values. The tap values are
        automatically normalized by the algorithm.

        Example:

        Perform a simple smoothing filter by defining three identical tap
        values::

            jim.ngbops.filter1d(tapz=[1, 1, 1])
        """
        kwargs.update({'filter': filter})
        kwargs.update({'pad': pad})

        if dz:
            kwargs.update({'dz': dz})
        if otype:
            kwargs.update({'otype': otype})

        self._jim_object._set(self._jim_object._jipjim.filter1d(kwargs))

    def filter2d(self, filter, dx=3, dy=3, pad='symmetric', otype=None,
                 **kwargs):
        """Subset raster dataset in spectral/temporal domain.

        This function is deprecated

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

            jim_filtered = jim.ngbops.filter2d('sobelxy')

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

            jim.ngbops.filter2d('dilate', dx=5, dy=5, circular=True)


        .. note::
            For a more comprehensive list of morphological operators, please
            refer to the corresponding methods, e.g.,
            :py:meth:`~._NgbOps.morphoDilate`

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
            `values <https://www.gnu.org/software/gsl/doc/html/interp.html>`_,
            removing the leading "gsl_interp").

        Example:

        Perform a median filter with kernel size of 3x3 pixels::

            jim.ngbops.filter2d('median', dx=5, dy=5)

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

            jim_multitemp.ngbops.filter2d('dwt_cut', threshold=20, otype=Int16)
            jim_multitemp[(jim_multitemp < 0) | (jim_multitemp > 255)] = 0
            jim_multitemp.convert(otype='Byte')
        """
        kwargs.update({'dx': dx, 'dy': dy, 'pad': pad})

        if otype is not None:
            kwargs.update({'otype': otype})

        if isinstance(filter, _np.ndarray):
            taps = kwargs.pop('filter', None)
            kwargs.update({'taps': taps})
            self._jim_object._set(self._jim_object._jipjim.filter2d(kwargs))
            self.firfilter2d(kwargs)
        else:
            kwargs.update({'filter': filter})
            self._jim_object._set(self._jim_object._jipjim.filter2d(kwargs))

    def firfilter1d(self, taps, pad='symmetric', **kwargs):
        """Compute the finite impulse response filter in time-spectral domain.

        :param taps: 1D array of filter taps
        :param pad: Padding method for filtering (how to handle edge effects).
            Choose between: symmetric, replicate, circular, zero (pad with 0)
        :return: filtered Jim object

        Example::

            jim = pj.Jim('/path/to/image.tif')

            jim.firfilter1d(jim, taps=[1, 2, 1], pad='symmetric')
        """
        if len(taps.shape) != 1:
            raise ValueError('Error: taps should be 1D array')

        taps = _np.array(taps).tolist()
        kwargs.update({'taps': taps})
        kwargs.update({'pad': pad})
        self._jim_object._set(self._jim_object._jipjim.firfilter1d(kwargs))

    def firfilter2d(self, taps, nodata=None, norm=None, **kwargs):
        """Compute the finite impulse response filter in spatial domain.

        :param taps: 2D array of filter taps
        :param nodata: list of no data values not to take into account when
            calculating the filter response
        :param norm: normalize tap values
        :return: filtered Jim object

        Example::

            jim = pj.Jim('/path/to/image.tif')

            jim.ngbops.firfilter2d(taps=[1, 2, 1], norm=True, pad='symmetric')
        """
        if len(taps.shape) != 2:
            raise ValueError('Error: taps should be 2D array')

        taps = _np.array(taps)
        kwargs.update({'taps': taps.flatten().tolist()})
        kwargs.update({'dimx': taps.shape[1]})
        kwargs.update({'dimy': taps.shape[0]})

        if nodata is not None:
            kwargs.update({'nodata': nodata})
        if norm is not None:
            kwargs.update({'norm': norm})

        self._jim_object._set(self._jim_object._jipjim.firfilter2d(kwargs))

    def getDissim(self, jimo=None, dissim_type=0):
        """Compute the dissimilarities.

        Compute the dissimilarities between horizontal and vertical pairs of
        adjacent pixels.

        :param jimo: a list of grey level Jim objects with the same
            definition domain. The dissimilarities are calculated for each
            image separately and composed using the point-wise maximum rule.
        :param dissim_type: integer value indicating the type of dissimilarity \
                                measure.
                            0 (default) for absolute difference.
                            1 for dissimilarity measure countering \
                              the chaining effect as described in \
                              :cite:`soille2011ismm`.
        """
        DIR_HORI = 0
        DIR_VERT = 1
        ABS_DIFF_op = 0
        MAX_op = 1
        MIN_op = 2

        if isinstance(jimo, _pj.Jim):
            jim_object_list = _pj.JimList([jimo])
        elif isinstance(jimo, list):
            jim_object_list = jimo
        else:
            jim_object_list = []

        if dissim_type == 0:
            # TODO: Check if looping through everything with edgeWeight() and
            #       then one call of supremum is not faster
            h_dissim = _pj.ngbops.edgeWeight(self._jim_object, DIR_HORI,
                                             ABS_DIFF_op)
            v_dissim = _pj.ngbops.edgeWeight(self._jim_object, DIR_VERT,
                                             ABS_DIFF_op)

            for im in jim_object_list:
                h_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_HORI,
                                                               ABS_DIFF_op))
                v_dissim.pixops.supremum(_pj.ngbops.edgeWeight(im, DIR_VERT,
                                                               ABS_DIFF_op))
        elif dissim_type == 1:
            mingraderograddil = _pj.pixops.infimum(
                _pj.ngbops.morphoGradientByDilationDiamond(
                    self._jim_object),
                _pj.ngbops.morphoGradientByErosionDiamond(self._jim_object))
            h_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_HORI,
                                             MAX_op)
            v_dissim = _pj.ngbops.edgeWeight(mingraderograddil, DIR_VERT,
                                             MAX_op)
            mingraderograddil = 0
            h_dissim.pixops.supremum(
                _pj.ngbops.edgeWeight(self._jim_object, DIR_HORI,
                                      ABS_DIFF_op))
            v_dissim.pixops.supremum(
                _pj.ngbops.edgeWeight(self._jim_object, DIR_VERT,
                                      ABS_DIFF_op))

            for im in jim_object_list:
                mingraderograddil = _pj.pixops.infimum(
                    _pj.ngbops.morphoGradientByDilationDiamond(im),
                    _pj.ngbops.morphoGradientByErosionDiamond(im))
                h_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil,
                                                     DIR_HORI,
                                                     MAX_op)
                v_dissim_crt = _pj.ngbops.edgeWeight(mingraderograddil,
                                                     DIR_VERT,
                                                     MAX_op)
                mingraderograddil = 0
                h_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(
                    im, DIR_HORI, ABS_DIFF_op))
                v_dissim_crt.pixops.supremum(_pj.ngbops.edgeWeight(
                    im, DIR_VERT, ABS_DIFF_op))

                h_dissim.pixops.supremum(h_dissim_crt)
                v_dissim.pixops.supremum(v_dissim_crt)
        else:
            raise ValueError('dissim_type {} not supported. Supported only '
                             'values 0 and 1'.format(dissim_type))

        return [h_dissim, v_dissim]


    def labelPixNgb(self, sec_jim_object, ox, oy, oz):
        """Label pix ngb.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate

        :return: labelled Jim object
        """
        self._jim_object._jipjim.d_labelPixNgb(sec_jim_object._jipjim,
                                               ox, oy, oz)

    def morphoDilate(self, sec_jim_object, ox, oy, oz, trFlag=0):
        """Output the dilation of im using the SE defined by imse.

        Its origin is set at coordinates (x,y,z). The reflection of the SE
        is considered if trflag equals 1 (no reflection by default). Points of
        the SE are points with a non zero value in imse.

        Modifies the instance on which the method was called.

        :param sec_jim_object: an image node for SE (UCHAR type)
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate
        :param trFlag: optional parameter (0 or 1)
        """
        self._jim_object._set(self._jim_object._jipjim.morphoDilate(
            sec_jim_object._jipjim, ox, oy, oz, trFlag))

    def morphoDilateDiamond(self):
        """Output the dilation of im using the elementary diamond shaped SE.

        Its origin is set at coordinates (x,y).

        Modifies the instance on which the method was called.
        """
        ox = 1
        oy = 1
        self._jim_object._jipjim.d_morphoDilateNgb4(ox, oy)

    def morphoDilateLine(self, dx, dy, k, o, type=0):
        """Output the dilation of im.

        Uses the line SE with slope dy/dx, length k, origin o, and line type
        (see details at  :cite:`soille-breen-jones96`).

        Modifies the instance on which the method was called.

        :param dx: integer for displacement along x-axis to set slope
        :param dy: integer for displacement along y-axis to set slope
        :param k: integer for number of pixels of line SE
        :param o: integer for origin of SE
        :param type: integer for line type (0 for plain and 1 for periodic).
            0 is the default value
        """
        self._jim_object._jipjim.d_morphoDilateLine(dx, dy, k, o, type)

    def morphoErode(self, sec_jim_object, ox, oy, oz, trFlag=0):
        """Output the erosion of im using the SE defined by imse.

        Its origin is set at coordinates (x,y,z). The reflection of the SE
        is considered if trflag equals 1 (no reflection by default). Points of
        the SE are points with a non zero value in imse.

        Modifies the instance on which the method was called.

        :param sec_jim_object: an image node for SE (UCHAR type)
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate
        :param trFlag: optional parameter (0 or 1)
        """
        self._jim_object._set(self._jim_object._jipjim.morphoErode(
            sec_jim_object._jipjim, ox, oy, oz, trFlag))

    def morphoErodeDiamond(self):
        """Output the erosion of im using the elementary diamond shaped SE.

        Its origin is set at coordinates (x,y).

        Modifies the instance on which the method was called.
        """
        ox = 1
        oy = 1
        self._jim_object._jipjim.d_morphoErodeNgb4(ox, oy)

    def morphoErodeLine(self, dx, dy, k, o, type=0):
        """Output the erosion of im using a line segment.

        Uses the line SE with slope dy/dx, length k, origin o, and line type.
        See details at  :cite:`soille-breen-jones96`.

        Modifies the instance on which the method was called.

        :param dx: integer for displacement along x-axis to set slope
        :param dy: integer for displacement along y-axis to set slope
        :param k: integer for number of pixels of line SE
            (must be odd; if not, it is extended by one pixel)
        :param o: integer for origin of SE
        :param type: integer for line type (0 for plain and 1 for periodic).
            0 is the default value
        """
        self._jim_object._jipjim.d_morphoErodeLine(dx, dy, k, o, type)

    def morphoGradientByDilationDiamond(self):
        """Output the gradient by dilation of Jim.

        Uses the elementary diamond shaped SE.

        Modifies the instance on which the method was called.
        """
        gradient = _pj.Jim(self._jim_object._jipjim.morphoDilateNgb4(1, 1)) - \
                   self._jim_object
        self._jim_object._set(gradient._jipjim)

    def morphoGradientByErosionDiamond(self):
        """Output the gradient by erosion of Jim.

        Uses the elementary diamond shaped SE.

        Modifies the instance on which the method was called.
        """
        self._jim_object -= _pj.Jim(
            self._jim_object._jipjim.morphoErodeNgb4(1, 1))

    def savgolay(self, **kwargs):
        """Compute the Savitzky-Golay filter in the time-spectral domain.

        :param nl: Number of leftward (past) data points used in Savitzky-
            Golay filter)
        :param nr: Number of rightward (future) data points used in Savitzky-
            Golay filter)
        :param ld: order of the derivative desired in Savitzky-Golay filter
            (e.g., ld=0 for smoothed function)
        :param m: order of the smoothing polynomial in Savitzky-Golay filter,
            also equal to the highest conserved moment; usual values are
            m=2 or m=4)

        Example:

        Perform a Savitzky-Golay filter to reconstruct a time series data set
        as in `J. Chen 2004 <https://doi.org/10.1016/j.rse.2004.03.014>`_::

            jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
            jim.pixops.convert('GDT_Float64')

            savgol = pj.ngbops.savgolay(jim, nl=7, nr=7, m=2, pad='replicate')
            for loop in range(0, 10):
                savgol.pixops.convert(otype=jim.properties.getDataType())
                savgol[savgol < jim] = jim
                savgol = pj.ngbops.savgolay(savgol, nl=4, nr=4, m=6,
                                            pad='replicate')
        """
        self._jim_object._set(self._jim_object._jipjim.savgolay(kwargs))

    def smoothNoData1d(self, nodata=0, **kwargs):
        """Smooth nodata in spectral/temporal domain.

        :param nodata: no data value to interpolate
        :param interpolationType: type of interpolation for spectral filtering
            (see https://www.gnu.org/software/gsl/doc/html/interp.html)

        Example::

            jim = pj.Jim('/path/to/multi-band/image.tif', band2plane=True)
            jim.pixops.convert('GDT_Float64')

            jim.ngbops.smoothNoData1d(0)
        """
        kwargs.update({'nodata': nodata})
        self._jim_object._set(self._jim_object._jipjim.smoothNoData1d(kwargs))

    # works but slower than directly in Numpy:  jim.np()[0]=jim.np().min(axis=0)
    # def stats1d(self, method, nodata=0, **kwargs):
    #     """Calculate statistics for each pixel in the spectral/temporal domain

    #     :param method: (list of) methods to calculate
    #     :param nodata: no data value to interpolate

    #     Example:

    #     jim=pj.Jim('/path/to/multi-band/image.tif',band2plane=True)

    #     jim.ngbops.stats1d()
    #     """

    #     kwargs.update({'method': method})
    #     kwargs.update({'nodata':nodata})
    #     self._jim_object._set(self._jim_object._jipjim.stats1d(kwargs))


class _NgbOpsList(_pj.modules.JimListModuleBase):
    """Define all NgbOps methods for JimLists."""

    pass


class _NgbOpsVect(_pj.modules.JimVectModuleBase):
    """Define all NgbOps methods for JimVects."""

    pass
