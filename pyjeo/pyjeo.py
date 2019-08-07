"""Basic file containing Jim, JimList and JimVect objects."""

from __future__ import division
import numpy
import gc as _gc
import warnings as _warnings
import os as _os

try:
    import jiplib as _jl
except ImportError:
    from jeodpp import jiplib as _jl

from modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, classify, demops, stats, all


del _jl.Jim.__del__


def np(aJim):
    """Return a pointer to numpy representation of values in a Jim object.

    The created pointer does not consume new memory.

    :param aJim: Jim object with values to which will the pointer point
    :return: a numpy representation of the Jim object
    """
    return _jl.np(aJim._jipjim)


def jim2np(aJim, band=0, copyData=True):
    """Return a numpy representation of a Jim object.

    :param aJim: Jim object to be converted
    :param band: band of Jim object to be converted
    :param copyData: Set to False if reference image is used as a template
        only, without copying actual pixel dat
    :return: a numpy representation of the Jim object
    """
    return _jl.jim2np(aJim._jipjim, band, copyData)


def np2jim(aNp):
    """Return a Jim representation of a numpy array.

    :param aNp: a numpy array
    :return: a Jim representation of a numpy array
    """
    return Jim(_jl.np2jim(aNp))


class _ParentJim(_jl.Jim):

    def __init__(self, image, kwargs):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        if kwargs:
            if image:
                if isinstance(image, Jim):
                    if 'copyData' in kwargs.keys():
                        super(_ParentJim, self).__init__(image._jipjim,
                                                         kwargs['copyData'])
                    else:
                        _warnings.warn(
                            'Not possible to create Jim image based on another'
                            ' one together with other kwargs than copyData. '
                            'kwargs ignored.', SyntaxWarning)
                        super(_ParentJim, self).__init__(image._jipjim)
                else:
                    kwargs.update({'filename': image})
                    super(_ParentJim, self).__init__(kwargs)
            else:
                super(_ParentJim, self).__init__(kwargs)
        else:
            if isinstance(image, Jim):
                super(_ParentJim, self).__init__(image._jipjim)
            else:
                super(_ParentJim, self).__init__(image)


class Jim():
    """Definition of Jim object."""

    def __init__(self, image=None, **kwargs):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        self._checkInitParamsSense(image, kwargs)

        # remove stdev and uniform from kwargs to use them in feed
        stdev = kwargs.pop('stdev', None)
        uniform = kwargs.pop('uniform', None)
        seed = kwargs.pop('seed', None)

        self._jipjim = _ParentJim(image, kwargs)

        self._all = all._All()
        self._ccops = ccops._CCOps()
        self._classify = classify._Classify()
        self._demops = demops._DEMOps()
        self._geometry = geometry._Geometry()
        self._io = io._IO()
        self._ngbops = ngbops._NgbOps()
        self._pixops = pixops._PixOps()
        self._properties = properties._Properties()
        self._stats = stats._Stats()

        if any(arg is not None for arg in [stdev, uniform, seed]):
            self._feed(stdev, uniform, seed, kwargs)

    @property
    def all(self):
        """Set up a caller and garbage cleaner for the module all."""
        self._all._set_caller(self)
        _gc.collect()
        return self._all

    @property
    def ccops(self):
        """Set up a caller and garbage cleaner for the module ccops."""
        self._ccops._set_caller(self)
        _gc.collect()
        return self._ccops

    @property
    def classify(self):
        """Set up a caller and garbage cleaner for the module classify."""
        self._classify._set_caller(self)
        _gc.collect()
        return self._classify

    @property
    def demops(self):
        """Set up a caller and garbage cleaner for the module demops."""
        self._demops._set_caller(self)
        _gc.collect()
        return self._demops

    @property
    def geometry(self):
        """Set up a caller and garbage cleaner for the module geometry."""
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    @property
    def io(self):
        """Set up a caller and garbage cleaner for the module io."""
        self._io._set_caller(self)
        _gc.collect()
        return self._io

    @property
    def ngbops(self):
        """Set up a caller and garbage cleaner for the module ngbops."""
        self._ngbops._set_caller(self)
        _gc.collect()
        return self._ngbops

    @property
    def pixops(self):
        """Set up a caller and garbage cleaner for the module pixops."""
        self._pixops._set_caller(self)
        _gc.collect()
        return self._pixops

    @property
    def properties(self):
        """Set up a caller and garbage cleaner for the module properties."""
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

    @property
    def stats(self):
        """Set up a caller and garbage cleaner for the module stats."""
        self._stats._set_caller(self)
        _gc.collect()
        return self._stats

    def getMethods(self, queried_module=None):
        """Print an overview of available methods in format module.method."""
        def treeStructure(module, queried_module):
            if queried_module and queried_module not in str(module):
                return ''

            module_methods = dir(module)
            for default_method in ['__init__', '__module__', '__doc__',
                                   '_set_caller']:
                module_methods.remove(default_method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:])] + \
                   module_methods

        methods = list()
        for module in [properties._Properties, io._IO, pixops._PixOps,
                       ngbops._NgbOps, geometry._Geometry, ccops._CCOps,
                       classify._Classify, demops._DEMOps, stats._Stats]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def np(self, band=0):
        """Return numpy array from Jim object.

        :param band: band index (starting from 0)
        :return: numpy array representation
        """
        try:
            self.properties.getDataType()
        except TypeError:
            raise ValueError(
                'Jim has to have a data type and dimensions to use Jim.np()')
        if band >= self.properties.nrOfBand():
            raise ValueError('Band out of bounds')
        return _jl.np(self._jipjim, band)

    def _checkInitParamsSense(self, image, kwargs):
        """Check if the combination of (kw)args for Jim init makes sense.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        keys = kwargs.keys()

        if not image:
            if keys:
                if ('seed' in keys or 'uniform' in keys or 'stdev' in keys) \
                        and ('ncol' not in keys or 'nrow' not in keys or
                             'otype' not in keys):
                    raise AttributeError('You cannot use any of parameters '
                                         '[seed, uniform, stdev] without '
                                         'specifying the geometry and otype '
                                         'of Jim.')

    def _checkNumberOfBands(self, another_jim):
        """Check if the number of bands matches for both Jim objects.

        Raise a Warning otherwise.

        :param another_jim: Another Jim object which nrOfBands is going to be
            checked
        """
        if self.properties.nrOfBand() > another_jim.properties.nrOfBand() > 1:
            raise IndexError('Jims used in conditions must have either the '
                             'same number of bands or one of them must have '
                             'number of bands == 1')
        elif 1 < self.properties.nrOfBand() < \
                another_jim.properties.nrOfBand():
            _warnings.warn('Left Jim has less bands than the right one and '
                           'more than one band. Only corresponding bands of '
                           'the right Jim will be taken in consideration.',
                           Warning)

    def _feed(self, stdev, uniform, seed, kwargs):
        """Feed the Jim object with either uniform or random seed of values.

        :param stdev: standard deviation
        :param uniform: if uniform distr, here are defined values [min, max]
        :param seed: numpy random seed
        """
        mean = kwargs.pop('mean', None)

        if seed is not None:
            numpy.random.seed(seed)

        if uniform:
            if isinstance(uniform, list):
                if len(uniform) != 2:
                    raise AttributeError(
                        'The list parsed as the uniform argument must be '
                        'in the form [min, max + 1]')
                self.np()[:] = numpy.random.uniform(uniform[0], uniform[1],
                                                    self.np().shape)
            else:
                self.np()[:] = numpy.random.uniform(0, uniform,
                                                    self.np().shape)
        else:
            if stdev is None:
                stdev = 1
            if mean is None:
                mean = 0

            self.np()[:] = numpy.random.normal(mean, stdev, self.np().shape)

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self._jipjim.__dict__.update(modified_object.__dict__)

    # *** unary operators *** #

    def __getitem__(self, item):
        """Evaluate operation of type self[key].

        :param item: key/index
        """
        if isinstance(item, JimVect):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Using a JimVect as an index not implemented '
                                 'for 3d Jim objects')
            nodata = self.properties.getNoDataVals()
            if nodata:
                nodata = nodata[0]
            else:
                nodata = 0
            return geometry.cropOgr(self, item,
                                    crop_to_cutline=True, nodata=nodata,
                                    align=True)
        elif isinstance(item, Jim):
            mask = item > 0
            return Jim(self * mask)
        else:
            npresult = numpy.array(self.np()[item], copy=True)
            # npresult=numpy.array(self.np()[item])
            if len(npresult.shape) == 3:
                nplane = npresult.shape[0]
                nrow = npresult.shape[1]
                ncol = npresult.shape[2]
            elif len(npresult.shape) == 2:
                nplane = 1
                nrow = npresult.shape[0]
                ncol = npresult.shape[1]
            elif len(npresult.shape) == 1:
                nplane = 1
                nrow = 1
                ncol = npresult.shape[0]
            elif len(npresult.shape) == 0:
                nplane = 1
                nrow = 1
                ncol = 1

            # [gs]item only supports single band image (use plane instead)
            nband = 1
            if self.properties.nrOfPlane() > 1:
                dim = 3
            else:
                dim = 2
            dx = self.properties.getDeltaX()
            dy = self.properties.getDeltaY()

            cropuli = 0
            # croplri=self.properties.nrOfCol()
            cropulj = 0
            # croplrj=self.properties.nrOfRow()
            if isinstance(item, tuple):
                # cols
                if len(item) > dim-1:
                    if isinstance(item[dim-1], slice):
                        if item[dim-1].start:
                            cropuli = item[dim-1].start
                        if item[dim-1].step:
                            dx *= item[dim-1].step
                        # croplri=item[dim-1].stop
                    else:
                        cropuli = item[dim-1]
                        # croplri=item[dim-1]+1
                # rows
                if len(item) > dim-2:
                    if isinstance(item[dim-2], slice):
                        if item[dim-2].start:
                            cropulj = item[dim-2].start
                        if item[dim-2].step:
                            dy *= item[dim-2].step
                        # croplrj=item[dim-2].stop
                    else:
                        cropulj = item[dim-2]
                        # croplrj=item[dim-2]+1

            upperLeft = self.geometry.image2geo(cropuli, cropulj)
            result = Jim(ncol=ncol, nrow=nrow, nband=nband, nplane=nplane,
                         otype=self.properties.getDataType())
            result.properties.setProjection(self.properties.getProjection())
            gt = self.properties.getGeoTransform()

            cropulx = upperLeft[0]-self.properties.getDeltaX()/2
            cropuly = upperLeft[1]+self.properties.getDeltaY()/2

            gt[0] = cropulx
            gt[1] = dx
            gt[2] = 0
            gt[3] = cropuly
            gt[4] = 0
            gt[5] = -dy
            result.properties.setGeoTransform(gt)
            result.np()[:] = npresult
            return result

    def __setitem__(self, item, value):
        """Evaluate operation of type self[key] = value.

        :param item: key/index
        :param value: value to save
        """
        if isinstance(item, JimVect):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __setitem__ with JimVect not '
                                 'implemented for 3d Jim objects')
            # TODO: decide on default behaviour of ALL_TOUCHED=TRUE
            # TODO: next lines should work, but problem with GML files when SRS
            #       is not defined as in S2 cloud masks

            # template=Jim(self)
            # template.geometry.rasterize(item,1.0)
            # self[template>0]=value

            if type(value) in (float, int):
                templateJim = Jim(self, copyData=False)
                templateJim = Jim(templateJim._jipjim.setMask(
                    item._jipjimvect, {'eo': ['ALL_TOUCHED=TRUE'],
                                       'nodata': 1}))
                self[templateJim > 0] = value
            elif isinstance(value, Jim):
                templateJim = Jim(self, copyData=False)
                templateJim = Jim(templateJim._jipjim.setMask(
                    item._jipjimvect, {'eo': ['ALL_TOUCHED=TRUE'],
                                       'nodata': 1}))
                self[templateJim > 0] = value
        elif isinstance(item, Jim):
            if value is not None:
                if isinstance(value, Jim):
                    self._jipjim.d_setMask(item._jipjim, value._jipjim)
                else:
                    self._jipjim.d_setMask(item._jipjim, value)
        elif isinstance(item, tuple):
            if isinstance(value, Jim):
                self.np()[item] = value.np()
            else:
                self.np()[item] = value
        else:
            raise ValueError('Error: __setitem__ only implemented for Vector, '
                             'Jim or tuples (dims of __getitem__ must '
                             'correspond with those of the tuple)')

    def __nonzero__(self):
        """Check if Jim contains data.

        :return: True if image contains data, False if image is empty
        """
        return self._jipjim.isInit()

    def __bool__(self):
        """Check if Jim contains data.

        :return: True if image contains data, False if image is empty
        """
        return self._jipjim.isInit()

    def __abs__(self):
        """Calculate the absolute value of Jim raster dataset.

        :return: Absolute value of Jim raster dataset
        """
        jim = Jim(self)
        for iband in range(0, self.properties.nrOfBand()):
            jim.np(iband)[:] = abs(jim.np(iband))
        return jim

    def __neg__(self):
        """Calculate the negation of Jim raster dataset.

        :return: Negation of Jim raster dataset (-dataset)
        """
        jim = Jim(self)
        for iband in range(0, self.properties.nrOfBand()):
            jim.np(iband)[:] = -(jim.np(iband))
        return jim

    def __invert__(self):
        """Calculate the complement of Jim raster dataset.

        :return: The complement of Jim raster dataset (~dataset)
        """
        jim = Jim(self)
        for iband in range(0, self.properties.nrOfBand()):
            jim.np(iband)[:] = ~(jim.np(iband))
        return jim

    # *** binary operators *** #
    def __eq__(self, right):
        """Pixel wise check for equality.

        :return: Jim object with pixels 1 if equal values, 0 otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) == right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) == right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) == right)
        return jim

    def __ne__(self, right):
        """Pixel wise check for non-equality.

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) != right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) != right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) != right)
        return jim

    def __lt__(self, right):
        """Pixel wise check for lesser than.

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) < right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) < right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) < right)
        return jim

    def __le__(self, right):
        """Pixel wise check for lesser or equal.

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) <= right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) <= right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) <= right)
        return jim

    def __gt__(self, right):
        """Pixel wise check for greater than.

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) > right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) > right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) > right)
        return jim

    def __ge__(self, right):
        """Pixel wise check for greater or equal.

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(), self.properties.nrOfRow(),
                          self.properties.nrOfBand(),
                          self.properties.nrOfPlane(), _jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] = (self.np(iband) >= right.np())
                else:
                    jim.np(iband)[:] = (self.np(iband) >= right.np(iband))
        else:
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] = (self.np(iband) >= right)
        return jim

    def __add__(self, right):
        """Pixel wise operation +."""
        result = Jim(self)
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] += right.np()
                else:
                    result.np(iband)[:] += right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] += right
        return result

    def __radd__(self, left):
        """Pixel wise operation + where self is the added value."""
        result = Jim(self)

        for iband in range(0, self.properties.nrOfBand()):
            result.np(iband)[:] += left

        return result

    def __iadd__(self, right):
        """Pixel wise operation +=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] += right.np()
                else:
                    self.np(iband)[:] += right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] += right
        return self

    def __sub__(self, right):
        """Pixel wise operation -."""
        result = Jim(self)
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] -= right.np()
                else:
                    result.np(iband)[:] -= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] -= right
        return result

    def __rsub__(self, left):
        """Pixel wise operation - where self is the one used to subtract."""
        result = -Jim(self)

        for iband in range(0, self.properties.nrOfBand()):
            result.np(iband)[:] += left

        return result

    def __isub__(self, right):
        """Pixel wise operation -=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] -= right.np()
                else:
                    self.np(iband)[:] -= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] -= right
        return self

    def __mul__(self, right):
        """Pixel wise operation *."""
        result = Jim(self)
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] *= right.np()
                else:
                    result.np(iband)[:] *= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] *= right
        return result
        # if isinstance(right, Jim):
        #     if right.properties.getDataType() == _jl.GDT_Float32 or \
        #        right.properties.getDataType() != _jl.GDT_Float64:
        #         if result.properties.getDataType() != _jl.GDT_Float32 and \
        #            result.properties.getDataType() != _jl.GDT_Float64:
        #             result.pixops.convert(otype=right.properties.getDataType())
        #     result.np()[:]*=right.np()
        # else:
        #     elif isinstance(right, int):
        #         result.np()[:]*=right
        #     elif isinstance(right, float):
        #         if result.properties.getDataType() != _jl.GDT_Float32 and \
        #            result.properties.getDataType() != _jl.GDT_Float64:
        #             result.pixops.convert(otype=right.properties.getDataType())
        # return self

    def __rmul__(self, left):
        """Pixel wise operation * where self is the multiplier."""
        result = Jim(self)

        for iband in range(0, self.properties.nrOfBand()):
            result.np(iband)[:] *= left

        return result

    def __imul__(self, right):
        """Pixel wise operation *=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] *= right.np()
                else:
                    self.np(iband)[:] *= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] *= right
        return self

    def __truediv__(self, right):
        """Pixel wise operation for true division."""
        result = Jim(self)
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] /= right.np()
                else:
                    result.np(iband)[:] /= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] /= right
        return result

    def __itruediv__(self, right):
        """Pixel wise operation for true division with /=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] /= right.np()
                else:
                    self.np(iband)[:] /= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] /= right
        return self

    def __div__(self, right):
        """Pixel wise operation /."""
        result = Jim(self)
        if 'int' in str(self.np().dtype):
            raise TypeError('You cannot divide a Jim object of int type')
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] /= right.np()
                else:
                    result.np(iband)[:] /= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] /= right
        return result

    def __idiv__(self, right):
        """Pixel wise operation /=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] /= right.np()
                else:
                    self.np(iband)[:] /= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] /= right
        return self

    def __mod__(self, right):
        """Pixel wise operation modulo."""
        result = Jim(self)
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    result.np(iband)[:] %= right.np()
                else:
                    result.np(iband)[:] %= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                result.np(iband)[:] %= right
        return result

    def __imod__(self, right):
        """Pixel wise operation modulo with %=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] %= right.np()
                else:
                    self.np(iband)[:] %= right.np(iband)
        else:
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] %= right
        return self

    def __pow__(self, right):
        """Pixel wise operation **."""
        result = Jim(self)
        for iband in range(0, self.properties.nrOfBand()):
            result.np(iband)[:] **= right
        return result

    def __ipow__(self, right):
        """Pixel wise operation **=."""
        for iband in range(0, self.properties.nrOfBand()):
            self.np(iband)[:] **= right
        return self

    def __lshift__(self, right):
        """Pixel wise operation <<."""
        if isinstance(right, int):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] <<= right
            return jim
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))

    def __ilshift__(self, right):
        """Pixel wise operation <<=."""
        if isinstance(right, int):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] <<= right
        else:
            raise TypeError('unsupported operand type for <<= : {}'.format(
                type(right)))
        return self

    def __rshift__(self, right):
        """Pixel wise operation >>."""
        if isinstance(right, int):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] >>= right
            return jim
        else:
            raise TypeError('unsupported operand type for >> : {}'.format(
                type(right)))

    def __irshift__(self, right):
        """Pixel wise operation >>=."""
        if isinstance(right, int):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] >>= right
        else:
            raise TypeError('unsupported operand type for >>= : {}'.format(
                type(right)))
        return self

    def __or__(self, right):
        """Pixel wise operation |."""
        if isinstance(right, Jim):
            jim = Jim(self)
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] |= right.np()
                else:
                    jim.np(iband)[:] |= right.np(iband)
            return jim
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] |= right
            return jim
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(right)))

    def __ior__(self, right):
        """Pixel wise operation |=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] |= right.np()
                else:
                    self.np(iband)[:] |= right.np(iband)
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] |= right
        else:
            raise TypeError('unsupported operand type for |= : {}'.format(
                type(right)))
        return self

    def __ror__(self, left):
        """Pixel wise operation | where self is the right object."""
        if isinstance(left, int) or isinstance(left, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] |= left
            return jim
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(left)))

    def __xor__(self, right):
        """Pixel wise operation ^."""
        if isinstance(right, Jim):
            jim = Jim(self)
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] ^= right.np()
                else:
                    jim.np(iband)[:] ^= right.np(iband)
            return jim
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] ^= right
            return jim
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(right)))

    def __ixor__(self, right):
        """Pixel wise operation ^=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] ^= right.np()
                else:
                    self.np(iband)[:] ^= right.np(iband)
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] ^= right
        else:
            raise TypeError('unsupported operand type for ^= : {}'.format(
                type(right)))
        return self

    def __rxor__(self, left):
        """Pixel wise operation ^ where self is the right object."""
        if isinstance(left, int) or isinstance(left, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] ^= left
            return jim
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(left)))

    def __and__(self, right):
        """Pixel wise operation &."""
        if isinstance(right, Jim):
            jim = Jim(self)
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    jim.np(iband)[:] &= right.np()
                else:
                    jim.np(iband)[:] &= right.np(iband)
            return jim
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] &= right
            return jim
        else:
            raise TypeError('unsupported operand type for & : {}'.format(
                type(right)))

    def __iand__(self, right):
        """Pixel wise operation &=."""
        if isinstance(right, Jim):
            self._checkNumberOfBands(right)

            for iband in range(0, self.properties.nrOfBand()):
                if right.properties.nrOfBand() == 1:
                    self.np(iband)[:] &= right.np()
                else:
                    self.np(iband)[:] &= right.np(iband)
        elif isinstance(right, int) or isinstance(right, numpy.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] &= right
        else:
            raise TypeError('unsupported operand type for &= : {}'.format(
                type(right)))
        return self

    def __rand__(self, left):
        """Pixel wise operation & where self is the right object."""
        if isinstance(left, int) or isinstance(left, numpy.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] &= left
            return jim
        else:
            raise TypeError('unsupported operand type for & : {}'.format(
                type(left)))


class _ParentList(_jl.JimList):

    def __init__(self, images_list, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        jiplib_images_list = [i._jipjim if isinstance(i, Jim) else i for i in
                              images_list]
        super(_ParentList, self).__init__(jiplib_images_list, *args)


class JimList(list):
    """Definition of JimList object."""

    def __init__(self, images_list=[], *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        super(JimList, self).__init__(images_list)
        self._jipjimlist = _ParentList(images_list, *args)

        self._all = all._AllList()
        self._ccops = ccops._CCOpsList()
        self._classify = classify._ClassifyList()
        self._demops = demops._DEMOpsList()
        self._geometry = geometry._GeometryList()
        self._io = io._IOList()
        self._ngbops = ngbops._NgbOpsList()
        self._pixops = pixops._PixOpsList()
        self._properties = properties._PropertiesList()
        self._stats = stats._StatsList()

    @property
    def geometry(self):
        """Set up a caller and garbage cleaner for the module geometry."""
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    @property
    def io(self):
        """Set up a caller and garbage cleaner for the module io."""
        self._io._set_caller(self)
        _gc.collect()
        return self._io

    @property
    def pixops(self):
        """Set up a caller and garbage cleaner for the module pixops."""
        self._pixops._set_caller(self)
        _gc.collect()
        return self._pixops

    @property
    def properties(self):
        """Set up a caller and garbage cleaner for the module properties."""
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

    @property
    def stats(self):
        """Set up a caller and garbage cleaner for the module stats."""
        self._stats._set_caller(self)
        _gc.collect()
        return self._stats

    def getMethods(self, queried_module=None):
        """Print an overview of available methods in format module.method."""
        def treeStructure(module, queried_module):
            if queried_module and queried_module not in str(module):
                return ''

            module_methods = dir(module)
            for default_method in ['__init__', '__module__', '__doc__',
                                   '_set_caller']:
                module_methods.remove(default_method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:-4] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:-4])] + \
                   module_methods

        methods = list()

        for module in [geometry._GeometryList, io._IOList, pixops._PixOpsList,
                       properties._PropertiesList, stats._StatsList]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def append(self, jim):
        """Add single element to the JimList."""
        assert isinstance(jim, Jim), \
            'Only Jim instances can be appended'
        super(JimList, self).append(jim)
        self._set(self, from_list=True)

    def count(self, jim):
        """Count the occurrences of element in the JimList.."""
        i = 0
        for jimobject in self:
            if jim.pixops.isEqual(jimobject):
                i += 1

        return i

    def extend(self, jim_list):
        """Add elements of a JimList to another JimList."""
        assert isinstance(jim_list, JimList), \
            'Only JimList instances can be used to extend JimList'
        super(JimList, self).extend(jim_list)
        self._set(self, from_list=True)

    def index(self, jim):
        """Return smallest index of element in the JimList."""
        for i in range(len(self)):
            if self[i].pixops.isEqual(jim):
                return i

    def insert(self, index, jim):
        """Insert elements to the JimList."""
        assert isinstance(jim, Jim), \
            'Only Jim instances can be inserted'
        super(JimList, self).insert(index, jim)
        self._set(self, from_list=True)

    def pop(self, index):
        """Remove and return element at given index."""
        popped = super(JimList, self).pop(index)
        self._set(self, from_list=True)
        return popped

    def remove(self, jim):
        """Remove the first occurence of an element from the JimList."""
        for i in range(len(self)):
            if self[i].pixops.isEqual(jim):
                del self[i]
                break
        self._set(self, from_list=True)

    def reverse(self):
        """Reverse the JimList."""
        super(JimList, self).reverse()
        self._set(self, from_list=True)

    def _set(self, modified_list, from_list=False):
        """Apply changes done in modified_list to the parent JimList instance.

        :param modified_object: modified JimList instance
        :param from_list: set True if changing function originates in list()
        """
        if not from_list:
            del self[:]
            for i in range(modified_list._jipjimlist.getSize()):
                im = modified_list._jipjimlist.getImage(i)
                self.append(Jim(im))
        else:
            for _ in range(self._jipjimlist.getSize()):
                self._jipjimlist.popImage()
            for image in modified_list:
                self._jipjimlist.pushImage(image._jipjim)

    def __dir__(self):
        """Change behaviour of the method whisperer to ignore jiplib methods.

        :return: a whispered module or method
        """
        return list(set(dir(JimList)) - {'sort'})


class _ParentVect(_jl.VectorOgr):

    def __init__(self, vector, kwargs):
        """Initialize the JimVect object and modules for methods.

        :param vector: path to a vector or another JimVect object as a basis
            for the JimVect object
        :param output: path to an output vector in case another JimVect object
            was provided (copy constructor)
        """
        self._checkInitParamsSense(vector, kwargs)

        if kwargs:
            if vector:
                if isinstance(vector, JimVect):
                    kwargs.update({'filename': kwargs.pop('output', None)})
                    super(_ParentVect, self).__init__(vector._jipjimvect,
                                                      kwargs)
                else:
                    kwargs.update({'filename': vector})
                    super(_ParentVect, self).__init__(kwargs)
            else:
                kwargs.update({'filename': kwargs.pop('output', None)})
                super(_ParentVect, self).__init__(kwargs)
        else:
            if vector:
                super(_ParentVect, self).__init__(vector)
            else:
                super(_ParentVect, self).__init__()

    def _checkInitParamsSense(self, vector, kwargs):
        """Check if the combination of (kw)args for JimVect init makes sense.

        :param vector: path to a vector or another JimVect object as a basis
            for the JimVect object
        """
        keys = kwargs.keys()

        if isinstance(vector, JimVect):
            if 'output' not in keys:
                raise AttributeError(
                    "Parameter output required for copy constructor")
        elif not vector:
            if 'output' in keys and not _os.path.isfile(kwargs['output']):
                raise AttributeError('Output path does not exist and the '
                                     'template vector is not specified')


class JimVect():
    """Definition of JimVect object."""

    def __init__(self, vector=None, **kwargs):
        """Create an empty VectorOgr object.

        Created object is an instance of the basis vector class of the pyJEO
        library

        :param vector: Path to a vector dataset
        :return: a JimVect object
        """
        self._jipjimvect = _ParentVect(vector, kwargs)

        self._all = all._AllVect()
        self._ccops = ccops._CCOpsVect()
        self._classify = classify._ClassifyVect()
        self._demops = demops._DEMOpsVect()
        self._geometry = geometry._GeometryVect()
        self._io = io._IOVect()
        self._ngbops = ngbops._NgbOpsVect()
        self._pixops = pixops._PixOpsVect()
        self._properties = properties._PropertiesVect()
        self._stats = stats._StatsVect()

    @property
    def classify(self):
        """Set up a caller and garbage cleaner for the module classify."""
        self._classify._set_caller(self)
        _gc.collect()
        return self._classify

    @property
    def io(self):
        """Set up a caller and garbage cleaner for the module io."""
        self._io._set_caller(self)
        _gc.collect()
        return self._io

    @property
    def properties(self):
        """Set up a caller and garbage cleaner for the module properties."""
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

    @property
    def geometry(self):
        """Set up a caller and garbage cleaner for the module geometry."""
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    def getMethods(self, queried_module=None):
        """Print an overview of available methods in format module.method."""
        def treeStructure(module, queried_module):
            if queried_module and queried_module not in str(module):
                return ''

            module_methods = dir(module)
            for default_method in ['__init__', '__module__', '__doc__',
                                   '_set_caller']:
                module_methods.remove(default_method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:-4] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:-4])] + \
                   module_methods

        methods = list()
        for module in [classify._ClassifyVect, io._IOVect,
                       properties._PropertiesVect]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def _set(self, modified_object):
        """Apply changes done in modified_object to parent VectorOgr instance.

        :param modified_object: modified VectorOgr instance
        """
        self._jipjimvect.__dict__.update(modified_object.__dict__)
