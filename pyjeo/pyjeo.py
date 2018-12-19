"""Basic file containing Jim and JimList objects."""

from __future__ import division
import numpy as _np
import gc as _gc

import jiplib as _jl

from modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, clssfy, demops, stats, all


del _jl.Jim.__del__


def jim2np(aJim, band=0):
    return _jl.jim2np(aJim._jipjim, band)


def np2jim(aNp):
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
                        import warnings
                        warnings.warn(
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
        self._jipjim = _ParentJim(image, kwargs)

        self._all = all._All()
        self._ccops = ccops._CCOps()
        self._clssfy = clssfy._Classify()
        self._demops = demops._DEMOps()
        self._geometry = geometry._Geometry()
        self._io = io._IO()
        self._ngbops = ngbops._NgbOps()
        self._pixops = pixops._PixOps()
        self._properties = properties._Properties()
        self._stats = stats._Stats()

    @property
    def all(self):
        self._all._set_caller(self)
        _gc.collect()
        return self._all

    @property
    def ccops(self):
        self._ccops._set_caller(self)
        _gc.collect()
        return self._ccops

    @property
    def clssfy(self):
        self._clssfy._set_caller(self)
        _gc.collect()
        return self._clssfy

    @property
    def demops(self):
        self._demops._set_caller(self)
        _gc.collect()
        return self._demops

    @property
    def geometry(self):
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    @property
    def io(self):
        self._io._set_caller(self)
        _gc.collect()
        return self._io

    @property
    def ngbops(self):
        self._ngbops._set_caller(self)
        _gc.collect()
        return self._ngbops

    @property
    def pixops(self):
        self._pixops._set_caller(self)
        _gc.collect()
        return self._pixops

    @property
    def properties(self):
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

    @property
    def stats(self):
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
                       clssfy._Classify, demops._DEMOps, stats._Stats]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self._jipjim.__dict__.update(modified_object.__dict__)

    # *** unary operators *** #

    def __getitem__(self, item):
        stridex = 1
        stridey = 1
        strideb = 1
        if isinstance(item, tuple):
            if isinstance(item[0], slice):
                if item[0].start is None:
                    minCol = 0
                else:
                    minCol =  item[0].start
                if item[0].stop is None:
                    maxCol = self.properties.nrOfCol()
                else:
                    maxCol = item[0].stop

                if item[0].step:
                    stridex = item[0].step

                if minCol < 0:
                    minCol = self.properties.nrOfCol()+minCol
                if maxCol < 0:
                    maxCol = self.properties.nrOfCol()+maxCol+1
            elif isinstance(item[0], int):
                if item[0] < 0:
                    minCol = self.properties.nrOfCol()+item[0]
                else:
                    minCol = item[0]
                maxCol = minCol + 1
            else:
                raise ValueError('column item must be slice or integer value')

            if isinstance(item[1], slice):
                if item[1].start is None:
                    minRow = 0
                else:
                    minRow =  item[1].start
                if item[1].stop is None:
                    maxRow = self.properties.nrOfRow()
                else:
                    maxRow = item[1].stop

                if item[1].step:
                    stridey = item[1].step

                if minRow < 0:
                    minRow = self.properties.nrOfRow()+minRow
                if maxRow < 0:
                    maxRow = self.properties.nrOfRow()+maxRow+1
            elif isinstance(item[1], int):
                if item[1] < 0:
                    minRow = self.properties.nrOfRow()+item[1]
                else:
                    minRow = item[1]
                maxRow = minRow+1
            else:
                raise ValueError('row item must be slice or integer value')

            bands = []
            if self.properties.nrOfPlane() > 1 and len(item) > 2:
                if isinstance(item[2], slice):
                    if item[2].start is None:
                        min_z = 0
                    else:
                        min_z = item[2].start
                    if item[2].stop is None:
                        max_z = self.properties.nrOfPlane()
                    else:
                        max_z = item[2].stop

                    if item[2].step:
                        stride_z = item[2].step

                    if min_z < 0:
                        min_z = self.properties.nrOfPlane()+min_z
                    if max_z < 0:
                        max_z = self.properties.nrOfPlane()+max_z+1
                elif isinstance(item[2], int):
                    if item[2] < 0:
                        min_z = self.properties.nrOfPlane()+item[2]
                    else:
                        min_z = item[2]
                    max_z = min_z+1
                else:
                    raise ValueError('Z index must be slice or integer value')

                if self.properties.nrOfBand() > 1:
                    if len(item) == 4:  # do slice x,y,z,band
                        if isinstance(item[3], slice):
                            if item[3].start < 0:
                                minBand = self.properties.nrOfBand()+item[
                                    3].start
                            else:
                                if item[3].start is None:
                                    minBand = 0
                                else:
                                    minBand = item[3].start
                            if item[3].stop < 0:
                                maxBand = self.properties.nrOfBand()+item[
                                    3].stop + 1
                            else:
                                if item[3].stop is None:
                                    maxBand = self.properties.nrOfBand()
                                else:
                                    maxBand = item[3].stop

                            if item[3].step:
                                strideb = item[3].step

                            bands = range(minBand, maxBand, strideb)
                        elif isinstance(item[3], tuple):
                            for band in item[3]:
                                if band < 0:
                                    bands.append(
                                        self.properties.nrOfBand()+band)
                                else:
                                    bands.append(band)
                        elif isinstance(item[3], int):
                            if item[3] < 0:
                                bands.append(
                                  self.properties.nrOfBand()+item[3])
                            else:
                                bands.append(item[3])
                        else:
                            raise ValueError('Error: band index must be '
                                             'slice, list or integer')
                        retJim = geometry.cropBand(self, band=bands)
                        if maxCol <= minCol or maxRow <= minRow or \
                              max_z <= min_z:
                            raise IndexError(
                                'Warning: index error, returning empty Jim')
                        retJim.geometry.crop(ulx=minCol, lrx=maxCol,
                                             uly=minRow, lry=maxRow,
                                             ulz=min_z, lrz=max_z,
                                             dx=stridex, dy=stridey,
                                             nogeo=True)
                        return retJim
                    else:
                        raise TypeError(
                            'Error: use 4 dimensions when slicing multiband '
                            '3-dim Jim object (x:y:z:band)')
                else:
                    if len(item) == 3:  # do slice x,y,z
                        if maxCol <= minCol or maxRow <= minRow or \
                                item[2].stop <= item[2].start:
                            raise IndexError('Warning: index error, '
                                             'returning empty Jim')
                        retJim = geometry.crop(self, ulx=minCol, uly=minRow,
                                               ulz=min_z, lrx=maxCol,
                                               lry=maxRow, lrz=max_z,
                                               dx=stridex, dy=stridey,
                                               nogeo=True)
                        return retJim
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing '
                                        '3-dim Jim object (x:y:z)')
            else:
                if self.properties.nrOfBand() > 1:
                    if len(item) == 3:  # do slice x,y,band
                        if isinstance(item[2], slice):
                            if item[2].start < 0:
                                minBand = self.properties.nrOfBand()+item[
                                    2].start
                            else:
                                if item[2].start is None:
                                    minBand = 0
                                else:
                                    minBand = item[2].start
                            if item[2].stop < 0:
                                maxBand = self.properties.nrOfBand()+item[
                                    2].stop
                            else:
                                if item[2].stop is None:
                                    maxBand = self.properties.nrOfBand()
                                else:
                                    maxBand = item[2].stop

                            if item[2].step:
                                strideb = item[2].step

                            bands = range(minBand, maxBand, strideb)
                        elif isinstance(item[2], tuple):
                            for band in item[2]:
                                if band < 0:
                                    bands.append(
                                        self.properties.nrOfBand()+band)
                                else:
                                    bands.append(band)
                        elif isinstance(item[2], int):
                            if item[2] < 0:
                                bands.append(
                                  self.properties.nrOfBand()+item[2])
                            else:
                                bands.append(item[2])
                        else:
                            raise ValueError('Error: band index must be '
                                             'slice, list or integer')
                        retJim = geometry.cropBand(self, band=bands)
                        if maxCol<=minCol or maxRow<=minRow:
                            raise IndexError(
                                'Warning: index error, returning empty Jim')
                        retJim.geometry.crop(ulx=minCol, uly=minRow, ulz=None,
                                             lrx=maxCol, lry=maxRow, lrz=None,
                                             dx=stridex, dy=stridey,
                                             nogeo=True)
                        return retJim
                    else:
                        raise TypeError(
                            'Error: use 3 dimensions when slicing multiband '
                            '2-dim Jim object (x:y:band)')
                else:
                    if len(item) == 2:  # do slice x,y
                        if maxCol<=minCol or maxRow<=minRow:
                            raise IndexError(
                                'Warning: index error, returning empty '
                                'Jim ({}, {}, {}, {})'.format(minCol,
                                                              minRow,
                                                              maxCol,
                                                              maxRow))
                        retJim = geometry.crop(self, ulx=minCol, uly=minRow,
                                               ulz=None, lrx=maxCol,
                                               lry=maxRow, lrz=None,
                                               dx=stridex, dy=stridey,
                                               nogeo=True)
                        return retJim
                    else:
                        raise TypeError('Error: use 2 dimensions when '
                                        'slicing 2-dim Jim object (x:y)')
        elif isinstance(item, Jim):
            mask = item>0
            return Jim(self*mask)
        elif isinstance(item, _jl.VectorOgr):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __getitem__ not implemented for 3d '
                                 'Jim objects')
            nodata = self.properties.getNoDataVals()
            if nodata:
                nodata = nodata[0]
            else:
                nodata = 0
            return geometry.cropOgr(self, item, crop_to_cutline=True,
                                    nodata=nodata, align=True)

    def __setitem__(self, item, value):
        if isinstance(item, Jim):  # or isinstance(value, Jim):
            if value is None:
                self._set(Jim(self, copyData=False))
            else:
                if isinstance(value, Jim):
                    self._jipjim.d_setMask(item._jipjim, value._jipjim)
                else:
                    self._jipjim.d_setMask(item._jipjim, value)
        elif isinstance(item, tuple):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __setitem__ not implemented for 3d '
                                 'Jim objects')
            else:
                bands = []
                if self.properties.nrOfBand() > 1:
                    if len(item) == 3:  # do slice x,y,band
                        stridex = 1
                        stridey = 1
                        strideb = 1
                        if isinstance(item[0], slice):
                            if item[0].start is None:
                                minCol = 0
                            else:
                                minCol =  item[0].start
                            if item[0].stop is None:
                                maxCol = self.properties.nrOfCol()
                            else:
                                maxCol = item[0].stop

                            if item[0].step:
                                stridex = item[0].step

                            if minCol < 0:
                                minCol = self.properties.nrOfCol()+minCol
                            if maxCol < 0:
                                maxCol = self.properties.nrOfCol()+maxCol+1
                        elif isinstance(item[0], int):
                            if item[0] < 0:
                                minCol = self.properties.nrOfCol()+item[0]
                            else:
                                minCol = item[0]
                            maxCol = minCol+1
                        else:
                            raise ValueError('column item must be slice or '
                                             'integer value')
                        if isinstance(item[1], slice):
                            if item[1].start is None:
                                minRow = 0
                            else:
                                minRow =  item[1].start
                            if item[1].stop is None:
                                maxRow = self.properties.nrOfRow()
                            else:
                                maxRow = item[1].stop

                            if item[1].step:
                                stridey = item[1].step

                            if minRow < 0:
                                minRow = self.properties.nrOfRow()+minRow
                            if maxRow < 0:
                                maxRow = self.properties.nrOfRow()+maxRow+1
                        elif isinstance(item[1], int):
                            if item[1] < 0:
                                minRow = self.properties.nrOfRow()+item[1]
                            else:
                                minRow = item[1]
                            maxRow = minRow+1
                        else:
                            raise ValueError('row item must be slice or '
                                             'integer value')
                        if isinstance(item[2], slice):
                              if item[2].step:
                                  strideb = item[2].step
                              if item[2].start < 0:
                                  minBand = self.properties.nrOfBand() + \
                                            item[2].start
                              else:
                                  if item[2].start is None:
                                      minBand = 0
                                  else:
                                      minBand = item[2].start
                              if item[2].stop < 0:
                                  maxBand = self.properties.nrOfBand()+\
                                            item[2].stop
                              else:
                                  if item[2].stop is None:
                                      maxBand = self.properties.nrOfBand()
                                  else:
                                      maxBand = item[2].stop
                              bands = range(minBand, maxBand, strideb)
                        elif isinstance(item[2], tuple):
                            for band in item[2]:
                                if band < 0:
                                    bands.append(
                                        self.properties.nrOfBand()+band)
                                else:
                                    bands.append(band)
                        elif isinstance(item[2], int):
                            if item[2] < 0:
                              bands.append(self.properties.nrOfBand()+item[2])
                            else:
                              bands.append(item[2])
                        if isinstance(value, float):
                            if self.properties.getDataType() != _jl.GDT_Float32\
                                    and self.properties.getDataType() != _jl.GDT_Float64:
                                self.pixops.convert(otype='GDT_Float32')
                        if type(value) in (float, int):
                            self.pixops.setData(value, ulx=minCol, uly=minRow,
                                                lrx=maxCol, lry=maxRow,
                                                bands=bands, dx=stridex,
                                                dy=stridey, nogeo=True)
                        else:
                            raise TypeError(
                                'Error: __setitem__ not implemented for value '
                                'type {}'.format(type(value)))
                    else:
                        raise TypeError(
                            'Error: use 3 dimensions when slicing multiband '
                            '2-dim Jim object (x:y:band)')
                else:
                    stridex = 1
                    stridey = 1
                    if len(item) == 2:  # do slice x,y
                        if isinstance(item[0], slice):
                            if item[0].start is None:
                                minCol = 0
                            else:
                                minCol =  item[0].start
                            if item[0].stop is None:
                                maxCol = self.properties.nrOfCol()
                            else:
                                maxCol = item[0].stop

                            if item[0].step:
                                stridex = item[0].step

                            if minCol<0:
                                minCol=self.properties.nrOfCol()+minCol
                            if maxCol<0:
                                maxCol=self.properties.nrOfCol()+maxCol+1
                        elif isinstance(item[0], int):
                            if item[0]<0:
                                minCol=self.properties.nrOfCol()+item[0]
                            else:
                                minCol = item[0]
                            maxCol = minCol+1
                        else:
                            raise ValueError('column item must be slice or '
                                             'integer value')
                        if isinstance(item[1], slice):
                            minRow = item[1].start
                            maxRow = item[1].stop
                            if item[1].step:
                                stridey = item[1].step
                            if minRow<0:
                                minRow=self.properties.nrOfRow()+minRow
                            if maxRow<0:
                                maxRow=self.properties.nrOfRow()+maxRow+1
                        elif isinstance(item[1], int):
                            if item[1]<0:
                                minRow=self.properties.nrOfRow()+item[1]
                            else:
                                minRow = item[1]
                            maxRow = minRow+1
                        else:
                            raise ValueError('row item must be slice or '
                                             'integer value')
                        if isinstance(value, float):
                            if self.properties.getDataType() != _jl.GDT_Float32\
                                    and self.properties.getDataType() != _jl.GDT_Float64:
                                self.pixops.convert(otype='GDT_Float32')
                        if type(value) in (float, int):
                            bands = [0]
                            self.pixops.setData(value, ulx=minCol, uly=minRow,
                                                lrx=maxCol, lry=maxRow,
                                                bands=bands, dx=stridex,
                                                dy=stridey, nogeo=True)
                        else:
                            raise TypeError(
                                'Error: __setitem__ not implemented for value '
                                'type {}'.format(type(value)))
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing '
                                        '2-dim Jim object (x:y)')
        elif isinstance(item, _jl.VectorOgr):
            # TODO: check if rasterize vector first and then set raster mask
            #       would be better
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __setitem__ not implemented for 3d '
                                 'Jim objects')
            # TODO: decide on default behaviour of ALL_TOUCHED=TRUE
            if type(value) in (float, int):
                nodataValues = self.properties.getNoDataVals()
                self._set(self.setMask(item, {'eo': ['ALL_TOUCHED=TRUE'],
                                              'nodata': value}))
                self.properties.clearNoData()
                self.properties.setNoDataVals(nodataValues)
            elif isinstance(value, Jim):
                templateJim = Jim(self, False)
                templateJim.setData(0)
                templateJim = Jim(templateJim.setMask(
                    item,
                    {'eo': ['ALL_TOUCHED=TRUE'], 'nodata': 1}))
                self[templateJim > 0] = value

    def __nonzero__(self):
        """Check if Jim contains data

        :return: True if image contains data, False if image is empty
        """
        return self._jipjim.isInit()

    def __bool__(self):
        """Check if Jim contains data

        :return: True if image contains data, False if image is empty
        """
        return self._jipjim.isInit()

    def __abs__(self):
        """Calculate the absolute value of Jim raster dataset

        :return: Absolute value of Jim raster dataset
        """
        jim = io.createJim(self)
        jim._jipjim.d_pointOpAbs()
        return Jim(jim)

    def __neg__(self):
        """Calculate the negation of Jim raster dataset

        :return: Negation of Jim raster dataset (-dataset)
        """
        if self.properties.getDataType() == _jl.GDT_Byte:
            print("Warning: converting data type to Int16")
            self.pixops.convert(otype='Int16')
        if self.properties.getDataType() == _jl.GDT_UInt16:
            print("Warning: converting data type to Int32")
            self.pixops.convert(otype='Int32')
        if self.properties.getDataType() == _jl.GDT_UInt32:
            print("Warning: converting data type to Int32, potential overflows"
                  " may occur!")
            self.pixops.convert(otype='Int32')
        return -1 * self

    # *** binary operators *** #
    def __eq__(self, right):
        """Pixel wise check for equality

        :return: Jim object with pixels 1 if equal values, 0 otherwise
        """
        if isinstance(right, Jim):
            return Jim(self._jipjim.eq(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.eq(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        # return Jim(self._jipjim.eq(right._jipjim))

    def __ne__(self, right):
        """Pixel wise check for non-equality

        :return: False if equal values, True otherwise
        """
        if isinstance(right, Jim):
            return Jim(self._jipjim.ne(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.ne(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))

    def __lt__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.lt(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.lt(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        # return Jim(self._jipjim.lt(right._jipjim))

    def __le__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.le(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.le(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        # return Jim(self._jipjim.le(right._jipjim))

    def __gt__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.gt(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.gt(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        # return Jim(self._jipjim.gt(right._jipjim))

    def __ge__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.ge(right._jipjim))
        elif type(right) in (int, float):
            return Jim(self._jipjim.ge(right))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        # return Jim(self._jipjim.ge(right._jipjim))

    def __add__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpArith(right._jipjim, _jl.ADD_op))
        elif type(right) in (int, float):
            return Jim(self._jipjim.pointOpArithCst(right, _jl.ADD_op))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(
                type(right)))

    def __radd__(self, left):
        if isinstance(left, Jim):
            return Jim(self._jipjim.pointOpArith(left._jipjim, _jl.ADD_op))
        elif type(left) in (int, float):
            return Jim(self._jipjim.pointOpArithCst(left, _jl.ADD_op))
        else:
            raise TypeError('unsupported operand type for + : {}'.format(
                type(left)))

    def __iadd__(self, right):
        if isinstance(right, Jim):
            self._jipjim.d_pointOpArith(right._jipjim, _jl.ADD_op)
        elif type(right) in (int, float):
            self._jipjim.d_pointOpArithCst(right, _jl.ADD_op)
        else:
            raise TypeError('unsupported operand type for + : {}'.format(
                type(right)))
        return self

    def __sub__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpArith(right._jipjim, _jl.SUB_op))
        elif type(right) in (int, float):
            return Jim(self._jipjim.pointOpArithCst(right, _jl.SUB_op))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(
                type(right)))

    def __rsub__(self, left):
        if isinstance(left, Jim):
            return -1 * Jim(self._jipjim.pointOpArith(left._jipjim,
                                                      _jl.SUB_op))
        elif type(left) in (int, float):
            return -1 * Jim(self._jipjim.pointOpArithCst(left, _jl.SUB_op))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(
                type(left)))

    def __isub__(self, right):
        if isinstance(right, Jim):
            self._jipjim.d_pointOpArith(right._jipjim, _jl.SUB_op)
        elif type(right) in (int, float):
            self._jipjim.d_pointOpArithCst(right, _jl.SUB_op)
        else:
            raise TypeError('unsupported operand type for - : {}'.format(
                type(right)))
        return self

    def __mul__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpArith(right._jipjim, _jl.MULT_op))
        elif isinstance(right, int):
            return Jim(self._jipjim.pointOpArithCst(right, _jl.MULT_op))
        elif isinstance(right, float):
            if self.properties.getDataType() != _jl.GDT_Float32 and \
                    self.properties.getDataType() != _jl.GDT_Float64:
                self_float = pixops.convert(self, otype='GDT_Float32')
                self_float._jipjim.d_pointOpArithCst(right, _jl.MULT_op)
                return Jim(self_float)
            else:
                return Jim(self._jipjim.pointOpArithCst(right, _jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(
                type(right)))

    def __rmul__(self, left):
        if isinstance(left, Jim):
            return Jim(self._jipjim.pointOpArith(left._jipjim, _jl.MULT_op))
        elif isinstance(left, int):
            return Jim(self._jipjim.pointOpArithCst(left, _jl.MULT_op))
        elif isinstance(left, float):
            if self.properties.getDataType() != _jl.GDT_Float32 and \
                    self.properties.getDataType() != _jl.GDT_Float64:
                self_float = pixops.convert(self, otype='GDT_Float32')
                self_float._jipjim.d_pointOpArithCst(left, _jl.MULT_op)
                return Jim(self_float)
            else:
                return Jim(self._jipjim.pointOpArithCst(left, _jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(
                type(left)))

    def __imul__(self, right):
        if isinstance(right, Jim):
            self._jipjim.d_pointOpArith(right._jipjim, _jl.MULT_op)
        elif type(right) in (int, float):
            self._jipjim.d_pointOpArithCst(right, _jl.MULT_op)
        else:
            raise TypeError('unsupported operand type for * : {}'.format(
                type(right)))
        return self

    def _trueDiv(self, right):
        # test
        # print("true division")
        right_float = Jim(None)
        self_float = Jim(None)
        if self.properties.getDataType() != _jl.GDT_Float32 and \
                self.properties.getDataType() != _jl.GDT_Float64:
            self_float = pixops.convert(self, otype='GDT_Float32')
        if isinstance(right, Jim):
            if right.properties.getDataType() != _jl.GDT_Float32 and \
                    right.properties.getDataType() != _jl.GDT_Float64:
                right_float = pixops.convert(right, otype='GDT_Float32')
                if self_float:
                    return Jim(self_float.pointOpArith(right_float,
                                                       _jl.DIV_op))
                else:
                    return Jim(self._jipjim.pointOpArith(right_float,
                                                         _jl.DIV_op))
            else:
                if self_float:
                    return Jim(self_float.pointOpArith(right._jipjim,
                                                       _jl.DIV_op))
                else:
                    return Jim(self._jipjim.pointOpArith(right._jipjim,
                                                         _jl.DIV_op))
        elif type(right) in (int, float):
            if self_float:
                return Jim(self_float.pointOpArithCst(right, _jl.DIV_op))
            else:
                return Jim(self._jipjim.pointOpArithCst(right, _jl.DIV_op))
        else:
            raise TypeError('unsupported operand type for / : {}'.format(
                type(right)))

    def _itrueDiv(self, right):
        if self.properties.getDataType() != _jl.GDT_Float32 and \
                self.properties.getDataType() != _jl.GDT_Float64:
            self.pixops.convert(otype='GDT_Float32')
        if isinstance(right, Jim):
            self._jipjim.d_pointOpArith(right._jipjim, _jl.DIV_op)
        elif type(right) in (int, float):
            self._jipjim.d_pointOpArithCst(right, _jl.DIV_op)
        else:
            raise TypeError('unsupported operand type for / : {}'.format(
                type(right)))
        return self

    def __div__(self, right):
        _trueDiv = False
        if self.properties.getDataType() == _jl.GDT_Float32 or \
                self.properties.getDataType() == _jl.GDT_Float64:
            _trueDiv = True
        else:
            if isinstance(right, Jim):
                if right.properties.getDataType() == _jl.GDT_Float32 or \
                        right.properties.getDataType() == _jl.GDT_Float64:
                    _trueDiv = True
            elif isinstance(right, float):
                    _trueDiv = True
        if _trueDiv:
            return(self._trueDiv(right))
        else:
            if isinstance(right, Jim):
                return Jim(self._jipjim.pointOpArith(right._jipjim, _jl.DIV_op))
            elif isinstance(right, int):
                return Jim(self._jipjim.pointOpArithCst(right, _jl.DIV_op))
            else:
                raise TypeError('unsupported operand type for / : {}'.format(
                    type(right)))

    def __idiv__(self, right):
        if self.properties.getDataType() != _jl.GDT_Float32 and \
                self.properties.getDataType() != _jl.GDT_Float64:
            self.pixops.convert(otype='GDT_Float32')
        if isinstance(right, Jim):
            self._jipjim.d_pointOpArith(right._jipjim, _jl.DIV_op)
        elif type(right) in (int, float):
            self._jipjim.d_pointOpArithCst(right, _jl.DIV_op)
        else:
            raise TypeError('unsupported operand type for / : {}'.format(
                type(right)))
        return self

    def __mod__(self, right):
        if isinstance(right, int):
            return Jim(self._jipjim.pointOpModulo(right))
        else:
            raise TypeError('unsupported operand type for % : {}'.format(
                type(right)))

    def __imod__(self, right):
        if isinstance(right, int):
            self._jipjim.d_pointOpModulo(right)
        else:
            raise TypeError('unsupported operand type for % : {}'.format(
                type(right)))
        return self

    def __lshift__(self, right):
        if isinstance(right, int):
            jim = io.createJim(self)
            for iband in range(0, self._jipjim.properties.nrOfBand()):
                jim.d_pointOpBitShift(-right, iband)
            return Jim(jim)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))

    def __ilshift__(self, right):
        if isinstance(right, int):
            for iband in range(0, self._jipjim.properties.nrOfBand()):
                self._jipjim.d_pointOpBitShift(-right, iband)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))
        return self

    def __rshift__(self, right):
        if isinstance(right, int):
            jim = io.createJim(self)
            for iband in range(0, self._jipjim.properties.nrOfBand()):
                jim.d_pointOpBitShift(right, iband)
            return Jim(jim)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))

    def __irshift__(self, right):
        if isinstance(right, int):
            for iband in range(0, self._jipjim.properties.nrOfBand()):
                self._jipjim.d_pointOpBitShift(right, iband)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))
        return self

    def __or__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpBitwise(right._jipjim, _jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(right)))

    def __ror__(self, left):
        if isinstance(left, Jim):
            return Jim(self._jipjim.pointOpBitwise(left._jipjim, _jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(left)))

    def __xor__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpBitwise(right._jipjim, _jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(right)))

    def __rxor__(self, left):
        if isinstance(left, Jim):
            return Jim(self._jipjim.pointOpBitwise(left._jipjim, _jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(left)))

    def __and__(self, right):
        if isinstance(right, Jim):
            return Jim(self._jipjim.pointOpBitwise(right._jipjim, _jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(
                type(right)))

    def __rand__(self, left):
        if isinstance(left, Jim):
            return Jim(self._jipjim.pointOpBitwise(left._jipjim, _jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(
                type(left)))


class _ParentList(_jl.JimList):

    def __init__(self, images_list, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        jiplib_images_list = [i._jipjim for i in images_list]
        super(_ParentList, self).__init__(jiplib_images_list, *args)


class JimList(list):
    """Definition of JimList object."""

    def __init__(self, images_list, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        super(JimList, self).__init__(images_list)
        self._jipjimlist = _ParentList(images_list, *args)

        self._all = all._AllList()
        self._ccops = ccops._CCOpsList()
        self._clssfy = clssfy._ClassifyList()
        self._demops = demops._DEMOpsList()
        self._geometry = geometry._GeometryList()
        self._io = io._IOList()
        self._ngbops = ngbops._NgbOpsList()
        self._pixops = pixops._PixOpsList()
        self._properties = properties._PropertiesList()
        self._stats = stats._StatsList()

    @property
    def geometry(self):
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    @property
    def pixops(self):
        self._pixops._set_caller(self)
        _gc.collect()
        return self._pixops

    @property
    def properties(self):
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

    @property
    def stats(self):
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
        for module in [geometry._GeometryList, pixops._PixOpsList,
                       properties._PropertiesList, stats._StatsList]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def append(self, jim_object):
        """Add single element to the JimList."""
        assert isinstance(jim_object, Jim), \
            'Only Jim instances can be appended'
        super(JimList, self).append(jim_object)
        self._set(self, from_list=True)

    def count(self, jim_object):
        """Count the occurrences of element in the JimList.."""
        i = 0
        for jim in self:
            if jim.pixops.isEqual(jim_object):
                i += 1

        return i

    def extend(self, jim_list):
        """Add elements of a JimList to another JimList."""
        assert isinstance(jim_list, JimList), \
            'Only JimList instances can be used to extend JimList'
        super(JimList, self).extend(jim_list)
        self._set(self, from_list=True)

    def index(self, jim_object):
        """Return smallest index of element in the JimList."""
        for i in range(len(self)):
            if self[i].pixops.isEqual(jim_object):
                return i

    def insert(self, index, jim_object):
        """Insert elements to the JimList."""
        assert isinstance(jim_object, Jim), \
            'Only Jim instances can be inserted'
        super(JimList, self).insert(index, jim_object)
        self._set(self, from_list=True)

    def pop(self, index):
        """Remove and return element at given index."""
        popped = super(JimList, self).pop(index)
        self._set(self, from_list=True)
        return popped

    def remove(self, jim_object):
        """Remove the first occurence of an element from the JimList."""
        for i in range(len(self)):
            if self[i].pixops.isEqual(jim_object):
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
                if isinstance(im, Jim):
                    self.append(im)
                else:
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

    def __init__(self, vector, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        super(_ParentVect, self).__init__(vector, *args)


class JimVect():
    """Definition of JimVect object."""

    def __init__(self, vector, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        self._jipjimvect = _ParentVect(vector, *args)

        self._all = all._AllVect()
        self._ccops = ccops._CCOpsVect()
        self._clssfy = clssfy._ClassifyVect()
        self._demops = demops._DEMOpsVect()
        self._geometry = geometry._GeometryVect()
        self._io = io._IOVect()
        self._ngbops = ngbops._NgbOpsVect()
        self._pixops = pixops._PixOpsVect()
        self._properties = properties._PropertiesVect()
        self._stats = stats._StatsVect()

    @property
    def clssfy(self):
        self._clssfy._set_caller(self)
        _gc.collect()
        return self._clssfy

    @property
    def io(self):
        self._io._set_caller(self)
        _gc.collect()
        return self._io

    @property
    def properties(self):
        self._properties._set_caller(self)
        _gc.collect()
        return self._properties

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
        for module in [clssfy._ClassifyVect, io._IOVect,
                       properties._PropertiesVect]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))
