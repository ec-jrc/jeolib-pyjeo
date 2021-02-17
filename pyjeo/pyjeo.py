"""Basic file containing Jim* objects and functions for numpy conversions."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
#
# Copyright (C) 2018-2020 European Union (Joint Research Centre)
#
# This file is part of pyjeo.
#
# pyjeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyjeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import division
import numpy as _np
import gc as _gc
import warnings as _warnings
import os as _os
import math
from pathlib import Path
from osgeo import ogr as _ogr

import jiplib as _jl

from .modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, classify, demops, stats, all
from . import exceptions
from .__init__ import _check_graph

if '__del__ ' in dir(_jl.Jim):
    del _jl.Jim.__del__


class _ParentJim(_jl.Jim):

    def __init__(self, image, kwargs):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        if kwargs:
            if image:
                if isinstance(image, Jim):
                    if 'copy_data' in kwargs.keys():
                        super(_ParentJim, self).__init__(image._jipjim,
                                                         kwargs['copy_data'])
                    else:
                        _warnings.warn(
                            'Not possible to create Jim image based on another'
                            ' one together with other kwargs than copy_data. '
                            'kwargs ignored.', SyntaxWarning)
                        super(_ParentJim, self).__init__(image._jipjim)
                else:
                    if 'tileindex' in kwargs.keys() and 'tiletotal' in kwargs.keys():
                        tileindex = kwargs.pop('tileindex')
                        tiletotal = kwargs.pop('tiletotal')
                        try:
                            overlap = kwargs.pop('overlap')
                        except KeyError:
                            overlap = 5
                        ajim = Jim(image, noread=True)
                        bbox = ajim.properties.getBBox()
                        ncol = int(round(math.sqrt(tiletotal)))
                        nrow = ncol
                        assert ncol * nrow == tiletotal, \
                            'Error: tiletotal must be squared integer'
                        assert tileindex < tiletotal, \
                            'Error: tileindex must be < ' + str(tiletotal)
                        assert len(bbox) == 4, \
                            'Error: bbox must be list of format ulx, uly, lrx, lry'
                        ncol = math.sqrt(tiletotal)
                        nrow = ncol
                        icol = tileindex % nrow
                        irow = tileindex // nrow
                        dx = (bbox[2] - bbox[0]) / ncol
                        dy = (bbox[1] - bbox[3]) / nrow
                        # overlap = dx*0.025
                        overlap = dx*overlap/200.0
                        ulx = bbox[0] + icol * dx - overlap
                        ulx = max(ulx, bbox[0])
                        uly = bbox[1] - irow * dy + overlap
                        uly = min(uly, bbox[1])
                        lrx = ulx + dx + overlap
                        lrx = min(lrx, bbox[2])
                        lry = uly - dy - overlap
                        lry = max(lry, bbox[3])
                        kwargs.update({'ulx':ulx})
                        kwargs.update({'uly':uly})
                        kwargs.update({'lrx':lrx})
                        kwargs.update({'lry':lry})
                    elif 'bbox' in kwargs:
                        bbox = kwargs.pop('bbox')
                        kwargs.update({'ulx':bbox[0]})
                        kwargs.update({'ulx':bbox[0]})
                        kwargs.update({'uly':bbox[1]})
                        kwargs.update({'lrx':bbox[2]})
                        kwargs.update({'lry':bbox[3]})
                    if isinstance(image, Path):
                        kwargs.update({'filename': str(image)})
                    else:
                        kwargs.update({'filename': image})
                    super(_ParentJim, self).__init__(kwargs)
            elif 'graph' in kwargs:
                graph = kwargs.pop('graph')
                _check_graph(graph, [4, 8])

                if graph == 4:
                    ngb = Jim(ncol=3, nrow=3, otype='Byte')
                    ngb[0, 1] = 1
                    ngb[1, 0] = 1
                    ngb[1, 2] = 1
                    ngb[2, 1] = 1
                else:
                    ngb = Jim(ncol=3, nrow=3, otype='Byte')
                    ngb.pixops.setData(1)
                    ngb[1, 1] = 0

                super(_ParentJim, self).__init__(ngb._jipjim)
            else:
                if 'bbox' in kwargs:
                    bbox = kwargs.pop('bbox')
                    kwargs.update({'ulx':bbox[0]})
                    kwargs.update({'ulx':bbox[0]})
                    kwargs.update({'uly':bbox[1]})
                    kwargs.update({'lrx':bbox[2]})
                    kwargs.update({'lry':bbox[3]})
                super(_ParentJim, self).__init__(kwargs)
        elif image:
            if isinstance(image, Jim):
                super(_ParentJim, self).__init__(image._jipjim)
            elif isinstance(image, _jl.Jim):
                super(_ParentJim, self).__init__(image)
            elif isinstance(image, Path):
                super(_ParentJim, self).__init__({'filename': str(image)})
            else:
                super(_ParentJim, self).__init__({'filename': image})
        else:
            if 'bbox' in kwargs:
                bbox = kwargs.pop('bbox')
                kwargs.update({'ulx':ulx})
                kwargs.update({'uly':uly})
                kwargs.update({'lrx':lrx})
                kwargs.update({'lry':lry})
            super(_ParentJim, self).__init__(image)


class Jim:
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

    @staticmethod
    def getMethods(queried_module: str = None):
        """Print an overview of available methods in format module.method."""
        def tree_structure(module, queried_module: str):
            """Change structure of the output to a tree one."""
            if queried_module and queried_module not in str(module):
                return ''

            mm = dir(module)
            module_methods = list(mm)
            for method in mm:
                if method[0] == '_':
                    module_methods.remove(method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:])] + \
                   module_methods

        methods = list()

        for module in [properties._Properties, io._IO, pixops._PixOps,
                       ngbops._NgbOps, geometry._Geometry, ccops._CCOps,
                       classify._Classify, demops._DEMOps, stats._Stats]:
            methods.extend(tree_structure(module, queried_module))

        print('\n'.join(methods))

    def np(self, band: int = 0):
        """Return numpy array from Jim object.

        The (contiguous) data array is organized
        as [planes][rows][columns] for the respective bands

        :param band: band index (starting from 0)
        :return: numpy array representation

        Some Numpy algorithms expect a 3D image that is organized as
        [rows][columns][planes] instead of [planes][rows][columns].
        Also to plot an RGB image using matplotlib, the axes of the Numpy array
        must be rolled::

            import matplotlib.pyplot as plt

            jim = pj.Jim('/path/to/rgb.tif')
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.imshow(np.rollaxis(jim.np(),0,3))
        """

        if not self:
            raise exceptions.JimEmptyError(
                'Jim has to have a data to use Jim.np()')
        try:
            self.properties.getDataType()
        except TypeError:
            raise exceptions.JimInnerParametersError(
                'Jim has to have a data type and dimensions to use Jim.np()')
        if band >= self.properties.nrOfBand():
            raise exceptions.JimBandsError('Band out of bounds')
        elif band < 0:
            band = self.properties.nrOfBand() + band
        return _jl.jim2np(self._jipjim, band, False)

    def xr(self):
        """Return xarray from Jim object.

        The (contiguous) data array is organized
        as [planes][rows][columns] for the respective bands

        :return: xarray representation
        """
        import xarray as _xr

        #jim is a multiband datacube (with multiple planes)
        planes = range(self.properties.nrOfPlane())
        bands = range(self.properties.nrOfBand())
        bbox = self.properties.getBBox()

        # xarray coordinates are pixel centered
        x = _np.arange(bbox[0]+self.properties.getDeltaX()/2,
                        bbox[2]+self.properties.getDeltaX()/2,
                        self.properties.getDeltaX())
        y = _np.arange(bbox[1]-self.properties.getDeltaY()/2,
                        bbox[3]-self.properties.getDeltaY()/2,
                        -self.properties.getDeltaY())


        # Build a xarray Dataset reference (without memory copy)
        # Do not alter shape or destroy x_dataset!
        if self.properties.nrOfPlane() > 1:
            x_dataset = _xr.Dataset({str(b):_xr.DataArray(self.np(b),
                                                        dims=['time', 'y', 'x'],
                                                        coords={'time': list(planes),
                                                        'x': x, 'y': y},
                                                        attrs={'_FillValue': 0})
                                    for b in bands})
        else:
            x_dataset = _xr.Dataset({str(b):_xr.DataArray(self.np(b),
                                                        dims=['y', 'x'],
                                                        coords={'x': x, 'y': y},
                                                        attrs={'_FillValue': 0})
                                    for b in bands})
        return x_dataset


    @staticmethod
    def _checkInitParamsSense(image, kwargs):
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
                    raise exceptions.JimIllegalArgumentError(
                        'You cannot use any of parameters '
                        '[seed, uniform, stdev] without specifying '
                        'the geometry and otype of Jim.')
        elif type(image) is str:
            non_standard_path = image[:4] == '/vsi' or ':' in image

            if not non_standard_path and not _os.path.isfile(image):
                raise ValueError('File does not exist')

    def _checkNumberOfBands(self, another_jim):
        """Check if the number of bands matches for both Jim objects.

        Raise a Warning otherwise.

        :param another_jim: Another Jim object which nrOfBands is going to be
            checked
        """
        if self.properties.nrOfBand() > another_jim.properties.nrOfBand() > 1:
            raise exceptions.JimBandsError(
                'Jims used in conditions must have either the same number of '
                'bands or one of them must have number of bands == 1')
        elif 1 < self.properties.nrOfBand() < \
                another_jim.properties.nrOfBand():
            _warnings.warn('Left Jim has less bands than the right one and '
                           'more than one band. Only corresponding bands of '
                           'the right Jim will be taken in consideration.',
                           Warning)

    def _feed(self, stdev: int, uniform, seed: int, kwargs):
        """Feed the Jim object with either uniform or random seed of values.

        :param stdev: standard deviation
        :param uniform: if uniform distr, here are defined values [min, max]
        :param seed: numpy random seed
        """
        mean = kwargs.pop('mean', None)

        if seed is not None:
            _np.random.seed(seed)

        if uniform:
            if isinstance(uniform, list):
                if len(uniform) != 2:
                    raise exceptions.JimIllegalArgumentError(
                        'The list parsed as the uniform argument must be '
                        'in the form [min, max + 1]')
                for band in range(0, self.properties.nrOfBand()):
                    self.np(band)[:] = _np.random.uniform(
                        uniform[0], uniform[1], self.np(band).shape)
            else:
                for band in range(0, self.properties.nrOfBand()):
                    self.np(band)[:] = _np.random.uniform(
                        0, uniform, self.np(band).shape)
        else:
            if stdev is None:
                stdev = 1
            if mean is None:
                mean = 0
            for band in range(0, self.properties.nrOfBand()):
                self.np(band)[:] = _np.random.normal(
                    mean, stdev, self.np(band).shape)

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
                raise exceptions.JimIllegalArgumentError(
                    'Using a JimVect as an index not implemented for Jim '
                    'objects with more than one plane')
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
            nband = self.properties.nrOfBand()
            for band in range(nband):
                npresult = _np.array(self.np(band)[item])
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
                else:
                    raise exceptions.JimIllegalArgumentError(
                        'Len of shape {} of the np object not supported. Supported'
                        ' lengths are [0, 1, 2, 3]'.format(len(npresult.shape)))

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

                ul = self.geometry.image2geo(cropuli, cropulj)
                if band == 0:
                    result = Jim(ncol=ncol, nrow=nrow, nband=nband, nplane=nplane,
                                otype=self.properties.getDataType())
                    result.properties.setProjection(self.properties.getProjection())
                    gt = self.properties.getGeoTransform()

                    cropulx = ul[0]-self.properties.getDeltaX()/2
                    cropuly = ul[1]+self.properties.getDeltaY()/2

                    gt[0] = cropulx
                    gt[1] = dx
                    gt[2] = 0
                    gt[3] = cropuly
                    gt[4] = 0
                    gt[5] = -dy
                    result.properties.setGeoTransform(gt)
                result.np(band)[:] = npresult
            return result

    def __setitem__(self, item, value):
        """Evaluate operation of type self[key] = value.

        :param item: key/index
        :param value: value to save
        """
        if isinstance(item, JimVect):
            if self.properties.nrOfPlane() > 1:
                raise exceptions.JimIllegalArgumentError(
                    '__setitem__ with JimVect not implemented for Jim objects '
                    'with more than one plane')
            # TODO: decide on default behaviour of ALL_TOUCHED=TRUE
            # TODO: next lines should work, but problem with GML files when SRS
            #       is not defined as in S2 cloud masks

            # template = Jim(self)
            # template.geometry.rasterize(item, 1.0)
            # self[template > 0] = value

            if type(value) in (float, int) or isinstance(value, Jim):
                template_jim = Jim(self, copy_data=False)
                template_jim._jipjim.d_setMask(
                    item._jipjimvect, {'eo': ['ALL_TOUCHED=TRUE'],
                                       'nodata': 1})
                self[template_jim > 0] = value
        elif isinstance(item, Jim):
            if value is not None:
                if isinstance(value, Jim):
                    self._jipjim.d_setMask(item._jipjim, value._jipjim)
                else:
                    self._jipjim.d_setMask(item._jipjim, value)
        elif isinstance(item, tuple):
            for band in range(self.properties.nrOfBand()):
                if isinstance(value, Jim):
                    self.np(band)[item] = value.np(band)
                else:
                    print(item)
                    print(self.np())
                    print("value: {}".format(value))
                    self.np(band)[item] = value
        else:
            raise exceptions.JimIllegalArgumentError(
                '__setitem__ only implemented for JimVector, Jim or tuples '
                '(dims of __getitem__ must correspond with those of the tuple)'
            )

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
            jim.np(iband)[:] = _np.abs(jim.np(iband))
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
            raise exceptions.JimIllegalArgumentError('You cannot divide a Jim '
                                                     'object of int type')
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
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for << : {}'.format(
                    type(right)))

    def __ilshift__(self, right):
        """Pixel wise operation <<=."""
        if isinstance(right, int):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] <<= right
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for <<= : {}'.format(
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
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for >> : {}'.format(
                    type(right)))

    def __irshift__(self, right):
        """Pixel wise operation >>=."""
        if isinstance(right, int):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] >>= right
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for >>= : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] |= right
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for | : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] |= right
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for |= : {}'.format(
                    type(right)))
        return self

    def __ror__(self, left):
        """Pixel wise operation | where self is the right object."""
        if isinstance(left, int) or isinstance(left, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] |= left
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for | : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] ^= right
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for ^ : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] ^= right
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for ^= : {}'.format(
                    type(right)))
        return self

    def __rxor__(self, left):
        """Pixel wise operation ^ where self is the right object."""
        if isinstance(left, int) or isinstance(left, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] ^= left
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for ^ : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] &= right
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for & : {}'.format(
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
        elif isinstance(right, int) or isinstance(right, _np.ndarray):
            for iband in range(0, self.properties.nrOfBand()):
                self.np(iband)[:] &= right
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for &= : {}'.format(
                    type(right)))
        return self

    def __rand__(self, left):
        """Pixel wise operation & where self is the right object."""
        if isinstance(left, int) or isinstance(left, _np.ndarray):
            jim = Jim(self)
            for iband in range(0, self.properties.nrOfBand()):
                jim.np(iband)[:] &= left
            return jim
        else:
            raise exceptions.JimIllegalArgumentError(
                'unsupported operand type for & : {}'.format(
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

    def __init__(self, images_list=None, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
            the Jim object
        """
        if images_list is None:
            images_list = []

        if isinstance(images_list, Jim):
            images_list = [images_list]
        elif not isinstance(images_list, list):
            raise exceptions.JimListIllegalArgumentError(
                'Argument images_list must be either list, JimList or a Jim '
                'object.')

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

    @staticmethod
    def getMethods(queried_module: str = None):
        """Print an overview of available methods in format module.method."""
        def tree_structure(module, queried_module: str):
            """Change structure of the output to a tree one."""
            if queried_module and queried_module not in str(module):
                return ''

            mm = dir(module)
            module_methods = list(mm)
            for method in mm:
                if method[0] == '_':
                    module_methods.remove(method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:-4] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:-4])] + \
                   module_methods

        methods = list()

        for module in [geometry._GeometryList, io._IOList, pixops._PixOpsList,
                       properties._PropertiesList, stats._StatsList]:
            methods.extend(tree_structure(module, queried_module))

        print('\n'.join(methods))

    def append(self, jim: Jim):
        """Add single element to the JimList."""
        assert isinstance(jim, Jim), \
            'Only Jim instances can be appended'
        super(JimList, self).append(jim)
        self._set(self, from_list=True)

    def count(self, jim: Jim):
        """Count the occurrences of element in the JimList.."""
        i = 0
        for jim_object in self:
            if jim.properties.isEqual(jim_object):
                i += 1

        return i

    def extend(self, jim_list):
        """Add elements of a JimList to another JimList."""
        assert isinstance(jim_list, JimList), \
            'Only JimList instances can be used to extend JimList'
        super(JimList, self).extend(jim_list)
        self._set(self, from_list=True)

    def index(self, jim: Jim, start: int = 0, end: int = None):
        """Return smallest index of element in the JimList."""
        if end is None:
            end = len(self)

        for i in range(start, end):
            if self[i].properties.isEqual(jim):
                return i

    def insert(self, index: int, jim: Jim):
        """Insert elements to the JimList."""
        assert isinstance(jim, Jim), \
            'Only Jim instances can be inserted'
        super(JimList, self).insert(index, jim)
        self._set(self, from_list=True)

    def pop(self,
            index: int = -1):
        """Remove and return element at given index."""
        popped = super(JimList, self).pop(index)
        self._set(self, from_list=True)
        return popped

    def remove(self, jim: Jim):
        """Remove the first occurrence of an element from the JimList."""
        for i in range(len(self)):
            if self[i].properties.isEqual(jim):
                del self[i]
                break
        self._set(self, from_list=True)

    def reverse(self):
        """Reverse the JimList."""
        super(JimList, self).reverse()
        self._set(self, from_list=True)

    def _set(self, modified_list, from_list: bool = False):
        """Apply changes done in modified_list to the parent JimList instance.

        :param modified_list: modified JimList instance
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
        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | output           | Path to an output vector in case another         |
        |                  | JimVect object was provided (copy constructor)   |
        +------------------+--------------------------------------------------+
        """
        self._checkInitParamsSense(vector, kwargs)

        if kwargs:
            if vector:
                if isinstance(vector, JimVect):
                    kwargs.update({'filename': kwargs.pop('output', None)})
                    super(_ParentVect, self).__init__(vector._jipjimvect,
                                                      kwargs)
                else:
                    if isinstance(vector, Path):
                        kwargs.update({'filename': str(vector)})
                    else:
                        kwargs.update({'filename': vector})
                    if 'bbox' in kwargs:
                        bbox = kwargs.pop('bbox')
                        kwargs.update({'ulx':bbox[0]})
                        kwargs.update({'uly':bbox[1]})
                        kwargs.update({'lrx':bbox[2]})
                        kwargs.update({'lry':bbox[3]})
                    super(_ParentVect, self).__init__(kwargs)
            elif 'wkt' in kwargs:
                geom = _ogr.CreateGeometryFromWkt(kwargs.pop('wkt'))
                super(_ParentVect, self).__init__(geom.ExportToJson())
            else:
                filename = kwargs.pop('output', None)
                if isinstance(filename, Path):
                    kwargs.update({'filename': str(filename)})
                else:
                    kwargs.update({'filename': filename})
                super(_ParentVect, self).__init__(kwargs)
        else:
            if vector:
                if isinstance(vector, Path):
                    super(_ParentVect, self).__init__(str(vector))
                else:
                    super(_ParentVect, self).__init__(vector)
            else:
                super(_ParentVect, self).__init__()

    @staticmethod
    def _checkInitParamsSense(vector, kwargs):
        """Check if the combination of (kw)args for JimVect init makes sense.

        :param vector: path to a vector or another JimVect object as a basis
            for the JimVect object
        """
        keys = kwargs.keys()

        if isinstance(vector, JimVect):
            if 'output' not in keys:
                raise exceptions.JimVectIllegalArgumentError(
                    "Parameter output required for copy constructor")
        elif not vector:
            if 'output' in keys and not _os.path.isfile(kwargs['output']):
                raise exceptions.JimVectIllegalArgumentError(
                    'Output path does not exist and the template vector is '
                    'not specified')


class JimVect:
    """Definition of JimVect object."""

    def __init__(self, vector=None, **kwargs):
        """Create an empty VectorOgr object.

        Created object is an instance of the basis vector class of the pyJEO
        library

        :param vector: Path to a vector dataset
        :return: a JimVect object
        """
        if isinstance(vector, Path):
            self._jipjimvect = _ParentVect(str(vector), kwargs)
        else:
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

    @staticmethod
    def getMethods(queried_module: str = None):
        """Print an overview of available methods in format module.method."""
        def tree_structure(module, queried_module: str):
            """Change structure of the output to a tree one."""
            if queried_module and queried_module not in str(module):
                return ''

            mm = dir(module)
            module_methods = list(mm)
            for method in mm:
                if method[0] == '_':
                    module_methods.remove(method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:-4] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:-4])] + \
                   module_methods

        methods = list()
        for module in [classify._ClassifyVect,
                       geometry._GeometryVect,
                       io._IOVect,
                       properties._PropertiesVect]:
            methods.extend(tree_structure(module, queried_module))

        print('\n'.join(methods))

    def _set(self, modified_object):
        """Apply changes done in modified_object to parent VectorOgr instance.

        :param modified_object: modified VectorOgr instance
        """
        self._jipjimvect.__dict__.update(modified_object.__dict__)

    def np(self, field: list = None, ln: int = 0):
        """Return numpy array from JimVect object.

        :param field: list of fields to consider
        :param ln: Layer to return
        :return: 2D numpy array representation of numerical fields of all features
        """
        if field is None:
            field = []
        if self.properties.getFeatureCount():
            return _jl.jimvect2np(self._jipjimvect, field, ln)
        else:
            return None

    def dict(self, field: list = None, ln: int = 0):
        """Return dictionary from JimVect object.

        :param field: list of fields to return
        :param ln: Layer to return
        :return: dictionary object representation of fields of all features

        Example: create a dictionary object and turn into pandas object

            import pandas as pd
            v = pj.Jim('/path/to/vector.sqlite')

            pob = pd.DataFrame(v.dict())
        """

        if field is None:
            field = []
        if self.properties.getFeatureCount():
            return _jl.dict(self._jipjimvect, field, ln)
        else:
            return None


def jim2np(jim_object: Jim,
           band: int = 0,
           copy_data: bool = True) -> _np.ndarray:
    """Return a numpy representation of a Jim object.

    :param jim_object: Jim object to be converted
    :param band: band of Jim object to be converted
    :param copy_data: Set to False if reference image is used as a template
        only, without copying actual pixel dat
    :return: a numpy representation of the Jim object
    """

    if isinstance(jim_object, Jim):
        if not jim_object:
            raise exceptions.JimEmptyError(
                'Jim has to have a data to use Jim.np()')
        return _jl.jim2np(jim_object._jipjim, band, copy_data)
    elif isinstance(jim_object, JimVect):
        raise exceptions.JimVectNotSupportedError(
            'input must be of Jim type, use jimvect2np for JimVect')
    else:
        raise exceptions.JimError(
            'input must be of Jim type')


def jim2xr(jim_object: Jim,
           copy_data: bool = True):
    """Return xarray from Jim object.

    The (contiguous) data array is organized
    as [planes][rows][columns] for the respective bands

    :return: xarray representation
    """
    import xarray as _xr

    #jim_object is a multiband datacube (with multiple planes)
    planes = range(jim_object.properties.nrOfPlane())
    bands = range(jim_object.properties.nrOfBand())
    bbox = jim_object.properties.getBBox()

    # xarray coordinates are pixel centered
    x = _np.arange(bbox[0]+jim_object.properties.getDeltaX()/2,
                    bbox[2]+jim_object.properties.getDeltaX()/2,
                    jim_object.properties.getDeltaX())
    y = _np.arange(bbox[1]-jim_object.properties.getDeltaY()/2,
                    bbox[3]-jim_object.properties.getDeltaY()/2,
                    -jim_object.properties.getDeltaY())
    # Build new copy of xarray Dataset (with memory copy):
    if jim_object.properties.nrOfPlane() > 1:
        new_dataset = _xr.Dataset({str(b):xr.DataArray(_pj.jim2np(jim_object,b),
                                                        dims=['time', 'y', 'x'],
                                                        coords={'time': ['t'+str(plane)
                                                        for plane in planes],
                                                        'x': x, 'y': y},
                                                        attrs={'_FillValue': 0})
                                for b in bands})
    else:
        x_dataset = _xr.Dataset({str(b):_xr.DataArray(_pj.jim2np(jim_object,b),
                                                    dims=['y', 'x'],
                                                    coords={'x': x, 'y': y},
                                                    attrs={'_FillValue': 0})
                                for b in bands})


def jimvect2np(jim_object: JimVect,
               field: list = None,
               ln: int = 0) -> _np.ndarray:
    """Return a numpy representation of a JimVect object.

    :param jim_object: JimVect object to be converted
        :param field: list of fields to consider
        :param ln: Layer to return
    :return: a numpy representation of the JimVect object
    """

    if isinstance(jim_object, JimVect):
        if field is None:
            field = []
        if jim_object.properties.getFeatureCount():
            return _jl.jimvect2np(jim_object._jipjimvect, field, ln)
        else:
            return None
    elif isinstance(jim_object, Jim):
        raise exceptions.JimNotSupportedError(
            'input must be of JimVect type, use jim2np for Jim')
    else:
        raise exceptions.JimVectError(
            'input must be of JimVect type')


def np2jim(np_object: _np.ndarray) -> Jim:
    """Return a Jim representation of a numpy array.

    The (contiguous) data array must be organized
    as [planes][rows][columns] for the respective bands

    :param np_object: a numpy array
    :return: a Jim representation of a numpy array
    """
    return Jim(_jl.np2jim(np_object))


def xr2jim(xr_object) -> Jim:
    """Return a Jim representation from an xarray.

    The (contiguous) data array must be organized
    as [planes][rows][columns] for the respective bands

    Returns a 3D Jim object with planes and bands
    xarray must be a cube with identical x and y coordinates for each dataset
    Notice the (contiguous) data array is organized as [planes][rows][columns]
    for the respective bands

    :param xr_object: an xarray
    :return: a Jim representation from  an xarray
    """
    import xarray as _xr

    jim = None
    projection = None
    for b in xr_object:
        if xr_object[b].attrs.get('spatial_ref') is not None:
            projection = xr_object[b].attrs.get('spatial_ref')
        elif jim is None:
            jim = np2jim(xr_object[b].values)
            gt = []
            dx = xr_object[b].coords['x'].values[1]-xr_object[b].coords['x'].values[0]
            dy = xr_object[b].coords['y'].values[0]-xr_object[b].coords['y'].values[1]
            ulx = xr_object[b].coords['x'].values[0]-dx/2.0
            uly = xr_object[b].coords['y'].values[0]+dy/2.0
            gt.append(ulx)
            gt.append(dx)
            gt.append(0)
            gt.append(uly)
            gt.append(0)
            gt.append(-dy)
            jim.properties.setGeoTransform(gt)
        else:
            jim.geometry.stackBand(np2jim(xr_object[b].values))
    return jim
