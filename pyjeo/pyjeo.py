"""Basic file containing Jim and JimList objects."""

from __future__ import division
import numpy
import gc as _gc

try:
    import jiplib as _jl
except ImportError:
    from jeodpp import jiplib as _jl

from modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, classify, demops, stats, all


del _jl.Jim.__del__

def np(aJim):
    return _jl.np(aJim._jipjim)


def jim2np(aJim, band=0, copyData=True):
    return _jl.jim2np(aJim._jipjim, band, copyData)


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
                    # kwargs.update({'band2plane':True})
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
        if 'stdev' in kwargs.keys() or 'uniform' in kwargs.keys():
            #remove stdev and uniform from kwargs
            stdev=kwargs.pop('stdev',None)
            uniform=kwargs.pop('uniform',None)
            seed=kwargs.pop('seed',None)
        self._jipjim = _ParentJim(image, kwargs)

        if stdev or uniform or seed:
            mean=kwargs.pop('mean',None)
            if seed:
                numpy.random.seed(seed)
            if len(uniform) == 2:
                min=uniform[0]
                max=uniform[1]
                scale=max-min
                offset=min
                self.np()[:]=numpy.random.rand(*(self.np().shape))
            else:
                if not stdev:
                    stdev=1
                if not mean:
                    mean=0
                scale=stdev
                offset=mean
            self.np()[:]=numpy.random.rand(*(self.np().shape))
            self*=scale
            self+=offset

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
    def classify(self):
        self._classify._set_caller(self)
        _gc.collect()
        return self._classify

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
                       classify._Classify, demops._DEMOps, stats._Stats]:
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def np(self):
        return _jl.np(self._jipjim)

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self._jipjim.__dict__.update(modified_object.__dict__)

    # *** unary operators *** #

    def __getitem__(self, item):
        if isinstance(item, JimVect):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __getitem__ not implemented for 3d '
                                 'Jim objects')
            nodata = self.properties.getNoDataVals()
            if nodata:
                nodata = nodata[0]
            else:
                nodata = 0
            return geometry.cropOgr(self, item._jipjimvect, crop_to_cutline=True, nodata=nodata, align=True)

        elif isinstance(item, Jim):
            mask = item>0
            return Jim(self*mask)
        else:
            npresult=numpy.array(self.np()[item],copy=True)
            # npresult=numpy.array(self.np()[item])
            if len(npresult.shape)==3:
                nplane=npresult.shape[0]
                nrow=npresult.shape[1]
                ncol=npresult.shape[2]
            elif len(npresult.shape)==2:
                nplane=1
                nrow=npresult.shape[0]
                ncol=npresult.shape[1]
            elif len(npresult.shape)==1:
                nplane=1
                nrow=1
                ncol=npresult.shape[0]
            elif len(npresult.shape)==0:
                nplane=1
                nrow=1
                ncol=1
            else:
                raise IndexError('Error: index in __getitem__ out of range')

            #[gs]item only supports single band image (use plane instead)
            nband=1
            if self.properties.nrOfPlane()>1:
                dim=3
            else:
                dim=2
            dx=self.properties.getDeltaX()
            dy=self.properties.getDeltaY()

            cropuli=0
            # croplri=self.properties.nrOfCol()
            cropulj=0
            # croplrj=self.properties.nrOfRow()
            if isinstance(item, tuple):
                #cols
                if len(item)>dim-1:
                    if isinstance(item[dim-1], slice):
                        if item[dim-1].start:
                            cropuli=item[dim-1].start
                        if item[dim-1].step:
                            dx*=item[dim-1].step
                        # croplri=item[dim-1].stop
                    else:
                        cropuli=item[dim-1]
                        # croplri=item[dim-1]+1
                else:
                    #test
                    print('we should not end up here for cols?')
                #rows
                if len(item)>dim-2:
                    if isinstance(item[dim-2], slice):
                        if item[dim-2].start:
                            cropulj=item[dim-2].start
                        if item[dim-2].step:
                            dy*=item[dim-2].step
                        # croplrj=item[dim-2].stop
                    else:
                        cropulj=item[dim-2]
                        # croplrj=item[dim-2]+1
                else:
                    #test
                    print('we should not end up here for rows?')

            upperLeft=self.geometry.image2geo(cropuli,cropulj);
            result=Jim(ncol=ncol,nrow=nrow,nband=nband,nplane=nplane,otype=self.properties.getDataType())
            result.properties.setProjection(self.properties.getProjection())
            gt=self.properties.getGeoTransform()

            cropulx=upperLeft[0]-self.properties.getDeltaX()/2;
            cropuly=upperLeft[1]+self.properties.getDeltaY()/2;

            gt[0]=cropulx;
            gt[1]=dx;
            gt[2]=0;
            gt[3]=cropuly;
            gt[4]=0;
            gt[5]=-dy;
            result.properties.setGeoTransform(gt)
            result.np()[:]=npresult
            return result

    def __setitem__(self, item, value):
        if isinstance(item, JimVect):
            if self.properties.nrOfPlane() > 1:
                raise ValueError('Error: __setitem__ with JimVect not implemented for 3d '
                                 'Jim objects')
            # TODO: decide on default behaviour of ALL_TOUCHED=TRUE
            # TODO: next lines should work, but problem with GML files when SRS is not defined as in S2 cloud masks
            # template=Jim(self)
            # print(item)
            # print("type of item: {}".format(type(item)))
            # template.geometry.rasterize(item,1.0)
            # self[template>0]=value
            if type(value) in (float, int):
                templateJim = Jim(self, copyData=False)
                templateJim = Jim(templateJim._jipjim.setMask(item._jipjimvect, {'eo': ['ALL_TOUCHED=TRUE'], 'nodata': 1}))
                self[templateJim>0]=value
            elif isinstance(value, Jim):
                templateJim = Jim(self, copyData=False)
                templateJim = Jim(templateJim._jipjim.setMask(item._jipjimvect, {'eo': ['ALL_TOUCHED=TRUE'], 'nodata': 1}))
                self[templateJim > 0] = value
        elif isinstance(item, Jim):  # or isinstance(value, Jim):
            if value is None:
                self._set(Jim(self, copyData=False))
            else:
                if isinstance(value, Jim):
                    self._jipjim.d_setMask(item._jipjim, value._jipjim)
                else:
                    self._jipjim.d_setMask(item._jipjim, value)
        elif isinstance(item, tuple):
            if isinstance(value, Jim):
                self.np()[item]=value.np()
            else:
                self.np()[item]=value
        else:
            raise ValueError('Error: __setitem__ only implemented for Vector, Jim or tuples')

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
        jim = Jim(self)
        jim.np()[:]=abs(jim.np())
        return jim

    def __neg__(self):
        """Calculate the negation of Jim raster dataset

        :return: Negation of Jim raster dataset (-dataset)
        """
        jim = Jim(self)
        jim.np()[:]=-(jim.np())
        return jim

    def __invert__(self):
        """Calculate the complement of Jim raster dataset

        :return: The complement of Jim raster dataset (~dataset)
        """
        jim = Jim(self)
        jim.np()[:]=~(jim.np())
        return jim

    # *** binary operators *** #
    def __eq__(self, right):
        """Pixel wise check for equality

        :return: Jim object with pixels 1 if equal values, 0 otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()==right.np())
        else:
            jim.np()[:]=(self.np()==right)
        return jim

    def __ne__(self, right):
        """Pixel wise check for non-equality

        :return: False if equal values, True otherwise
        """
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()!=right.np())
        else:
            jim.np()[:]=(self.np()!=right)
        return jim

    def __lt__(self, right):
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()<right.np())
        else:
            jim.np()[:]=(self.np()<right)
        return jim

    def __le__(self, right):
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()<=right.np())
        else:
            jim.np()[:]=(self.np()<=right)
        return jim

    def __gt__(self, right):
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()>right.np())
        else:
            jim.np()[:]=(self.np()>right)
        return jim

    def __ge__(self, right):
        jim = Jim(_jl.Jim(self.properties.nrOfCol(),self.properties.nrOfRow(),self.properties.nrOfBand(),self.properties.nrOfPlane(),_jl.GDT_Byte))
        jim.properties.setGeoTransform(self.properties.getGeoTransform())
        jim.properties.setProjection(self.properties.getProjection())
        if isinstance(right, Jim):
            jim.np()[:]=(self.np()>=right.np())
        else:
            jim.np()[:]=(self.np()>=right)
        return jim

    def __add__(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]+=right.np()
        else:
            result.np()[:]+=right
        return result

    def __radd__(self, left):
        result=Jim(self)
        if isinstance(left, Jim):
            result.np()[:]+=left.np()
        else:
            result.np()[:]+=left
        return result

    def __iadd__(self, right):
        if isinstance(right, Jim):
            self.np()[:]+=right.np()
        else:
            self.np()[:]+=right
        return self

    def __sub__(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]-=right.np()
        else:
            result.np()[:]-=right
        return result

    def __rsub__(self, left):
        result=Jim(self)
        if isinstance(left, Jim):
            result.np()[:]-=left.np()
        else:
            result.np()[:]-=left
        return result

    def __isub__(self, right):
        if isinstance(right, Jim):
            self.np()[:]-=right.np()
        else:
            self.np()[:]-=right
        return self

    def __mul__(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]*=right.np()
        else:
            result.np()[:]*=right
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
        result=Jim(self)
        if isinstance(left, Jim):
            result.np()[:]*=left.np()
        else:
            result.np()[:]*=left
        return result

    def __imul__(self, right):
        if isinstance(right, Jim):
            self.np()[:]*=right.np()
        else:
            self.np()[:]*=right
        return self

    def _trueDiv(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]/=right.np()
        else:
            result.np()[:]/=right
        return result

    def _itrueDiv(self, right):
        if isinstance(right, Jim):
            self.np()[:]/=right.np()
        else:
            self.np()[:]/=right
        return self

    def __div__(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]/=right.np()
        else:
            result.np()[:]/=right
        return result

    def __idiv__(self, right):
        if isinstance(right, Jim):
            self.np()[:]/=right.np()
        else:
            self.np()[:]/=right
        return self

    def __mod__(self, right):
        result=Jim(self)
        if isinstance(right, Jim):
            result.np()[:]%=right.np()
        else:
            result.np()[:]%=right
        return result

    def __imod__(self, right):
        if isinstance(right, Jim):
            self.np()[:]%=right.np()
        else:
            self.np()[:]%=right
        return self

    def __lshift__(self, right):
        if isinstance(right, int):
            jim = Jim(self)
            jim.np()[:]<<=right
            return jim
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))

    def __ilshift__(self, right):
        if isinstance(right, int):
            self.np()[:]<<=right
        else:
            raise TypeError('unsupported operand type for << : {}'.format(
                type(right)))
        return self

    def __rshift__(self, right):
        if isinstance(right, int):
            jim = Jim(self)
            jim.np()[:]>>=right
            return jim
        else:
            raise TypeError('unsupported operand type for >> : {}'.format(
                type(right)))

    def __irshift__(self, right):
        if isinstance(right, int):
            self.np()[:]>>=right
        else:
            raise TypeError('unsupported operand type for >> : {}'.format(
                type(right)))
        return self

    def __or__(self, right):
        if isinstance(right, Jim):
            jim = Jim(self)
            jim.np()[:]|=right.np()
            return jim
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(right)))

    def __ior__(self, right):
        if isinstance(right, Jim):
            self.np()[:]|=right.np()
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(right)))
        return self

    def __ror__(self, left):
        if isinstance(left, Jim):
            jim = Jim(self)
            jim.np()[:]|=left.np()
            return jim
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(left)))

    def __xor__(self, right):
        if isinstance(right, Jim):
            jim = Jim(self)
            jim.np()[:]^=right.np()
            return jim
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(right)))

    def __ixor__(self, right):
        if isinstance(right, Jim):
            self.np()[:]^=right.np()
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(right)))
        return self

    def __rxor__(self, left):
        if isinstance(left, Jim):
            jim = Jim(self)
            jim.np()[:]^=left.np()
            return jim
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(
                type(left)))

    def __and__(self, right):
        if isinstance(right, Jim):
            jim = Jim(self)
            jim.np()[:]&=right.np()
            return jim
        else:
            raise TypeError('unsupported operand type for & : {}'.format(
                type(right)))

    def __iand__(self, right):
        if isinstance(right, Jim):
            self.np()[:]&=right.np()
        else:
            raise TypeError('unsupported operand type for | : {}'.format(
                type(right)))
        return self

    def __rand__(self, left):
        if isinstance(left, Jim):
            jim = Jim(self)
            jim.np()[:]&=left.np()
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
        self._geometry._set_caller(self)
        _gc.collect()
        return self._geometry

    @property
    def io(self):
        self._io._set_caller(self)
        _gc.collect()
        return self._io

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

        for module in [geometry._GeometryList, io._IOList, pixops._PixOpsList,
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

    # def __init__(self, vector, *args):
    def __init__(self, vector, kwargs):
        """Initialize the JimVect object and modules for methods.

        :param vector: path to a vector or another JimVect object as a basis for
            the JimVect object
        """
        # super(_ParentVect, self).__init__(vector, *args)
        if kwargs:
            if vector:
                if isinstance(vector, JimVect):
                    super(_ParentVect, self).__init__(vector._jipjimvect,kwargs)
                else:
                    kwargs.update({'filename': vector})
                    super(_ParentVect, self).__init__(kwargs)
            else:
                super(_ParentVect, self).__init__(kwargs)
        else:
            if isinstance(vector, JimVect):
                super(_ParentVect, self).__init__(vector._jipjimvect)
            else:
                if vector:
                    super(_ParentVect, self).__init__(vector)
                else:
                    super(_ParentVect, self).__init__()


class JimVect():
    """Definition of JimVect object."""

    # def __init__(self, vector, *args):
    def __init__(self, vector=None, **kwargs):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        # self._jipjimvect = _ParentVect(vector, *args)
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
        self._classify._set_caller(self)
        _gc.collect()
        return self._classify

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

    @property
    def geometry(self):
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
        """Apply changes done in modified_object to the parent VectorOgr instance.

        :param modified_object: modified VectorOgr instance
        """
        self._jipjimvect.__dict__.update(modified_object.__dict__)
