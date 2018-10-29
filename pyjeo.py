from jiplib import Jim as _jipJim
import numpy as np
import jiplib as _jl
from modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, clssfy, demops, all


class Jim(_jipJim):
    def __init__(self, image):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        super(Jim, self).__init__(image)

        self.properties = properties._Properties(self)
        self.io = io._IO(self)
        self.pixops = pixops._PixOps(self)
        self.ngbops = ngbops._NgbOps(self)
        self.geometry = geometry._Geometry(self)
        self.ccops = ccops._CCOps(self)
        self.clssfy = clssfy._Classify(self)
        self.demops = demops._DEMOps(self)
        self.all = all._All(self)

    def __dir__(self):
        """Change behaviour of the method whisperer to ignore jiplib methods.

        :return: a whispered module or method
        """
        pyjeo_Jim_methods = list(set(dir(Jim)) - set(dir(_jipJim)))
        return [i for i in self.__dict__.keys() if i != 'this'] + \
               pyjeo_Jim_methods

    def getMethods(self):
        """Print an overview of available methods in format module.method."""
        def treeStructure(module):
            module_methods = dir(module)
            for default_method in ['__init__', '__module__', '__doc__']:
                module_methods.remove(default_method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:] + '.' + \
                                    module_methods[i]

            return ['\nmodule {}:'.format(module.__name__.lower()[1:])] + \
                   module_methods

        methods = list()
        for module in [properties._Properties, io._IO, pixops._PixOps,
                       ngbops._NgbOps, geometry._Geometry, ccops._CCOps,
                       clssfy._Classify, demops._DEMOps]:
            methods.extend(treeStructure(module))

        print('\n'.join(methods))

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self.__dict__.update(modified_object.__dict__)

    ### unary operators ###

    def __set__(self, value):
        if value is None:
            self=Jim()
        elif type(value) in (int):
            if self.getDataType() == _jl.GDT_Byte:
                self.d_pointOpBlank(value)
        elif isinstance(value, Jim):
            self=Jim(value)


    def __getitem__(self, item):
        if isinstance(item, tuple):
            if self.nrOfPlane()>1:
                if self.nrOfBand()>1:
                    if len(item) == 4:#do slice x,y,z,band
                        if isinstance(item[3], slice):
                            if item[3].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            retJim=self.cropBand({'startband':item[3].start,'endband':item[3].stop-1})
                        elif type(item[3]) in (int,tuple):
                            retJim=self.cropBand({'band':item[3]})
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        if isinstance(item[0], slice) and isinstance(item[1], slice) and isinstance(item[2], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None) or item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,retJim.nrOfBand()):
                                retJim.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, item[2].start, item[2].stop],band)
                            return retJim
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 4 dimensions when slicing multiband 3-dim Jim object (x:y:z:band)')
                else:
                    if len(item) == 3:#do slice x,y,z
                        if isinstance(item[0], slice) and isinstance(item[1], slice) and isinstance(item[2], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None) or item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,retJim.nrOfBand()):
                                retJim.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, item[2].start, item[2].stop],band)
                            return retJim
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing 3-dim Jim object (x:y:z)')
            else:
                if self.nrOfBand()>1:
                    if len(item) == 3:#do slice x,y,band
                        if isinstance(item[2], slice):
                            if item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            retJim=self.cropBand({'startband':item[2].start,'endband':item[2].stop-1})
                        elif type(item[2]) in (int,tuple):
                            retJim=self.cropBand({'band':item[2]})
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        if isinstance(item[0], slice) and isinstance(item[1], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,retJim.nrOfBand()):
                                retJim.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, 0, 0],band)
                            return retJim
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing multiband 2-dim Jim object (x:y:band)')
                else:
                    if len(item) == 2:#do slice x,y
                        if isinstance(item[0], slice) and isinstance(item[1], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            retJim.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, 0, 0])
                            return retJim
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing 2-dim Jim object (x:y)')

    def __setitem__(self, item, value):
        if isinstance(item, tuple):
            if self.nrOfPlane()>1:
                if self.nrOfBand()>1:
                    if len(item) == 4:#do slice x,y,z,band
                        if isinstance(item[3], slice):
                            if item[3].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            self.d_cropBand({'startband':item[3].start,'endband':item[3].stop-1})
                        elif type(item[3]) in (int,tuple):
                            self.d_cropBand({'band':item[3]})
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        if isinstance(item[0], slice) and isinstance(item[1], slice) and isinstance(item[2], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None) or item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,self.nrOfBand()):
                                self.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, item[2].start, item[2].stop],band)
                            if isinstance(value,Jim):
                                if self.nrOfBand() != value.nrOfBand():
                                    raise ValueError('Error: number of bands do not match')
                                if self.nrOfCol() != value.nrOfCol():
                                    raise ValueError('Error: number of cols do not match')
                                if self.nrOfRow() != value.nrOfRow():
                                    raise ValueError('Error: number of rows do not match')
                                for band in range(0,self.nrOfBand()):
                                    value.copyData(self.getDataPointer(band))
                            elif type(value) in (int,float):
                                for band in range(0,self.nrOfBand()):
                                    self.setData(value,band)
                            return self
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 4 dimensions when slicing multiband 3-dim Jim object (x:y:z:band)')
                else:
                    if len(item) == 3:#do slice x,y,z
                        if isinstance(item[0], slice) and isinstance(item[1], slice) and isinstance(item[2], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None) or item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,self.nrOfBand()):
                                self.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, item[2].start, item[2].stop],band)
                            if isinstance(value,Jim):
                                if self.nrOfBand() != value.nrOfBand():
                                    raise ValueError('Error: number of bands do not match')
                                if self.nrOfCol() != value.nrOfCol():
                                    raise ValueError('Error: number of cols do not match')
                                if self.nrOfRow() != value.nrOfRow():
                                    raise ValueError('Error: number of rows do not match')
                                for band in range(0,self.nrOfBand()):
                                    value.copyData(self.getDataPointer(band))
                            elif type(value) in (int,float):
                                for band in range(0,self.nrOfBand()):
                                    self.setData(value,band)
                            return self
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing 3-dim Jim object (x:y:z)')
            else:
                if self.nrOfBand()>1:
                    if len(item) == 3:#do slice x,y,band
                        if isinstance(item[2], slice):
                            if item[2].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            self.d_cropBand({'startband':item[2].start,'endband':item[2].stop-1})
                        elif type(item[2]) in (int,tuple):
                            self.cropBand({'band':item[2]})
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        if isinstance(item[0], slice) and isinstance(item[1], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            for band in range(0,self.nrOfBand()):
                                self.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, 0, 0],band)
                            if isinstance(value,Jim):
                                if self.nrOfBand() != value.nrOfBand():
                                    raise ValueError('Error: number of bands do not match')
                                if self.nrOfCol() != value.nrOfCol():
                                    raise ValueError('Error: number of cols do not match')
                                if self.nrOfRow() != value.nrOfRow():
                                    raise ValueError('Error: number of rows do not match')
                                for band in range(0,self.nrOfBand()):
                                    value.copyData(self.getDataPointer(band))
                            elif type(value) in (int,float):
                                for band in range(0,self.nrOfBand()):
                                    self.setData(value,band)
                            return self
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing multiband 2-dim Jim object (x:y:band)')
                else:
                    if len(item) == 2:#do slice x,y
                        if isinstance(item[0], slice) and isinstance(item[1], slice):
                            if item[0].step not in (1, None) or item[1].step not in (1, None):
                                raise ValueError('Error: only step=1 supported')
                            self.geometry.crop([item[0].start, item[0].stop, item[1].start, item[1].stop, 0, 0])
                            if isinstance(value,Jim):
                                if self.nrOfBand() != value.nrOfBand():
                                    raise ValueError('Error: number of bands do not match')
                                if self.nrOfCol() != value.nrOfCol():
                                    raise ValueError('Error: number of cols do not match')
                                if self.nrOfRow() != value.nrOfRow():
                                    raise ValueError('Error: number of rows do not match')
                                for band in range(0,self.nrOfBand()):
                                    value.copyData(self.getDataPointer(band))
                            elif type(value) in (int,float):
                                for band in range(0,self.nrOfBand()):
                                    self.setData(value,band)
                            return self
                        else:
                            raise TypeError('items must be slice')
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing 2-dim Jim object (x:y)')

    # def __setitem__(self, item, value):
    #     if value is None:
    #         raise AttributeError("can't set item of Jim")
    #     if isinstance(item, slice):
    #         #todo: implement slicing to replace cropband and MIA window imageFrameSet/Add/imageInsert/imageFrameSet/imageFrameAdd 
    #         raise typerror('slicing not supported')
    #         # if item.step not in (1, None):
    #         #     raise ValueError('only step=1 supported')
    #     if isinstance(value, Jim):
    #         projection=self.properties.getProjection()
    #         gt=self.properties.getGeoTransform()
    #         selfnp=_jl.jim2np(self)
    #         valuenp=_jl.jim2np(value)
    #         itemnp=_jl.jim2np(item)
    #         itemnp=itemnp>0
    #         selfnp[itemnp]=valuenp[itemnp]
    #         jim=Jim(_jl.np2jim(selfnp))
    #         jim.properties.setProjection(projection)
    #         jim.properties.setGeoTransform(gt)
    #         return jim
    #     else:
    #         projection=self.properties.getProjection()
    #         gt=self.properties.getGeoTransform()
    #         selfnp=_jl.jim2np(self)
    #         itemnp=_jl.jim2np(item)
    #         itemnp=itemnp>0
    #         selfnp[itemnp]=value
    #         jim=Jim(_jl.np2jim(selfnp))
    #         jim.properties.setProjection(projection)
    #         jim.properties.setGeoTransform(gt)
    #         return jim

    def __nonzero__(self):
        """Check if Jim contains data

        :return: True if image contains data, False if image is empty
        """
        return self.isInit()

    def __bool__(self):
        """Check if Jim contains data

        :return: True if image contains data, False if image is empty
        """
        return self.isInit()

    def __abs__(self):
        """Calculate the absolute value of Jim raster dataset

        :return: Absolute value of Jim raster dataset
        """
        self.d_pointOpAbs()
        return self

    ### binary operators ###
    def __eq__(self, aJim):
        """Change behaviour of == to check values, not memory alloc pointer.

        :return: True if equal values, False otherwise
        """
        projection = self.properties.getProjection()
        gt = self.properties.getGeoTransform()
        selfnp = _jl.jim2np(self)
        anp = _jl.jim2np(aJim)
        selfnp = np.uint8(1) * (selfnp == anp)
        jim = Jim(_jl.np2jim(selfnp))
        jim.properties.setProjection(projection)
        jim.properties.setGeoTransform(gt)
        return jim

    def __neq__(self, other):
        """Change behaviour of != to check values, not memory alloc pointer.

        :return: False if equal values, True otherwise
        """
        if isinstance(other, Jim):
            return not self.isEqual(other)
        else:
            return False

    def __lt__(self, aJim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            anp=_jl.jim2np(aJim)
            selfnp=np.uint8(1)*(selfnp<anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __le__(self, aJim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            anp=_jl.jim2np(aJim)
            selfnp=np.uint8(1)*(selfnp<=anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __gt__(self, aJim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            anp=_jl.jim2np(aJim)
            selfnp=np.uint8(1)*(selfnp>anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __ge__(self, aJim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            anp=_jl.jim2np(aJim)
            selfnp=np.uint8(1)*(selfnp>=anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim

    def __add__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpArith(right,_jl.ADD_op))
        elif type(right) in (int,float):
            return Jim(self.pointOpArithCst(right,_jl.ADD_op))
        else:
            raise TypeError('unsupported operand type for + : {}'.format(type(right)))
    def __radd__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpArith(left,_jl.ADD_op))
        elif type(left) in (int,float):
            return Jim(self.pointOpArithCst(left,_jl.ADD_op))
        else:
            raise TypeError('unsupported operand type for + : {}'.format(type(right)))
    def __iadd__(self, right):
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.ADD_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.ADD_op)
        else:
            raise TypeError('unsupported operand type for + : {}'.format(type(right)))
        return self

    def __sub__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpArith(right,_jl.SUB_op))
        elif type(right) in (int,float):
            return Jim(self.pointOpArithCst(right,_jl.SUB_op))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
    def __rsub__(self, left):
        if isinstance(left, Jim):
            return -1*Jim(self.pointOpArith(left,_jl.SUB_op))
        elif type(left) in (int,float):
            return -1*Jim(self.pointOpArithCst(left,_jl.SUB_op))
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
    def __isub__(self, right):
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.SUB_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.SUB_op)
        else:
            raise TypeError('unsupported operand type for - : {}'.format(type(right)))
        return self

    def __mul__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpArith(right,_jl.MULT_op))
        elif type(right) in (int,float):
            return Jim(self.pointOpArithCst(right,_jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(right)))
    def __rmul__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpArith(left,_jl.MULT_op))
        elif type(left) in (int,float):
            return Jim(self.pointOpArithCst(left,_jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(right)))
    def __imul__(self, right):
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.MULT_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.MULT_op)
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(right)))
        return self

    def __truediv__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpArith(right,_jl.DIV_op))
        elif type(right) in (int,float):
            return Jim(self.pointOpArithCst(right,_jl.DIV_op))
        else:
            raise TypeError('unsupported operand type for / : {}'.format(type(right)))
    # def __rtruediv__(self, left):
    #     if isinstance(left, Jim):
    #         return Jim(self.pointOpArith(left,_jl.DIV_op))
    #     elif type(left) in (int,float):
    #         return Jim(self.pointOpArithCst(left,_jl.DIV_op))
    #     else:
    #         raise TypeError('unsupported operand type for / : {}'.format(type(right)))
    def __itruediv__(self, right):
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.DIV_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.DIV_op)
        else:
            raise TypeError('unsupported operand type for / : {}'.format(type(right)))
        return self

    def __mod__(self, right):
        if type(right) in (int):
            return Jim(self.pointOpModulo(right))
        else:
            raise TypeError('unsupported operand type for % : {}'.format(type(right)))
    def __imod__(self, right):
        if type(right) in (int):
            self.d_pointOpModulo(right)
        else:
            raise TypeError('unsupported operand type for % : {}'.format(type(right)))
        return self

    def __lshift__(self, right):
        if type(right) in (int):
            return Jim(self.pointOpBitShift(-right))
        else:
            raise TypeError('unsupported operand type for << : {}'.format(type(right)))
    def __ilshift__(self, right):
        if type(right) in (int):
            self.d_pointOpBitShift(-right)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(type(right)))
        return self

    def __rshift__(self, right):
        if type(right) in (int):
            return Jim(self.pointOpBitShift(right))
        else:
            raise TypeError('unsupported operand type for << : {}'.format(type(right)))
    def __rlshift__(self, right):
        if type(right) in (int):
            self.d_pointOpBitShift(right)
        else:
            raise TypeError('unsupported operand type for << : {}'.format(type(right)))
        return self

    def __or__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpBitWise(right,_jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(type(right)))
    def __ror__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitWise(left,_jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(type(right)))
    # def __ior__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitWise(right,_jl.OR_op)
    #     else:
    #         raise TypeError('unsupported operand type for | : {}'.format(type(right)))
    #     return self

    def __xor__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpBitWise(right,_jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(type(right)))
    def __rxor__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitWise(left,_jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(type(right)))
    # def __ixor__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitWise(right,_jl.XOR_op)
    #     else:
    #         raise TypeError('unsupported operand type for ^ : {}'.format(type(right)))
    #     return self

    def __and__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpBitWise(right,_jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(type(right)))
    def __rand__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitWise(left,_jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(type(right)))
    # def __iand__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitWise(right,_jl.AND_op)
    #     else:
    #         raise TypeError('unsupported operand type for & : {}'.format(type(right)))
    #     return self
