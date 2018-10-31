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

    def getMethods(self, queried_module=None):
        """Print an overview of available methods in format module.method."""
        def treeStructure(module, queried_module):
            if queried_module and queried_module not in str(module):
                return ''

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
            methods.extend(treeStructure(module, queried_module))

        print('\n'.join(methods))

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self.__dict__.update(modified_object.__dict__)

    ### unary operators ###

    def __getitem__(self, item):
        stridex=1
        stridey=1
        strideb=1
        if isinstance(item, tuple):
            if isinstance(item[0],slice):
                minCol=item[0].start
                maxCol=item[0].stop-1
                if item[0].step:
                    stridex=item[0].step
            elif isinstance(item[0],int):
                minCol=item[0]
                maxCol=item[0]+1
            else:
                raise ValueError('column item must be slice or integer value')
            if isinstance(item[1],slice):
                minRow=item[1].start
                maxRow=item[1].stop-1
                if item[1].step:
                    stridey=item[1].step
            elif isinstance(item[1],int):
                minRow=item[1]
                maxRow=item[1]+1
            else:
                raise ValueError('row item must be slice or integer value')
            bands=[]
            if self.nrOfPlane()>1:
                if self.properties.nrOfBand()>1:
                    if len(item) == 4:#do slice x,y,z,band
                        if isinstance(item[3], slice):
                            if item[2].step:
                                strideb=item[2].step
                            bands=range(item[3].start,item[3].stop,strideb)
                        if isinstance(item[3], tuple):
                            bands=item[3]
                        elif isinstance(item[3], int):
                            bands.append(item[3])
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        retJim=geometry.cropBand(self,band=bands)
                        retJim.geometry.crop(ulx=minCol, lrx=maxCol, uly=minRow, lry=maxRow, ulz=item[2].start, lrz=item[2].stop, dx=stridex, dy=stridey, geo=False)
                        # retJim.geometry.crop(ulx=ulx, lrx=lrx, uly=uly, lry=lry, ulz=item[2].start, lrz=item[2].stop,band=band)
                        return retJim
                    else:
                        raise TypeError('Error: use 4 dimensions when slicing multiband 3-dim Jim object (x:y:z:band)')
                else:
                    if len(item) == 3:#do slice x,y,z
                        retJim=geometry.crop(self,ulx=minCol, uly=minRow, ulz=item[2].start, lrx=maxCol, lry=maxRow, lrz=item[2].stop, dx=stridex, dy=stridey, geo=False)
                        # retJim=geometry.crop(self,ulx=ulx, lrx=lrx, uly=uly, lry=lry, ulz=item[2].start, lrz=item[2].stop)
                        return retJim
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing 3-dim Jim object (x:y:z)')
            else:
                if self.properties.nrOfBand()>1:
                    if len(item) == 3:#do slice x,y,band
                        if isinstance(item[2], slice):
                            if item[2].step:
                                strideb=item[2].step
                            bands=range(item[2].start,item[2].stop,strideb)
                        if isinstance(item[2], tuple):
                            bands=item[2]
                        elif isinstance(item[2], int):
                            bands.append(item[2])
                        else:
                            raise ValueError('Error: band index must be slice, list or integer')
                        retJim=geometry.cropBand(self,band=bands)
                        retJim.geometry.crop(ulx=minCol, uly=minRow, ulz=None, lrx=maxCol, lry=maxRow, lrz=None, dx=stridex, dy=stridey, geo=False)
                        # retJim.geometry.crop(ulx=ulx, lrx=lrx, uly=uly, lry=lry,band=band)
                        return retJim
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing multiband 2-dim Jim object (x:y:band)')
                else:
                    if len(item) == 2:#do slice x,y
                        retJim=geometry.crop(self,ulx=minCol, uly=minRow, ulz=None, lrx=maxCol, lry=maxRow, lrz=None, dx=stridex, dy=stridey, geo=False)
                        # retJim=geometry.crop(self,ulx=ulx, lrx=lrx, uly=uly, lry=lry,band=0)
                        return retJim
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing 2-dim Jim object (x:y)')

    def __setitem__(self, item, value):
        if isinstance(item, Jim) or isinstance(value, Jim):
            if value is None:
                #todo set empty Jim?
                raise AttributeError("can't set item of Jim")
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            valuenp=_jl.jim2np(value)
            itemnp=_jl.jim2np(item)
            itemnp=itemnp>0
            if isinstance(value, Jim):
                selfnp[itemnp]=valuenp[itemnp]
            else:
                selfnp[itemnp]=value
            self._set(_jl.np2jim(selfnp))
            self.properties.setProjection(projection)
            self.properties.setGeoTransform(gt)
        elif isinstance(item, tuple):
            if self.nrOfPlane()>1:
                raise ValueError('Error: __setitem__ not implemented for 3d Jim objects')
            else:
                if self.properties.nrOfBand()>1:
                    if len(item) == 3:#do slice x,y,band
                        stridex=1
                        stridey=1
                        strideb=1
                        if isinstance(item[0],slice):
                            minCol=item[0].start
                            maxCol=item[0].stop-1
                            if item[0].step:
                                stridex=item[0].step
                        elif isinstance(item[0],int):
                            minCol=item[0]
                            maxCol=item[0]
                        else:
                            raise ValueError('column item must be slice or integer value')
                        if isinstance(item[1],slice):
                            minRow=item[1].start
                            maxRow=item[1].stop-1
                            if item[1].step:
                                stridey=item[1].step
                        elif isinstance(item[1],int):
                            minRow=item[1]
                            maxRow=item[1]
                        else:
                            raise ValueError('row item must be slice or integer value')
                        if type(value) in (int,float):
                            if isinstance(item[2],slice):
                                if item[2].step:
                                    strideb=item[2].step
                                bands=range(item[2].start,item[2].stop,strideb)
                            else:
                                bands=[item[2]]
                            self.pixops.setData(value,ulx=minCol,uly=minRow,lrx=maxCol,lry=maxRow,bands=bands,dx=stridex,dy=stridey,geo=False)
                        else:
                            raise TypeError('Error: __setitem__ not implemented for value type {}'.format(type(value)))
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing multiband 2-dim Jim object (x:y:band)')
                else:
                    stridex=1
                    stridey=1
                    if len(item) == 2:#do slice x,y
                        if isinstance(item[0],slice):
                            minCol=item[0].start
                            maxCol=item[0].stop-1
                            if item[0].step:
                                stridex=item[0].step
                        elif isinstance(item[0],int):
                            minCol=item[0]
                            maxCol=item[0]
                        else:
                            raise ValueError('column item must be slice or integer value')
                        if isinstance(item[1],slice):
                            minRow=item[1].start
                            maxRow=item[1].stop-1
                            if item[1].step:
                                stridey=item[1].step
                        elif isinstance(item[1],int):
                            minRow=item[1]
                            maxRow=item[1]
                        else:
                            raise ValueError('row item must be slice or integer value')
                        if type(value) in (int,float):
                            bands=[0]
                            self.pixops.setData(value,ulx=minCol,uly=minRow,lrx=maxCol,lry=maxRow,bands=bands,dx=stridex,dy=stridey,geo=False)
                        else:
                            raise TypeError('Error: __setitem__ not implemented for value type {}'.format(type(value)))
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
