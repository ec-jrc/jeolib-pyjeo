from __future__ import division
from jiplib import Jim as _jipJim
import numpy as np
import jiplib as _jl
from modules import pjio as io, properties, pixops, ngbops, geometry, \
    ccops, clssfy, demops, all


class Jim(_jipJim):
    def __init__(self, image, *args):
        """Initialize the Jim object and modules for methods.

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        super(Jim, self).__init__(image, *args)

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
                        retJim.geometry.crop(ulx=minCol, lrx=maxCol, uly=minRow, lry=maxRow, ulz=item[2].start, lrz=item[2].stop, dx=stridex, dy=stridey, nogeo=True)
                        # retJim.geometry.crop(ulx=ulx, lrx=lrx, uly=uly, lry=lry, ulz=item[2].start, lrz=item[2].stop,band=band)
                        return retJim
                    else:
                        raise TypeError('Error: use 4 dimensions when slicing multiband 3-dim Jim object (x:y:z:band)')
                else:
                    if len(item) == 3:#do slice x,y,z
                        retJim=geometry.crop(self,ulx=minCol, uly=minRow, ulz=item[2].start, lrx=maxCol, lry=maxRow, lrz=item[2].stop, dx=stridex, dy=stridey, nogeo=True)
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
                        retJim.geometry.crop(ulx=minCol, uly=minRow, ulz=None, lrx=maxCol, lry=maxRow, lrz=None, dx=stridex, dy=stridey, nogeo=True)
                        # retJim.geometry.crop(ulx=ulx, lrx=lrx, uly=uly, lry=lry,band=band)
                        return retJim
                    else:
                        raise TypeError('Error: use 3 dimensions when slicing multiband 2-dim Jim object (x:y:band)')
                else:
                    if len(item) == 2:#do slice x,y
                        retJim=geometry.crop(self,ulx=minCol, uly=minRow, ulz=None, lrx=maxCol, lry=maxRow, lrz=None, dx=stridex, dy=stridey, nogeo=True)
                        # retJim=geometry.crop(self,ulx=ulx, lrx=lrx, uly=uly, lry=lry,band=0)
                        return retJim
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing 2-dim Jim object (x:y)')
        elif isinstance(item, Jim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            itemnp=_jl.jim2np(item)
            # itemnp=itemnp==0
            selfnp[itemnp==0]=0
            retJim=Jim(_jl.np2jim(selfnp))
            retJim.properties.setProjection(projection)
            retJim.properties.setGeoTransform(gt)
            return retJim
        elif isinstance(item, _jl.VectorOgr):
            if self.nrOfPlane()>1:
                raise ValueError('Error: __getitem__ not implemented for 3d Jim objects')
            nodata=self.properties.getNoDataVals()
            if nodata:
                nodata=nodata[0]
            else:
                nodata=0
            return geometry.cropOgr(self,item,crop_to_cutline=True,nodata=nodata,align=True)

    def __setitem__(self, item, value):
        #todo: checkif rasterize vector first and then set raster mask would be better
        if isinstance(item, Jim):# or isinstance(value, Jim):
            if value is None:
                #todo set empty Jim?
                raise AttributeError("can't set item of Jim")
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            itemnp=_jl.jim2np(item)
            itemnp=itemnp>0
            if isinstance(value, Jim):
                valuenp=_jl.jim2np(value)
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
                        if isinstance(value,float):
                            if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
                                self.pixops.convert(otype='GDT_Float32')
                        if type(value) in (float,int):
                            bands=[0]
                            if isinstance(item[2],slice):
                                if item[2].step:
                                    strideb=item[2].step
                                bands=range(item[2].start,item[2].stop,strideb)
                            else:
                                bands=[item[2]]
                            print("setData to {}".format(value))
                            self.pixops.setData(value,ulx=minCol,uly=minRow,lrx=maxCol,lry=maxRow,bands=bands,dx=stridex,dy=stridey,nogeo=True)
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
                        if isinstance(value,float):
                            if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
                                self.pixops.convert(otype='GDT_Float32')
                        if type(value) in (float,int):
                            bands=[0]
                            print("setData to {}".format(value))
                            self.pixops.setData(value,ulx=minCol,uly=minRow,lrx=maxCol,lry=maxRow,bands=bands,dx=stridex,dy=stridey,nogeo=True)
                        else:
                            raise TypeError('Error: __setitem__ not implemented for value type {}'.format(type(value)))
                    else:
                        raise TypeError('Error: use 2 dimensions when slicing 2-dim Jim object (x:y)')
        elif isinstance(item, _jl.VectorOgr):
            if self.nrOfPlane()>1:
                raise ValueError('Error: __setitem__ not implemented for 3d Jim objects')
            #todo: decide on default behaviour of ALL_TOUCHED=TRUE
            if type(value) in (float,int):
                nodataValues=self.properties.getNoDataVals()
                self._set(self.setMask(item,{'eo':['ALL_TOUCHED=TRUE'],'nodata':value}))
                self.properties.clearNoData()
                self.properties.setNoDataVals(nodataValues)
            elif isinstance(value, Jim):
                templateJim=Jim(self,False)
                templateJim.setData(0)
                templateJim=Jim(templateJim.setMask(item,{'eo':['ALL_TOUCHED=TRUE'],'nodata':1}))
                self[templateJim>0]=value

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
        return Jim(self.pointOpAbs())

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
            print("Warning: converting data type to Int32, potential overflows may occur!")
            self.pixops.convert(otype='Int32')
        return -1*self

    ### binary operators ###
    def __eq__(self, aJim):
        """Change behaviour of == to check values, not memory alloc pointer.

        :return: Jim object with pixels 1 if equal values, 0 otherwise
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

    def __lt__(self, right):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            if isinstance(right, Jim):
                anp=_jl.jim2np(right)
            else:
                anp=right
            selfnp=np.uint8(1)*(selfnp<anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __le__(self, right):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            if isinstance(right, Jim):
                anp=_jl.jim2np(right)
            else:
                anp=right
            selfnp=np.uint8(1)*(selfnp<=anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __gt__(self, right):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            if isinstance(right, Jim):
                anp=_jl.jim2np(right)
            else:
                anp=right
            selfnp=np.uint8(1)*(selfnp>anp)
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
    def __ge__(self, right):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            if isinstance(right, Jim):
                anp=_jl.jim2np(right)
            else:
                anp=right
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
            raise TypeError('unsupported operand type for + : {}'.format(type(left)))
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
            raise TypeError('unsupported operand type for - : {}'.format(type(left)))
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
        elif isinstance(right, int):
            return Jim(self.pointOpArithCst(right,_jl.MULT_op))
        elif isinstance(right, float):
            if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
                self_float=pixops.convert(self,otype='GDT_Float32')
                self_float.d_pointOpArithCst(right,_jl.MULT_op)
                return Jim(self_float)
            else:
                return Jim(self.pointOpArithCst(right,_jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(right)))
    def __rmul__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpArith(left,_jl.MULT_op))
        elif isinstance(left, int):
            return Jim(self.pointOpArithCst(left,_jl.MULT_op))
        elif isinstance(left, float):
            if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
                self_float=pixops.convert(self,otype='GDT_Float32')
                self_float.d_pointOpArithCst(left,_jl.MULT_op)
                return Jim(self_float)
            else:
                return Jim(self.pointOpArithCst(left,_jl.MULT_op))
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(left)))
    def __imul__(self, right):
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.MULT_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.MULT_op)
        else:
            raise TypeError('unsupported operand type for * : {}'.format(type(right)))
        return self

    def trueDiv(self, right):
        #test
        print("true division")
        right_float=Jim(None)
        self_float=Jim(None)
        if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
            self_float=pixops.convert(self,otype='GDT_Float32')
        if isinstance(right, Jim):
            if right.properties.getDataType() != _jl.GDT_Float32 and right.properties.getDataType() != _jl.GDT_Float64:
                right_float=pixops.convert(right,otype='GDT_Float32')
                if self_float:
                    return Jim(self_float.pointOpArith(right_float,_jl.DIV_op))
                else:
                    return Jim(self.pointOpArith(right_float,_jl.DIV_op))
            else:
                if self_float:
                    return Jim(self_float.pointOpArith(right,_jl.DIV_op))
                else:
                    return Jim(self.pointOpArith(right,_jl.DIV_op))
        elif type(right) in (int,float):
            if self_float:
                return Jim(self_float.pointOpArithCst(right,_jl.DIV_op))
            else:
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
    def itrueDiv(self, right):
        #test
        print("true division")
        if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
            self.pixops.convert(otype='GDT_Float32')
        if isinstance(right, Jim):
            self.d_pointOpArith(right,_jl.DIV_op)
        elif type(right) in (int,float):
            self.d_pointOpArithCst(right,_jl.DIV_op)
        else:
            raise TypeError('unsupported operand type for / : {}'.format(type(right)))
        return self

    def __div__(self, right):
        #test
        print("division")
        truediv=False
        if self.properties.getDataType() == _jl.GDT_Float32 or self.properties.getDataType() == _jl.GDT_Float64:
            truediv=True
        else:
            if isinstance(right, Jim):
                if right.properties.getDataType() == _jl.GDT_Float32 or right.properties.getDataType() == _jl.GDT_Float64:
                    truediv=True
            elif isinstance(right,float):
                    truediv=True
        if truediv:
            return(self.trueDiv(right))
        else:
            if isinstance(right, Jim):
                return Jim(self.pointOpArith(right,_jl.DIV_op))
            elif type(right) in (int):
                return Jim(self.pointOpArithCst(right,_jl.DIV_op))
            else:
                raise TypeError('unsupported operand type for / : {}'.format(type(right)))

    def __idiv__(self, right):
        if self.properties.getDataType() != _jl.GDT_Float32 and self.properties.getDataType() != _jl.GDT_Float64:
            self.pixops.convert(otype='GDT_Float32')
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
            return Jim(self.pointOpBitwise(right,_jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(type(right)))
    def __ror__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitwise(left,_jl.OR_op))
        else:
            raise TypeError('unsupported operand type for | : {}'.format(type(left)))
    # def __ior__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitwise(right,_jl.OR_op)
    #     else:
    #         raise TypeError('unsupported operand type for | : {}'.format(type(right)))
    #     return self

    def __xor__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpBitwise(right,_jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(type(right)))
    def __rxor__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitwise(left,_jl.XOR_op))
        else:
            raise TypeError('unsupported operand type for ^ : {}'.format(type(left)))
    # def __ixor__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitwise(right,_jl.XOR_op)
    #     else:
    #         raise TypeError('unsupported operand type for ^ : {}'.format(type(right)))
    #     return self

    def __and__(self, right):
        if isinstance(right, Jim):
            return Jim(self.pointOpBitwise(right,_jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(type(right)))
    def __rand__(self, left):
        if isinstance(left, Jim):
            return Jim(self.pointOpBitwise(left,_jl.AND_op))
        else:
            raise TypeError('unsupported operand type for & : {}'.format(type(left)))
    # def __iand__(self, right):
    #     if isinstance(right, Jim):
    #         self.d_pointOpBitwise(right,_jl.AND_op)
    #     else:
    #         raise TypeError('unsupported operand type for & : {}'.format(type(right)))
    #     return self
