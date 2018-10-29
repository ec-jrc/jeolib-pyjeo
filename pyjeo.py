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

    def __setitem__(self, item, value):
        if isinstance(item, slice):
            raise typerror('slicing not supported')
            # if item.step not in (1, None):
            #     raise ValueError('only step=1 supported')
        if isinstance(value, Jim):
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            valuenp=_jl.jim2np(value)
            itemnp=_jl.jim2np(item)
            itemnp=itemnp>0
            selfnp[itemnp]=valuenp[itemnp]
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim
        else:
            projection=self.properties.getProjection()
            gt=self.properties.getGeoTransform()
            selfnp=_jl.jim2np(self)
            itemnp=_jl.jim2np(item)
            itemnp=itemnp>0
            selfnp[itemnp]=value
            jim=Jim(_jl.np2jim(selfnp))
            jim.properties.setProjection(projection)
            jim.properties.setGeoTransform(gt)
            return jim

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
