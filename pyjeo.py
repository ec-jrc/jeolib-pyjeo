from jiplib import Jim as _jipJim
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

    def __eq__(self, other):
        """Change behaviour of == to check values, not memory alloc pointer.

        :return: True if equal values, False otherwise
        """
        if isinstance(other, Jim):
            return self.isEqual(other)
        else:
            return False

    def _set(self, modified_object):
        """Apply changes done in modified_object to the parent Jim instance.

        :param modified_object: modified Jim instance
        """
        self.__dict__.update(modified_object.__dict__)

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
