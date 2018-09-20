from jiplib import Jim as _jipJim
from modules import properties, io, pixops, ngbops, geometry, \
    ccops, clssfy, demops, all


class Jim(_jipJim):
    def __init__(self, image):
        """initialize the Jim object and modules for methods

        :param image: path to a raster or another Jim object as a basis for
        the Jim object
        """
        super(Jim, self).__init__(image)
        print('you created the Jim Object')

        self.dx = 5
        self.dy = 10
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
        """change behaviour of the method whisperer to ignore jiplib methods

        :return: a whispered package or method
        """
        pyjeo_Jim_methods = list(set(dir(Jim)) - set(dir(_jipJim)))
        return [i for i in self.__dict__.keys() if i != 'this'] + \
               pyjeo_Jim_methods

    def _set(self, modified_object):
        """apply changes done in the modified_object to the parent Jim instance

        :param modified_object: modified Jim instance
        """
        try:
            self.__dict__.update(modified_object.__dict__)
        except AttributeError:
            from inspect import stack
            raise ValueError(
                'A problem happened in function {}. Are you sure that you used'
                ' it the right way?'.format(stack()[1][3]))

    def get_methods(self):
        """print an overview of available methods in format module.method"""
        def tree_structure(module):
            module_methods = dir(module)
            for default_method in ['__init__', '__module__', '__doc__']:
                module_methods.remove(default_method)

            for i in range(len(module_methods)):
                module_methods[i] = module.__name__.lower()[1:] + '.' + \
                                    module_methods[i]

            return module_methods

        methods = list()
        for module in [properties._Properties, io._IO, pixops._PixOps,
                       ngbops._NgbOps, geometry._Geometry, ccops._CCOps,
                       clssfy._Classify, demops._DEMOps]:
            methods.extend(tree_structure(module))

        print(', '.join(methods))
