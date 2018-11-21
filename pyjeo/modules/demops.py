try:
    import pyjeo as _pj
except:
    from jeodpp import pyjeo as _pj


def demFlowDirectionD8(jim_object):
    """Compute the D8 steepest slope direction of each pixel.

    The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
    E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
    to 0.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object.demFlowDirectionD8())


class _DEMOps():
    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def demFlowDirectionD8(self):
        """Compute the D8 steepest slope direction of each pixel.

        The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
        E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
        to 0.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(self._jim_object.demFlowDirectionD8())
