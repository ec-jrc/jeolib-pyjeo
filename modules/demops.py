class _DEMOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def demFlowDirectionD8(self):
        """Compute the D8 steepest slope direction of each pixel.

        The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
        E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
        to 0.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(self._jim_object.demFlowDirectionD8())
