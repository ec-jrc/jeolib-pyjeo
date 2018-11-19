import jiplib as _jl


class _NgbOps():
    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def morphoErode(self, sec_jim_object, ox, oy, oz, trFlag=0):
        """Output the dilation of im using the SE defined by imse.

        Its origin is set at coordinates (x,y,z). The reflection of the SE
        is considered if trflag equals 1 (no reflection by default). Points of
        the SE are points with a non zero value in imse.

        :param sec_jim_object: an image node for SE (UCHAR type)
        :param ox: x coordinate
        :param oy: y coordinate
        :param oz: z coordinate
        :param trFlag: optional parameter (0 or 1)
        """
        self._jim_object.morphoErode(sec_jim_object, ox, oy, oz, trFlag)
