import jiplib as _jl


def morphoErode(im, imSE, ox, oy, oz, trFlag=0):
    """outputs the dilation of im using the SE defined by imse

    Its origin is set at coordinates (x,y,z). The reflection of the SE
    is considered if trflag equals 1 (no reflection by default). Points of
    the SE are points with a non zero value in imse.

    :param im: an image node
    :param imSE: an image node for SE (UCHAR type)
    :param ox: x coordinate
    :param oy: y coordinate
    :param oz: z coordinate
    :param trFlag: optional parameter (0 or 1)
    :return:
    """
    return _jl.erode(im, imSE, ox, oy, oz, trFlag)


class _NgbOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object
