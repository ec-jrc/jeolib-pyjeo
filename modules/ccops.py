import jiplib as _jl
import pyjeo as _pj


def labelImagePixels(jim_object):
    """Label each non-zero pixel of im with a unique label.

    Labels unless label overflows occurs.

    :param jim_object: a Jim object
    :return: labeled Jim object
    """
    return _pj.Jim(_jl.Jim.labelPix(jim_object))


def distance2dEuclideanSquared(jim_object, iband=0):
    """Compute the squared Euclidean distance transform of im.

    im must be a 2-D binary image. Original algorihtm proposed by Saito
    and Toriwaki (1994) and then optimised independently by (Hirata,
    1996) and (Meijster et al., 2000). See also *edt for the actual
    Euclidean distance transform. Note that a temporary buffer of type
    UINT16 is used for sums along/lines and columns so that uncontrolled
    results will occur if an object shows more than 16 2 /2 foreground
    pixels along a given line or column.

    :param jim_object: a Jim object
    :param iband: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object.distance2dEuclideanSquared(iband))


class _CCOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def labelImagePixels(self):
        """Label each non-zero pixel of im with a unique label.

        Labels unless label overflows occurs.

        Modifies the instance on which the method was called.
        """
        self._jim_object.d_labelPix()

    def distance2dEuclideanSquared(self, iband=0):
        """Compute the squared Euclidean distance transform of im.

        im must be a 2-D binary image. Original algorihtm proposed by Saito
        and Toriwaki (1994) and then optimised independently by (Hirata,
        1996) and (Meijster et al., 2000). See also *edt for the actual
        Euclidean distance transform. Note that a temporary buffer of type
        UINT16 is used for sums along/lines and columns so that uncontrolled
        results will occur if an object shows more than 16 2 /2 foreground
        pixels along a given line or column.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(self._jim_object.distance2dEuclideanSquared(
            iband))
