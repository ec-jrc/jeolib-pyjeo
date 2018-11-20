import jiplib as _jl
import pyjeo as _pj


def stretch(jim_object, **kwargs):
    return _pj.Jim(jim_object.stretch(kwargs))


class _Stats():
    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def getStats(self, **kwargs):
        return self._jim_object.getStats(kwargs)

    def getStatProfile(self, **kwargs):
        return _pj.Jim(self._jim_object.statProfile(kwargs))

    def getHisto1d(self):
        """Compute the frequency distribution of the grey levels of im.

        :return: a Jim object
        """
        return _pj.Jim(self._jim_object.histo1d())

    def getHisto2d(self, im2):
        """Compute the frequency distribution of the grey levels pairs.

        :param im2: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object.histo2d(im2))

    def getHisto3d(self, im2, im3):
        """Compute the frequency distribution of the grey levels pairs.

        :param im2: a Jim object
        :param im3: a Jim object
        :return: a Jim object
        """
        # TODO: Doesn't work
        return _pj.Jim(self._jim_object.histo3d(im2, im3))

    def stretch(self, **kwargs):
        self._jim_object._set(self._jim_object.stretch(kwargs))
