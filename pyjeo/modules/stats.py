"""Module for statistical functions and interpolations."""

import pyjeo as _pj


def getStats(jim_object, **kwargs):
    return jim_object._jipjim.getStats(kwargs)


def getStatProfile(jim_object, function, **kwargs):
    kwargs.update({'function': function})
    return _pj.Jim(jim_object._jipjim.statProfile(kwargs))


def getHisto1d(jim_object):
    """Compute the frequency distribution of the grey levels of im.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo1d())


def getHisto2d(jim_object, jim_object2):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo2d(jim_object2._jipjim))


def getHisto3d(jim_object, jim_object2, jim_object3):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :param jim_object3: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.histo3d(jim_object2._jipjim,
                                              jim_object3._jipjim))


def stretch(jim_object, **kwargs):
    return _pj.Jim(jim_object._jipjim.stretch(kwargs))


class _Stats():
    """Define all statistical methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def getHisto1d(self):
        """Compute the frequency distribution of the grey levels of im.

        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo1d())

    def getHisto2d(self, jim_object2):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo2d(jim_object2._jipjim))

    def getHisto3d(self, jim_object2, jim_object3):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :param jim_object3: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object._jipjim.histo3d(jim_object2._jipjim,
                                                        jim_object3._jipjim))

    def getHistoCumulative(self):
        """Compute the cumulative frequency distribution of the grey levels.

        :return: a Jim object
        """
        if self._jim_object.properties.getDataType() not in [4, 5]:
            raise TypeError('Object must be of type UInt32 or Int32')
        return _pj.Jim(self._jim_object._jipjim.histo1dCumulative())

    def getStats(self, **kwargs):
        return self._jim_object._jipjim.getStats(kwargs)

    def statProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        self._jim_object._set(self._jim_object._jipjim.statProfile(kwargs))

    def stretch(self, **kwargs):
        self._jim_object._set(self._jim_object._jipjim.stretch(kwargs))


class _StatsList():
    """Define all statistical methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def getStatProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        return _pj.Jim(self._jim_list._jipjimlist.statProfile(kwargs))

    def getStats(self, **kwargs):
        return self._jim_list._jipjimlist.getStats(kwargs)


class _StatsVect():
    """Define all statistical methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
