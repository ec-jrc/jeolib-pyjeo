"""Module for statistical functions and interpolations."""

try:
    import pyjeo as _pj
except ImportError:
    try:
        from jeodpp import pyjeo as _pj
    except ImportError:
        import jeodpp.pyjeo as _pj


def getStats(jim_object, **kwargs):
    return jim_object.getStats(kwargs)


def getStatProfile(jim_object, **kwargs):
    return _pj.Jim(jim_object.statProfile(kwargs))


def getHisto1d(jim_object):
    """Compute the frequency distribution of the grey levels of im.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object.histo1d())


def getHisto2d(jim_object, jim_object2):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object.histo2d(jim_object2))


def getHisto3d(jim_object, jim_object2, jim_object3):
    """Compute the frequency distribution of the grey levels pairs.

    :param jim_object: a Jim object
    :param jim_object2: a Jim object
    :param jim_object3: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object.histo3d(jim_object2, jim_object3))


def stretch(jim_object, **kwargs):
    return _pj.Jim(jim_object.stretch(kwargs))


class _Stats():
    """Define all statistical methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def getStats(self, **kwargs):
        return self._jim_object.getStats(kwargs)

    def getStatProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        return _pj.Jim(self._jim_object.statProfile(kwargs))

    def getHisto1d(self):
        """Compute the frequency distribution of the grey levels of im.

        :return: a Jim object
        """
        return _pj.Jim(self._jim_object.histo1d())

    def getHisto2d(self, jim_object2):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object.histo2d(jim_object2))

    def getHisto3d(self, jim_object2, jim_object3):
        """Compute the frequency distribution of the grey levels pairs.

        :param jim_object2: a Jim object
        :param jim_object3: a Jim object
        :return: a Jim object
        """
        return _pj.Jim(self._jim_object.histo3d(jim_object2, jim_object3))

    def stretch(self, **kwargs):
        self._jim_object._set(self._jim_object.stretch(kwargs))


class _StatsList():
    """Define all statistical methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def getStats(self, **kwargs):
        return self._jim_list.getStats(kwargs)

    def getStatProfile(self, function, **kwargs):
        kwargs.update({'function': function})
        return _pj.Jim(self._jim_list.statProfile(kwargs))


class _StatsVect():
    """Define all statistical methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
