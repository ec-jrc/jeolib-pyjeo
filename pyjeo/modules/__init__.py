"""Define a Base class for all modules."""


class JimModuleBase:
    """Base class for Jim modules."""

    def __init__(self):
        pass

    def _set_caller(self, caller):
        self._jim_object = caller


class JimListModuleBase:
    """Base class for JimList modules."""

    def __init__(self):
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class JimVectModuleBase:
    """Base class for JimVect modules."""

    def __init__(self):
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
