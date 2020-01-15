"""Define a Base class for all modules."""


class JimModuleBase:
    """Base class for Jim modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original Jim object."""
        self._jim_object = caller


class JimListModuleBase:
    """Base class for JimList modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original JimList object."""
        self._jim_list = caller


class JimVectModuleBase:
    """Base class for JimVect modules."""

    def __init__(self):
        """Initialize the object."""
        pass

    def _set_caller(self, caller):
        """Set the reference to the original JimVect object."""
        self._jim_vect = caller
