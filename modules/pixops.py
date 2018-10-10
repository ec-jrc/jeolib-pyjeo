import jiplib as _jl


class _PixOps():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def pointOpBitWise(self, sec_jim_object, operation_code):
        """Bitwise operation between two images.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param operation_code: 10 or AND op, 11 or OR op, and 12 or XOR op
        """
        self._jim_object._set(self._jim_object.pointOpBitwise(
            self._jim_object, sec_jim_object, operation_code))

    def pointOpBlank(self, value):
        """Set all pixels of image to value.

        Modifies the instance on which the method was called.

        :param value: new value for pixels of Jim object
        """
        self._jim_object.d_pointOpBlank(value)
