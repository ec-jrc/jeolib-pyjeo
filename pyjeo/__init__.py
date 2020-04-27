from .pyjeo import *

import os as _os
import tempfile as _tempfile
import random as _random
import string as _string


__version__ = '0.5.1'


def _check_graph(graph, allowed_values):
    """Check whether values used as graph are allowed.

    :param graph: an integer holding for the graph connectivity
    :param allowed_values: values allowed for the graph parameter
    """
    if graph not in allowed_values:
        raise ValueError('Value {} not allow as a graph parameter. Only values'
                         ' {} are allowed.'.format(graph, allowed_values))


def _get_random_path():
    """Return path of non-existing file in the temp directory.

    Needed for intermediate products for JimVect destructive methods.
    """
    random_string = ''.join(_random.sample(_string.ascii_letters, 5))
    temp_dir = _tempfile.gettempdir()
    non_existing_path = _os.path.join(temp_dir, random_string)
    while _os.path.isfile(non_existing_path):
        random_string = ''.join(_random.sample(_string.ascii_letters, 5))
        non_existing_path = _os.path.join(temp_dir, random_string)

    return non_existing_path
