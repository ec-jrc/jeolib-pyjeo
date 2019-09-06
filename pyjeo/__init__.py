from .pyjeo import *


__version__ = '0.5.0'


def _check_graph(graph, allowed_values):
    """Check whether values used as graph are allowed.

    :param graph: an integer holding for the graph connectivity
    :param allowed_values: values allowed for the graph parameter
    """
    if graph not in allowed_values:
        raise ValueError('Value {} not allow as a graph parameter. Only values'
                         ' {} are allowed.'.format(graph, allowed_values))
