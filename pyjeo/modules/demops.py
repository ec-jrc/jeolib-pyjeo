"""Module for operations connected to digital elevation models."""

import pyjeo as _pj


def catchmentBasinConfluence(jim_object, d8):
    """Compute the catchment basin confluence.

    :param jim_object: an image node holding labelled outlet pixels with
        value 1 and river pixels with value 2
    :param d8: an image node holding d8 flow directions
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demCatchmenBasinOutlet(d8._jipjim))


def catchmentBasinOutlet(jim_object, d8):
    """Compute the catchment basin outlet.

    :param jim_object: an image node holding labelled outlets
    :param d8: an image node holding d8 flow directions
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demCatchmentBasinOutlet(d8._jipjim))


def contribDrainArea(jim_object, graph):
    """Output the contributing drainage areas of a DEM.

    Outputs the contributing drainage areas of a DEM given its
    graph-connected drainage directions coded as follows: NW=5, N=3, NE=7,
    W=1, E=2, SW=6, S=4, SE=8.

    :param jim_object: an image node with D8 drainage directions (UCHAR)
    :param graph: integer for number of possible flow directions (either 4 or 8)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demContributingDrainageArea(graph))


def contribDrainAreaInf(jim_object):
    """Output the contributing drainage areas of a DEM.

    Outputs the contributing drainage areas of a DEM given its dinf drainage
    directions.

    :param jim_object: a Jim object with Dinf drainage directions
        (t FLOAT, -1.0 for undefined direction)
    :param graph: integer for number of nearest neighbours to consider
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demContributingDrainageAreaDInf())


def contribDrainAreaStrat(cda, threshold, dir):
    """Extract river networks.

    Do it by flagging the downstreams of all points whose contributing
    drainage areas exceed those given by the threshold image. The dowstreams
    are detected by following the drainage directions stored in the image dir.

    :param cda: an image node (INT32) for contributing drainage area
    :param threshold: an image node (USHORT) for cda threshold levels
    :param dir: an image node (UCHAR) for flow directions
    :return:
    """
    return _pj.Jim(
        cda._jipjim.demContributingDrainageAreaStratify(threshold._jipjim,
                                                        dir._jipjim))


def floodDir(jim_object, graph):
    """Compute the local flow directions.

    Compute them as the inverse of the flood wave
    direction occurring during an immersion simulation (i.e., flooding starting
    from the lowest elevations). The codes for each direction are as follows:
    NW=5, N=3, NE=7, W=1, E=2, SW=6, S=4, SE=8. When a pixel has no lower
    neighbour, it is set to 0.

    :param jim_object: a Jim object
    :param graph: integer for number of nearest neighbours to consider
        (either 4 or 8)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFloodDirection(graph))


def flow(jim_object, graph):
    """Compute the contributing drainage areas using D8 drainage directions.

    :param jim_object: a Jim object
    :param graph: integer for number of nearest neighbours to consider
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFlow(graph))


def flowDirectionD8(jim_object):
    """Compute the D8 steepest slope direction of each pixel.

    The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
    E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
    to 0.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFlowDirectionD8())


def flowDirectionDInf(jim_object):
    """Compute the dinf steepest slope direction of each pixel.

    computes the dinf steepest slope direction of each pixel according
    to (Tarboton, 1997). Slope directions are measured counter-clockwise from
    east, i.e., range equals (0,2pi)values are in the range (0,2pi),
    while pixels having no dowslope (plateaus and pits) are set to -1.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFlowDirectionDInf())


def flowDirectionFlat(jim_object, dem_jim, graph):
    """See publication (Soille, 2002).

    Flat regions (i.e., no flow direction) must be of type USHORT (with flat
    regions set to 65533) or INT32 (with flat regions set to INT32 MAX-2).

    :param jim_object: an image node for flat regions (USHORT or INT32)
    :param dem_jim: an image node for corresponding DEM (USHORT)
    :param graph: integer for number of nearest neighbours to consider
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFlowDirectionFlat(dem_jim._jipjim,
                                                           graph))


def flowDirectionFlatGeodesic(jim_object, dem_jim, graph):
    """Inverse geodesic distance Away From Ascending Border.

    :param jim_object: a Jim object
    :param dem_jim: a Jim object containing DEM
    :param graph: integer for number of nearest neighbours to consider
    :return: a Jim object
    """
    return _pj.Jim(
        jim_object._jipjim.demFlowDirectionFlatGeodesic(dem_jim._jipjim,
                                                        graph))


def flowNew(jim_object, drain_image, graph=8):
    """Compute the contributing drainage area of each pixel.

    Computes the contributing drainage area of each pixel of im given the
    graph-connected drainage directions stored in imdir.

    :param jim_object: a Jim object
    :param drain_image: the d8 drainage directions for each pixel of im
    :param graph: integer for connectivity (must be 8)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demFlowNew(drain_image._jipjim, graph))


def pitRemovalCarve(labeled_jim, grey_jim, graph, maxfl):
    """Use for carving, algorithm description in Soille et al. 2003.

    :param labeled_jim: an image node with labelled relevant minima
    :param grey_jim: an image node with grey tone image
    :param graph: an integer for connectivity
    :param maxfl: an integer for highest flooding level
    :return: a Jim object
    """
    return _pj.Jim(labeled_jim._jipjim.demPitRemovalCarve(grey_jim._jipjim,
                                                          graph, maxfl))


def pitRemovalOptimal(labeled_jim, grey_jim, graph, maxfl, flag):
    """Optimal removal of spurious pits in grid digital elevation models.

    Note that irrelevant minima must have all an intensity greater than that
    of the lowest minimum! The actual carved image is stored in imr.

    :param labeled_jim: an image node with labelled relevant minima
    :param grey_jim: an image node with grey tone image
    :param graph: an integer for connectivity
    :param maxfl: an integer for highest flooding level
    :param flag: 0 (default) for energy based, area based otherwise
    :return: a Jim object
    """
    return _pj.Jim(labeled_jim._jipjim.demPitRemovalOptimal(grey_jim._jipjim,
                                                          graph, maxfl, flag))


def slopeD8(jim_object):
    """Compute the steepest slope within a 3x3 neighbourhood for each pixel.

    It corresponds to the slope along the D8 direction.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demSlopeD8())


def slopeDInf(jim_object):
    """Output the slope along the dinf drainage directions.

    :param jim_object: a Jim object
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demSlopeDInf())


def strahler(d8):
    """Compute the Strahler thing.

    :param d8: an image node holding d8 directions on river networks,
        0 elsewhere
    :return: a Jim object
    """
    return _pj.Jim(d8._jipjim.demStrahlerOrder())


class _DEMOps():
    """Define all DEMOps methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def catchmentBasinConfluence(self, d8):
        """Compute the catchment basin confluence.

        :param jim_object: an image node holding labelled outlet pixels with
            value 1 and river pixels with value 2
        :param d8: an image node holding d8 flow directions
        :return: a Jim object
        """
        self._jim_object._jipjim.d_demCatchmenBasinConfluence(d8._jipjim)

    def catchmentBasinOutlet(self, d8):
        """Compute the catchment basin outlet.

        :param jim_object: an image node holding labelled outlets
        :param d8: an image node holding d8 flow directions
        :return: a Jim object
        """
        self._jim_object._jipjim.d_demCatchmentBasinOutlet(d8._jipjim)

    def contribDrainArea(self, graph):
        """Output the contributing drainage areas of a DEM.

        Outputs the contributing drainage areas of a DEM given its
        graph-connected drainage directions coded as follows: NW=5, N=3, NE=7,
        W=1, E=2, SW=6, S=4, SE=8.

        :param jim_object: an image node with D8 drainage directions (UCHAR)
        :param graph: integer for number of possible flow directions
            (either 4 or 8)
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demContributingDrainageArea(graph))

    def contribDrainAreaInf(self):
        """Output the contributing drainage areas of a DEM.

        Outputs the contributing drainage areas of a DEM given its dinf
        drainage directions. Jim object must be with Dinf drainage directions
        (t FLOAT, -1.0 for undefined direction).

        :param graph: integer for number of nearest neighbours to consider
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demContributingDrainageAreaDInf())

    def contribDrainAreaStrat(self, threshold, dir):
        """Extract river networks.

        Do it by flagging the downstreams of all points whose contributing
        drainage areas exceed those given by the threshold image. The
        dowstreams are detected by following the drainage directions stored
        in the image dir. Jim must be an image node (INT32) for contributing
        drainage area.

        :param threshold: an image node (USHORT) for cda threshold levels
        :param dir: an image node (UCHAR) for flow directions
        :return:
        """
        self._jim_object._set(
            self._jim_object._jipjim.demContributingDrainageAreaStratify(
                threshold._jipjim, dir._jipjim))

    def floodDir(self, graph):
        """Compute the local flow directions.

        Compute them as the inverse of the flood wave
        direction occurring during an immersion simulation (i.e., flooding starting
        from the lowest elevations). The codes for each direction are as follows:
        NW=5, N=3, NE=7, W=1, E=2, SW=6, S=4, SE=8. When a pixel has no lower
        neighbour, it is set to 0.

        :param jim_object: a Jim object
        :param graph: integer for number of nearest neighbours to consider
            (either 4 or 8)
        :return: a Jim object
        """
        self._jim_object._jipjim.d_demFloodDirection(graph)

    def flow(self, graph):
        """Compute the contributing drainage areas using D8 drainage directions.

        :param jim_object: a Jim object
        :param graph: integer for number of nearest neighbours to consider
        :return: a Jim object
        """
        self._jim_object._set(self._jim_object._jipjim.demFlow(graph))

    def flowDirectionD8(self):
        """Compute the D8 steepest slope direction of each pixel.

        The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
        E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
        to 0.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(self._jim_object._jipjim.demFlowDirectionD8())

    def flowDirectionDInf(self):
        """Compute the dinf steepest slope direction of each pixel.

        computes the dinf steepest slope direction of each pixel according
        to (Tarboton, 1997). Slope directions are measured counter-clockwise from
        east, i.e., range equals (0,2pi)values are in the range (0,2pi),
        while pixels having no dowslope (plateaus and pits) are set to -1.

        :param jim_object: a Jim object
        :return: a Jim object
        """
        self._jim_object._set(self._jim_object._jipjim.demFlowDirectionDInf())

    def flowDirectionFlat(self, dem_jim, graph):
        """See publication (Soille, 2002).

        Flat regions (i.e., no flow direction) must be of type USHORT (with flat
        regions set to 65533) or INT32 (with flat regions set to INT32 MAX-2).

        :param jim_object: an image node for flat regions (USHORT or INT32)
        :param dem_jim: an image node for corresponding DEM (USHORT)
        :param graph: integer for number of nearest neighbours to consider
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demFlowDirectionFlat(dem_jim._jipjim,
                                                          graph))

    def flowDirectionFlatGeodesic(self, dem_jim, graph):
        """Inverse geodesic distance Away From Ascending Border.

        :param jim_object: a Jim object
        :param dem_jim: a Jim object containing DEM
        :param graph: integer for number of nearest neighbours to consider
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demFlowDirectionFlatGeodesic(
                dem_jim._jipjim, graph))

    def flowNew(self, drain_image, graph=8):
        """Compute the contributing drainage area of each pixel.

        Computes the contributing drainage area of each pixel of im given the
        graph-connected drainage directions stored in imdir.

        :param jim_object: a Jim object
        :param drain_image: the d8 drainage directions for each pixel of im
        :param graph: integer for connectivity (must be 8)
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demFlowNew(drain_image._jipjim,
                                                graph))

    def pitRemovalCarve(self, grey_jim, graph, maxfl):
        """Use for carving, algorithm description in Soille et al. 2003.

        :param labeled_jim: an image node with labelled relevant minima
        :param grey_jim: an image node with grey tone image
        :param graph: an integer for connectivity
        :param maxfl: an integer for highest flooding level
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demPitRemovalCarve(grey_jim._jipjim,
                                                        graph, maxfl))

    def pitRemovalOptimal(self, grey_jim, graph, maxfl, flag):
        """Optimal removal of spurious pits in grid digital elevation models.

        Note that irrelevant minima must have all an intensity greater than that
        of the lowest minimum! The actual carved image is stored in imr.

        :param labeled_jim: an image node with labelled relevant minima
        :param grey_jim: an image node with grey tone image
        :param graph: an integer for connectivity
        :param maxfl: an integer for highest flooding level
        :param flag: 0 (default) for energy based, area based otherwise
        :return: a Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.demPitRemovalOptimal(grey_jim._jipjim,
                                                        graph, maxfl, flag))

    def slopeD8(self):
        """Compute the steepest slope within a 3x3 neighbourhood for each pixel.

        It corresponds to the slope along the D8 direction.

        :param jim_object: a Jim object
        :return: a Jim object
        """
        self._jim_object._set(self._jim_object._jipjim.demSlopeD8())

    def slopeDInf(self):
        """Output the slope along the dinf drainage directions.

        :param jim_object: a Jim object
        :return: a Jim object
        """
        self._jim_object._set(self._jim_object._jipjim.demSlopeDInf())

    def strahler(self):
        """Compute the Strahler thing.

        :param d8: an image node holding d8 directions on river networks,
            0 elsewhere
        :return: a Jim object
        """
        self._jim_object._jipjim.d_demStrahlerOrder()


class _DEMOpsList():
    """Define all DEMOps methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class _DEMOpsVect():
    """Define all DEMOps methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
