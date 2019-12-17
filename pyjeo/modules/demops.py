"""Module for operations connected to digital elevation models."""

import numpy

import pyjeo as _pj
from . import JimModuleBase as _JimModuleBase
from . import JimListModuleBase as _JimListModuleBase
from . import JimVectModuleBase as _JimVectModuleBase


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
    :param graph: integer for number of possible flow directions
        (either 4 or 8)
    :return: a Jim object
    """
    _pj._check_graph(graph, [4, 8])

    return _pj.Jim(jim_object._jipjim.demContributingDrainageArea(graph))


def contribDrainAreaInf(jim_object):
    """Output the contributing drainage areas of a DEM.

    Outputs the contributing drainage areas of a DEM given its dinf drainage
    directions.

    :param jim_object: a Jim object with Dinf drainage directions
        (t FLOAT, -1.0 for undefined direction)
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
    _pj._check_graph(graph, [4, 8])

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
    """See publication :cite:`soille2002dgci`).

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
    _pj._check_graph(graph, [8])

    return _pj.Jim(jim_object._jipjim.demFlowNew(drain_image._jipjim, graph))


def pitRemovalCarve(labeled_jim, grey_jim, graph, maxfl):
    """Use for carving.

    Algorithm description in :cite:`soille-vogt-colombo2003wrr` and
    :cite:`soille2004prl`

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
    :cite:`soille2004wrr`

    :param labeled_jim: an image node with labelled relevant minima
    :param grey_jim: an image node with grey tone image
    :param graph: an integer for connectivity
    :param maxfl: an integer for highest flooding level
    :param flag: 0 (default) for energy based, area based otherwise
    :return: a Jim object
    """
    return _pj.Jim(labeled_jim._jipjim.demPitRemovalOptimal(grey_jim._jipjim,
                                                            graph, maxfl,
                                                            flag))


def slope(jim_object, scale=1.0, zscale=1.0, percent=False):
    """Compute the slope of a Jim object.

    :param jim_object: Jim
    :param scale: horizontal scale
    :param zscale: vertical scale
    :param percent: if True, return value in percents, degrees otherwise
    :return: a Jim object representing the slope
    """
    tapsdx = numpy.array(
        [[-1.0, 0.0, 1.0], [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]])
    tapsdy = numpy.array(
        [[-1.0, -2.0, -1.0], [0.0, 0.0, 0.0], [1.0, 2.0, 1.0]])
    tapsdx *= zscale
    tapsdy *= zscale
    jimdx = _pj.Jim(jim_object)
    jimdy = _pj.Jim(jim_object)
    if jim_object.properties.getDataType() != 'Float32' and \
       jim_object.properties.getDataType() != 'Float64':
        jimdx.pixops.convert(otype="Float32")
        jimdy.pixops.convert(otype="Float32")
    jimdx.ngbops.firfilter2d(
        tapsdx, nodata=jim_object.properties.getNoDataVals(), norm=True)
    jimdx = abs(jimdx)
    jimdx /= jimdx.properties.getDeltaX() * scale
    jimdx *= jimdx
    jimdy.ngbops.firfilter2d(
        tapsdy, nodata=jim_object.properties.getNoDataVals(), norm=True)
    jimdy = abs(jimdy)
    jimdy /= jimdy.properties.getDeltaX()*scale
    jimdy *= jimdy
    rad2deg = 180.0 / numpy.pi
    jimdx += jimdy
    jimdx.np()[:] = numpy.sqrt(jimdx.np())
    if percent:
        jimdx *= 100
    else:
        jimdx.np()[:] = numpy.arctan(jimdx.np())
        jimdx *= rad2deg
    return jimdx


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


# this numpy implementation works and is faster but does not ignore nodata
# values in the image during the convolution process
# def slopenp(jim_object, scale=1.0, zscale=1.0, percent=False, nodata=None):
#     if jim_object.properties.getNoDataVals() and not nodata:
#         nodata=jim_object.properties.getNoDataVals()[0]
#     tapsdx=numpy.array([[-1.0,0.0,1.0],[-2.0,0.0,2.0],[-1.0,0.0,1.0]])
#     tapsdy=numpy.array([[-1.0,-2.0,-1.0],[0.0,0.0,0.0],[1.0,2.0,1.0]])
#     tapsdx*=zscale
#     tapsdy*=zscale
#     jimdx=_pj.Jim(jim_object)
#     jimdy=_pj.Jim(jim_object)
#     if jim_object.properties.getDataType() != 'Float32' and \
#        jim_object.properties.getDataType() != 'Float64':
#         jimdx.pixops.convert(otype="Float32")
#         jimdy.pixops.convert(otype="Float32")
#     if jim_object.properties.getNoDataVals():
#         for ndval in jim_object.properties.getNoDataVals():
#             jimdx[jimdx==ndval]=nodata
#         jimdx.np()[jimdx.np()!=nodata]=signal.convolve2d(jim_object.np(),tapsdx,boundary='symm',mode='same')
#     else:
#         jimdx.np()[:]=signal.convolve2d(jim_object.np(),tapsdx,boundary='symm',mode='same')
#     jimdx/=8.0*jimdx.properties.getDeltaX()*scale
#     jimdx*=jimdx
#     if jim_object.properties.getNoDataVals():
#         for ndval in jim_object.properties.getNoDataVals():
#             jimdy[jimdy==ndval]=nodata
#         jimdy.np()[jimdy.np()!=nodata]=signal.convolve2d(jim_object.np(),tapsdy,boundary='symm',mode='same')
#     else:
#         jimdy.np()[:]=signal.convolve2d(jim_object.np(),tapsdy,boundary='symm',mode='same')
#     jimdy/=8.0*jimdy.properties.getDeltaX()*scale
#     jimdy*=jimdy
#     rad2deg=180.0/numpy.pi
#     jimdx+=jimdy
#     jimdx.np()[:]=numpy.sqrt(jimdx.np())
#     if percent:
#         jimdx*=100
#     else:
#         jimdx.np()[:]=numpy.arctan(jimdx.np())
#         jimdx*=rad2deg
#         # jimdx=90-jimdx
#     return jimdx


def strahler(jim_object):
    """Compute the Strahler thing.

    :param jim_object: an image node holding d8 directions on river networks,
        0 elsewhere
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.demStrahlerOrder())


class _DEMOps(_JimModuleBase):
    """Define all DEMOps methods."""

    def catchmentBasinConfluence(self, d8):
        """Compute the catchment basin confluence.

        :param jim_object: an image node holding labelled outlet pixels with
            value 1 and river pixels with value 2
        :param d8: an image node holding d8 flow directions

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_demCatchmenBasinConfluence(d8._jipjim)

    def catchmentBasinOutlet(self, d8):
        """Compute the catchment basin outlet.

        :param jim_object: an image node holding labelled outlets
        :param d8: an image node holding d8 flow directions

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_demCatchmentBasinOutlet(d8._jipjim)

    def contribDrainArea(self, graph):
        """Output the contributing drainage areas of a DEM.

        Outputs the contributing drainage areas of a DEM given its
        graph-connected drainage directions coded as follows: NW=5, N=3, NE=7,
        W=1, E=2, SW=6, S=4, SE=8.

        Modifies the instance on which the method was called.

        :param jim_object: an image node with D8 drainage directions (UCHAR)
        :param graph: integer for number of possible flow directions
            (either 4 or 8)

        Example: Compute the contributing drain area of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          # input for contribDrainArea() must be a D8 drainage object
          jim.demops.flowDirectionD8()
          jim.demops.contribDrainArea(8)
        """
        _pj._check_graph(graph, [4, 8])

        self._jim_object._set(
            self._jim_object._jipjim.demContributingDrainageArea(graph))

    def contribDrainAreaInf(self):
        """Output the contributing drainage areas of a DEM.

        Outputs the contributing drainage areas of a DEM given its dinf
        drainage directions. Jim object must be with Dinf drainage directions
        (t FLOAT, -1.0 for undefined direction).

        Modifies the instance on which the method was called.

        :param graph: integer for number of nearest neighbours to consider
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

        Modifies the instance on which the method was called.

        :param threshold: an image node (USHORT) for cda threshold levels
        :param dir: an image node (UCHAR) for flow directions

        Example: Compute the contributing drain area strat of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          # we need a d8 object
          d8 = pj.demops.flowDirectionD8(jim)
          # we need a cda object
          cda = pj.demops.contribDrainArea(d8, 8)
          # we need a threshold Jim
          thresh = pj.pixops.setData(jim, 5)
          # now we can compute the contribDrainAreaStrat()
          cda.demops.contribDrainAreaStrat(thresh, d8)
        """
        self._jim_object._set(
            self._jim_object._jipjim.demContributingDrainageAreaStratify(
                threshold._jipjim, dir._jipjim))

    def floodDir(self, graph):
        """Compute the local flow directions.

        Compute them as the inverse of the flood wave
        direction occurring during an immersion simulation (i.e., flooding
        starting from the lowest elevations). The codes for each direction
        are as follows: NW=5, N=3, NE=7, W=1, E=2, SW=6, S=4, SE=8. When a
        pixel has no lower neighbour, it is set to 0.

        Modifies the instance on which the method was called.

        :param graph: integer for number of nearest neighbours to consider
            (either 4 or 8)

        Example: Compute the local flow directions of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          jim.demops.floodDir()
        """
        _pj._check_graph(graph, [4, 8])

        self._jim_object._jipjim.d_demFloodDirection(graph)

    def flow(self, graph):
        """
        Compute the contributing drainage areas using D8 drainage directions.

        Modifies the instance on which the method was called.

        :param graph: integer for number of nearest neighbours to consider

        Example: Compute the contributing areas of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          # we need a flow direction d8 object
          jim.demops.flowDirectionD8()
          jim.demops.flow(8)
        """
        self._jim_object._set(self._jim_object._jipjim.demFlow(graph))

    def flowDirectionD8(self):
        """Compute the D8 steepest slope direction of each pixel.

        The codes for each direction are as follows: NW=5, N=3, NE=7, W=1,
        E=2, SW=6, S=4, SE=8. When a pixel has no lower neighbour, it is set
        to 0.

        Modifies the instance on which the method was called.

        Example: Compute the D8 steepest slope direction of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          jim.demops.flowDirectionD8()
        """
        self._jim_object._set(self._jim_object._jipjim.demFlowDirectionD8())

    def flowDirectionDInf(self):
        """Compute the dinf steepest slope direction of each pixel.

        computes the dinf steepest slope direction of each pixel according
        to (Tarboton, 1997). Slope directions are measured counter-clockwise
        from east, i.e., range equals (0,2pi)values are in the range (0,2pi),
        while pixels having no dowslope (plateaus and pits) are set to -1.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(self._jim_object._jipjim.demFlowDirectionDInf())

    def flowDirectionFlat(self, dem_jim, graph):
        """See publication (Soille, 2002).

        Flat regions (i.e., no flow direction) must be of type USHORT (with
        flat regions set to 65533) or INT32 (with flat regions set to INT32
        MAX-2).

        Modifies the instance on which the method was called.

        :param jim_object: an image node for flat regions (USHORT or INT32)
        :param dem_jim: an image node for corresponding DEM (USHORT)
        :param graph: integer for number of nearest neighbours to consider
        """
        self._jim_object._set(
            self._jim_object._jipjim.demFlowDirectionFlat(dem_jim._jipjim,
                                                          graph))

    def flowDirectionFlatGeodesic(self, dem_jim, graph):
        """Inverse geodesic distance Away From Ascending Border.

        Modifies the instance on which the method was called.

        :param dem_jim: a Jim object containing DEM
        :param graph: integer for number of nearest neighbours to consider
        """
        self._jim_object._set(
            self._jim_object._jipjim.demFlowDirectionFlatGeodesic(
                dem_jim._jipjim, graph))

    def flowNew(self, drain_image, graph=8):
        """Compute the contributing drainage area of each pixel.

        Computes the contributing drainage area of each pixel of im given the
        graph-connected drainage directions stored in imdir.

        Modifies the instance on which the method was called.

        :param drain_image: the d8 drainage directions for each pixel of im
        :param graph: integer for connectivity (must be 8)

        Example: Compute the contributing drainage area of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          # we need a flow direction d8 object
          flow = pj.demops.flowDirectionD8(jim)
          jim.demops.flowNew(flow)
        """
        _pj._check_graph(graph, [8])

        self._jim_object._set(
            self._jim_object._jipjim.demFlowNew(drain_image._jipjim,
                                                graph))

    def pitRemovalCarve(self, grey_jim, graph, maxfl):
        """Use for carving, algorithm description in Soille et al. 2003.

        Modifies the instance on which the method was called.

        :param labeled_jim: an image node with labelled relevant minima
        :param grey_jim: an image node with grey tone image
        :param graph: an integer for connectivity
        :param maxfl: an integer for highest flooding level

        Example: Compute the pit removal carve of a Jim::

          jim0 = pj.Jim('/path/to/raster.tif')
          # we need a labeled object
          jim1 = pj.ccops.labelImagePixels(jim0)
          jim1.demops.pitRemovalCarve(jim0, 8, 212)
        """
        self._jim_object._set(
            self._jim_object._jipjim.demPitRemovalCarve(grey_jim._jipjim,
                                                        graph, maxfl))

    def pitRemovalOptimal(self, grey_jim, graph, maxfl, flag):
        """Optimal removal of spurious pits in grid digital elevation models.

        Note that irrelevant minima must have all an intensity greater than
        that of the lowest minimum! The actual carved image is stored in imr.

        Modifies the instance on which the method was called.

        :param labeled_jim: an image node with labelled relevant minima
        :param grey_jim: an image node with grey tone image
        :param graph: an integer for connectivity
        :param maxfl: an integer for highest flooding level
        :param flag: 0 (default) for energy based, area based otherwise

        Example: Compute the optimal pit removal carve of a Jim::

          jim0 = pj.Jim('/path/to/raster.tif')
          # we need a labeled object
          jim1 = pj.ccops.labelImagePixels(jim0)
          jim1.demops.pitRemovalOptimal(jim0, 8, 212, 0)
        """
        self._jim_object._set(
            self._jim_object._jipjim.demPitRemovalOptimal(grey_jim._jipjim,
                                                          graph, maxfl, flag))

    def slopeD8(self):
        """
        Compute the steepest slope within a 3x3 neighbourhood for each pixel.

        It corresponds to the slope along the D8 direction.

        Modifies the instance on which the method was called.

        Example: Compute the steepest slope of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          jim.demops.slopeD8()
        """
        self._jim_object._set(self._jim_object._jipjim.demSlopeD8())

    def slopeDInf(self):
        """Output the slope along the dinf drainage directions."""
        self._jim_object._set(self._jim_object._jipjim.demSlopeDInf())

    def strahler(self):
        """Compute the Strahler order.

        Computes the Strahler order on a Jim object holding d8 directions on
        river networks, 0 elsewhere.

        Modifies the instance on which the method was called.

        Example: Compute the Strahler order of a Jim::

          jim = pj.Jim('/path/to/raster.tif')
          # we need a flow direction d8 object
          jim.demops.flowDirectionD8()
          jim.demops.strahler()
        """
        self._jim_object._jipjim.d_demStrahlerOrder()


class _DEMOpsList(_JimListModuleBase):
    """Define all DEMOps methods for JimLists."""

    pass


class _DEMOpsVect(_JimVectModuleBase):
    """Define all DEMOps methods for JimVects."""

    pass
