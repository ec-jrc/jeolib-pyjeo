"""Module for connected-component operations."""
# Author(s): Pieter.Kempeneers@ec.europa.eu,
#            Ondrej Pesek,
#            Pierre.Soille@ec.europa.eu
# Copyright (C) 2018-2022 European Union (Joint Research Centre)
#
# This file is part of pyjeo.
#
# pyjeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyjeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyjeo.  If not, see <https://www.gnu.org/licenses/>.

import pyjeo as _pj


# TODO: Test
def alphaTree(atree,
              attr0cc,
              type: int):
    """Compute the alpha tree.

    :param atree: a Jim object
    :param attr0cc: a Jim object
    :param type: an integer
    :return: Jim object with the alpha tree
    """
    return _pj.Jim(atree._jipjim.alphaTree(attr0cc._jipjim, type))


def alphaTreeDissim(dissimh,
                    dissimv,
                    alpha: int):
    """Create Jim holding the tree.

    This function does not have the member function alternative.

    :param dissimh: Jim for horizontal edge dissimilarities
    :param dissimv: Jim for vertical edge dissimilarities
    :param alpha: integer for dissimilarity threshold value
    :return: the tree stored in an JimList
    """
    return _pj.JimList(dissimh._jipjim.alphaTreeDissimGet(dissimv._jipjim,
                                                          alpha))


def alphaTreeNextLevel(atree,
                       label_jim,
                       alpha: int):
    """Compute the alpha tree next level.

    :param atree: a Jim object
    :param label_jim: a Jim object
    :param alpha: integer
    :return: Jim object representing the next level of alpha tree
    """
    return _pj.Jim(atree._jipjim.alphaTreeNextLevel(label_jim._jipjim, alpha))


def alphaTreePersistenceLut(atree):
    """Compute the alpha tree persistence lut.

    :param atree: a Jim object
    :return: Jim object representing the alpha tree persistence lut
    """
    return _pj.Jim(atree._jipjim.alphaTreeGetPersistenceLut())


def alphaTreeToCCs(atree,
                   label_jim,
                   lut: bool,
                   rule: int):
    """Convert an alpha tree to connected components.

    :param atree: a Jim object
    :param label_jim: a Jim object
    :param lut: flag
    :param rule: integer
    :return: Jim object representing alpha tree converted to CCs
    """
    return _pj.Jim(atree._jipjim.alphaTreeToCCs(label_jim._jipjim, lut, rule))


def convertRgbToHsx(jim,
                    x_type: str):
    """Convert RGB to HSX.

    Returns the hue, saturation, and value, lightness, or intensity channels
    of an input RGB colour image. The hue component is identical for all 3
    models. The luminance is equal to max(R,G,B) for HSV, (max-min)/2 for HSL
    and (R+G+B)/3 for HSI. See specific formulae for the saturation at
    `<http://en.wikipedia.org/wiki/HSL_and_HSV>`__`.

    :param jim: multi-band Jim with three bands
        representing red, green and blue channels
    :param x_type: string with key ('V' (default) for Value, 'L' for Lightness,
        and 'I' for Intensity)
    :return: Jim with three bands containing the HSX channels
    """
    assert jim.properties.nrOfBand() == 3, \
        'Error: input jim must be multi-band image with three bands (r, g, b)'
    jimr = _pj.geometry.cropBand(jim, 0)
    jimg = _pj.geometry.cropBand(jim, 1)
    jimb = _pj.geometry.cropBand(jim, 2)

    jimlist = _pj.JimList(jimr._jipjim.convertRgbToHsx(jimg, jimb, x_type))
    return jimlist.geometry.stackBand()


def convertHsiToRgb(jim):
    """Convert HSI to RGB.

    Takes the hue, saturation, and intensity channels of a colour image
    and returns an image node containing a colour RGB image.

    :param jim: multi-band Jim with three bands representing hue, saturation,
        and intensity channels
    :return: Jim with three bands containing the RGB channels
    """
    assert jim.properties.nrOfBand() == 3, \
        'Error: input jim must be multi-band image with three bands (h, s, i)'
    jimh = _pj.geometry.cropBand(jim, 0)
    jims = _pj.geometry.cropBand(jim, 1)
    jimi = _pj.geometry.cropBand(jim, 2)

    return _pj.Jim(jimh._jipjim.convertRgbToHsx(jims, jimi))


def convertHlsToRgb(jim):
    """Convert HLS to RGB.

    Takes the hue, lightness, and saturation channels of a colour image
    and returns an image node containing a colour RGB image.

    :param jim: multi-band Jim with three bands representing hue, lightness,
        and saturation channels
    :return: Jim with three bands containing the RGB channels
    """
    assert jim.properties.nrOfBand() == 3, \
        'Error: input jim must be multi-band image with three bands (h, s, i)'
    jimh = _pj.geometry.cropBand(jim, 0)
    jiml = _pj.geometry.cropBand(jim, 1)
    jims = _pj.geometry.cropBand(jim, 2)

    return _pj.Jim(jimh._jipjim.convertRgbToHsx(jiml, jims))


def dbscan(jim_dissim,
           eps: float,
           min_pts: int):
    """Compute dbscan.

    :param jim_dissim: a Jim object containing the dissimilarity matrix
    :param eps: float
    :param min_pts: integer
    :return: a Jim object
    """
    return _pj.Jim(jim_dissim._jipjim.dbscan(eps, min_pts))


def dissim(jim_map,
           jim_mask,
           nc: int,
           type: int):
    """Compute the dissimilarity matrix.

    :param jim_map: a Jim object
    :param jim_mask: a Jim object
    :param nc: integer
    :param type: integer
    :return: Jim object representing the dissimilarity matrix
    """
    return _pj.Jim(jim_map._jipjim.dissimilarityMatrix(jim_mask._jipjim,
                                                       nc, type))


def dissimToAlphaCCs(dissimh,
                     dissimv,
                     alpha: int):
    """Create Jim holding the labelled alpha-connected component.

    :param dissimh: Jim for horizontal edge dissimilarities
    :param dissimv: Jim for vertical edge dissimilarities
    :param alpha: integer for dissimilarity threshold value
    :return: Jim object
    """
    return _pj.Jim(dissimh._jipjim.dissimToAlphaCCs(dissimv._jipjim,
                                                    alpha))


def distance2d4(jim,
                band: int = 0):
    """Compute the 2-dimensional 4-connected distance function of a Jim.

    :param jim: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distance2d4(band))


def distance2dChamfer(jim,
                      type: int,
                      band: int = 0):
    """Compute the chamfer distance function of Jim.

    Note that the type of Jim must be large enough to hold the largest
    distance otherwise overflow may occur and generate artefacts in the
    resulting distance function.

    :param jim: a Jim object
    :param type: Integer for type of chamfer distance {1, 11, 34, 57, 5711}
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distance2dChamfer(type, band))


def distance2dEuclideanConstrained(marker,
                                   mask,
                                   band: int = 0):
    """Compute the Euclidean geodesic distance from the marker set.

    Defined by the image marker and within the geodesic mask defined by
    the image mask. The algorithm is described in (:cite:`soille91`).

    Modifies the instance on which the method was called.

    :param marker: Image marker
    :param mask: binary UCHAR image with geodesic mask pixels set to 1
    :param band: List of band indices to crop (index is 0 based)
    """
    return _pj.Jim(marker._jipjim.distance2dEuclideanConstrained(mask._jipjim,
                                                                 band))


def distance2dEuclideanSquared(jim,
                               band: int = 0):
    """Compute the squared Euclidean distance transform in 2-D.

    jim must be a binary image. Multi-plane images will be processed
    plane by plane as 2-D images. The original algorithm was proposed
    by :cite:`saito-toriwaki94` and then optimised independently by
    :cite:`hirata96` and :cite:`meijster-roerdink-hesselink2000`.
    Based on the Euclidean distance transform. Note that a temporary buffer
    of type UINT16 is used for sums along/lines and columns so that
    uncontrolled results will occur if an object shows more than 16 2 /2
    foreground pixels along a given line or column.

    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """

    if jim.properties.nrOfPlane() < 2:
        return _pj.Jim(jim._jipjim.distance2dEuclideanSquared(band))
    else:
        return_jim = _pj.Jim(jim)
        return_jim.ccops.distance2dEuclideanSquared(band)
        return return_jim


def distanceGeodesic(mask,
                     marker,
                     graph: int,
                     band: int = 0):
    """Compute geodesic distance function from marker within mask.

    Using graph connectivity.

    :param mask: an image node for geodesic mask
    :param marker: an image node for marker image
    :param graph: integer for connectivity
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(mask._jipjim.jldistanceGeodesic(
        marker._jipjim, graph, band))


def distanceInfluenceZones2dEuclidean(jim,
                                      band: int = 0):
    """Output error-free influence zones of the labelled regions in Jim.

    Jim must be a 2-D image, its last bit being reserved (e.g., for a UCHAR
    images, authorized label values range from 1 to 127). Algorithm based
    on the Euclidean distance transform.

    :param jim: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distanceInfluenceZones2dEuclidean(band))


def flatZonesSeeded(jim1,
                    jim2,
                    jim3,
                    ox: float,
                    oy: float,
                    oz: float):
    """Compute seeded flat zones.

    :param jim1: a Jim object
    :param jim2: a Jim object
    :param jim3: a Jim object
    :param ox: x coordinate of origin
    :param oy: y coordinate of origin
    :param oz: z coordinate of origin
    :return: a Jim object
    """
    return _pj.Jim(jim1._jipjim.flatZonesSeeded(jim2, jim3, ox, oy, oz))


def getRegionalMinima(jim,
                      graph: int):
    """Compute regional minima of the input image.

    The pixels belonging to a regional minimum are set to 1, all other pixels
    are set to 0.

    :param jim: a Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object of type unsigned char containing regional
        minima of the input Jim object
    """
    return _pj.Jim(jim._jipjim.getRegionalMinima(graph))


def label(jim,
          ngb):
    """Label each connected component (non zero values) with a unique integer value using graph-connectivity.

    The labels start from value 2.

    :param jim: a Jim object holding a binary image,
        data must be of unsigned integer type.
    :param ngb: Jim object for neighbourhood, e.g., create with:

         - pj.Jim(graph=4): horizontally and vertically connected in 2-D
         - pj.Jim(graph=6): horizontally and vertically connected in 3-D
         - pj.Jim(graph=8): horizontally, vertically, and diagonally connected in 2-D
         - pj.Jim(graph=26): horizontally, vertically, and diagonally connected in 3-D
    :return: labeled Jim object with unique integer values (starting from value 2) for each connected component.

    Example:

    Label horizontally and vertically connected pixels in 2-D image with unique labels::

        jim = pj.Jim(ncol = 5, nrow = 5, otype = 'GDT_UInt16', uniform = [0,2])
        labeled = pj.ccops.label(jim, pj.Jim(graph=4))

    binary input image:

    .. image:: figures/ccops_label2d_input.png
       :width: 65 %

    labeled output image with unique integer values for connected components (here represented by different gray scales):

    .. image:: figures/ccops_label2d_output.png
       :width: 65 %

    Label horizontally and vertically connected pixels in 3-D image with unique labels::

        jim = pj.Jim(ncol = 3, nrow = 3, nplane = 3, otype = 'GDT_Byte', uniform = [0,2])
        labeled = pj.ccops.label(jim, pj.Jim(graph=6))

    binary input 3-D image:

    .. image:: figures/ccops_label3d_input.png
       :width: 65 %


    labeled output image with unique integer values for connected components (here represented by different gray scales):

    .. image:: figures/ccops_label3d_output.png
       :width: 65 %


    """

    dim = 0
    if jim.properties.nrOfPlane() > 1:
        dim = 1;
    jim.geometry.imageFrameAdd(1, 1, 1, 1, dim, dim, 0)
    ret_jim = _pj.Jim(jim._jipjim.labelBinary(ngb._jipjim, 1, 1, dim))
    jim.geometry.imageFrameSubtract(1, 1, 1, 1, dim, dim)
    ret_jim.geometry.imageFrameSubtract(1, 1, 1, 1, dim, dim)
    return ret_jim


def labelConstrainedCCs(jim,
                        local_range: int,
                        global_range: int,
                        ngb):
    """Label each alpha-omega connected component.

    Label with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim: a Jim or Jim list of grey level images having all the same
        definition domain and data type.
    :param local_range: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param global_range: integer value indicating maximum global difference
        (difference between the maximum and minimum values of each resulting
        connected component)
    :param ngb: Jim object for neighbourhood, e.g., create with pj.Jim(graph=4)
    :return: labeled Jim object
    """
    if isinstance(jim, _pj.Jim):
        if jim.properties.nrOfBand() == 1:
            return _pj.Jim(jim._jipjim.labelConstrainedCCs(
                ngb._jipjim, 1, 1, 0, global_range, local_range))
        else:
            return _pj.Jim(jim._jipjim.labelConstrainedCCsMultiband(
                ngb._jipjim, 1, 1, 0, global_range, local_range))
    elif isinstance(jim, _pj.JimList):
        return _pj.Jim(jim._jipjimlist.labelConstrainedCCsMultiband(
            ngb._jipjim, 1, 1, 0, global_range, local_range))
    else:
        raise _pj.exceptions.JimIllegalArgumentError(
            "input must be Jim or JimList object")


def labelConstrainedCCsAttr(jim,
                            graph: int,
                            rg: int,
                            rl: int):
    """Label image, in development.

    :param jim: a Jim object to label
    :param graph: integer for graph connectivity
    :param rg: integer
    :param rl: integer for range parameter lambda l under the strongly
        connected assumption
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelConstrainedCCsAttr(
        graph, rg, rl))


def labelConstrainedCCsCi(jim,
                          ngb,
                          ox: float,
                          oy: float,
                          oz: float,
                          rl: int):
    """Label image, in development.

    :param jim: a Jim object to label
    :param ngb: a Jim object for neighbourhood
    :param ox: x coordinate of origin of ngb Jim
    :param oy: y coordinate of origin of ngb Jim
    :param oz: z coordinate of origin of ngb Jim
    :param rl: integer for range parameter lambda l under the strongly
        connected assumption
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelConstrainedCCsCi(
        ngb._jipjim, ox, oy, oz, rl))


def labelConstrainedCCsDissim(jim,
                              local_range: int,
                              global_range: int,
                              dissim_type: int = 0):
    """Label each alpha-omega connected components with a unique label.

     Label using graph-connectivity and the dissimilarity measure countering
     the chaining effect as described in :cite:`soille2011ismm`

    :param jim: a Jim or a Jim list of grey level images having all the same
        definition domain and data type.
    :param local_range: integer value indicating maximum absolute local
        difference between 2 adjacent pixels along alpha-connected paths
    :param global_range: integer value indicating maximum global difference
        (difference between the maximum and minimum values of each resulting
        connected component)
    :param dissim_type: int value indicating type of dissimilarity measure. \
                        0 (default) for absolute difference. \
                        1 for dissimilarity measure countering the chaining \
                          effect as described in :cite:`soille2011ismm`
    :return: labeled Jim object
    """
    diss = _pj.ngbops.getDissim(jim, dissim_type)

    if isinstance(jim, _pj.Jim):
        if jim.properties.nrOfBand() == 1:
            return _pj.Jim(jim._jipjim.labelConstrainedCCsDissim(
                diss[0]._jipjim, diss[1]._jipjim,
                global_range, local_range))
        else:
            return _pj.Jim(jim._jipjim.labelConstrainedCCsMultibandDissim(
                diss[0]._jipjim, diss[1]._jipjim,
                global_range, local_range))
    elif isinstance(jim, _pj.JimList):
        return _pj.Jim(jim._jipjimlist.labelConstrainedCCsMultibandDissim(
            diss[0]._jipjim, diss[1]._jipjim, global_range, local_range))
    else:
        raise _pj.exceptions.JimIllegalArgumentError(
            "input must be Jim or JimList object")


def labelConstrainedCCsMi(jim,
                          jim_mi,
                          jim_se,
                          ox: float,
                          oy: float,
                          oz: float,
                          rg: int,
                          rl: int):
    """Label image, in development.

    :param jim: a Jim object to label
    :param jim_mi: a Jim object
    :param jim_se: a Jim object
    :param ox: x coordinate of origin of
    :param oy: y coordinate of origin of
    :param oz: z coordinate of origin of
    :param rg: integer
    :param rl: integer for range parameter lambda l under the strongly
        connected assumption
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelConstrainedCCsMi(
        jim_mi._jipjim, jim_se._jipjim, ox, oy, oz, rg, rl))


def labelConstrainedCCsVariance(jim,
                                ox: float,
                                oy: float,
                                oz: float,
                                rg: int,
                                rl: int,
                                varmax: float,
                                ngb):
    """Label image.

    :param jim: a Jim object
    :param ox: x coordinate of origin of ngb Jim
    :param oy: y coordinate of origin of ngb Jim
    :param oz: z coordinate of origin of ngb Jim
    :param rg: integer for range parameter lambda g
    :param rl: integer for range parameter lambda l
    :param varmax: float for maximum variance of cc
    :param ngb: Jim object for neighbourhood, e.g., create with pj.Jim(graph=4)
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelConstrainedCCsVariance(
        ngb._jipjim, ox, oy, oz, rg, rl, varmax))


def labelErode(jim,
               graph: int = 4):
    """Label to contain the union of the erosions of the iso-intensity CCs.

    :param jim: a Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    return _pj.Jim(jim._jipjim.labelErode(graph))


def labelFlatZones(jim,
                   ngb):
    """Label each image flat zone with a unique label using graph-connectivity.

    :param jim: a Jim object holding a grey level image
    :param ngb: Jim object for neighbourhood, e.g., create with pj.Jim(graph=4)
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelFlatZones(ngb._jipjim, 1, 1, 0))


def labelFlatZonesSeeded(jim,
                         jim_ngb,
                         jim_seeds,
                         ox: float,
                         oy: float,
                         oz: float):
    """Label the flat regions using the ngb Jim and its origin.

    :param jim: a Jim object
    :param jim_ngb: a Jim objects for neighbourhood
    :param jim_seeds: A Jim object for seeds (must be of type UCHAR)
    :param ox: x coordinate of origin of jim_ngb
    :param oy: y coordinate of origin of jim_ngb
    :param oz: z coordinate of origin of jim_ngb
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelFlatZonesSeeded(
        jim_ngb._jipjim, jim_seeds._jipjim, ox, oy, oz))


def labelPixels(jim):
    """Label each non-zero pixel of im with a unique label.

    Labels unless label overflows occurs.

    :param jim: a Jim object
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelPix())


def labelsOuterContour(jim_label,
                       graph: int):
    """Get outer contour of labels.

    :param jim_label: a Jim object
    :param graph: integer for connectivity
    :return: a Jim object containing the outher edge
    """
    return _pj.Jim(jim_label._jipjim.labelsGetOuterContour(
        graph))


def labelsOuterEdge(jim_label,
                    graph: int):
    """Get outer edge of labels.

    :param jim_label: a Jim object containing labels
    :param graph: integer for connectivity
    :return: a Jim object containing the outher edge
    """
    return _pj.Jim(jim_label._jipjim.labelsGetOuterEdge(
        graph))


def labelsOuterEdgeLut(jim_label,
                       jim_edge_label):
    """Get outer edge lut of labels.

    :param jim_label: a Jim object containing labels
    :param jim_edge_label: a Jim object containing the labels outer edge
    :return: a Jim object containing the outher edge lut
    """
    return _pj.Jim(jim_label._jipjim.labelsGetOuterEdgeLut(
        jim_edge_label._jipjim))


def labelsSet(label_jim,
              ival_jim,
              indic: int):
    """Set labels to regions.

    :param label_jim: Jim object with labels
    :param ival_jim: a Jim object
    :param indic: an integer
    :return: a Jim object with set region labels
    """
    return _pj.Jim(label_jim._jipjim.labelsSet(ival_jim, indic))


def labelsSetArea(jim):
    """Set area to regions based on Tessel surface.

    :param jim: Jim object with labels (data must be of unsigned integer type)
    :return: a Jim object with area set to region labels

    Example:

    Label horizontally and vertically connected pixels in 2-D image with unique labels and set area::

        jim = pj.Jim(ncol = 5, nrow = 5, otype = 'GDT_UInt16', uniform = [0,2])
        labeled = pj.ccops.label(jim, pj.Jim(graph=4))
        area = pj.ccops.labelsSetArea(labeled)

    binary input image:

    .. image:: figures/ccops_label2d_input.png
       :width: 65 %

    labeled output image with unique integer values for connected components (here represented by different gray scales):

    .. image:: figures/ccops_label2d_output.png
       :width: 65 %

    output image with area for connected components:

    .. image:: figures/ccops_labelsSetArea2d_output.png
       :width: 80 %


    """
    return _pj.Jim(jim._jipjim.labelsSet())


def labelsSetGraph(label_jim,
                   ival_jim,
                   indic: int,
                   graph: int):
    """Set labels to regions.

    :param label_jim: Jim object with labels
    :param ival_jim: a Jim object
    :param indic: an integer
    :param graph: an integer holding for the graph connectivity
    :return: a Jim object with set region labels
    """
    return _pj.Jim(label_jim._jipjim.labelsSetGraph(ival_jim, indic, graph))


def labelStronglyCCs(jim,
                     local_range: int,
                     ngb):
    """Label each strongly alpha-connected component.

    Label with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim: a Jim or a Jim list of grey level images having all the same
        definition domain and data type.
    :param local_range: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param ngb: Jim object for neighbourhood, e.g., create with pj.Jim(graph=4)
    :return: labeled Jim object
    """
    if isinstance(jim, _pj.Jim):
        if jim.properties.nrOfBand() == 1:
            return _pj.Jim(jim._jipjim.labelStronglyCCs(
                ngb._jipjim, 1, 1, 0, local_range))
        else:
            return _pj.Jim(jim._jipjim.labelStronglyCCsMultiband(
                ngb._jipjim, 1, 1, 0, local_range))
    elif isinstance(jim, _pj.JimList):
        return _pj.Jim(jim._jipjimlist.labelStronglyCCsMultiband(
            ngb._jipjim, 1, 1, 0, local_range))
    else:
        raise _pj.exceptions.JimIllegalArgumentError(
            "Input must be Jim or JimList object")


def morphoFillHoles(jim_object,
                    graph: int,
                    border_flag: int = 1):
    """Remove the not border-connected regional minima of the image.

    Uses graph connectivity (originally proposed for removing pits in
    digital elevation models. See :cite:`soille-ansoult90` and
    :cite:`soille-gratin94` for a fast implementation)

    :param jim_object: input Jim object
    :param border_flag:
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    maxval = jim_object.stats.getStats('max')['max']
    marker = _pj.pixops.setData(jim_object, maxval)
    if border_flag:
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
    else:
        marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
    marker.pixops.supremum(jim_object)
    marker.ccops.morphoGeodesicReconstructionByErosion(jim_object, graph,
                                                       border_flag)
    return marker


def morphoGeodesicReconstructionByDilation(jim_object_mark,
                                           jim_object_mask,
                                           graph: int,
                                           border_flag: int = 1):
    """Compute the morphological reconstruction by dilation of mask image.

    Mask image is from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param border_flag: integer (border values are set to PIX MIN in BOTH
        images if flag equals 0, otherwise the image are internally processed
        by adding a border which is then removed at the end of the
        processing). Default value is 1.
    :return: jim_object_mark containing the result of the morphological
        reconstruction by dilation
    """
    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByDilation(
        jim_object_mask._jipjim, graph, border_flag))


def morphoGeodesicReconstructionByErosion(jim_object_mark,
                                          jim_object_mask,
                                          graph: int,
                                          border_flag: int = 1):
    """Compute the morphological reconstruction by erosion of mask image.

    Mask image is from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param border_flag: integer (border values are set to PIX MIN in BOTH
        images if flag equals 0, otherwise the image are internally processed
        by adding a border which is then removed at the end of the processing).
        Default value is 1.
    :return: jim_object_mark containing the result of the morphological
        reconstruction by erosion
    """
    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByErosion(
        jim_object_mask._jipjim, graph, border_flag))


def morphoRemoveBorder(jim_object,
                       graph: int):
    """Remove the border-connected components of an image.

    Uses graph connectivity.

    :param jim_object: input Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    minmax = jim_object.stats.getStats(function=['min', 'max'])
    minval = minmax['min']
    maxval = minmax['max']
    marker = _pj.pixops.setData(jim_object, minval)
    marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
    marker.pixops.infimum(jim_object)
    marker.ccops.morphoGeodesicReconstructionByDilation(jim_object, graph)
    marker.pixops.simpleArithOp(jim_object, 17)  # 17 for SUBSWAP_op_ovf
    return marker


def partitionSimilarity(jim1,
                        jim2,
                        graph: int):
    """Create a list of 4 1-D images.

    They contain the following information: correspondence table between the
    labels of im1 and im2, similarity measure between these labels,
    correspondence table between the labels of im2 and im1, similarity measure
    between these labels. Create Jim holding the tree.

    :param jim1: first image
    :param jim2: second image
    :param graph: an integer for connectivity
    :return: a list of images
    """
    assert jim1.properties.getDataType() == jim2.properties.getDataType(),\
        'Error: inputs must have same data type'

    return _pj.JimList(jim1._jipjim.partitionSimilarity(jim2._jipjim, graph))


def propagate(label_jim,
              dst_jim,
              imap_jim,
              nc: int,
              graph: int):
    """Perform propagation.

    :param label_jim: a Jim object with labels
    :param dst_jim: a Jim object
    :param imap_jim: a Jim object
    :param nc: an integer
    :param graph: an integer for connectivity
    :return: propagated Jim object
    """
    return _pj.Jim(label_jim._jipjim.propagate(dst_jim, imap_jim, nc, graph))


def relabel(jim_label1,
            jim_label2,
            jim_area):
    """Perform propagation.

    :param jim_label1: a Jim object
    :param jim_label2: a Jim object
    :param jim_area: a Jim object
    :return: relabeled Jim object
    """
    return _pj.Jim(jim_label1._jipjim.labelRelabel(jim_label2, jim_area))


def seededRegionGrowing(jim,
                        seeds,
                        ngb):
    """Calculate the seeded region growing.

    Seeded region growing :cite:`adams-bischof94` including adaptations
    presented in :cite:`mehnert-jackway97`.

    A seeded region algorithm whereby labelled seeds (seeds) are grown in
    a multi-channel image using graph-connectivity. The growth is driven by
    the spectral distances (L2 norm) are calculated between pixels along
    the external boundary of the already grown regions and the corresponding
    pixels along the internal boundary of the seeds. Both jim and seeds are
    modified by this function. The image of seeds is modified by expanding
    the corresponding initial values of the seeds.

    :param jim: a Jim or Jim list of grey level images having all the same
        definition domain and data type.
    :param seeds: a Jim image for labelled seeds (UINT32 type)
    :param ngb: Jim object for neighbourhood, e.g., create with pj.Jim(graph=4)
    :return: labeled Jim object
    """
    if isinstance(jim, _pj.Jim):
        jim_object_list = _pj.JimList([jim])
    else:
        jim_object_list = jim

    return _pj.Jim(
        jim_object_list._jipjimlist.segmentationSeededRegionGrowingMultiband(
            seeds._jipjim, ngb._jipjim, 1, 1, 0))


def segmentBinaryPatterns(jim_object,
                          graph: int = 8, size: float = 1.0, transition: int = 1, internal: int = 1):
    """Morphological segmentation of binary patterns.

    Described in :cite:`soille-vogt2009`, see also :cite:`soille-vogt2022foss4g`

    :param jim_object: input Jim object with pixels of type unsigned char and with foreground pixels set to 2, background pixels set to 1, and no data pixels set to 0.
    :param size: a float value >=1.0 indicating the width of the edges;
    :param graph: an integer  (4 or 8) holding for the graph connectivity (default 8)
    :param transition: a Boolean integer value (0 or 1) indicating how transitions should be processed (default 1)
    :param internal: a Boolean integer value (0 or 1) indicating how embedded components should be processed with 0 for no special treatment or 1 for assigning special values to pixels belonging to embedded components like core components fully surrounded by a larger core component (default 1)
    :return: a new Jim object holding the morphological segmentation of the input
        binary image
    """
    return _pj.Jim(jim_object._jipjim.segmentBinaryPatterns(size, graph, transition, internal))


def segmentImageMultiband(jimlist,
                          local_range: int,
                          region_size: int,
                          contrast: int = 0,
                          version: int = 0,
                          graph: int = 4,
                          filename_prefix: str = ""):
    """Do multiband image segmentation.

    Based on the method described in :cite:`brunner-soille2007`

    The contrast threshold value is used for merging the regions with similar
    contrast as follows: < 0 (do not perform region merge), 0 (determine best
    contrast value automatically), and > 0 (use this value as threshold value).
    Authorised version values are: 0 (compare to whole region), 1 (compare to
    original seeds), and 2 (compare to pixel neighbours).  If the optional
    string filename_prefix is given, data files to use with gnuplot are
    stored in filename_prefix_xxx.dat, otherwise data files are not
    generated (default).

    :param jimlist: a Jim list of grey level images having all the same
        definition domain and data type.
    :param local_range: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param region_size: integer value for minimum size of iso-intensity region
        in output image (must be >= 2 pixels)
    :param contrast: (default is 0)
    :param version:  (default is 0)
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :param filename_prefix: Prefix for filenames
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if isinstance(jimlist, _pj.Jim):
        return _pj.Jim(jimlist._jipjim.segmentImageMultiband(
            graph, local_range, region_size, contrast, version,
            filename_prefix))
    elif isinstance(jimlist, _pj.JimList):
        return _pj.Jim(jimlist._jipjimlist.segmentImageMultiband(
            graph, local_range, region_size, contrast, version,
            filename_prefix))
    else:
        raise _pj.exceptions.JimIllegalArgumentError(
            "Input must be Jim or JimList object")


def vertexConnectedness(jim_object,
                        alpha: int,
                        graph: int = 8,
                        deg: int = None):
    """Label Jim based on vertex connectedness.

    :param jim_object: input Jim object
    :param alpha: value
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param deg: integer
    :return: vertex connectedness labeled Jim
    """
    return _pj.Jim(jim_object._jipjim.vertexDegreeAlpha(alpha, graph, deg))


def vertexDegreeAlpha(jim_object,
                      alpha: int,
                      graph: int = 8):
    """Label Jim based on vertex degree alpha.

    :param jim_object: input Jim object
    :param alpha: value
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: vertex degree alpha labeled Jim
    """
    return _pj.Jim(jim_object._jipjim.vertexDegreeAlpha(alpha, graph))


def vertexSeparation(jim_object,
                     graph: int = 8,
                     type: int = None):
    """Label Jim based on vertex separation.

    :param jim_object: input Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param type: integer
    :return: vertex separated labeled Jim
    """
    return _pj.Jim(jim_object._jipjim.vertexSeparation(graph, type))


def watershed(jim_object,
              graph: int = 8):
    """Watershed segmentation based on immersion simulation.

    Described in :cite:`soille-vincent90`, see also :cite:`vincent-soille91`

    :param jim_object: input Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    return _pj.Jim(jim_object._jipjim.segmentationWatershed(graph))


class _CCOps(_pj.modules.JimModuleBase):
    """Define all CCOps methods."""

    def alphaTree(self,
                  attr0cc,
                  type: int):
        """Compute the alpha tree.

        Modifies the instance on which the method was called.

        :param attr0cc: a Jim object
        :param type: an integer
        """
        self._jim_object._set(self._jim_object._jipjim.alphaTree(
            attr0cc._jipjim, type))

    def alphaTreeNextLevel(self,
                           label_jim,
                           alpha: int):
        """Convert an alpha tree to connected components.

        Modifies the instance on which the method was called.

        :param label_jim: a Jim object
        :param alpha: integer
        """
        self._jim_object._set(self._jim_object._jipjim.alphaTreeNextLevel(
            label_jim._jipjim, alpha))

    def alphaTreePersistenceLut(self):
        """Compute the alpha tree persistence lut.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(
            self._jim_object._jipjim.alphaTreeGetPersistenceLut())

    def alphaTreeToCCs(self,
                       label_jim,
                       lut: bool,
                       rule: int):
        """Convert an alpha tree to connected components.

        Modifies the instance on which the method was called.

        :param label_jim: a Jim object
        :param lut: flag
        :param rule: integer
        """
        self._jim_object._set(self._jim_object._jipjim.alphaTreeToCCs(
            label_jim._jipjim, lut, rule))

    def convertRgbToHsx(self,
                        x_type: str):
        """Convert RGB to HSX.

        Returns the hue, saturation, and value, lightness, or intensity
        channels of an input RGB colour image. The hue component is identical
        for all 3 models. The luminance is equal to max(R,G,B) for HSV,
        (max-min)/2 for HSL and (R+G+B)/3 for HSI. See specific formulae for
        the saturation at `<http://en.wikipedia.org/wiki/HSL_and_HSV>`__.

        Modifies the instance on which the method was called.

        :param x_type: string with key ('V' (default) for Value,
            'L' for Lightness, and 'I' for Intensity)
        """
        assert self._jim_object.properties.nrOfBand() == 3, \
            'Error: input jim must be multi-band image with three bands ' \
            '(r, g, b)'
        jimr = _pj.geometry.cropBand(self._jim_object, 0)
        jimg = _pj.geometry.cropBand(self._jim_object, 1)
        jimb = _pj.geometry.cropBand(self._jim_object, 2)

        jimlist = _pj.JimList(jimr._jipjim.convertRgbToHsx(jimg, jimb,
                                                           x_type))
        self._jim_object._set(jimlist.geometry.stackBand())

    def convertHsiToRgb(self):
        """Convert HSI to RGB.

        Takes the hue, saturation, and intensity channels of a colour image
        and returns an image node containing a colour RGB image.

        Modifies the instance on which the method was called.
        """
        assert self._jim_object.properties.nrOfBand() == 3, \
            'Input jim must be multi-band image with three bands ' \
            '(h, s, i)'
        jimh = _pj.geometry.cropBand(self._jim_object, 0)
        jims = _pj.geometry.cropBand(self._jim_object, 1)
        jimi = _pj.geometry.cropBand(self._jim_object, 2)

        self._jim_object._set(jimh._jipjim.convertHsiToRgb(jims, jimi))

    def convertHlsToRgb(self):
        """Convert HLS to RGB.

        Takes the hue, lightness, and saturation channels of a colour image
        and returns an image node containing a colour RGB image.

        Modifies the instance on which the method was called.
        """
        assert self._jim_object.properties.nrOfBand() == 3, \
            'Input jim must be multi-band image with three bands ' \
            '(h, s, i)'
        jimh = _pj.geometry.cropBand(self._jim_object, 0)
        jiml = _pj.geometry.cropBand(self._jim_object, 1)
        jims = _pj.geometry.cropBand(self._jim_object, 2)

        self._jim_object._set(jimh._jipjim.convertHlsToRgb(jiml, jims))

    def dbscan(self,
               eps: float,
               min_pts: int):
        """Compute dbscan.

        This method should be called on an object which contains
        a dissimilarity matrix.

        Modifies the instance on which the method was called.

        :param eps: float
        :param min_pts: integer
        """
        self._jim_object._set(self._jim_object._jipjim.dbscan(eps, min_pts))

    def dissim(self,
               jim_mask,
               nc: int,
               type: int):
        """Compute the dissimilarity matrix.

        Modifies the instance on which the method was called.

        :param jim_mask: a Jim object
        :param nc: integer
        :param type: integer
        """
        self._jim_object._set(self._jim_object._jipjim.dissimilarityMatrix(
            jim_mask._jipjim, nc, type))

    def distance2d4(self,
                    band: int = 0):
        """Compute the 2-dimensional 4-connected distance function of a Jim.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_distance2d4(band)

    def distance2dChamfer(self,
                          type: int,
                          band: int = 0):
        """Compute the chamfer distance function of Jim.

        Note that the type of Jim must be large enough to hold the largest
        distance otherwise overflow may occur and generate artefacts in the
        resulting distance function.

        Modifies the instance on which the method was called.

        :param type: Integer for type of chamfer distance {1, 11, 34, 57, 5711}
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_distance2dChamfer(type, band)

    def distance2dEuclideanConstrained(self,
                                       mask,
                                       band: int = 0):
        """Compute the Euclidean geodesic distance from the marker set.

        Defined by the image marker and within the geodesic mask defined by
        the image mask. The algorithm is described in (:cite:`soille91`).

        Modifies the instance on which the method was called.

        :param mask: binary UCHAR image with geodesic mask pixels set to 1
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._set(
            self._jim_object._jipjim.distance2dEuclideanConstrained(
                mask._jipjim, band))

    def distance2dEuclideanSquared(self,
                                   band: int = 0):
        """Compute the squared Euclidean distance transform in 2-D.

        The instance must be a binary image. Multi-plane images will be processed
        plane by plane as 2-D images. The original algorithm was proposed
        by :cite:`saito-toriwaki94` and then optimised independently by
        :cite:`hirata96` and :cite:`meijster-roerdink-hesselink2000`.
        Based on the Euclidean distance transform. Note that a temporary buffer
        of type UINT16 is used for sums along/lines and columns so that
        uncontrolled results will occur if an object shows more than 2^16 /2
        foreground pixels along a given line or column.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)

        Example:

        Create a binary image that is 1 at the center and border and 0 elsewhere::

            jim = pj.Jim(nrow=30, ncol=30, otype='Byte')

            jim[15-1:15+1,15-1:15+1] = 1
            jim[15-1:15+1,15-1:15+1] = 1
            jim[0,:] = 1
            jim[-1,:] = 1
            jim[:,0] = 1
            jim[:,-1] = 1

        binary input image:

        .. image:: figures/distance_input.png
           :width: 65 %

        To Euclidean distance is calculated from the background pixels (value 0) to the foreground pixels (value 1). If we want to calculate the distance from the foreground to the background (e.g., distance from binary cloud mask to nearest clear pixel), the image must be negated first::

            jim = jim != 1
            jim.ccops.distance2dEuclideanSquared()

        Squared Euclidean distance image:

        .. image:: figures/distance_output.png
           :width: 65 %

        The cloud mask can then be extended by thresholding the distance image::


            jim = jim < 4

        .. image:: figures/distance_threshold.png
           :width: 65 %
        """

        jim = _pj.Jim(_pj.geometry.cropPlane(
            self._jim_object, 0))
        jim._set(
            jim._jipjim.distance2dEuclideanSquared(band))

        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.Jim(_pj.geometry.cropPlane(
                self._jim_object, iplane))
            jimplane._set(
                jimplane._jipjim.distance2dEuclideanSquared(band))
            jim.geometry.stackPlane(jimplane)
        self._jim_object._set(jim._jipjim)

        # self._jim_object._set(
        #     self._jim_object._jipjim.distance2dEuclideanSquared(band))

    def distanceGeodesic(self,
                         marker,
                         graph: int,
                         band: int = 0):
        """Compute geodesic distance function from marker within mask.

        Using graph connectivity.

        Modifies the instance on which the method was called.

        :param marker: an image node for marker image
        :param graph: integer for connectivity
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_jldistanceGeodesic(
            marker._jipjim, graph, band)

    def distanceInfluenceZones2dEuclidean(self,
                                          band: int = 0):
        """Output error-free influence zones of the labelled regions in Jim.

        Jim must be a 2-D image, its last bit being reserved (e.g., for a UCHAR
        images, authorized label values range from 1 to 127).

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._set(
            self._jim_object._jipjim.distanceInfluenceZones2dEuclidean(band))

    def flatZonesSeeded(self,
                        jim2,
                        jim3,
                        ox: float,
                        oy: float,
                        oz: float):
        """Compute seeded flat zones.

        Modifies the instance on which the method was called.

        :param jim2: a Jim object
        :param jim3: a Jim object
        :param ox: x coordinate of origin
        :param oy: y coordinate of origin
        :param oz: z coordinate of origin
        """
        self._jim_object._jipjim.d_flatZonesSeeded(jim2, jim3, ox, oy, oz)

    def label(self,
              ngb):
        """Label each connected component (non zero values) with a unique integer value using graph-connectivity.

        The labels start from value 2.
        Modifies the instance on which the method was called.

        :param ngb: Jim object for neighbourhood, e.g., create with:

            - pj.Jim(graph=4): horizontally and vertically connected in 2-D
            - pj.Jim(graph=6): horizontally and vertically connected in 3-D
            - pj.Jim(graph=8): horizontally, vertically, and diagonally connected in 2-D
            - pj.Jim(graph=26): horizontally, vertically, and diagonally connected in 3-D

        Example:

        Label horizontally and vertically connected pixels in 2-D image with unique labels::

            jim = pj.Jim(ncol = 5, nrow = 5, otype = 'GDT_UInt16', uniform = [0,2])
            jim.ccops.label(pj.Jim(graph=4))

        binary input image:

        .. image:: figures/ccops_label2d_input.png
           :width: 65 %

        labeled output image with unique integer values for connected components (here represented by different gray scales):

        .. image:: figures/ccops_label2d_output.png
           :width: 65 %

        Label horizontally and vertically connected pixels in 3-D image with unique labels::

            jim = pj.Jim(ncol = 3, nrow = 3, nplane = 3, otype = 'GDT_Byte', uniform = [0,2])
            jim.ccops.label(pj.Jim(graph=6))

        binary input 3-D image:

        .. image:: figures/ccops_label3d_input.png
           :width: 65 %


        labeled output image with unique integer values for connected components (here represented by different gray scales):

        .. image:: figures/ccops_label3d_output.png
           :width: 65 %
        """

        self._jim_object.geometry.imageFrameAdd(1, 1, 1, 1, 0, 0, 0)
        self._jim_object._jipjim.d_labelBinary(ngb._jipjim, 1, 1, 0)
        self._jim_object.geometry.imageFrameSubtract(1, 1, 1, 1, 0, 0)

    def labelConstrainedCCs(self,
                            local_range: int,
                            global_range: int,
                            ngb):
        """Label each alpha-omega connected component.

        Label with a unique label using graph-connectivity
        :cite:`soille2008pami`

        :param local_range: integer value indicating maximum absolute local
            difference between 2 adjacent pixels
        :param global_range: integer value indicating maximum global difference
            (difference between the maximum and minimum values of each
            resulting connected component)
        :param ngb: Jim object for neighbourhood, e.g., create with
            pj.Jim(graph=4)
        """
        if self._jim_object.properties.nrOfBand() == 1:
            self._jim_object._set(
                self._jim_object._jipjim.labelConstrainedCCs(
                    ngb._jipjim, 1, 1, 0, global_range, local_range))
        else:
            self._jim_object._set(
                self._jim_object._jipjim.labelConstrainedCCsMultiband(
                    ngb._jipjim, 1, 1, 0, global_range, local_range))

    def labelConstrainedCCsAttr(self,
                                graph: int,
                                rg: int,
                                rl: int):
        """Label image, in development.

        :param graph: integer for graph connectivity
        :param rg: integer
        :param rl: integer for range parameter lambda l under the strongly
            connected assumption
        :return: labeled Jim object
        """
        self._jim_object._set(
            self._jim_object._jipjim.labelConstrainedCCsAttr(
                graph, rg, rl))

    def labelConstrainedCCsCi(self,
                              ngb,
                              ox: float,
                              oy: float,
                              oz: float,
                              rl: int):
        """Label image, in development.

        Modifies the instance on which the method was called.

        :param ngb: a Jim object for neighbourhood
        :param ox: x coordinate of origin of ngb Jim
        :param oy: y coordinate of origin of ngb Jim
        :param oz: z coordinate of origin of ngb Jim
        :param rl: integer for range parameter lambda l
            under the strongly connected assumption
        """
        self._jim_object._set(
            self._jim_object._jipjim.labelConstrainedCCsCi(
                ngb._jipjim, ox, oy, oz, rl))

    def labelConstrainedCCsMi(self,
                              jim_mi,
                              jim_se,
                              ox: float,
                              oy: float,
                              oz: float,
                              rg: int,
                              rl: int):
        """Label image, in development.

        Modifies the instance on which the method was called.

        :param jim_mi: a Jim object
        :param jim_se: a Jim object
        :param ox: x coordinate of origin of ngb Jim
        :param oy: y coordinate of origin of ngb Jim
        :param oz: z coordinate of origin of ngb Jim
        :param rg: integer
        :param rl: integer for range parameter lambda l
            under the strongly connected assumption
        """
        self._jim_object._set(
            self._jim_object._jipjim.labelConstrainedCCsMi(
                jim_mi._jipjim, jim_se._jipjim, ox, oy, oz, rg, rl))

    def labelConstrainedCCsVariance(self,
                                    ox: float,
                                    oy: float,
                                    oz: float,
                                    rg: int,
                                    rl: int,
                                    varmax: float,
                                    ngb):
        """Label image.

        Modifies the instance on which the method was called.

        :param ox: x coordinate of origin of ngb Jim
        :param oy: y coordinate of origin of ngb Jim
        :param oz: z coordinate of origin of ngb Jim
        :param rg: integer for range parameter lambda g
        :param rl: integer for range parameter lambda l
        :param varmax: float for maximum variance of cc
        :param ngb: Jim object for neighbourhood, e.g., create with
            pj.Jim(graph=4)
        """
        self._jim_object._set(
            self._jim_object._jipjim.labelConstrainedCCsVariance(
                ngb._jipjim, ox, oy, oz, rg, rl, varmax))

    def labelErode(self,
                   graph: int = 4):
        """Label to contain the union of the erosions of the iso-intensity CCs.

        Modifies the instance on which the method was called.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images, default is 4)
        """
        _pj._check_graph(graph, [4, 8])

        self._jim_object._set(self._jim_object._jipjim.labelErode(graph))

    def labelFlatZonesSeeded(self,
                             jim_ngb,
                             jim_seeds,
                             ox: float,
                             oy: float,
                             oz: float):
        """Label the flat regions using the ngb Jim and its origin.

        :param jim_ngb: a Jim objects for neighbourhood
        :param jim_seeds: A Jim object for seeds (must be of type UCHAR)
        :param ox: x coordinate of origin of jim_ngb
        :param oy: y coordinate of origin of jim_ngb
        :param oz: z coordinate of origin of jim_ngb
        """
        self._jim_object._jipjim.d_labelFlatZonesSeeded(
            jim_ngb._jipjim, jim_seeds._jipjim, ox, oy, oz)

    def labelPixels(self):
        """Label each non-zero pixel of im with a unique label.

        Labels unless label overflows occurs.

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_labelPix()

    def labelsOuterContour(self,
                           graph: int):
        """Compute the outer contour of labels.

        This method should be called on a Jim object which contains labels.

        Modifies the instance on which the method was called.

        :param graph: integer for connectivity
        """
        self._jim_object._set(self._jim_object._jipjim.labelsGetOuterContour(
            graph))

    def labelsOuterEdge(self,
                        graph: int):
        """Compute the outer edge of labels.

        This method should be called on a Jim object which contains labels.

        Modifies the instance on which the method was called.

        :param graph: integer for connectivity
        """
        self._jim_object._set(self._jim_object._jipjim.labelsGetOuterEdge(
            graph))

    def labelsOuterEdgeLut(self,
                           jim_edge_label):
        """Compute the outer edge lut of labels.

        This method should be called on a Jim object which contains labels.

        Modifies the instance on which the method was called.

        :param jim_edge_label: a Jim object containing the labels outer edge
        """
        self._jim_object._set(self._jim_object._jipjim.labelsGetOuterEdgeLut(
            jim_edge_label._jipjim))

    def labelsSet(self,
                  ival_jim,
                  indic: int):
        """Set labels to regions.

        Modifies the instance on which the method was called.

        :param ival_jim: a Jim object
        :param indic: an integer
        :return: a Jim object with set region labels
        """
        self._jim_object._jipjim.d_labelsSet(ival_jim, indic)

    def labelsSetArea(self):
        """Set area to regions based on Tessel surface.

        Example:

        Label horizontally and vertically connected pixels in 2-D image with unique labels and set area::

            jim = pj.Jim(ncol = 5, nrow = 5, otype = 'GDT_UInt16', uniform = [0,2])
            jim.ccops.label(pj.Jim(graph=4))
            jim.ccops.labelsSetArea()

        binary input image:

        .. image:: figures/ccops_label2d_input.png
           :width: 65 %

        labeled output image with unique integer values for connected components (here represented by different gray scales):

        .. image:: figures/ccops_label2d_output.png
           :width: 65 %

        output image with area for connected components:

        .. image:: figures/ccops_labelsSetArea2d_output.png
           :width: 80 %

        """

        self._jim_object._jipjim.d_labelsSetArea()

    def labelsSetGraph(self,
                       ival_jim,
                       indic: int,
                       graph: int):
        """Set labels to regions.

        Modifies the instance on which the method was called.

        :param ival_jim: a Jim object
        :param indic: an integer
        :param graph: an integer holding for the graph connectivity
        :return: a Jim object with set region labels
        """
        self._jim_object._jipjim.d_labelsSetGraph(ival_jim, indic, graph)

        #todo: not working (self is not modified)
    def morphoFillHoles(self,
                        graph: int,
                        border_flag: int = 1):
        """Remove the image-border-connected regional minima.

        Uses graph connectivity.

        Modifies the instance on which the method was called.

        :param border_flag:
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        """
        # jim_object=_pj.Jim(self._jim_object._jipjim)
        maxval = self._jim_object.stats.getStats('max')['max']
        marker = _pj.pixops.setData(self._jim_object, maxval)
        if border_flag:
            marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
        else:
            marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
        marker.pixops.supremum(self._jim_object)

        marker.ccops.morphoGeodesicReconstructionByErosion(self._jim_object,
                                                           graph, border_flag)
        self._jim_object._set(marker._jipjim)

    def morphoGeodesicReconstructionByDilation(self,
                                               jim_object_mask,
                                               graph: int,
                                               flag: int = 1):
        """Compute the morphological reconstruction by dilation.

        Dilation of the current object is from mark image using graph
        connectivity.

        Modifies the instance on which the method was called.

        :param jim_object_mask: a Jim object for mask image
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :param flag: integer (border values are set to PIX MIN in BOTH images
            if flag equals 0, otherwise the image are internally processed by
            adding a border which is then removed at the end of
            the processing). Default value is 1.
        """
        self._jim_object._jipjim.d_geodesicReconstructionByDilation(
            jim_object_mask._jipjim, graph, flag)

    def morphoGeodesicReconstructionByErosion(self,
                                              jim_object_mask,
                                              graph: int,
                                              flag: int = 1):
        """Compute the morphological reconstruction by erosion.

        Erosion of the current object is from mark image using graph
        connectivity.

        Modifies the instance on which the method was called.

        :param jim_object_mask: a Jim object for mask image
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :param flag: integer (border values are set to PIX MIN in BOTH images
            if flag equals 0, otherwise the image are internally processed by
            adding a border which is then removed at the end of the
            processing). Default value is 1.
        """
        self._jim_object._jipjim.d_geodesicReconstructionByErosion(
            jim_object_mask._jipjim, graph, flag)

    def morphoRemoveBorder(self,
                           graph: int):
        """Remove the border-connected components using graph connectivity.

        Modifies the instance on which the method was called.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        """
        jim_object = self._jim_object._jipjim
        minmax = jim_object.stats.getStats(function=['min', 'max'])
        minval = minmax['min']
        maxval = minmax['max']
        marker = _pj.pixops.setData(jim_object, minval)
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
        marker.pixops.infimum(jim_object)
        marker.ccops.morphoGeodesicReconstructionByDilation(jim_object, graph)
        marker.pixops.simpleArithOp(jim_object, 17)  # 17 for SUBSWAP_op_ovf
        self._jim_object._set(marker)

    def propagate(self,
                  dst_jim,
                  imap_jim,
                  nc: int,
                  graph: int):
        """Perform propagation.

        :param dst_jim: a Jim object
        :param imap_jim: a Jim object
        :param nc: an integer
        :param graph: an integer for connectivity
        """
        self._jim_object._jipjim.d_propagate(dst_jim, imap_jim, nc, graph)

    def relabel(self,
                jim_label2,
                jim_area):
        """Perform propagation.

        :param jim_label2: a Jim object
        :param jim_area: a Jim object
        """
        self._jim_object._jipjim.d_labelRelabel(jim_label2, jim_area)

    def vertexConnectedness(self,
                            alpha: float,
                            graph: int = 8,
                            deg: int = None):
        """Label Jim based on vertex connectedness.

        :param alpha: value
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :param deg: integer
        :return: vertex connectedness labeled Jim
        """
        self._jim_object._set(self._jim_object._jipjim.vertexConnectedness(
            alpha, graph, deg))

    def vertexDegreeAlpha(self,
                          alpha: float,
                          graph: int = 8):
        """Label Jim based on vertex degree alpha.

        :param alpha: value
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :return: vertex degree alpha labeled Jim
        """
        self._jim_object._set(self._jim_object._jipjim.vertexDegreeAlpha(
            alpha, graph))

    def vertexSeparation(self,
                         graph: int = 8,
                         type: int = None):
        """Label Jim based on vertex degree alpha.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :param type: integer
        :return: vertex separated labeled Jim
        """
        self._jim_object._set(self._jim_object._jipjim.vertexSeparation(
            graph, type))


class _CCOpsList(_pj.modules.JimListModuleBase):
    """Define all CCOps methods for JimLists."""

    def labelConstrainedCCs(self,
                            local_range: int,
                            global_range: int,
                            ngb):
        """Label each alpha-omega connected component.

        Label everyting with a unique label using graph-connectivity
        :cite:`soille2008pami`.

        :param local_range: integer value indicating maximum absolute local
            difference between 2 adjacent pixels
        :param global_range: integer value indicating maximum global difference
            (difference between the maximum and minimum values of each
            resulting connected component)
        :param ngb: Jim object for neighbourhood, e.g., create with
            pj.Jim(graph=4)
        :return: labeled Jim object
        """
        return _pj.Jim(self._jim_list._jipjimlist.labelConstrainedCCsMultiband(
            ngb._jipjim, 1, 1, 0, global_range, local_range))


class _CCOpsVect(_pj.modules.JimVectModuleBase):
    """Define all CCOps methods for JimVects."""

    pass
