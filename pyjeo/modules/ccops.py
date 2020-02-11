"""Module for connected-component operations."""

import pyjeo as _pj
from . import JimModuleBase as _JimModuleBase
from . import JimListModuleBase as _JimListModuleBase
from . import JimVectModuleBase as _JimVectModuleBase


# TODO: Test
def alphaTreeDissim(dissimh, dissimv, alpha):
    """Create Jim holding the tree.

    :param dissimh: Jim for horizontal edge dissimilarities
    :param dissimv: Jim for vertical edge dissimilarities
    :param alpha: integer for dissimilarity threshold value
    :return: the tree stored in an JimList
    """
    return _pj.JimList(dissimh._jipjim.alphaTreeDissimGet(dissimv._jipjim,
                                                          alpha))


def convertRgbToHsx(jim, theType):
    """Convert RGB to HSX.

    Returns the hue, saturation, and value, lightness, or intensity channels
    of an input RGB colour image. The hue component is identical for all 3
    models. The luminance is equal to max(R,G,B) for HSV, (max-min)/2 for HSL
    and (R+G+B)/3 for HSI. See specific formulae for the saturation at
    http://en.wikipedia.org/wiki/HSL_and_HSV.

    :param jim: multi-band Jim with three bands
        representing red, green and blue channels
    :param type: string with key ('V' (default) for Value, 'L' for Lightness,
        and 'I' for Intensity)
    :return: Jim with three bands containing the HSX channels
    """
    assert jim.properties.nrOfBand() == 3, \
        'Error: input jim must be multi-band image with three bands (r, g, b)'
    jimr =_pj.geometry.cropBand(jim, 0)
    jimg =_pj.geometry.cropBand(jim, 1)
    jimb =_pj.geometry.cropBand(jim, 2)

    jimlist =_pj.JimList(jimr._jipjim.convertRgbToHsx(jimg, jimb, theType))
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
    jimh =_pj.geometry.cropBand(jim, 0)
    jims =_pj.geometry.cropBand(jim, 1)
    jimi =_pj.geometry.cropBand(jim, 2)

    return pj.Jim(jimh._jipjim.convertRgbToHsx(jims, jimi))


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
    jimh =_pj.geometry.cropBand(jim, 0)
    jiml =_pj.geometry.cropBand(jim, 1)
    jims =_pj.geometry.cropBand(jim, 2)

    return pj.Jim(jimh._jipjim.convertRgbToHsx(jiml, jims))


def dissimToAlphaCCs(dissimh, dissimv, alpha):
    """Create Jim holding the labelled alpha-connected component.

    :param dissimh: Jim for horizontal edge dissimilarities
    :param dissimv: Jim for vertical edge dissimilarities
    :param alpha: integer for dissimilarity threshold value
    :return: Jim object
    """
    return _pj.Jim(dissimh._jipjim.dissimToAlphaCCs(dissimv._jipjim,
                                                    alpha))


def distance2d4(jim, band=0):
    """Compute the 2-dimensional 4-connected distance function of a Jim.

    :param jim: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distance2d4(band))


def distance2dChamfer(jim, type, band=0):
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


def distance2dEuclideanConstrained(marker, mask, band=0):
    """Compute the Euclidean geodesic distance from the marker set.

    Defined by the image marker and within the geodesic mask defined by
    the image mask. The algorithm is described in (Soille, 1991).

    Modifies the instance on which the method was called.

    :param mask: binary UCHAR image with geodesic mask pixels set to 1
    :param band: List of band indices to crop (index is 0 based)
    """
    return _pj.Jim(marker._jipjim.distance2dEuclideanConstrained(mask._jipjim,
                                                                 band))


def distance2dEuclideanSquared(jim, band=0):
    """Compute the squared Euclidean distance transform of im.

    Jim must be a 2-D binary image. Original algorithm proposed by Saito
    and Toriwaki (1994) and then optimised independently by (Hirata,
    1996) and (Meijster et al., 2000). Based on the
    Euclidean distance transform. Note that a temporary buffer of type
    UINT16 is used for sums along/lines and columns so that uncontrolled
    results will occur if an object shows more than 16 2 /2 foreground
    pixels along a given line or column.

    :param jim: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distance2dEuclideanSquared(band))


def distanceGeodesic(mask, marker, graph, band=0):
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


def distanceInfluenceZones2dEuclidean(jim, band=0):
    """Output error-free influence zones of the labelled regions in Jim.

    Jim must be a 2-D image, its last bit being reserved (e.g., for a UCHAR
    images, authorized label values range from 1 to 127). Algorithm based
    on the Euclidean distance transform.

    :param jim: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim._jipjim.distanceInfluenceZones2dEuclidean(band))


def getRegionalMinima(jim, graph):
    """Compute the regional minima of the input image.

    The pixels belonging to a regional minimum are set to 1, all other pixels
    are set to 0.

    :param jim: a Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object of type unsigned char containing the regional
        minima of the input Jim object
    """
    return _pj.Jim(jim._jipjim.getRegionalMinima(graph))


def labelConstrainedCCs(jim, localRange, globalRange, graph=4):
    """Label each alpha-omega connected component.

    Label with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim: a Jim or Jim list of grey level images having all the same
        definition domain and data type.
    :param localRange: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param globalRange: integer value indicating maximum global difference
        (difference between the maximum and minimum values of each resulting
        connected component)
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    if isinstance(jim, _pj.Jim):
        return _pj.Jim(jim._jipjim.labelConstrainedCCs(
            ngb._jipjim, 1, 1, 0, globalRange, localRange))
    else:
        return _pj.Jim(jim._jipjimlist.labelConstrainedCCsMultiband(
            ngb._jipjim, 1, 1, 0, globalRange, localRange))


def labelConstrainedCCsDissim(jim, localRange, globalRange, dissimType=0):
    """Label each alpha-omega connected components with a unique label.

     Label using graph-connectivity and the dissimilarity measure countering
     the chaining effect as described in :cite:`soille2011ismm`

    :param jim: a Jim or a Jim list of grey level images having all the same
        definition domain and data type.
    :param localRange: integer value indicating maximum absolute local
        difference between 2 adjacent pixels along alpha-connected paths
    :param globalRange: integer value indicating maximum global difference
        (difference between the maximum and minimum values of each resulting
        connected component)
    :param dissimType: int value indicating type of dissimilarity measure. \
                       0 (default) for absolute difference. \
                       1 for dissimilarity measure countering the chaining \
                         effect as described in :cite:`soille2011ismm`
    :return: labeled Jim object
    """
    dissim = _pj.ngbops.getDissim(jim, dissimType)

    if isinstance(jim, _pj.Jim):
        return _pj.Jim(jim._jipjim.labelConstrainedCCsDissim(
            dissim[0]._jipjim, dissim[1]._jipjim, globalRange, localRange))
    else:
        return _pj.Jim(jim._jipjimlist.labelConstrainedCCsMultibandDissim(
            dissim[0]._jipjim, dissim[1]._jipjim, globalRange, localRange))


def labelConstrainedCCsVariance(jim, ox, oy, oz, rg, rl, varmax, graph=4):
    """Label image.

    :param jim: a Jim object
    :param ox: x coordinate of origin of imngb
    :param oy: y coordinate of origin of imngb
    :param oz: z coordinate of origin of imngb
    :param rg: integer for range parameter lambda g
    :param rl: integer for range parameter lambda l
    :param varmax: float for maximum variance of cc
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    return _pj.Jim(jim._jipjim.labelConstrainedCCsVariance(
        ngb._jipjim, ox, oy, oz, rg, rl, varmax))


def labelErode(jim, graph=4):
    """Label to contain the union of the erosions of the iso-intensity CCs.

    :param jim: a Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    return _pj.Jim(jim._jipjim.labelErode(graph))


def labelFlatZonesGraph(jim, graph=4):
    """Label each image flat zone with a unique label using graph-connectivity.

    :param jim: a Jim object holding a grey level image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    return _pj.Jim(jim._jipjim.labelFlatZones(ngb._jipjim, 1, 1, 0))


def labelFlatZonesSeeded(jim, jim_ngb, jim_seeds, ox, oy, oz):
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


def labelGraph(jim, graph=4):
    """Label each non-zero connected component with a unique label.

    Label using graph-connectivity.

    :param jim: a Jim object holding a binary image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    return _pj.Jim(jim._jipjim.labelBinary(ngb._jipjim, 1, 1, 0))


def labelImagePixels(jim):
    """Label each non-zero pixel of im with a unique label.

    Labels unless label overflows occurs.

    :param jim: a Jim object
    :return: labeled Jim object
    """
    return _pj.Jim(jim._jipjim.labelPix())


def labelStronglyCCs(jim, localRange, graph=4):
    """Label each strongly alpha-connected component.

    Label with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim: a Jim or a Jim list of grey level images having all the same
        definition domain and data type.
    :param localRange: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    if isinstance(jim, _pj.Jim):
        return _pj.Jim(jim._jipjim.labelConstrainedCCsCi(
            ngb._jipjim, 1, 1, 0, localRange))
    else:
        return _pj.Jim(jim._jipjimlist.labelStronglyCCsMultiband(
            ngb._jipjim, 1, 1, 0, localRange))


def morphoFillHoles(ajim, graph, borderFlag=1):
    """Remove the not border-connected regional minima of the image.

    Uses graph connectivity (originally proposed for removing pits in
    digital elevation models. See :cite:`soille-ansoult90` and
    :cite:`soille-gratin94` for a fast implementation)

    :param ajim: input Jim object
    :param borderFlag:
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    maxval = ajim.stats.getStats('max')['max']
    marker = _pj.pixops.setData(ajim, maxval)
    if borderFlag:
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
    else:
        marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
    marker.pixops.supremum(ajim)
    marker.ccops.morphoGeodesicReconstructionByErosion(ajim, graph, borderFlag)
    return marker


def morphoGeodesicReconstructionByDilation(jim_object_mark, jim_object_mask,
                                           graph, borderFlag=1):
    """Compute the morphological reconstruction by dilation of mask image.

    Mask image is from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param borderFlag: integer (border values are set to PIX MIN in BOTH
        images if flag equals 0, otherwise the image are internally processed
        by adding a border which is then removed at the end of the
        processing). Default value is 1.
    :return: jim_object_mark containing the result of the morphological
        reconstruction by dilation
    """
    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByDilation(
        jim_object_mask._jipjim, graph, borderFlag))


def morphoGeodesicReconstructionByErosion(jim_object_mark, jim_object_mask,
                                          graph, borderFlag=1):
    """Compute the morphological reconstruction by erosion of mask image.

    Mask image is from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :param borderFlag: integer (border values are set to PIX MIN in BOTH
        images if flag equals 0, otherwise the image are internally processed
        by adding a border which is then removed at the end of the processing).
        Default value is 1.
    :return: jim_object_mark containing the result of the morphological
        reconstruction by erosion
    """
    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByErosion(
        jim_object_mask._jipjim, graph, borderFlag))


def morphoRemoveBorder(ajim, graph):
    """Remove the border-connected components of an image.

    Uses graph connectivity.

    :param ajim: input Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    minmax = ajim.stats.getStats(function=['min', 'max'])
    minval = minmax['min']
    maxval = minmax['max']
    marker = _pj.pixops.setData(ajim, minval)
    marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
    marker.pixops.infimum(ajim)
    marker.ccops.morphoGeodesicReconstructionByDilation(ajim, graph)
    marker.pixops.simpleArithOp(ajim, 17)  # 17 for SUBSWAP_op_ovf
    return marker


def partitionSimilarity(jim1, jim2, graph):
    """Create a list of 4 1-D images.

    They contain the following information: correspondence table between the
    labels of im1 and im2, similarity measure between these labels,
    correspondence table between the labels of im2 and im1, similarity measure
    between these labels. Create Jim holding the tree.

    :param jim1: first image
    :param jim2: seconf image
    :param graph: an integer for connectivity
    :return: a list of images
    """
    assert jim1.properties.getDataType() == jim2.properties.getDataType(),\
        'Error: inputs must have same data type'

    return _pj.JimList(jim1._jipjim.partitionSimilarity(jim2._jipjim, graph))


def seededRegionGrowing(jim, seeds, graph=4):
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
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    if isinstance(jim, _pj.Jim):
        jim_object_list = _pj.JimList([jim])
    else:
        jim_object_list = jim

    return _pj.Jim(
        jim_object_list._jipjimlist.segmentationSeededRegionGrowingMultiband(
            seeds._jipjim, ngb._jipjim, 1, 1, 0))


def segmentImageMultiband(jimList, localRange, regionSize, contrast=0,
                          version=0, graph=4, dataFileNamePrefix=""):
    """Do multiband image segmentation.

    Based on the method described in :cite:`brunner-soille2007`

    The contrast threshold value is used for merging the regions with similar
    contrast as follows: < 0 (do not perform region merge), 0 (determine best
    contrast value automatically), and > 0 (use this value as threshold value).
    Authorised version values are: 0 (compare to whole region), 1 (compare to
    original seeds), and 2 (compare to pixel neighbours).  If the optional
    string dataFileNamePrefix is given, data files to use with gnuplot are
    stored in dataFileNamePrefix_xxx.dat, otherwise data files are not
    generated (default).

    :param jimList: a Jim list of grey level images having all the same
        definition domain and data type.
    :param localRange: integer value indicating maximum absolute local
        difference between 2 adjacent pixels
    :param regionSize: integer value for minimum size of iso-intensity region
        in output image (must be >= 2 pixels)
    :param contrast: (default is 0)
    :param version:  (default is 0)
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images, default is 4)
    :param dataFileName:
    :return: labeled Jim object
    """
    _pj._check_graph(graph, [4, 8])

    if graph == 4:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0, 1] = 1
        ngb[1, 0] = 1
        ngb[1, 2] = 1
        ngb[2, 1] = 1
    else:
        ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1, 1] = 0

    return _pj.Jim(jimList._jipjimlist.segmentImageMultiband(
        graph, localRange, regionSize, contrast, version, dataFileNamePrefix))


def watershed(ajim, graph=8):
    """Watershed segmentation based on immersion simulation.

    Described in :cite:`soille-vincent90`, see also :cite:`vincent-soille91`

    :param ajim: input Jim object
    :param graph: an integer holding for the graph connectivity
        (4 or 8 for 2-D images)
    :return: a new Jim object with the connected component of the input
        object removed
    """
    return _pj.Jim(ajim._jipjim.segmentationWatershed(ajim._jipjim, graph))


class _CCOps(_JimModuleBase):
    """Define all CCOps methods."""

    def convertRgbToHsx(self, theType, band):
        """Convert RGB to HSX.

        Returns the hue, saturation, and value, lightness, or intensity
        channels of an input RGB colour image. The hue component is identical
        for all 3 models. The luminance is equal to max(R,G,B) for HSV,
        (max-min)/2 for HSL and (R+G+B)/3 for HSI. See specific formulae for
        the saturation at http://en.wikipedia.org/wiki/HSL_and_HSV.

        Modifies the instance on which the method was called.

        :param theType: string with key ('V' (default) for Value,
            'L' for Lightness, and 'I' for Intensity)
        """
        assert self._jim_object.properties.nrOfBand() == 3, \
            'Error: input jim must be multi-band image with three bands ' \
            '(r, g, b)'
        jimr =_pj.geometry.cropBand(self._jim_object, 0)
        jimg =_pj.geometry.cropBand(self._jim_object, 1)
        jimb =_pj.geometry.cropBand(self._jim_object, 2)

        jimlist =_pj.JimList(jimr._jipjim.convertRgbToHsx(jimg, jimb, theType))
        self._jim_object._set(jimlist.geometry.stackBand())

    def convertHsiToRgb(self):
        """Convert HSI to RGB.

        Takes the hue, saturation, and intensity channels of a colour image
        and returns an image node containing a colour RGB image.

        Modifies the instance on which the method was called.

        :param jim: multi-band Jim with three bands representing hue,
            saturation, and intensity channels
        :return: Jim with three bands containing the RGB channels
        """
        assert jim.properties.nrOfBand() == 3, \
            'Error: input jim must be multi-band image with three bands ' \
            '(h, s, i)'
        jimh =_pj.geometry.cropBand(self._jim_object, 0)
        jims =_pj.geometry.cropBand(self._jim_object, 1)
        jimi =_pj.geometry.cropBand(self._jim_object, 2)

        self._jim_object._set(jimh._jipjim.convertHsiToRgb(jims, jimi))

    def convertHlsToRgb(self):
        """Convert HLS to RGB.

        Takes the hue, lightness, and saturation channels of a colour image
        and returns an image node containing a colour RGB image.

        :param jim: multi-band Jim with three bands representing hue, saturation,
        and intensity channels
        :return: Jim with three bands containing the RGB channels
        """
        assert jim.properties.nrOfBand() == 3, \
            'Error: input jim must be multi-band image with three bands ' \
            '(h, s, i)'
        jimh =_pj.geometry.cropBand(jim, 0)
        jiml =_pj.geometry.cropBand(jim, 1)
        jims =_pj.geometry.cropBand(jim, 2)

        return pj.Jim(jimh._jipjim.convertHlsToRgb(jiml, jims))

    def distance2d4(self, band=0):
        """Compute the 2-dimensional 4-connected distance function of a Jim.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_distance2d4(band)

    def distance2dChamfer(self, type, band=0):
        """Compute the chamfer distance function of Jim.

        Note that the type of Jim must be large enough to hold the largest
        distance otherwise overflow may occur and generate artefacts in the
        resulting distance function.

        Modifies the instance on which the method was called.

        :param type: Integer for type of chamfer distance {1, 11, 34, 57, 5711}
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_distance2dChamfer(type, band)

    def distance2dEuclideanConstrained(self, mask, band=0):
        """Compute the Euclidean geodesic distance from the marker set.

        Defined by the image marker and within the geodesic mask defined by
        the image mask. The algorithm is described in (Soille, 1991).

        Modifies the instance on which the method was called.

        :param mask: binary UCHAR image with geodesic mask pixels set to 1
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._set(
            self._jim_object._jipjim.distance2dEuclideanConstrained(
                mask._jipjim, band))

    def distance2dEuclideanSquared(self, band=0):
        """Compute the squared Euclidean distance transform.

        im must be a 2-D binary image. Original algorithm proposed by Saito
        and Toriwaki (1994) and then optimised independently by (Hirata,
        1996) and (Meijster et al., 2000). Based on the
        Euclidean distance transform. Note that a temporary buffer of type
        UINT16 is used for sums along/lines and columns so that uncontrolled
        results will occur if an object shows more than 16 2 /2 foreground
        pixels along a given line or column.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._set(
            self._jim_object._jipjim.distance2dEuclideanSquared(band))

    def distanceGeodesic(self, marker, graph, band=0):
        """Compute geodesic distance function from marker within mask.

        Using graph connectivity.

        Modifies the instance on which the method was called.

        :param marker: an image node for marker image
        :param graph: integer for connectivity
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_jldistanceGeodesic(
            marker._jipjim, graph, band)

    def distanceInfluenceZones2dEuclidean(self, band=0):
        """Output error-free influence zones of the labelled regions in Jim.

        Jim must be a 2-D image, its last bit being reserved (e.g., for a UCHAR
        images, authorized label values range from 1 to 127).

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._set(
            self._jim_object._jipjim.distanceInfluenceZones2dEuclidean(band))

    def labelConstrainedCCsVariance(self, ox, oy, oz, rg, rl, varmax, graph=4):
        """Label image.

        Modifies the instance on which the method was called.

        :param jim: a Jim object
        :param ox: x coordinate of origin of imngb
        :param oy: y coordinate of origin of imngb
        :param oz: z coordinate of origin of imngb
        :param rg: integer for range parameter lambda g
        :param rl: integer for range parameter lambda l
        :param varmax: float for maximum variance of cc
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images, default is 4)
        :return: labeled Jim object
        """
        _pj._check_graph(graph, [4, 8])

        if graph == 4:
            ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb[0, 1] = 1
            ngb[1, 0] = 1
            ngb[1, 2] = 1
            ngb[2, 1] = 1
        else:
            ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb.pixops.setData(1)
            ngb[1, 1] = 0

        self._jim_object._set(
            self._jim_object._jipjim.labelConstrainedCCsVariance(
                ngb._jipjim, ox, oy, oz, rg, rl, varmax))

    def labelErode(self, graph=4):
        """Label to contain the union of the erosions of the iso-intensity CCs.

        Modifies the instance on which the method was called.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images, default is 4)
        """
        _pj._check_graph(graph, [4, 8])

        self._jim_object._set(self._jim_object._jipjim.labelErode(graph))

    def labelFlatZonesSeeded(self, jim_ngb, jim_seeds, ox, oy, oz):
        """Label the flat regions using the ngb Jim and its origin.

        :param jim_ngb: a Jim objects for neighbourhood
        :param jim_seeds: A Jim object for seeds (must be of type UCHAR)
        :param ox: x coordinate of origin of jim_ngb
        :param oy: y coordinate of origin of jim_ngb
        :param oz: z coordinate of origin of jim_ngb
        """
        self._jim_object._jipjim.d_labelFlatZonesSeeded(
            jim_ngb._jipjim, jim_seeds._jipjim, ox, oy, oz)

    def labelGraph(self, graph=8):
        """Label each non-zero connected component with a unique label.

        Uses graph-connectivity

        Modifies the instance on which the method was called.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images, default is 8)
        :return: labeled Jim object
        """
        _pj._check_graph(graph, [4, 8])

        if graph == 4:
            ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb[0, 1] = 1
            ngb[1, 0] = 1
            ngb[1, 2] = 1
            ngb[2, 1] = 1
        else:
            ngb = _pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb.pixops.setData(1)
            ngb[1, 1] = 0

        self._jim_object._jipjim.d_labelBinary(ngb._jipjim, 1, 1, 0)

    def labelImagePixels(self):
        """Label each non-zero pixel of im with a unique label.

        Labels unless label overflows occurs.

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_labelPix()

        #todo: not working (self is not modified)
    def morphoFillHoles(self, graph, borderFlag=1):
        """Remove the image-border-connected regional minima.

        Uses graph connectivity.

        Modifies the instance on which the method was called.

        :param borderFlag:
        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        :return: a new Jim object with the connected component of the input
            object removed
        """
        # ajim=_pj.Jim(self._jim_object._jipjim)
        maxval = self._jim_object.stats.getStats('max')['max']
        marker = _pj.pixops.setData(self._jim_object, maxval)
        if borderFlag:
            marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
        else:
            marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
        marker.pixops.supremum(self._jim_object)

        marker.ccops.morphoGeodesicReconstructionByErosion(self._jim_object,
                                                           graph, borderFlag)
        self._jim_object._set(marker._jipjim)

    def morphoGeodesicReconstructionByDilation(self, jim_object_mask, graph,
                                               flag=1):
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
        :return: jim_object_mark containing the result of the morphological
            reconstruction by dilation
        """
        self._jim_object._jipjim.d_geodesicReconstructionByDilation(
            jim_object_mask._jipjim, graph, flag)

    def morphoGeodesicReconstructionByErosion(self, jim_object_mask, graph,
                                              flag=1):
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
        :return: jim_object_mark containing the result of the morphological
            reconstruction by erosion
        """
        self._jim_object._jipjim.d_geodesicReconstructionByErosion(
            jim_object_mask._jipjim, graph, flag)

    def morphoRemoveBorder(self, graph):
        """Remove the border-connected components using graph connectivity.

        Modifies the instance on which the method was called.

        :param graph: an integer holding for the graph connectivity
            (4 or 8 for 2-D images)
        """
        ajim = self._jim_object._jipjim
        minmax = ajim.stats.getStats(function=['min', 'max'])
        minval = minmax['min']
        maxval = minmax['max']
        marker = _pj.pixops.setData(ajim, minval)
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
        marker.pixops.infimum(ajim)
        marker.ccops.morphoGeodesicReconstructionByDilation(ajim, graph)
        marker.pixops.simpleArithOp(ajim, 17)  # 17 for SUBSWAP_op_ovf
        self._jim_object._set(marker)


class _CCOpsList(_JimListModuleBase):
    """Define all CCOps methods for JimLists."""

    pass


class _CCOpsVect(_JimVectModuleBase):
    """Define all CCOps methods for JimVects."""

    pass
