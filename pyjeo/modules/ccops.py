"""Module for connected-component operations."""

import pyjeo as _pj


def labelImagePixels(jim_object):
    """Label each non-zero pixel of im with a unique label.

    Labels unless label overflows occurs.

    :param jim_object: a Jim object
    :return: labeled Jim object
    """
    return _pj.Jim(jim_object._jipjim.labelPix())


def labelGraph(jim_object, graph=8):
    """Label each non-zero connected component with a unique label using graph-connectivity

    :param jim_object: a Jim object holding a binary image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 8)
    :return: labeled Jim object
    """
    if (graph==4):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0,1]=1
        ngb[1,0]=1
        ngb[1,2]=1
        ngb[2,1]=1
    elif (graph==8):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1,1]=0
    else:
        print("graph must be equal to 4 or 8")
        raise ValueError('graph must be equal to 4 or 8')

    return _pj.Jim(jim_object._jipjim.labelBinary(ngb._jipjim, 1, 1, 0))


def labelFlatZonesGraph(jim_object, graph=8):
    """Label each image flat zone with a unique label using graph-connectivity

    :param jim_object: a Jim object holding a grey level image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 8)
    :return: labeled Jim object
    """
    if (graph==4):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0,1]=1
        ngb[1,0]=1
        ngb[1,2]=1
        ngb[2,1]=1
    elif (graph==8):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1,1]=0
    else:
        print("graph must be equal to 4 or 8")
        raise ValueError('graph must be equal to 4 or 8')

    return _pj.Jim(jim_object._jipjim.labelFlatZones(ngb._jipjim, 1, 1, 0))


def labelAlphaCCsGraph(jim_object, localRange, globalRange, graph=8):
    """Label each alpha-connected component with a unique label using graph-connectivity

    :param jim_object: a Jim object holding a grey level image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 8)
    :return: labeled Jim object
    """
    if (graph==4):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb[0,1]=1
        ngb[1,0]=1
        ngb[1,2]=1
        ngb[2,1]=1
    elif (graph==8):
        ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
        ngb.pixops.setData(1)
        ngb[1,1]=0
    else:
        print("graph must be equal to 4 or 8")
        raise ValueError('graph must be equal to 4 or 8')

    return _pj.Jim(jim_object._jipjim.labelAlphaCCs(ngb._jipjim, 1, 1, 0, localRange, globalRange))


def labelConstrainedCCsMultiband(jim_object, se, ox, oy, oz, r1, r2):
    """TODO

    :param jim_object: a multi-band Jim object
    :param se:
    :param ox:
    :param oy:
    :param oz:
    :param r1:
    :param r2:
    :return: labeled Jim object
    """
    if jim_object.properties.nrOfBand() < 2:
        print("Jim image must be a multi-band image")
        raise ValueError('Jim image must be a multi-band image')

    return _pj.Jim(jim_object._jipjim.labelConstrainedCCsMultiband(se._jipjim, ox, oy, oz, r1, r2))


def distance2dEuclideanSquared(jim_object, band=0):
    """Compute the squared Euclidean distance transform of im.

    im must be a 2-D binary image. Original algorihtm proposed by Saito
    and Toriwaki (1994) and then optimised independently by (Hirata,
    1996) and (Meijster et al., 2000). See also *edt for the actual
    Euclidean distance transform. Note that a temporary buffer of type
    UINT16 is used for sums along/lines and columns so that uncontrolled
    results will occur if an object shows more than 16 2 /2 foreground
    pixels along a given line or column.

    :param jim_object: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.distance2dEuclideanSquared(band))


def getRegionalMinima(jim_object, graph):
    """Compute the regional minima of the input image.  The pixels belonging to a regional minimum are set to 1, all other pixels are set to 0.

    :param jim_object: a Jim object
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)
    :return: a new Jim object of type unsigned char containing the regional minima of the input Jim object 
    """
    return _pj.Jim(jim_object._jipjim.getRegionalMinima(graph))


def morphoGeodesicReconstructionByDilation(jim_object_mark, jim_object_mask, graph, borderFlag=1):
    """Compute the morphological reconstruction by dilation of mask image from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)
    :param borderFlag: integer (border values are set to PIX MIN in BOTH images if flag equals
                 0, otherwise the image are internally processed by adding a border
                 which is then removed at the end of the processing). Default value is 1.
    :return: jim_object_mark containing the result of the morphological reconstruction by dilation
    """

    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByDilation(jim_object_mask._jipjim, graph, borderFlag))


def morphoGeodesicReconstructionByErosion(jim_object_mark, jim_object_mask, graph, borderFlag=1):
    """Compute the morphological reconstruction by erosion of mask image from mark image using graph connectivity.

    :param jim_object_mark: a Jim object for marker image
    :param jim_object_mask: a Jim object for mask image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)
    :param borderFlag: integer (border values are set to PIX MIN in BOTH images if flag equals
                 0, otherwise the image are internally processed by adding a border
                 which is then removed at the end of the processing). Default value is 1.
    :return: jim_object_mark containing the result of the morphological reconstruction by erosion
    """

    return _pj.Jim(jim_object_mark._jipjim.geodesicReconstructionByErosion(jim_object_mask._jipjim, graph, borderFlag))


def morphoRemoveBorder(ajim, graph):
    """
    Remove the connected components of an image that are connected to the image border using graph connectivity
    :param ajim: input Jim object
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)

    :return: a new Jim object with the connected component of the imput object removed
    """
    minmax=ajim.stats.getStats(function=['min','max'])
    minval=minmax['min']
    maxval=minmax['max']
    # marker=_pj.pixops.blank(ajim, minval)
    marker=_pj.pixops.setData(ajim, minval)
    marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
    marker.pixops.infimum(ajim)
    marker.ccops.morphoGeodesicReconstructionByDilation(ajim, graph)
    marker.pixops.simpleArithOp(17, ajim) # 17 for SUBSWAP_op_ovf
    return marker


def morphoFillHoles(ajim, graph, borderFlag=1):
    """
    Remove the regional minima of the image that are not connected to the image border using graph connectivity
    :param ajim: input Jim object
    :param borderFlag:
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)

    :return: a new Jim object with the connected component of the imput object removed
    """
    maxval=ajim.stats.getStats('max')['max']
    marker=_pj.pixops.setData(ajim, maxval)
    if borderFlag:
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
    else:
        marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
    marker.pixops.supremum(ajim)
    marker.ccops.morphoGeodesicReconstructionByErosion(ajim, graph, borderFlag)
    return marker

class _CCOps():
    """Define all CCOps methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def labelImagePixels(self):
        """Label each non-zero pixel of im with a unique label.

        Labels unless label overflows occurs.

        Modifies the instance on which the method was called.
        """
        self._jim_object._jipjim.d_labelPix()

    def labelGraph(self, graph=8):
        """Label each non-zero connected component with a unique label using graph-connectivity

           :param jim_object: a Jim object holding a binary image
           :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 8)
           :return: labeled Jim object
         """
        if (graph==4):
            ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb[0,1]=1
            ngb[1,0]=1
            ngb[1,2]=1
            ngb[2,1]=1
        elif (graph==8):
            ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
            ngb.pixops.setData(1)
            ngb[1,1]=0
        else:
            print("graph must be equal to 4 or 8")
            raise ValueError('graph must be equal to 4 or 8')

        self._jim_object._jipjim.d_labelBinary(ngb._jipjim, 1, 1, 0)

    def distance2dEuclideanSquared(self, band=0):
        """Compute the squared Euclidean distance transform

        im must be a 2-D binary image. Original algorihtm proposed by Saito
        and Toriwaki (1994) and then optimised independently by (Hirata,
        1996) and (Meijster et al., 2000). See also *edt for the actual
        Euclidean distance transform. Note that a temporary buffer of type
        UINT16 is used for sums along/lines and columns so that uncontrolled
        results will occur if an object shows more than 16 2 /2 foreground
        pixels along a given line or column.

        Modifies the instance on which the method was called.
        """
        self._jim_object._set(
            self._jim_object._jipjim.distance2dEuclideanSquared(band))

    def morphoGeodesicReconstructionByDilation(self, jim_object_mask, graph, flag=1):
        """Compute the morphological reconstruction by dilation of the current object from mark image using graph connectivity.

        :param jim_object_mask: a Jim object for mask image
        :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)
        :param flag: integer (border values are set to PIX MIN in BOTH images if flag equals
                 0, otherwise the image are internally processed by adding a border
                 which is then removed at the end of the processing). Default value is 1.
        :return: jim_object_mark containing the result of the morphological reconstruction by dilation
        """
        self._jim_object._jipjim.d_geodesicReconstructionByDilation(jim_object_mask._jipjim, graph, flag)

    def morphoGeodesicReconstructionByErosion(self, jim_object_mask, graph, flag=1):
        """Compute the morphological reconstruction by erosion of the current object from mark image using graph connectivity.

        :param jim_object_mask: a Jim object for mask image
        :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)
        :param flag: integer (border values are set to PIX MIN in BOTH images if flag equals
                 0, otherwise the image are internally processed by adding a border
                 which is then removed at the end of the processing). Default value is 1.
        :return: jim_object_mark containing the result of the morphological reconstruction by erosion
        """
        self._jim_object._jipjim.d_geodesicReconstructionByErosion(jim_object_mask._jipjim, graph, flag)

    def morphoRemoveBorder(self, graph):
        """
        Remove the connected components that are connected to the image border using graph connectivity

        :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)

        """
        ajim=self._jim_object._jipjim
        minmax=ajim.stats.getStats(function=['min','max'])
        minval=minmax['min']
        maxval=minmax['max']
        # marker=_pj.pixops.blank(ajim, minval)
        marker=_pj.pixops.setData(ajim, minval)
        marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, maxval)
        marker.pixops.infimum(ajim)
        marker.ccops.morphoGeodesicReconstructionByDilation(ajim, graph)
        marker.pixops.simpleArithOp(17, ajim) # 17 for SUBSWAP_op_ovf
        self._jim_object._set(marker)

        #todo: not working (self is not modified)
    def morphoFillHoles(self, graph, borderFlag=1):
        """
        Remove the regional minima of the image that are not connected to the image border using graph connectivity
        :param borderFlag:
        :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images)

        :return: a new Jim object with the connected component of the imput object removed
        """
        # ajim=_pj.Jim(self._jim_object._jipjim)
        maxval=self._jim_object.stats.getStats('max')['max']
        marker=_pj.pixops.setData(self._jim_object, maxval)
        if borderFlag:
            marker.geometry.imageFrameSet(1, 1, 1, 1, 0, 0, 0)
        else:
            marker.geometry.imageFrameSet(2, 2, 2, 2, 0, 0, 0)
        marker.pixops.supremum(self._jim_object)

        marker.ccops.morphoGeodesicReconstructionByErosion(self._jim_object, graph, borderFlag)
        self._jim_object._set(marker._jipjim)


class _CCOpsList():
    """Define all CCOps methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller


class _CCOpsVect():
    """Define all CCOps methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller
