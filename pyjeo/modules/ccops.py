"""Module for connected-component operations."""

import pyjeo as _pj


def labelImagePixels(jim_object):
    """Label each non-zero pixel of im with a unique label.

    Labels unless label overflows occurs.

    :param jim_object: a Jim object
    :return: labeled Jim object
    """
    return _pj.Jim(jim_object._jipjim.labelPix())


def labelGraph(jim_object, graph=4):
    """Label each non-zero connected component with a unique label using graph-connectivity

    :param jim_object: a Jim object holding a binary image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
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


def labelFlatZonesGraph(jim_object, graph=4):
    """Label each image flat zone with a unique label using graph-connectivity

    :param jim_object: a Jim object holding a grey level image
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
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


def labelConstrainedCCsGraph(jim_object, localRange, globalRange, graph=4):
    """Label each alpha-omega connected components with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim_object: a Jim object holding a grey level image
    :param localRange: integer value indicating maximum absolute local difference between 2 adjacent pixels
    :param globalRange: integer value indicating maximum global difference (difference between the maximum and minimum values of each resulting connected component)
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
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

    return _pj.Jim(jim_object._jipjim.labelConstrainedCCs(ngb._jipjim, 1, 1, 0, localRange, globalRange))


def labelConstrainedCCsMultibandGraph(jim_object, localRange, globalRange, graph=4):
    """Label each alpha-omega connected component with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim_object: a Jim object holding a multi-band image
    :param localRange: integer value indicating maximum absolute local difference between 2 adjacent pixels
    :param globalRange: integer value indicating maximum global difference (difference between the maximum and minimum values of each resulting connected component)
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    # if jim_object.properties.nrOfBand() < 2:
    #     print("Jim image must be a multi-band image")
    #     raise ValueError('Jim image must be a multi-band image')

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

    return _pj.Jim(jim_object._jipjim.labelConstrainedCCsMultiband(ngb._jipjim, 1, 1, 0, localRange, globalRange))


def labelStronglyCCsMultibandGraph(jim_object, localRange, graph=4):
    """Label each strongly alpha-connected component with a unique label using graph-connectivity :cite:`soille2008pami`

    :param jim_object: a Jim object holding a multi-band image
    :param localRange: integer value indicating maximum absolute local difference between 2 adjacent pixels
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
    :return: labeled Jim object
    """
    # if jim_object.properties.nrOfBand() < 2:
    #     print("Jim image must be a multi-band image")
    #     raise ValueError('Jim image must be a multi-band image')

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

    return _pj.Jim(jim_object._jipjim.labelStronglyCCsMultiband(ngb._jipjim, 1, 1, 0, localRange))


def segmentImageMultiband(jim_object, localRange, regionSize, contrast=0, version=0, graph=4, dataFileNamePrefix=""):
    """Multiband image segmentation based on the method described in :cite:`brunner-soille2007`

The contrast threshold value is used for merging the regions with similar contrast as follows: < 0 (do not perform region merge), 0 (determine best contrast value automatically), and > 0 (use this value as threshold value).  Authorised version values are: 0 (compare to whole region), 1 (compare to original seeds), and 2 (compare to pixel neighbours).  If the optional string dataFileNamePrefix is given, data files to use with gnuplot are stored in dataFileNamePrefix_xxx.dat, otherwise data files are not generated (default).

    :param jim_object: a Jim object holding a multi-band image
    :param localRange: integer value indicating maximum absolute local difference between 2 adjacent pixels
    :param regionSize: integer value for minimum size of iso-intensity region in output image (must be >= 2 pixels)
    :param contrast: (default is 0)
    :param version:  (default is 0)
    :param graph: an integer holding for the graph connecvity (4 or 8 for 2-D images, default is 4)
    :param dataFileName:  
    :return: labeled Jim object
    """
    # if jim_object.properties.nrOfBand() < 2:
    #     print("Jim image must be a multi-band image")
    #     raise ValueError('Jim image must be a multi-band image')

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

    return _pj.Jim(jim_object._jipjim.segmentImageMultiband(graph, localRange, regionSize, contrast, version, dataFileNamePrefix))


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
    """Remove the connected components of an image that are connected to the image border using graph connectivity
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
    """Remove the regional minima of the image that are not connected to the image border using graph connectivity (originally proposed for removing pits in digital elevation models, see :cite:`soille-ansoult90` and  :cite:`soille-gratin94` for a fast implementation)
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



# (defun *labelccdissim_ismm2011 (im omega alpha &key (diamondwidth 2))
# ; \lspfunction{*}{labelccdissim_ismm2011}{im omega alpha &key (diamondwidth 2)}
# ; \param{im}{}
# ; \param{omega}{}
# ; \param{alpha}{}
# ; \param{(diamondwidth}{}
# ; \return{}
# ; \desc{see \cite{soille2011ismm}}
# ; \myseealso{}
# ; \lspfile{\crtlspfile}
# ; \example{}{}
#   (let* (
# 	 (mingraderograddil
# 	  (@inf (*graderographframe im diamondwidth 4)
# 		(*graddilgraphframe im diamondwidth 4)
# 		)
# 	  )
# 	 (DIR-HORI 0)
# 	 (DIR-VERT 1)
# 	 (ABS-DIFF_op 0)
# 	 (MAX_op 1)
# 	 (MIN_op 2)
# 	 ; (h_dissim) (v_dissim)
# 	 )
#     (setq h_dissim (*edgeweight mingraderograddil DIR-HORI MAX_op))
#     (setq v_dissim (*edgeweight mingraderograddil DIR-VERT MAX_op))
#     (*imfree mingraderograddil)
#     (@sup h_dissim
# 	  (*edgeweight im DIR-HORI ABS-DIFF_op)
# 	  )
#     (@sup v_dissim
# 	  (*edgeweight im DIR-VERT ABS-DIFF_op))

#     (@subframebox
#      (*labelccdissim (*addframebox im 1 1 1 1 0 0 255)
# 		     (*addframebox h_dissim 1 1 1 1 0 0 255)
# 		     (*addframebox v_dissim 1 1 1 1 0 0 255)
# 		     omega
# 		     alpha)
#      1 1 1 1 0 0)
#     )
#   )



def labelConstrainedCCsGraphDissim(jim_object, localRange, globalRange):
    """Label each alpha-omega connected components with a unique label using graph-connectivity and the dissimilarity measure countering the chaining effect as described in :cite:`soille2011ismm`

    :param jim_object: a Jim object holding a grey level image
    :param localRange: integer value indicating maximum absolute local difference between 2 adjacent pixels
    :param globalRange: integer value indicating maximum global difference (difference between the maximum and minimum values of each resulting connected component)
    :return: labeled Jim object
    """

    # Graph is 4 for all dissimilarity based on edge weights
    ngb=_pj.Jim(ncol=3, nrow=3, otype='Byte')
    ngb[0,1]=1
    ngb[1,0]=1
    ngb[1,2]=1
    ngb[2,1]=1

    mingraderograddil = _pj.pixops.INF(_pj.ngb.morphoErode

    return _pj.Jim(jim_object._jipjim.labelConstrainedCCs(ngb._jipjim, 1, 1, 0, localRange, globalRange))




#    ajim.ngb.edgeWeight()
#    out=_pj.ngb.edgeWeight(ajim)

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
