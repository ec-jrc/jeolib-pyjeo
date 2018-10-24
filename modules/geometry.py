import pyjeo as _pj


def crop3D(jim_object, x1, y1, z1, x2, y2, z2, iband=0):
    """Extract a subimage from Jim object by cutting it from given coords.

    :param jim_object: a Jim object
    :param x1: x coordinate of 1st pixel
    :param y1: y coordinate of 1st pixel
    :param z1: z coordinate of 1st pixel
    :param x2: x coordinate of 2nd pixel
    :param y2: y coordinate of 2nd pixel
    :param z2: z coordinate of 2nd pixel
    :param iband: List of band indices to crop (index is 0 based)
    :return: Extracted subimage as Jim instance
    """
    return _pj.Jim(jim_object.imageCut(x1, y1, z1, x2, y2, z2, iband))


def plotLine(jim_object, x1, y1, x2, y2, val):
    """Draw a line from [x1, y1] to [x2, y2] by setting pixels of Jim to val

    :param jim_object: a Jim object
    :param x1: an integer for x-coordinate of 1st point
    :param y1: an integer for y-coordinate of 1st point
    :param x2: an integer for x-coordinate of 2nd point
    :param y2: an integer for y-coordinate of 2nd point
    :return:
    """
    return _pj.Jim(jim_object.plotLine(x1, y1, x2, y2, val))


class _Geometry():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def crop3D(self, leftSize, rightSize, topSize, belowSize, upSize,
               downSize):
        """Peel off a frame of im.

        Modifies the instance on which the method was called.

        :param leftSize: width of left frame
        :param rightSize: width if right frame
        :param topSize: width of top frame
        :param belowSize: width if bottom frame
        :param upSize: width of upper frame
        :param downSize: width if lower frame
        """
        self._jim_object._set(self._jim_object.d_imageFrameSubstract(
            [leftSize, rightSize, topSize, belowSize, upSize, downSize]))

    def crop(self, ulx=0, uly=0, lrx=0, lry=0, **kwargs):
        """Subset raster dataset.

        Subset raster dataset according in spatial (subset region) or
        spectral/temporal domain (subset bands).

        Modifies the instance on which the method was called.

        :param ulx: Upper left x value of bounding box to crop
        :param uly: Upper left y value of bounding box to crop
        :param lrx: Lower right x value of bounding box to crop
        :param lry: Lower right y value of bounding box to crop
        :param kwargs: See table below

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | dx               | Output resolution in x (default: keep original resolution)                      |
        +------------------+---------------------------------------------------------------------------------+
        | dy               | Output resolution in y (default: keep original resolution)                      |
        +------------------+---------------------------------------------------------------------------------+
        | nodata           | Nodata value to put in image if out of bounds                                   |
        +------------------+---------------------------------------------------------------------------------+
        | align            | Align output bounding box to input image                                        |
        +------------------+---------------------------------------------------------------------------------+
        """
        kwargs.update({'ulx': ulx})
        kwargs.update({'uly': uly})
        kwargs.update({'lrx': lrx})
        kwargs.update({'lry': lry})
        self._jim_object._set(self._jim_object.crop(kwargs))

    def cropOgr(self, extent, **kwargs):
        """Subset raster dataset.

        Subset raster dataset in spatial domain defined by a vector dataset.

        Modifies the instance on which the method was called.

        :param extent: Get boundary from extent from polygons in vector file
        :param kwargs: See table below

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | ln               | Layer name of extent to crop                                                    |
        +------------------+---------------------------------------------------------------------------------+
        | eo               | Special extent options controlling rasterization                                |
        +------------------+---------------------------------------------------------------------------------+
        | crop_to_cutline  | True will crop the extent of the target dataset to the extent of the cutline    |
        |                  | The outside area will be set to no data (the value defined by the key 'nodata') |
        +------------------+---------------------------------------------------------------------------------+
        | crop_in_cutline  | True: inverse operation to crop_to_cutline                                      |
        |                  | The inside area will be set to no data (the value defined by the key 'nodata')  |
        +------------------+---------------------------------------------------------------------------------+
        | dx               | Output resolution in x (default: keep original resolution)                      |
        +------------------+---------------------------------------------------------------------------------+
        | dy               | Output resolution in y (default: keep original resolution)                      |
        +------------------+---------------------------------------------------------------------------------+
        | nodata           | Nodata value to put in image if out of bounds                                   |
        +------------------+---------------------------------------------------------------------------------+
        | align            | Align output bounding box to input image                                        |
        +------------------+---------------------------------------------------------------------------------+

        .. note::
           Possible values for the key 'eo' are: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG. For instance you can use 'eo':'ATTRIBUTE=fieldname'
        """
        self._jim_object._set(self._jim_object.cropOgr(extent,kwargs))

    def cropBand(self, band):
        """Subset raster dataset.

        Subset raster dataset in spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0=jl.io.createJim('/path/to/raster0.tif')
            jim0.cropBand(band=[0,1,2])

        """
        self._jim_object.d_cropBand({'band': band})

    def cropBandRange(self, startband, endband):
        """Subset raster dataset.

        Subset raster dataset in range of spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param startband: Start band sequence number (index is 0 based)
        :param endband: End band sequence number (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0=jl.io.createJim('/path/to/raster0.tif')
            jim0.cropBandRange(startband=0,startBand=2)

        """
        self._jim_object.d_cropBand({'startband': startband,'endband': endband})

    def stackBand(self, band):
        """Stack bands from another raster dataset to current raster dataset.

        Modifies the instance on which the method was called.

        :param band: List of band indices to stack (index is 0 based)

        Example:

        Append all the bands of raster dataset jim1 to the current image jim0::

            jim0=jl.io.createJim('/path/to/raster0.tif')
            jim1=jl.io.createJim('/path/to/raster1.tif')
            jim0.stackBand(jim1)

        Append the first three bands of raster dataset jim1 to the current image jim0::

            jim0=jl.io.createJim('/path/to/raster0.tif')
            jim1=jl.io.createJim('/path/to/raster1.tif')
            jim0.stackBand(jim1,band=[0,1,2])
        """
        self._jim_object.d_stackBand({'band': band})

    def stackBandRange(self, startband, endband):
        """Subset raster dataset.

        Stack range of bands from another raster dataset to current raster dataset.

        Modifies the instance on which the method was called.

        :param startband: Start band sequence number (index is 0 based)
        :param endband: End band sequence number (index is 0 based)

        Example:

        Append the first three bands of raster dataset jim1 to the current image jim0::

            jim0=jl.io.createJim('/path/to/raster0.tif')
            jim1=jl.io.createJim('/path/to/raster1.tif')
            jim0.stackBandRange(jim1,startband=0,endband=2])
        """
        self._jim_object.d_stackBand({'startband': startband,'endband': endband})

    def extractOgr(self, jim_ref, **kwargs):
        """Extract pixel values from raster image based on a vector dataset sample.

        :param jim_ref: reference Jim instance
        :param kwargs: See table below
        :return: A VectorOgr with the same geometry as the sample vector
            dataset and an extra field for each of the calculated raster value
            (zonal) statistics. The same layer name(s) of the sample will be
            used for the output vector dataset


        +------------------+--------------------------------------------------------------------------------------------------------+
        | key              | value                                                                                                  |
        +==================+========================================================================================================+
        | rule             | Rule how to calculate zonal statistics per feature (see list of :ref:`supported rules <extract_rules>`)|
        +------------------+--------------------------------------------------------------------------------------------------------+
        | copy             | Copy these fields from the sample vector dataset (default is to copy all fields)                       |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | label            | Create extra field named 'label' with this value                                                       |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field identifier (sequence of features)                       |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is to use extract all bands                              |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands to extract                                            |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)                                                                 |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | endband          | End band sequence number (0 indexed)                                                                   |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | output           | Name of the output vector dataset in which the zonal statistics are saved                              |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | oformat          | Output vector dataset format                                                                           |
        +------------------+--------------------------------------------------------------------------------------------------------+
        | co               | Creation option for output vector dataset                                                              |
        +------------------+--------------------------------------------------------------------------------------------------------+

        .. _extract_rules:

        :Supported rules to extract:

        +------------------+---------------------------------------------------------------------------------------------------+
        | rule             | description                                                                                       |
        +==================+===================================================================================================+
        | point            | extract a single pixel within the polygon or on each point feature                                |
        +------------------+---------------------------------------------------------------------------------------------------+
        | allpoints        | Extract all pixel values covered by the polygon                                                   |
        +------------------+---------------------------------------------------------------------------------------------------+
        | centroid         | Extract pixel value at the centroid of the polygon                                                |
        +------------------+---------------------------------------------------------------------------------------------------+
        | mean             | Extract average of all pixel values within the polygon                                            |
        +------------------+---------------------------------------------------------------------------------------------------+
        | stdev            | Extract standard deviation of all pixel values within the polygon                                 |
        +------------------+---------------------------------------------------------------------------------------------------+
        | median           | Extract median of all pixel values within the polygon                                             |
        +------------------+---------------------------------------------------------------------------------------------------+
        | min              | Extract minimum value of all pixels within the polygon                                            |
        +------------------+---------------------------------------------------------------------------------------------------+
        | max              | Extract maximum value of all pixels within the polygon                                            |
        +------------------+---------------------------------------------------------------------------------------------------+
        | sum              | Extract sum of the values of all pixels within the polygon                                        |
        +------------------+---------------------------------------------------------------------------------------------------+
        | mode             | Extract the mode of classes within the polygon (classes must be set with the option class)        |
        +------------------+---------------------------------------------------------------------------------------------------+
        | proportion       | Extract proportion of class(es) within the polygon (classes must be set with the option class)    |
        +------------------+---------------------------------------------------------------------------------------------------+
        | count            | Extract count of class(es) within the polygon (classes must be set with the option class)         |
        +------------------+---------------------------------------------------------------------------------------------------+
        | percentile       | Extract percentile as defined by option perc (e.g, 95th percentile of values covered by polygon)  |
        +------------------+---------------------------------------------------------------------------------------------------+


        .. note::
            To ignore some pixels from the extraction process, see list of :ref:`mask <extract_mask>` key values:

        .. _extract_mask:

        :Supported key values to mask pixels that must be ignored in the extraction process:

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | srcnodata        | List of nodata values not to extract                                            |
        +------------------+---------------------------------------------------------------------------------+
        | bndnodata        | List of band in input image to check if pixel is valid (used for srcnodata)     |
        +------------------+---------------------------------------------------------------------------------+
        | mask             | Use the the specified file as a validity mask                                   |
        +------------------+---------------------------------------------------------------------------------+
        | mskband          | Use the the specified band of the mask file defined                             |
        +------------------+---------------------------------------------------------------------------------+
        | msknodata        | List of mask values not to extract                                              |
        +------------------+---------------------------------------------------------------------------------+
        | threshold        | Maximum number of features to extract. Use percentage value as string           |
        |                  | (e.g., '10%') or integer value for absolute threshold                           |
        +------------------+---------------------------------------------------------------------------------+

        Example:

        Extract a random sample of 100 points, calculating the mean value based on a 3x3 window (buffer value of 1 pixel neighborhood) in a vector dataset in memory::

            reference=jl.io.createVector('/path/to/reference.sqlite')
            jim0=jl.io.createJim('/path/to/raster.tif')
            v=jim0.extractOgr(reference,rule=['mean'],output='/path/to/output.sqlite',oformat='SQLite')
            v.write()
        """
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold']=float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold']=-kwargs['threshold']
        return self._jim_object.extractOgr(jim_ref, kwargs)

    def extractSample(self, **kwargs):
        """Extract a random or grid sample from a raster dataset.

        :param kwargs: See table below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | random           | Extract a random sample with a size equal to the defined value                  |
        +------------------+---------------------------------------------------------------------------------+
        | grid             | Extract a grid sample with a grid size equal to the defined value               |
        +------------------+---------------------------------------------------------------------------------+
        | rule             | Rule how to calculate zonal statistics per feature                              |
        |                  | (see list of :ref:`supported rules <extract_rules>`)                            |
        +------------------+---------------------------------------------------------------------------------+
        | buffer           | Buffer for calculating statistics for point features (in number of pixels)      |
        +------------------+---------------------------------------------------------------------------------+
        | label            | Create extra field named 'label' with this value                                |
        +------------------+---------------------------------------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field identifier (sequence of features)|
        +------------------+---------------------------------------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is to use extract all bands       |
        +------------------+---------------------------------------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands to extract                     |
        +------------------+---------------------------------------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)                                          |
        +------------------+---------------------------------------------------------------------------------+
        | endband          | End band sequence number (0 indexed)                                            |
        +------------------+---------------------------------------------------------------------------------+
        | output           | Name of the output vector dataset in which the zonal statistics are saved       |
        +------------------+---------------------------------------------------------------------------------+
        | ln               | Layer name of output vector dataset                                             |
        +------------------+---------------------------------------------------------------------------------+
        | oformat          | Output vector dataset format                                                    |
        +------------------+---------------------------------------------------------------------------------+
        | co               | Creation option for output vector dataset                                       |
        +------------------+---------------------------------------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list of :ref:`mask <extract_mask>` key values:

        Example:

        Extract a random sample of 100 points, calculating the mean value based on a 3x3 window (buffer value of 1 pixel neighborhood) in a vector dataset in memory::

            v01=jim0.extractSample(random=100,buffer=1,rule=['mean'],output='mem01',oformat='Memory')

        Extract a sample of 100 points using a regular grid sampling scheme. For each grid point, calculate the median value based on a 3x3 window (buffer value of 1 pixel neighborhood). Write the result in a SQLite vector dataset on disk::

            outputfn='/path/to/output.sqlite'
            npoint=100
            gridsize=int(jim.nrOfCol()*jim.getDeltaX()/math.sqrt(npoint))
            v=jim.extractSample(grid=gridsize,buffer=1,rule=['median'],output=outputfn,oformat='SQLite')
            v.write()

        """
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold']=float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold']=-kwargs['threshold']
        return self._jim_object.extractSample(kwargs)

    def extractImg(self, reference, **kwargs):
        """Extract pixel values from an input based on a raster sample dataset.

        :reference: thematic raster dataset with integer values, typically a land cover map
        :param kwargs: See table below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | rule             | Rule how to calculate zonal statistics per feature                              |
        |                  | (see list of :ref:`supported rules <extract_rules>`)                            |
        +------------------+---------------------------------------------------------------------------------+
        | class            | List of classes to extract from the raster sample dataset.                      |
        |                  | Leave empty to extract all valid data pixels from thee sample                   |
        +------------------+---------------------------------------------------------------------------------+
        | cname            | Name of the class label in the output vector dataset (default is 'label')       |
        +------------------+---------------------------------------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field identifier (sequence of features)|
        +------------------+---------------------------------------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is to use extract all bands       |
        +------------------+---------------------------------------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands to extract                     |
        +------------------+---------------------------------------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)                                          |
        +------------------+---------------------------------------------------------------------------------+
        | endband          | End band sequence number (0 indexed)                                            |
        +------------------+---------------------------------------------------------------------------------+
        | down             | Down sampling factor to extract a subset of the sample based on a grid          |
        +------------------+---------------------------------------------------------------------------------+
        | output           | Name of the output vector dataset in which the zonal statistics are saved       |
        +------------------+---------------------------------------------------------------------------------+
        | ln               | Layer name of output vector dataset                                             |
        +------------------+---------------------------------------------------------------------------------+
        | oformat          | Output vector dataset format                                                    |
        +------------------+---------------------------------------------------------------------------------+
        | co               | Creation option for output vector dataset                                       |
        +------------------+---------------------------------------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list of :ref:`nodata <extract_nodata>` key values:

        .. _extract_nodata:

        :Supported key values to mask pixels that must be ignored in the extraction process:

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | srcnodata        | List of nodata values not to extract                                            |
        +------------------+---------------------------------------------------------------------------------+
        | bndnodata        | List of band in input image to check if pixel is valid (used for srcnodata)     |
        +------------------+---------------------------------------------------------------------------------+
        | threshold        | Maximum number of features to extract. Use percentage value as string           |
        |                  | (e.g., '10%') or integer value for absolute threshold.                          |
        |                  | You can provide a list of threshold values, one for each class.                 |
        +------------------+---------------------------------------------------------------------------------+
        Example:

        Open a raster sample dataset based on land cover map (e.g., Corine) and use it to extract a stratified sample of 100 points from an input raster dataset with four spectral bands ('B02', 'B03', 'B04', 'B08'). Only sample classes 2 (urban), 12 (agriculture), 25 (forest), 41 (water) and an aggregated (rest) class 50::

            jim_ref=jl.createJim('/path/to/landcovermap.tif')

            classes=[2,12,25,41,50]
            thresholds=['20%','25%','25%','10%','5%']

            jim_ref=jl.createJim('/path/to/multiband.tif','dx'=jim.getDeltaX(),'dy'=jim.getDeltaY(),'ulx'=jim.getUlx(),'uly'=jim.getUly(),'lrx'=jim.getLrx(),'lry'=jim.getLry())

            outputfn='/path/to/output.sqlite'
            sample=jim.extractImg(jim_ref,srcnodata=[0],output=outputfn,class=classes,threshold=thresholds,bandname=['B02','B03','B04','B08'],band=[0,1,2,3])
        """
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold']=float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold']=-kwargs['threshold']
        return self._jim_object.extractImg(reference, kwargs)

    def warp(self, t_srs, **kwargs):
        """
        Warp a raster dataset to a target spatial reference system

        :param t_srs: Target spatial reference system
        :param kwargs: See table below

        +------------------+---------------------------------------------------------------------------------+
        | key              | value                                                                           |
        +==================+=================================================================================+
        | s_srs            | Source spatial reference system (default is to read from input)                 |
        +------------------+---------------------------------------------------------------------------------+
        | resample         | Resample algorithm used for reading pixel data in case of interpolation         |
        |                  | (default: GRIORA_NearestNeighbour).                                             |
        |                  | Check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a         |
        |                  | or available options.                                                           |
        +------------------+---------------------------------------------------------------------------------+
        | nodata           | Nodata value to put in image if out of bounds                                   |
        +------------------+---------------------------------------------------------------------------------+

        Example:

        Read a raster dataset from disk and warp to the target spatial reference system::

            jim=jl.createJim('/path/to/file.tif')
            jim.warp('epsg:3035')

        Read a raster dataset from disk that is in lat lon (epsg:4326), select a bounding box in a different spatial reference system (epsg:3035). Notice the raster dataset read is still in the original projection (epsg:4326). Then warp the raster dataset to the target spatial reference system (epsg:3035)::

            jim=jl.createJim('/path/to/file.tif',t_srs='epsg:3035',ulx=1000000,uly=4000000,lrx=1500000,lry=3500000)
            jim.warp('epsg:3035',s_srs='epsg:4326')

        """
        kwargs.update({'t_srs': t_srs})
        self._jim_object._set(self._jim_object.warp(kwargs))

    def imageInsert(self, sec_jim_object, x, y, z, iband=0):
        """Merge Jim instance with values of sec_jim_object in given coords.

        Modifies the instance on which the method was called.

        :param jim_object: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param iband: List of band indices to crop (index is 0 based)
        """
        # TODO: Doesn't work
        self._jim_object.d_imageInsert(sec_jim_object, x, y, z, iband)

    def imageInsertCompose(self, imlbl, im2, x, y, z, val, iband=0):
        """Merge Jim instance with values of im2 if val of imlbl == val.

        Modifies the instance on which the method was called.

        :param imRaster_imlbl: a Jim object
        :param imRaster_im2: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param val: integer for label value
        :param iband: List of band indices to crop (index is 0 based)
        """
        self._jim_object.d_imageInsertCompose(imlbl, im2, x, y, z, val, iband)

    def imageFrameSet(self, l, r, t, b, u, d, val, iband=0):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        :param iband: List of band indices to crop (index is 0 based)
        """
        self._jim_object.d_imageFrameSet([l, r, t, b, u, d], val, iband)

    def imageFrameAdd(self, l, r, t, b, u, d, val, iband=0):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        :param iband: List of band indices to crop (index is 0 based)
        """
        self._jim_object.d_imageFrameAdd([l, r, t, b, u, d], val, iband)

    def magnify(self, n):
        """Magnify the image.

        Modifies the instance on which the method was called.

        :param n: a positive integer for magnifying size by pixel replication
        """
        self._jim_object._set(self._jim_object.imageMagnify(n))

    def plotLine(self, x1, y1, x2, y2, val):
        """Draw a line from [x1, y1] to [x2, y2] by setting pixels to value val

        Modifies the instance on which the method was called.

        :param x1: an integer for x-coordinate of 1st point
        :param y1: an integer for y-coordinate of 1st point
        :param x2: an integer for x-coordinate of 2nd point
        :param y2: an integer for y-coordinate of 2nd point
        """
        self._jim_object.d_plotLine(x1, y1, x2, y2, val)

    # TODO: getboundingbox, Histograms?, how to work with
    # compose (jim1.compose(jim2, jim3, jim4, 2)), how to ovlmatrix, how to
    # grid (jim1.gridding(jim1, jim3, jim4, 1)), how and where nni,
    # where addframeboxelem, where subframeboxelem,
>>>>>>> 94c7b09cdc7f901c1ba408900b22b3ce73da5779
