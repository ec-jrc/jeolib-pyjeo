class _Geometry():
    def __init__(self, jim_object):
        """Initialize the module.

        :param jim_object: parent Jim object to have access to its attributes
        """
        self._jim_object = jim_object

    def imageFrameSubstract(self, leftSize, rightSize, topSize, belowSize,
                            upSize, downSize):
        """peels off a frame of im

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

    def crop(self, ulx, uly, lrx, lry, **kwargs):
        """Subset raster dataset

        Subset raster dataset according in spatial (subset region) or
        spectral/temporal domain (subset bands).

        Modifies the instance on which the method was called.

        :param ulx: Upper left x value of bounding box to crop
        :param uly: Upper left y value of bounding box to crop
        :param lrx: Lower right x value of bounding box to crop
        :param lry: Lower right y value of bounding box to crop
        :param kwargs: See below

        :keyword arguments:
            :ln: Layer name of extent to crop
            :eo: Special extent options controlling rasterization
            :crop_to_cutline: True will crop the extent of the target dataset
                to the extent of the cutline The outside area will be set to no
                data (the value defined by the key 'nodata')
            :crop_in_cutline: True: inverse operation to crop_to_cutline The
                inside area will be set to no data (the value defined by
                the key 'nodata')
            :dx: Output resolution in x (default: keep original resolution)
            :dy: Output resolution in y (default: keep original resolution)
            :nodata: Nodata value to put in image if out of bounds
            :align: Align output bounding box to input image
        """
        kwargs.update({'ulx': ulx})
        kwargs.update({'uly': uly})
        kwargs.update({'lrx': lrx})
        kwargs.update({'lry': lry})
        self._jim_object._set(self._jim_object.crop(kwargs))

    def cropOgr(self, extent, **kwargs):
        """Subset raster dataset

        Subset raster dataset according in spatial.

        Modifies the instance on which the method was called.

        :param extent: Get boundary from extent from polygons in vector file
        :param kwargs: See below

        :keyword arguments:
            :ln: Layer name of extent to crop
            :eo: Special extent options controlling rasterization
            :crop_to_cutline: True will crop the extent of the target dataset
                to the extent of the cutline The outside area will be set to no
                data (the value defined by the key 'nodata')
            :crop_in_cutline: True: inverse operation to crop_to_cutline The
                inside area will be set to no data (the value defined by
                the key 'nodata')
            :dx: Output resolution in x (default: keep original resolution)
            :dy: Output resolution in y (default: keep original resolution)
            :nodata: Nodata value to put in image if out of bounds
            :align: Align output bounding box to input image
        """
        kwargs.update({'extent': extent})
        self._jim_object._set(self._jim_object.crop(kwargs))

    def cropBand(self, band):
        """Subset raster dataset

        Subset raster dataset according in spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)
        :param kwargs: See below

        :keyword arguments:
            :ln: Layer name of extent to crop
            :eo: Special extent options controlling rasterization
            :crop_to_cutline: True will crop the extent of the target dataset
                to the extent of the cutline The outside area will be set to no
                data (the value defined by the key 'nodata')
            :crop_in_cutline: True: inverse operation to crop_to_cutline The
                inside area will be set to no data (the value defined by
                the key 'nodata')
            :dx: Output resolution in x (default: keep original resolution)
            :dy: Output resolution in y (default: keep original resolution)
            :nodata: Nodata value to put in image if out of bounds
            :align: Align output bounding box to input image
        """
        self._jim_object._set(self._jim_object.crop({'band': band}))

    def cropBandRange(self, startband, endband):
        """Subset raster dataset

        Subset raster dataset according in spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param startband: Start band sequence number (index is 0 based)
        :param endband: End band sequence number (index is 0 based)
        :param kwargs: See below

        :keyword arguments:
            :ln: Layer name of extent to crop
            :eo: Special extent options controlling rasterization
            :crop_to_cutline: True will crop the extent of the target dataset
                to the extent of the cutline The outside area will be set to no
                data (the value defined by the key 'nodata')
            :crop_in_cutline: True: inverse operation to crop_to_cutline The
                inside area will be set to no data (the value defined by
                the key 'nodata')
            :dx: Output resolution in x (default: keep original resolution)
            :dy: Output resolution in y (default: keep original resolution)
            :nodata: Nodata value to put in image if out of bounds
            :align: Align output bounding box to input image
        """
        self._jim_object._set(self._jim_object.crop({'startband': startband,
                                                     'endband': endband}))

    def extractOgr(self, jim_ref, **kwargs):
        """Extract pixel values from raster image using a vector dataset sample

        :param jim_ref: reference Jim instance
        :param kwargs: See below
        :return: A VectorOgr with the same geometry as the sample vector
            dataset and an extra field for each of the calculated raster value
            (zonal) statistics. The same layer name(s) of the sample will be
            used for the output vector dataset

        :keyword arguments:
            :rule: Rule how to calculate zonal statistics per feature
                (point, allpoints, centroid, mean, stdev, median, min, max,
                sum, mode, proportion, count, percentile)
            :copy: Copy these fields from the sample vector dataset (default
                is to copy all fields)
            :label:	Create extra field named 'label' with this value
            :fid: Create extra field named 'fid' with this field identifier
                (sequence of features)
            :band: List of bands to extract (0 indexed). Default is to use
                extract all bands
            :bandname: List of band name corresponding to list of bands to
                extract
            :startband: Start band sequence number (0 indexed)
            :endband: End band sequence number (0 indexed)
            :output: Name of the output vector dataset in which the zonal
                statistics are saved
            :oformat: Output vector dataset format
            :co: Creation option for output vector dataset
            :srcnodata: List of nodata values not to extract
            :bndnodata: List of band in input image to check if pixel is
                valid (used for srcnodata)
            :mask: Use the the specified file as a validity mask
            :mskband: Use the the specified band of the mask file defined
            :msknodata: List of mask values not to extract
            :threshold: Maximum number of features to extract (use positive
                values for percentage value and negative value for absolute
                threshold)
        """
        return self._jim_object.extractOgr(jim_ref, kwargs)

    def extractSample(self, **kwargs):
        """Extract a random or grid sample from raster image

        :param kwargs: See below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        :keyword arguments:
            :rule: Rule how to calculate zonal statistics per feature
                (point, allpoints, centroid, mean, stdev, median, min, max,
                sum, mode, proportion, count, percentile)
            :buffer: Buffer for calculating statistics for point features (in
                number of pixels)
            :label:	Create extra field named 'label' with this value
            :fid: Create extra field named 'fid' with this field identifier
                (sequence of features)
            :band: List of bands to extract (0 indexed). Default is to use
                extract all bands
            :bandname: List of band name corresponding to list of bands to
                extract
            :startband: Start band sequence number (0 indexed)
            :endband: End band sequence number (0 indexed)
            :output: Name of the output vector dataset in which the zonal
                statistics are saved
            :ln: Layer name of output vector dataset
            :oformat: Output vector dataset format
            :co: Creation option for output vector dataset
            :srcnodata: List of nodata values not to extract
            :bndnodata: List of band in input image to check if pixel is
                valid (used for srcnodata)
            :mask: Use the the specified file as a validity mask
            :mskband: Use the the specified band of the mask file defined
            :msknodata: List of mask values not to extract
            :threshold: Maximum number of features to extract (use positive
                values for percentage value and negative value for absolute
                threshold)
        """
        return self._jim_object.extractSample(kwargs)

    def extractImg(self, **kwargs):
        """Extract pixel values from an input based on a raster sample dataset

        :param kwargs: See below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        :keyword arguments:
            :rule: Rule how to calculate zonal statistics per feature
                (point, allpoints, centroid, mean, stdev, median, min, max,
                sum, mode, proportion, count, percentile)
            :class: List of classes to extract from the raster sample dataset.
                Leave empty to extract all valid data pixels from thee sample
            :cname:	Name of the class label in the output vector dataset
                (default is 'label')
            :fid: Create extra field named 'fid' with this field identifier
                (sequence of features)
            :band: List of bands to extract (0 indexed). Default is to use
                extract all bands
            :bandname: List of band name corresponding to list of bands to
                extract
            :startband: Start band sequence number (0 indexed)
            :endband: End band sequence number (0 indexed)
            :down: Down sampling factor to extract a subset of the sample based
                on a grid
            :output: Name of the output vector dataset in which the zonal
                statistics are saved
            :ln: Layer name of output vector dataset
            :oformat: Output vector dataset format
            :co: Creation option for output vector dataset
            :srcnodata: List of nodata values not to extract
            :bndnodata: List of band in input image to check if pixel is
                valid (used for srcnodata)
            :mask: Use the the specified file as a validity mask
            :mskband: Use the the specified band of the mask file defined
            :msknodata: List of mask values not to extract
            :threshold: Maximum number of features to extract (use positive
                values for percentage value and negative value for absolute
                threshold)
        """
        return self._jim_object.extractImg(jim_ref, kwargs)
