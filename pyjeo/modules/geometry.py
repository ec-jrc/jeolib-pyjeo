"""Module for operations working with the geometry of the Jim objects."""

import pyjeo as _pj


def image2geo(jim_object, i, j):
    """ Convert image coordinates to georeferenced.

    Convert image coordinates (column and row) to georeferenced
    coordinates (x and y)

    :param i: image column number (starting from 0)
    :param j: image row number (starting from 0)
    :return: georeferenced coordinates according to the object spatial
        reference system
    """
    return jim_object._jipjim.image2geo(i, j)


def geo2image(jim_object, x, y):
    """ Convert georeferenced coordinates to image.

    Convert georeferenced coordinates (column and row) to image
    coordinates (x and y).

    :param x: georeferenced coordinate in x according to the object spatial
        reference system
    :param y: georeferenced coordinate in y according to the object spatial
        reference system
    :return: image coordinates (row and column, starting from 0)
    """
    coords = jim_object._jipjim.geo2image(x, y)
    return [int(coords[0]), int(coords[1])]


def crop(jim_object, ulx=None, uly=None, ulz=None, lrx=None, lry=None,
         lrz=None, dx=None, dy=None, nogeo=False, **kwargs):
    """Subset raster dataset.

    Subset raster dataset according in spatial (subset region) domain

    :param ulx: Upper left x value of bounding box to crop
    :param uly: Upper left y value of bounding box to crop
    :param ulz: Upper left z value of bounding box to crop
    :param lrx: Lower right x value of bounding box to crop
    :param lry: Lower right y value of bounding box to crop
    :param lrz: Lower right y value of bounding box to crop
    :param dx: spatial resolution in x to crop (stride if nogeo is True)
    :param dy: spatial resolution in y to crop (stride if nogeo is True)
    :param nogeo: use image coordinates if True, default is spatial reference
        system coordinates
    """
    if ulz is not None or lrz is not None:
        assert len(kwargs) == 0, 'It is not supported to use both z coords ' \
                                 'and special cropping parameters'
        jim = _pj.Jim(jim_object)
        jim.geometry.crop(ulx, uly, ulz, lrx, lry, lrz, dx, dy, nogeo,
                          **kwargs)
        return jim
    elif len(kwargs) == 0:
        if nogeo:
            if ulx is None:
                ulx = 0
            if uly is None:
                uly = 0
            if lrx is None:
                lrx = jim_object.properties.nrOfCol()-1
            if lry is None:
                lry = 0
            if dx is None:
                dx = 1
        else:
            if ulx is None:
                ulx = jim_object.properties.getUlx()
            if uly is None:
                uly = jim_object.properties.getUly()
            if lrx is None:
                lrx = jim_object.properties.getLrx()
            if lry is None:
                lry = jim_object.properties.getLry()
            if dx is None:
                dx = jim_object.properties.getDeltaX()
            if dy is None:
                dy = jim_object.properties.getDeltaY()

        kwargs.update({'ulx': ulx})
        kwargs.update({'uly': uly})
        kwargs.update({'lrx': lrx})
        kwargs.update({'lry': lry})
        kwargs.update({'dx': dx})
        kwargs.update({'dy': dy})
        kwargs.update({'nogeo': nogeo})
        return _pj.Jim(jim_object._jipjim.crop(kwargs))

    else:
        if nogeo:
            uli = ulx
            ulj = uly
            lri = lrx
            lrj = lry
            upperLeft = jim_object._jipjim.image2geo(ulx, uly)
            lowerRight = jim_object._jipjim.image2geo(lrx, lry)
            ulx = upperLeft[0]
            uly = upperLeft[1]
            lrx = lowerRight[0]+jim_object.properties.getDeltaX()/2.0
            lry = lowerRight[1]-jim_object.properties.getDeltaY()/2.0
        kwargs.update({'ulx': ulx})
        kwargs.update({'uly': uly})
        kwargs.update({'lrx': lrx})
        kwargs.update({'lry': lry})
        return _pj.Jim(jim_object._jipjim.crop(kwargs))


def cropOgr(jim_object, extent, **kwargs):
    """Subset raster dataset.

    Subset raster dataset in spatial domain defined by a vector dataset.

    :param extent: Get boundary from extent from polygons in vector file
    :param kwargs: See table below

    +------------------+------------------------------------------------------+
    | key              | value                                                |
    +==================+======================================================+
    | ln               | Layer name of extent to crop                         |
    +------------------+------------------------------------------------------+
    | eo               | Special extent options controlling rasterization     |
    +------------------+------------------------------------------------------+
    | crop_to_cutline  | The outside area will be set to no data (the value   |
    |                  | defined by the key 'nodata')                         |
    +------------------+------------------------------------------------------+
    | crop_in_cutline  | The inside area will be set to no data (the value    |
    |                  | defined by the key 'nodata')                         |
    +------------------+------------------------------------------------------+
    | dx               | Output resolution in x (default: keep original       |
    |                  | resolution)                                          |
    +------------------+------------------------------------------------------+
    | dy               | Output resolution in y (default: keep original       |
    |                  | resolution)                                          |
    +------------------+------------------------------------------------------+
    | nodata           | Nodata value to put in image if out of bounds        |
    +------------------+------------------------------------------------------+
    | align            | Align output bounding box to input image             |
    +------------------+------------------------------------------------------+

    .. note::
       Possible values for the key 'eo' are:

       ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

       For instance you can use 'eo':'ATTRIBUTE=fieldname'
    """
    return _pj.Jim(jim_object._jipjim.cropOgr(extent, kwargs))


def cropBand(jim_object, band):
    """Subset raster dataset.

    Subset raster dataset in spectral/temporal domain.

    :param jim_object: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three bands from raster dataset jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim0.cropBand(band=[0,1,2])

    """
    return _pj.Jim(jim_object._jipjim.cropBand({'band': band}))


def cropBandRange(jim_object, startband, endband):
    """Subset raster dataset.

    Subset raster dataset in range of spectral/temporal domain.

    :param jim_object: a Jim object
    :param startband: Start band sequence number (index is 0 based)
    :param endband: End band sequence number (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three bands from raster dataset jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim0.cropBandRange(startband=0,startBand=2)

    """
    return _pj.Jim(jim_object._jipjim.cropBand({'startband': startband,
                                                'endband': endband}))


def stackBand(jim_object, jim_other=None, band=None, **kwargs):
    """Stack bands from another raster dataset to current raster dataset.

    :param jim_object: a Jim object to stack bands
    :param jim_other: a Jim object from which to copy bands
    :param band: List of band indices to stack (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Append all the bands of raster dataset jim1 to image jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        pj.geometry.stackBand(jim0,jim1)

    Append the first three bands of raster dataset jim1 to the image jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        pj.geometry.stackBand(jim0,jim1,band=[0,1,2])
    """
    if isinstance(jim_object, _pj.JimList):
        if band:
            kwargs.update({'band': band})
        return _pj.Jim(jim_object._jipjim.stackBand(kwargs))
    else:
        if band:
            return _pj.Jim(jim_object._jipjim.stackBand(jim_other._jipjim,
                                                        {'band': band}))
        else:
            return _pj.Jim(jim_object._jipjim.stackBand(jim_other._jipjim))


def stackBandRange(jim_object, jim_other, startband, endband):
    """Subset raster dataset.

    Stack range of bands from another raster dataset to current raster dataset.

    :param jim_object: a Jim object to stack bands
    :param jim_other: a Jim object from which to copy bands
    :param startband: Start band sequence number (index is 0 based)
    :param endband: End band sequence number (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Append the first three bands of raster dataset jim1 to image jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        pj.geometry.stackBand(jim0,jim1,startband=0,endband=0)
    """
    return _pj.Jim(jim_object._jipjim.stackBand(jim_other._jipjim,
                                                {'startband': startband,
                                                 'endband': endband}))


def warp(jim_object, t_srs, **kwargs):
    """
    Warp a raster dataset to a target spatial reference system

    :param jim_object: a Jim object
    :param t_srs: Target spatial reference system
    :param kwargs: See table below
    :return: Cropped subimage as Jim instance

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

    Read a raster dataset from disk and warp to the target spatial reference
    system::

        jim = pj.Jim('/path/to/file.tif')
        jim.warp('epsg:3035')

    Read a raster dataset from disk that is in lat lon (epsg:4326), select
    a bounding box in a different spatial reference system (epsg:3035).
    Notice the raster dataset read is still in the original
    projection (epsg:4326). Then warp the raster dataset to the target
    spatial reference system (epsg:3035)::

        jim = pj.Jim('/path/to/file.tif', t_srs='epsg:3035', ulx=1000000, uly=4000000, lrx=1500000, lry=3500000)
        jim.warp('epsg:3035', s_srs='epsg:4326')

    """
    kwargs.update({'t_srs': t_srs})
    # return _pj.Jim(jim_object._set(self._jim_object._jipjim.warp(kwargs)))
    return _pj.Jim(jim_object._jipjim.warp(kwargs))


def imageInsert(jim_object, sec_jim_object, x, y, z, band=0):
    """Merge Jim instance with values of sec_jim_object in given coords.

    :param jim_object: a Jim object
    :param jim_object: a Jim object
    :param x: x coordinate of 1st pixel
    :param y: y coordinate of 1st pixel
    :param z: z coordinate of 1st pixel
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.imageInsert(sec_jim_object._jipjim,
                                                  x, y, z, band))


def imageInsertCompose(jim_object, imlbl, im2, x, y, z, val, band=0):
    """Merge Jim instance with values of im2 if val of imlbl == val.

    :param jim_object: a Jim object
    :param imRaster_imlbl: a Jim object
    :param imRaster_im2: a Jim object
    :param x: x coordinate of 1st pixel
    :param y: y coordinate of 1st pixel
    :param z: z coordinate of 1st pixel
    :param val: integer for label value
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.imageInsertCompose(imlbl._jipjim,
                                                         im2._jipjim,
                                                         x, y, z, val, band))


def imageFrameSet(jim_object, l, r, t, b, u, d, val, band=0):
    """Set the values of the image frame to value val.

    :param jim_object: a Jim object
    :param l: width of left frame
    :param r: width of right frame
    :param t: width of top frame
    :param b: width of bottom frame
    :param u: width of upper frame
    :param d: width of lower frame
    :param val: value of frame
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.imageFrameSet([l, r, t, b, u, d], val,
                                                    band))


def imageFrameAdd(jim_object, l, r, t, b, u, d, val, band=0):
    """Set the values of the image frame to value val.

    :param jim_object: a Jim object
    :param l: width of left frame
    :param r: width of right frame
    :param t: width of top frame
    :param b: width of bottom frame
    :param u: width of upper frame
    :param d: width of lower frame
    :param val: value of frame
    :param band: List of band indices to crop (index is 0 based)
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.imageFrameAdd([l, r, t, b, u, d], val,
                                                    band))


def magnify(jim_object, n):
    """Magnify the image.

    :param jim_object: a Jim object
    :param n: a positive integer for magnifying size by pixel replication
    :return: a Jim object containing the magnified image
    """
    return _pj.Jim(jim_object._jipjim.imageMagnify(n))


def plotLine(jim_object, x1, y1, x2, y2, val):
    """Draw a line from [x1, y1] to [x2, y2] by setting pixels of Jim to val

    :param jim_object: a Jim object
    :param x1: an integer for x-coordinate of 1st point
    :param y1: an integer for y-coordinate of 1st point
    :param x2: an integer for x-coordinate of 2nd point
    :param y2: an integer for y-coordinate of 2nd point
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.plotLine(x1, y1, x2, y2, val))


def polygonize(jim_object, output, **kwargs):
    """Polygonize Jim object based on GDALPolygonize


    :param jim_object: a Jim object
    :param output: output filename of JimVect object that is returned. Use  /vsimem for in memory vectors
    :param kwargs: See table below
    :return: JimVect object with polygons

    +------------------+------------------------------------------------------+
    | key              | value                                                |
    +==================+======================================================+
    | ln               | Output layer name                                    |
    +------------------+------------------------------------------------------+
    | oformat          | Output vector dataset format                         |
    +------------------+------------------------------------------------------+
    | co               | Creation option for output vector dataset            |
    +------------------+------------------------------------------------------+
    | name             | Field name of the output layer (default is DN)       |
    +------------------+------------------------------------------------------+
    | nodata           | Disgard this nodata value when creating polygons     |
    +------------------+------------------------------------------------------+
    """
    kwargs.update({'output': output})
    if isinstance(jim_object, _pj.Jim):
        avect=jim_object._jipjim.polygonize(kwargs)
        pjvect=_pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only intersect with Jim object')


def rasterize(jim_object, jim_vect, burnValue=1,eo=['ALL_TOUCHED'],ln=None):
    """Rasterize Jim object based on GDALRasterizeLayersBuf

    :param jim_object: a template Jim object
    :param jim_vect: JimVect object that needs to be polygonized
    :param burnValue: burn value
    :param eo: option (default is ALL_TOUCHED)
    :param ln: layer names (optional)
    :return: rasterized Jim object

    .. note::
       Possible values for the key 'eo' are:

       ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

       For instance you can use 'eo':'ATTRIBUTE=fieldname'
    """
    if not isinstance(jim_object, _pj.Jim):
        raise TypeError('Error: template must be a Jim object')
    if not isinstance(jim_vect, _pj.JimVect):
        raise TypeError('Error: can only rasterize a JimVect')

    ajim=_pj.Jim(jim_object)
    ajim._jipjim.d_rasterizeBuf(item._jipjimvect,burnValue,eo,ln)
    return ajim

class _Geometry():
    """Define all Geometry methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def band2plane(self):
        """Convert 2-dimensional multi-band image to 3-dimensional single band multi-plane image)

        """
        self._jim_object._jipjim.band2plane()

    def crop(self, ulx=None, uly=None, ulz=None, lrx=None, lry=None,
             lrz=None, dx=None, dy=None, nogeo=False, **kwargs):
        """Subset raster dataset.

        Subset raster dataset according in spatial (subset region) domain

        :param ulx: Upper left x value of bounding box to crop
        :param uly: Upper left y value of bounding box to crop
        :param ulz: Upper left z value of bounding box to crop
        :param lrx: Lower right x value of bounding box to crop
        :param lry: Lower right y value of bounding box to crop
        :param lrz: Lower right y value of bounding box to crop
        :param dx: spatial resolution in x to crop (stride if geo is False)
        :param dy: spatial resolution in y to crop (stride if geo is False)
        :param nogeo: use image coordinates if True, default is spatial
            reference system coordinates
        """
        if ulz is not None or lrz is not None:
            assert len(kwargs) == 0, \
                'It is not supported to use both z coords and special ' \
                'cropping parameters'
            gt = self._jim_object.properties.getGeoTransform()
            nr_of_cols = self._jim_object.properties.nrOfCol()
            nr_of_rows = self._jim_object.properties.nrOfRow()
            if nogeo:
                uli = ulx
                ulj = uly
                lri = lrx
                lrj = lry
                upperLeft = self._jim_object.geometry.image2geo(ulx, uly)
                lowerRight = self._jim_object.geometry.image2geo(lrx, lry)
                ulx = upperLeft[0]
                uly = upperLeft[1]
                # lrx=lowerRight[0]+self._jim_object.properties.getDeltaX()/2.0
                # lry=lowerRight[1]-self._jim_object.properties.getDeltaY()/2.0
            else:
                upperLeftImage = self._jim_object.geometry.geo2image(ulx, uly)
                uli = upperLeftImage[0]
                ulj = upperLeftImage[1]
                lowerRightImage = self._jim_object.geometry.geo2image(lrx, lry)
                lri = lowerRightImage[0]
                lrj = lowerRightImage[1]
            for iband in range(0, self._jim_object.properties.nrOfBand()):
                self._jim_object._jipjim.d_imageFrameSubstract([
                    uli, nr_of_cols - lri,
                    ulj, nr_of_rows - lrj,
                    ulz, self._jim_object.properties.nrOfPlane() - lrz],
                    iband)
            gt[0] = ulx
            gt[3] = uly
            self._jim_object.properties.setGeoTransform(gt)
        elif len(kwargs) == 0:
            if nogeo:
                if ulx is None:
                    ulx = 0
                if uly is None:
                    uly = 0
                if lrx is None:
                    lrx = self._jim_object.properties.nrOfCol()-1
                if lry is None:
                    lry = 0
                if dx is None:
                    dx = 1
            else:
                if ulx is None:
                    ulx = self._jim_object.properties.getUlx()
                if uly is None:
                    uly = self._jim_object.properties.getUly()
                if lrx is None:
                    lrx = self._jim_object.properties.getLrx()
                if lry is None:
                    lry = self._jim_object.properties.getLry()
                if dx is None:
                    dx = self._jim_object.properties.getDeltaX()
                if dy is None:
                    dy = self._jim_object.properties.getDeltaY()
            kwargs.update({'ulx': ulx})
            kwargs.update({'uly': uly})
            kwargs.update({'lrx': lrx})
            kwargs.update({'lry': lry})
            kwargs.update({'dx': dx})
            kwargs.update({'dy': dy})
            kwargs.update({'nogeo': nogeo})
            self._jim_object._set(self._jim_object._jipjim.crop(kwargs))
            # self._jim_object._set(self._jim_object.crop(ulx, uly, lrx, lry,
            #                                             dx, dy, geo))
        else:
            if nogeo:
                uli = ulx
                ulj = uly
                lri = lrx
                lrj = lry
                upperLeft = self._jim_object.geometry.image2geo(ulx, uly)
                lowerRight = self._jim_object.geometry.image2geo(lrx, lry)
                ulx = upperLeft[0]
                uly = upperLeft[1]
                lrx = lowerRight[0]+self._jim_object.properties.getDeltaX()/2.0
                lry = lowerRight[1]-self._jim_object.properties.getDeltaY()/2.0
            kwargs.update({'ulx': ulx})
            kwargs.update({'uly': uly})
            kwargs.update({'lrx': lrx})
            kwargs.update({'lry': lry})
            self._jim_object._set(self._jim_object._jipjim.crop(kwargs))
            # return _pj.Jim(self._jim_object.crop(kwargs))

    def cropOgr(self, extent, **kwargs):
        """Subset raster dataset.

        Subset raster dataset in spatial domain defined by a vector dataset.

        Modifies the instance on which the method was called.

        :param extent: Get boundary from extent from polygons in vector file
        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | ln               | Layer name of extent to crop                     |
        +------------------+--------------------------------------------------+
        | eo               | Special extent options controlling rasterization |
        +------------------+--------------------------------------------------+
        | crop_to_cutline  | The outside area will be set to no data (the     |
        |                  | value defined by the key 'nodata')               |
        +------------------+--------------------------------------------------+
        | crop_in_cutline  | The inside area will be set to no data (the      |
        |                  | value defined by the key 'nodata')               |
        +------------------+--------------------------------------------------+
        | dx               | Output resolution in x (default: keep original   |
        |                  | resolution)                                      |
        +------------------+--------------------------------------------------+
        | dy               | Output resolution in y (default: keep original   |
        |                  | resolution)                                      |
        +------------------+--------------------------------------------------+
        | nodata           | Nodata value to put in image if out of bounds    |
        +------------------+--------------------------------------------------+
        | align            | Align output bounding box to input image         |
        +------------------+--------------------------------------------------+

        .. note::
           Possible values for the key 'eo' are:

           ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

           For instance you can use 'eo':'ATTRIBUTE=fieldname'
        """
        self._jim_object._set(self._jim_object._jipjim.cropOgr(extent, kwargs))

    def cropBand(self, band):
        """Subset raster dataset.

        Subset raster dataset in spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim0.cropBand(band=[0,1,2])

        """
        self._jim_object._jipjim.d_cropBand({'band': band})

    def cropBandRange(self, startband, endband):
        """Subset raster dataset.

        Subset raster dataset in range of spectral/temporal domain.

        Modifies the instance on which the method was called.

        :param startband: Start band sequence number (index is 0 based)
        :param endband: End band sequence number (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim0.cropBandRange(startband=0,startBand=2)

        """
        self._jim_object._jipjim.d_cropBand({'startband': startband,
                                             'endband': endband})

    def stackBand(self, jim_other, band=None):
        """Stack bands from another raster dataset to current raster dataset.

        Modifies the instance on which the method was called.

        :param jim_other: a Jim object from which to copy bands
        :param band: List of band indices to stack (index is 0 based)

        Example:

        Append all the bands of raster dataset jim1 to the current image jim0::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim0.stackBand(jim1)

        Append the first three bands of raster dataset jim1 to the current
        image jim0::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim0.stackBand(jim1, band=[0, 1, 2])
        """
        if band:
            self._jim_object._jipjim.d_stackBand(jim_other._jipjim,
                                                 {'band': band})
        else:
            self._jim_object._jipjim.d_stackBand(jim_other._jipjim)

    def stackBandRange(self, jim_other, startband, endband):
        """Subset raster dataset.

        Stack range of bands from another raster dataset to current raster
        dataset.

        Modifies the instance on which the method was called.

        :param jim_other: a Jim object from which to copy bands
        :param startband: Start band sequence number (index is 0 based)
        :param endband: End band sequence number (index is 0 based)

        Example:

        Append the first three bands of raster dataset jim1 to the current
        image jim0::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim0.stackBandRange(jim1, startband=0, endband=2)
        """
        self._jim_object._jipjim.d_stackBand(jim_other._jipjim,
                                             {'startband': startband,
                                              'endband': endband})

    def extractOgr(self, jim_ref, rule, output, **kwargs):
        """Extract pixel values from raster image based on a vector dataset.

        :param jim_ref: reference JimVect instance
        :rule: Rule how to calculate zonal statistics per feature (see list of :ref:`supported rules <extract_rules>`)
        :output: Name of the output vector dataset in which the zonal
            statistics will be saved
        :param kwargs: See table below
        :return: A VectorOgr with the same geometry as the sample vector
            dataset and an extra field for each of the calculated raster value
            (zonal) statistics. The same layer name(s) of the sample will be
            used for the output vector dataset


        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | copy             | Copy these fields from the sample vector dataset |
        |                  | (default is to copy all fields)                  |
        +------------------+--------------------------------------------------+
        | label            | Create extra field named 'label' with this value |
        +------------------+--------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field   |
        |                  | identifier (sequence of features)                |
        +------------------+--------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is |
        |                  | to use extract all bands                         |
        +------------------+--------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands |
        |                  | to extract                                       |
        +------------------+--------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)           |
        +------------------+--------------------------------------------------+
        | endband          | End band sequence number (0 indexed)             |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+

        .. _extract_rules:

        :Supported rules to extract:

        +------------------+--------------------------------------------------+
        | rule             | description                                      |
        +==================+==================================================+
        | point            | extract a single pixel within the polygon or on  |
        |                  | each point feature                               |
        +------------------+--------------------------------------------------+
        | allpoints        | Extract all pixel values covered by the polygon  |
        +------------------+--------------------------------------------------+
        | centroid         | Extract pixel value at the centroid of           |
        |                  | the polygon                                      |
        +------------------+--------------------------------------------------+
        | mean             | Extract average of all pixel values within the   |
        |                  | polygon                                          |
        +------------------+--------------------------------------------------+
        | stdev            | Extract standard deviation of all pixel values   |
        |                  | within the polygon                               |
        +------------------+--------------------------------------------------+
        | median           | Extract median of all pixel values within        |
        |                  | the polygon                                      |
        +------------------+--------------------------------------------------+
        | min              | Extract minimum value of all pixels within       |
        |                  | the polygon                                      |
        +------------------+--------------------------------------------------+
        | max              | Extract maximum value of all pixels within       |
        |                  | the polygon                                      |
        +------------------+--------------------------------------------------+
        | sum              | Extract sum of the values of all pixels within   |
        |                  | the polygon                                      |
        +------------------+--------------------------------------------------+
        | mode             | Extract the mode of classes within the polygon   |
        |                  | (classes must be set with the option class)      |
        +------------------+--------------------------------------------------+
        | proportion       | Extract proportion of class(es) within           |
        |                  | the polygon                                      |
        |                  | (classes must be set with the option class)      |
        +------------------+--------------------------------------------------+
        | count            | Extract count of class(es) within the polygon    |
        |                  | (classes must be set with the option class)      |
        +------------------+--------------------------------------------------+
        | percentile       | Extract percentile as defined by option perc     |
        |                  | (e.g, 95th percentile of values covered by       |
        |                  | polygon)                                         |
        +------------------+--------------------------------------------------+


        .. note::
            To ignore some pixels from the extraction process, see list
            of :ref:`mask <extract_mask>` key values:

        .. _extract_mask:

        :Supported key values to mask pixels that must be ignored in the extraction process:

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | srcnodata        | List of nodata values not to extract             |
        +------------------+--------------------------------------------------+
        | bndnodata        | List of band in input image to check if pixel is |
        |                  | valid (used for srcnodata)                       |
        +------------------+--------------------------------------------------+
        | mask             | Use the the specified file as a validity mask    |
        +------------------+--------------------------------------------------+
        | mskband          | Use the the specified band of the mask file      |
        |                  | defined                                          |
        +------------------+------------------------------------------------  +
        | msknodata        | List of mask values not to extract               |
        +------------------+--------------------------------------------------+
        | threshold        | Maximum number of features to extract. Use       |
        |                  | percentage value as string                       |
        |                  | (e.g., '10%') or integer value for absolute      |
        |                  | threshold                                        |
        +------------------+--------------------------------------------------+

        Example:

        Extract a random sample of 100 points, calculating the mean value based
        on a 3x3 window (buffer value of 1 pixel neighborhood) in a vector
        dataset in memory::

            reference = pj.JimVect('/path/to/reference.sqlite')
            jim0 = pj.Jim('/path/to/raster.tif')
            v = jim0.extractOgr(reference, rule=['mean'], output='/path/to/output.sqlite', oformat='SQLite')
            v.write()

        """

        kwargs.update({'output': output})
        kwargs.update({'rule': rule})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']

        return _pj.JimVect(self._jim_object._jipjim.extractOgr(jim_ref._jim_vect, kwargs))

    def extractSample(self, output, **kwargs):
        """Extract a random or grid sample from a raster dataset.

        :output: Name of the output vector dataset in which the zonal
            statistics will be saved
        :param kwargs: See table below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | random           | Extract a random sample with a size equal to     |
        |                  | the defined value                                |
        +------------------+--------------------------------------------------+
        | grid             | Extract a grid sample with a grid size equal to  |
        |                  | the defined value                                |
        +------------------+--------------------------------------------------+
        | rule             | Rule how to calculate zonal statistics per       |
        |                  | feature                                          |
        |                  | (see list of                                     |
        |                  | :ref:`supported rules <extract_rules>`)          |
        +------------------+--------------------------------------------------+
        | buffer           | Buffer for calculating statistics for point      |
        |                  | features (in number of pixels)                   |
        +------------------+--------------------------------------------------+
        | label            | Create extra field named 'label' with this value |
        +------------------+--------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field   |
        |                  | identifier (sequence of features)                |
        +------------------+--------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is |
        |                  | to use extract all bands                         |
        +------------------+--------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands |
        |                  | to extract                                       |
        +------------------+--------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)           |
        +------------------+--------------------------------------------------+
        | endband          | End band sequence number (0 indexed)             |
        +------------------+--------------------------------------------------+
        | ln               | Layer name of output vector dataset              |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list
            of :ref:`mask <extract_mask>` key values:

        Example:

        Extract a random sample of 100 points, calculating the mean value based
        on a 3x3 window (buffer value of 1 pixel neighborhood) in a vector
        dataset in memory::

            v01 = jim0.extractSample(random=100, buffer=1, rule=['mean'], output='mem01', oformat='Memory')

        Extract a sample of 100 points using a regular grid sampling scheme.
        For each grid point, calculate the median value based on a 3x3 window
        (buffer value of 1 pixel neighborhood). Write the result in a SQLite
        vector dataset on disk::

            outputfn = '/path/to/output.sqlite'
            npoint = 100
            gridsize = int(jim.nrOfCol()*jim.getDeltaX()/math.sqrt(npoint))
            v = jim.extractSample(grid=gridsize, buffer=1, rule=['median'], output=outputfn, oformat='SQLite')
            v.write()

        """
        kwargs.update({'output': output})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']
        return self._jim_object._jipjim.extractSample(kwargs)

    def extractImg(self, reference, output, **kwargs):
        """Extract pixel values from an input based on a raster sample dataset.

        :param reference: thematic raster dataset with integer values,
            typically a land cover map
        :param kwargs: See table below
        :return: A VectorOgr with fields for each of the calculated raster
            value (zonal) statistics

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | rule             | Rule how to calculate zonal statistics per       |
        |                  | feature                                          |
        |                  | (see list of                                     |
        |                  | :ref:`supported rules <extract_rules>`)          |
        +------------------+--------------------------------------------------+
        | class            | List of classes to extract from the raster sample|
        |                  | dataset.                                         |
        |                  | Leave empty to extract all valid data pixels from|
        |                  | thee sample                                      |
        +------------------+--------------------------------------------------+
        | cname            | Name of the class label in the output vector     |
        |                  | dataset (default is 'label')                     |
        +------------------+--------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field   |
        |                  | identifier (sequence of features)                |
        +------------------+--------------------------------------------------+
        | band             | List of bands to extract (0 indexed). Default is |
        |                  | to use extract all bands                         |
        +------------------+--------------------------------------------------+
        | bandname         | List of band name corresponding to list of bands |
        |                  | to extract                                       |
        +------------------+--------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)           |
        +------------------+--------------------------------------------------+
        | endband          | End band sequence number (0 indexed)             |
        +------------------+--------------------------------------------------+
        | down             | Down sampling factor to extract a subset of      |
        |                  | the sample based on a grid                       |
        +------------------+--------------------------------------------------+
        | ln               | Layer name of output vector dataset              |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+

        .. note::
            To ignore some pixels from the extraction process, see list
            of :ref:`nodata <extract_nodata>` key values:

        .. _extract_nodata:

        :Supported key values to mask pixels that must be ignored in the extraction process:

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | srcnodata        | List of nodata values not to extract             |
        +------------------+--------------------------------------------------+
        | bndnodata        | List of band in input image to check if pixel is |
        |                  | valid (used for srcnodata)                       |
        +------------------+--------------------------------------------------+
        | threshold        | Maximum number of features to extract. Use       |
        |                  | percentage value as string                       |
        |                  | (e.g., '10%') or integer value for absolute      |
        |                  | threshold.                                       |
        |                  | You can provide a list of threshold values, one  |
        |                  | for each class.                                  |
        +------------------+--------------------------------------------------+

        Example:

        Open a raster sample dataset based on land cover map (e.g., Corine) and
        use it to extract a stratified sample of 100 points from an input
        raster dataset with four spectral bands ('B02', 'B03', 'B04', 'B08').
        Only sample classes 2 (urban), 12 (agriculture), 25 (forest),
        41 (water) and an aggregated (rest) class 50::

            jim_ref=pj.Jim('/path/to/landcovermap.tif')

            classes=[2,12,25,41,50]
            thresholds=['20%','25%','25%','10%','5%']

            jim_ref=pj.Jim('/path/to/multiband.tif','dx'=jim.getDeltaX(),'dy'=jim.getDeltaY(),'ulx'=jim.getUlx(),'uly'=jim.getUly(),'lrx'=jim.getLrx(),'lry'=jim.getLry())

            outputfn='/path/to/output.sqlite'
            sample=jim.extractImg(jim_ref,srcnodata=[0],output=outputfn,class=classes,threshold=thresholds,bandname=['B02','B03','B04','B08'],band=[0,1,2,3])
        """
        kwargs.update({'output': output})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']
        return self._jim_object._jipjim.extractImg(reference._jipjim, kwargs)

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

        Read a raster dataset from disk and warp to the target spatial
        reference system::

            jim=pj.Jim('/path/to/file.tif')
            jim.warp('epsg:3035')

        Read a raster dataset from disk that is in lat lon (epsg:4326), select
        a bounding box in a different spatial reference system (epsg:3035).
        Notice the raster dataset read is still in the original projection
        (epsg:4326). Then warp the raster dataset to the target spatial
        reference system (epsg:3035)::

            jim=pj.Jim('/path/to/file.tif',t_srs='epsg:3035',ulx=1000000,uly=4000000,lrx=1500000,lry=3500000)
            jim.warp('epsg:3035',s_srs='epsg:4326')

        """
        kwargs.update({'t_srs': t_srs})
        self._jim_object._set(self._jim_object._jipjim.warp(kwargs))

    def imageInsert(self, sec_jim_object, x, y, z, band=0):
        """Merge Jim instance with values of sec_jim_object in given coords.

        Modifies the instance on which the method was called.

        :param jim_object: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_imageInsert(sec_jim_object._jipjim, x, y, z,
                                               band)

    def imageInsertCompose(self, imlbl, im2, x, y, z, val, band=0):
        """Merge Jim instance with values of im2 if val of imlbl == val.

        Modifies the instance on which the method was called.

        :param imRaster_imlbl: a Jim object
        :param imRaster_im2: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param val: integer for label value
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_imageInsertCompose(imlbl._jipjim,
                                                      im2._jipjim,
                                                      x, y, z, val, band)

    def imageFrameSet(self, l, r, t, b, u, d, val, band=0):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_imageFrameSet([l, r, t, b, u, d], val, band)

    def imageFrameAdd(self, l, r, t, b, u, d, val, band=0):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        :param band: List of band indices to crop (index is 0 based)
        """
        self._jim_object._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val, band)

    def magnify(self, n):
        """Magnify the image.

        Modifies the instance on which the method was called.

        :param n: a positive integer for magnifying size by pixel replication
        """
        self._jim_object._set(self._jim_object._jipjim.imageMagnify(n))

    def plotLine(self, x1, y1, x2, y2, val):
        """Draw a line from [x1, y1] to [x2, y2] by setting pixels to value val

        Modifies the instance on which the method was called.

        :param x1: an integer for x-coordinate of 1st point
        :param y1: an integer for y-coordinate of 1st point
        :param x2: an integer for x-coordinate of 2nd point
        :param y2: an integer for y-coordinate of 2nd point
        """
        self._jim_object._jipjim.d_plotLine(x1, y1, x2, y2, val)

    def image2geo(self, i, j):
        """Convert image coordinates to georeferenced coordinates.

        Convert image coordinates (column and row) to georeferenced
        coordinates (x and y)

        :param i: image column number (starting from 0)
        :param j: image row number (starting from 0)
        :return: georeferenced coordinates according to the object spatial
            reference system
        """
        return self._jim_object._jipjim.image2geo(i, j)

    def geo2image(self, x, y):
        """ Convert georeferenced coordinates to image coordinates.

        Convert georeferenced coordinates (x and y) to image
        coordinates (column and row).

        :param x: georeferenced coordinate in x according to the object spatial
            reference system
        :param y: georeferenced coordinate in y according to the object spatial
            reference system
        :return: image coordinates (row and column, starting from 0)
        """
        coord = self._jim_object._jipjim.geo2image(x, y)
        return [int(coord[0]), int(coord[1])]

    def polygonize(self, output, **kwargs):
        """Polygonize Jim object based on GDALPolygonize


        :param output: output filename of JimVect object that is returned. Use  /vsimem for in memory vectors
        :param kwargs: See table below
        :return: JimVect object with polygons

        +------------------+------------------------------------------------------+
        | key              | value                                                |
        +==================+======================================================+
        | ln               | Output layer name                                    |
        +------------------+------------------------------------------------------+
        | oformat          | Output vector dataset format                         |
        +------------------+------------------------------------------------------+
        | co               | Creation option for output vector dataset            |
        +------------------+------------------------------------------------------+
        | name             | Field name of the output layer (default is DN)       |
        +------------------+------------------------------------------------------+
        | nodata           | Disgard this nodata value when creating polygons     |
        +------------------+------------------------------------------------------+
        """
        kwargs.update({'output': output})
        avect=self._jim_object._jipjim.polygonize(kwargs)
        pjvect=_pj.JimVect()
        pjvect._set(avect)
        return pjvect

    def rasterize(self, jim_vect, burnValue=1,eo=['ALL_TOUCHED'],ln=None):
        """Rasterize Jim object based on GDALRasterizeLayersBuf

    CPLErr Jim::rasterizeBuf(VectorOgr& ogrReader, double burnValue, const std::vector<std::string>& eoption, const std::vector<std::string>& layernames ){

        :param jim_vect: JimVect object that needs to be polygonized
        :param burnValue: burn value
        :param eo: option (default is ALL_TOUCHED)
        :param ln: layer names (optional)

        .. note::
          Possible values for the key 'eo' are:

          ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

          For instance you can use 'eo':'ATTRIBUTE=fieldname'
        """

        # print("type of jim_vect {}:".format(type(jim_vect)))
        # if not isinstance(jim_vect, _pj.JimVect):
        #     raise TypeError('Error: can only rasterize JimVect')

        # self._jim_object._jipjim.d_rasterizeBuf(jim_vect._jipjimvect,float(burnValue),eo,ln)
        self._jim_object._jipjim.d_rasterizeBuf(jim_vect._jipjimvect,float(burnValue))

class _GeometryList():
    """Define all Geometry methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def stackBand(self, **kwargs):
        return _pj.Jim(self._jim_list._jipjimlist.stackBand(kwargs))

    def extractOgr(self, sample, rule, output, **kwargs):
        if isinstance(sample, _pj.JimVect):
            kwargs.update({'rule': rule})
            kwargs.update({'output': output})
            if 'threshold' in kwargs:
                if '%' in kwargs['threshold']:
                    kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
                else:
                    kwargs['threshold'] = -kwargs['threshold']
            avect=self._jim_list._jipjimlist.extractOgr(sample._jipjimvect,kwargs)
            pjvect=_pj.JimVect()
            pjvect._set(avect)
            return pjvect
            # return _pj.JimVect(self._jim_list._jipjimlist.extractOgr(sample._jipjimvect,kwargs))
        else:
            raise TypeError('Error: extractOgr must operate on vector sample of type JimVect')

class _GeometryVect():
    """Define all Geometry methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    def intersect(self, jim, output, **kwargs):
        """Intersect JimVect object with Jim object and return only those features with an intersect

        :param jim: Jim object to intersect
        :param output: output filename of JimVect object that is returned. Use  /vsimem for in memory vectors
        :param kwargs: See table below

        +------------------+------------------------------------------------------+
        | key              | value                                                |
        +==================+======================================================+
        | oformat          | Output vector dataset format                         |
        +------------------+------------------------------------------------------+
        | co               | Creation option for output vector dataset            |
        +------------------+------------------------------------------------------+
        """
        kwargs.update({'output': output})
        if isinstance(jim, _pj.Jim):
            avect=self._jim_vect._jipjimvect.intersect(jim._jipjim,kwargs)
            pjvect=_pj.JimVect()
            pjvect._set(avect)
            #todo: do not return but overwrite self
            return pjvect
            # return _pj.JimVect(self._jim_vect._jipjimvect.intersect(jim._jipjim,kwargs))
        else:
            raise TypeError('Error: can only intersect with Jim object')

    def join(self, jvec, output, **kwargs):
        """Join JimVect object with another JimVect object.

        :param jim: Jim object to intersect
        :param output: output filename of JimVect object that is returned. Use  /vsimem for in memory vectors
        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | key              | Key(s) used to join (default is fid)             |
        +------------------+--------------------------------------------------+
        | method           | Join method: "INNER","OUTER_LEFT","OUTER_RIGHT", |
        |                  | "OUTER_FULL". (default is INNER)                 |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+
        """
        kwargs.update({'output': output})
        if isinstance(jvec, _pj.JimVect):
            avect=self._jim_vect._jipjimvect.join(jvec._jipjimvect,kwargs)
            pjvect=_pj.JimVect()
            pjvect._set(avect)
            #todo: do not return but overwrite self
            return pjvect
            # return _pj.JimVect(self._jim_vect._jipjimvect.join(jvec._jipjimvect,kwargs))
        else:
            raise TypeError('Error: can only join with JimVect object')

    def convexHull(self, output, **kwargs):
        """Create the convex hull on a JimVect object.

        :param output: output filename of JimVect object that is returned. Use  /vsimem for in memory vectors
        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+
        """
        kwargs.update({'output': output})
        print("kwargs: {}".format(kwargs))
        avect=self._jim_vect._jipjimvect.convexHull(kwargs)
        pjvect=_pj.JimVect()
        pjvect._set(avect)
        #todo: do not return but overwrite self
        return pjvect
