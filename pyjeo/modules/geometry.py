"""Module for operations working with the geometry of the Jim objects."""

import pyjeo as _pj
import os
import numpy

from collections import Iterable


# def append(jvec1, jvec2, output, **kwargs):
#     """Append JimVect object with another JimVect object.

#     :param jvec1: first JimVect
#     :param jvec2: second JimVect to append
#     :param output: output filename of JimVect object that is returned.
#         Use /vsimem for in memory vectors
#     :param oformat: output vector dataset format
#     """
#     kwargs.update({'output': output})
#     if isinstance(jvec2, _pj.JimVect):
#         avect = _pj.JimVect(jvec1, **kwargs)
#         avect._jipjimvect.append(jvec2._jipjimvect)
#         return _pj.JimVect(avect)
#     else:
#         raise TypeError('Error: can only join with JimVect object')


def band2plane(jim):
    """Convert 2-dimensional multi-band object to a 3-dimensional.

    The result will be a single band multi-plane object

    Example: convert a multi-band object with 12 bands to a 3-dimensional
    single band object with 12 planes::

        jim2d=pj.Jim('/path/to/multi/band/image.tif')
        jim2d.properties.nrOfBand()
        12
        jim2d.properties.nrOfPlane()
        1
        jim3d=pj.geometry.band2plane(jim2d)
        jim3d.properties.nrOfPlane()
        12
        jim3d.properties.nrOfBand()
        1

    Notice that a multi-band image can also be read directly as
    a multi-plane object::

        jim3d=pj.Jim('/path/to/multi/band/image.tif', band2plane=True)
        jim3d.properties.nrOfBand()
        1
        jim3d.properties.nrOfPlane()
        12

    """
    result = _pj.Jim(jim)
    result.geometry.band2plane()
    return result


def convexHull(jim_vect, output, **kwargs):
    """Create the convex hull on a JimVect object.

    :param jim_vect: JimVect object to be used for the hull generation
    :param output: output filename of JimVect object that is returned.
        Use /vsimem for in memory vectors
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

    avect = jim_vect._jipjimvect.convexHull(kwargs)

    pjvect = _pj.JimVect()
    pjvect._set(avect)

    return pjvect


def crop(jim_object, ulx=None, uly=None, ulz=None, lrx=None, lry=None,
         lrz=None, dx=None, dy=None, nogeo=False, **kwargs):
    """Subset raster dataset.

    Subset raster dataset according in spatial (subset region) domain

    :param jim_object: a Jim object
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

    see :py:meth:`~_Geometry.crop` for an example how to use this function
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
                lry = jim_object.properties.nrOfRow()-1
            if dx is None:
                dx = 1
            if dy is None:
                dy = 1
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

        jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.crop(
            kwargs))
        for iplane in range(1, jim_object.properties.nrOfPlane()):
            jimplane = _pj.Jim(_pj.geometry.cropPlane(
                jim_object, iplane)._jipjim.crop(kwargs))
            jim.geometry.stackPlane(jimplane)
        return jim

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
            lrx = lowerRight[0] + jim_object.properties.getDeltaX() / 2.0
            lry = lowerRight[1] - jim_object.properties.getDeltaY() / 2.0
            if dx is None:
                dx = 1
            if dy is None:
                dy = 1
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

        jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.crop(
            kwargs))
        for iplane in range(1, jim_object.properties.nrOfPlane()):
            jimplane = _pj.Jim(_pj.geometry.cropPlane(
                jim_object, iplane)._jipjim.crop(kwargs))
            jim.geometry.stackPlane(jimplane)
        return jim
        # return _pj.Jim(jim_object._jipjim.crop(kwargs))


def cropBand(jim_object, band):
    """Subset raster dataset.

    Subset raster dataset in spectral domain.

    :param jim_object: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three bands from raster dataset jim::

        jim=pj.Jim('/path/to/raster.tif')
        jim3=pj.geometry.cropBand(jim,band=[0,1,2])

    """
    return _pj.Jim(jim_object._jipjim.cropBand({'band': band}))


def cropOgr(jim_object, extent, **kwargs):
    """Subset raster dataset.

    Subset raster dataset in spatial domain defined by a vector dataset.

    :param jim_object: a Jim object
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
    jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.cropOgr(
        extent._jipjimvect, kwargs))
    for iplane in range(1, jim_object.properties.nrOfPlane()):
        jimplane = _pj.Jim(_pj.geometry.cropPlane(
            jim_object, iplane)._jipjim.cropOgr(extent._jipjimvect, kwargs))
        jim.geometry.stackPlane(jimplane)
    return jim


def cropPlane(jim_object, plane):
    """Subset raster dataset.

    Subset raster dataset in temporal domain.

    :param jim_object: a Jim object
    :param plane: List of plane indices to crop (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three planes from raster dataset jim::

        jim=pj.Jim('/path/to/raster.tif')
        jim3=pj.geometry.cropPlane(jim,plane=[0,1,2])

    """
    return _pj.Jim(jim_object._jipjim.cropPlane({'plane': plane}))


def geo2image(jim_object, x, y):
    """Convert georeferenced coordinates (column and row) to image coordinates.

    :param jim_object: a Jim object
    :param x: georeferenced coordinate in x according to the object spatial
        reference system
    :param y: georeferenced coordinate in y according to the object spatial
        reference system
    :return: image coordinates (row and column, starting from 0)

    Get column and row index (0 based) of some georeferenced coordinates
    x and y (in this case first pixel: 0, 0)::

      jim=pj.Jim('/path/to/raster.tif')
      x=jim.properties.getUlx()
      y=jim.properties.getUly()
      pj.geometry.geo2image(jim,x,y)
    """
    coords = jim_object._jipjim.geo2image(x, y)
    return [int(coords[0]), int(coords[1])]


def image2geo(jim_object, i, j):
    """Convert image coordinates (col and row) to georeferenced coordinates.

    :param jim_object: a Jim object
    :param i: image column number (starting from 0)
    :param j: image row number (starting from 0)
    :return: georeferenced coordinates according to the object spatial
        reference system

    Get upper left corner in georeferenced coordinates
    (in SRS of the Jim object)::

      jim=pj.Jim('/path/to/raster.tif')
      pj.geometry.image2geo(jim,0,0)

    """
    return jim_object._jipjim.image2geo(i, j)


def imageFrameAdd(jim_object, l=0, r=0, t=0, b=0, u=0, d=0, val=0):
    """Add an image frame and set its values to value val.

    :param jim_object: a Jim object
    :param l: width of left frame
    :param r: width of right frame
    :param t: width of top frame
    :param b: width of bottom frame
    :param u: width of upper frame
    :param d: width of lower frame
    :param val: value of frame
    :return: a Jim object
    """
    if jim_object.properties.nrOfBand() > 1:
        returnJim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if returnJim:
                jimband = _pj.geometry.cropBand(jim_object, band=band)
                jimband._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
                returnJim.geometry.stackBand(jimband)
            else:
                returnJim = _pj.geometry.cropBand(jim_object, band=band)
                returnJim._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
        return _pj.Jim(returnJim)
    else:
        return _pj.Jim(jim_object._jipjim.imageFrameAdd(
            [l, r, t, b, u, d], val))


def imageFrameSet(jim_object, l=0, r=0, t=0, b=0, u=0, d=0, val=0, band=None):
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
    if band is None:
        return _pj.Jim(jim_object._jipjim.imageFrameSet(
            [l, r, t, b, u, d], val, -1))
    else:
        return _pj.Jim(jim_object._jipjim.imageFrameSet(
            [l, r, t, b, u, d], val, band))


def imageFrameSubtract(jim_object, l=0, r=0, t=0, b=0, u=0, d=0):
    """Subtract an image frame.

    :param jim_object: a Jim object
    :param l: width of left frame
    :param r: width of right frame
    :param t: width of top frame
    :param b: width of bottom frame
    :param u: width of upper frame
    :param d: width of lower frame
    :return: a Jim object
    """
    if jim_object.properties.nrOfBand() > 1:
        returnJim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if returnJim:
                jimband = _pj.geometry.cropBand(jim_object,
                                                band=band)
                jimband._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
                returnJim.geometry.stackBand(jimband)
            else:
                returnJim = _pj.geometry.cropBand(jim_object,
                                                  band=band)
                returnJim._jipjim.d_imageFrameSubtract(
                    [l, r, t, b, u, d])
        return _pj.Jim(returnJim)
    else:
        return _pj.Jim(jim_object._jipjim.imageFrameSubtract(
            [l, r, t, b, u, d]))


def imageInsert(jim_object, sec_jim_object, x, y, z, band=None):
    """Merge Jim instance with values of sec_jim_object in given coords.

    :param jim_object: a Jim object
    :param sec_jim_object: a Jim object
    :param x: x coordinate of 1st pixel
    :param y: y coordinate of 1st pixel
    :param z: z coordinate of 1st pixel
    :param band: List of band indices to insert (index is 0 based)
    :return: a Jim object
    """
    bands = []
    if band is None:
        bands = range(0, jim_object.properties.nrOfBand())
    else:
        bands.extend(band)

    returnJim = None

    for band in bands:
        if returnJim:
            jimband = jim_object._jipjim.imageInsert(
                sec_jim_object._jipjim, x, y, z, band)
            returnJim.geometry.stackBand(jimband)
        else:
            returnJim = jim_object._jipjim.imageInsert(
                sec_jim_object._jipjim, x, y, z, band)

    return _pj.Jim(returnJim)


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
    bands = []
    if band is None:
        bands = range(0, jim_object.properties.nrOfBand())
    else:
        try:
            bands.extend(band)
        except TypeError:
            bands.append(band)

    returnJim = None

    for band in bands:
        if returnJim:
            jimband = jim_object._jipjim.imageInsertCompose(
                imlbl._jipjim, im2._jipjim, x, y, z, val, band)
            returnJim.geometry.stackBand(jimband)
        else:
            returnJim = jim_object._jipjim.imageInsertCompose(
                imlbl._jipjim, im2._jipjim, x, y, z, val, band)

    return _pj.Jim(returnJim)


def intersect(jvec, jim, output, **kwargs):
    """Intersect JimVect object with Jim object.

    Return only those features with an intersect

    :param jvec: JimVect object to intersect
    :param jim: Jim object with which to intersect
    :param output: output filename of JimVect object that is returned.
        Use /vsimem for in memory vectors
    :param kwargs: See table below
    :return: intersected JimVect object

    +------------------+------------------------------------------------------+
    | key              | value                                                |
    +==================+======================================================+
    | oformat          | Output vector dataset format                         |
    +------------------+------------------------------------------------------+
    | co               | Creation option for output vector dataset            |
    +------------------+------------------------------------------------------+

    Example: intersect a sample with a Jim object::

      jim = pj.Jim('/path/to/raster.tif')
      v = pj.JimVect('/path/to/vector.sqlite')
      sampleintersect = pj.geometry.intersect(
          v, jim, output='/vsimem/intersect', oformat='SQLite',
          co=['OVERWRITE=YES'])
      sampleintersect.io.write('/path/to/output.sqlite')

    """
    kwargs.update({'output': output})
    if isinstance(jim, _pj.Jim):
        avect = jvec._jipjimvect.intersect(jim._jipjim, kwargs)
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only intersect with Jim object')


def join(jvec1, jvec2, output, **kwargs):
    """Join JimVect object with another JimVect object.

    A key field is used to find corresponding features in both objects.

    :param jvec: first JimVect object to join
    :param jvec: second JimVect object to join
    :param output: output filename of JimVect object that is returned.
        Use /vsimem for in memory vectors
    :param kwargs: See table below
    :return: joined JimVect object

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

    .. |inner| image:: figures/join_inner.png
         :width: 20 %
    .. |outer_left| image:: figures/join_outer_left.png
         :width: 20 %
    .. |outer_right| image:: figures/join_outer_right.png
         :width: 20 %
    .. |outer_full| image:: figures/join_outer_full.png
         :width: 20 %

    The join methods currently supported are:

        :INNER |inner|: join two JimVect objects, keeping only those \
                        features for which identical keys in both objects \
                        are found
        :OUTER_LEFT |outer_left|: join two JimVect objects, keeping all \
                                  features from first object
        :OUTER_RIGHT |outer_right|: join two JimVect objects, keeping all \
                                    features from second object
        :OUTER_FULL |outer_full|: join two JimVect objects, keeping all \
                                  features from both objects

    Example: join two vectors, based on the key 'id', which is a common
    field shared between v1 and v2. Use OUTER_FULL as the join method::

      v1 = pj.JimVect('/path/to/vector1.sqlite')
      v2 = pj.JimVect('/path/to/vector2.sqlite')
      v3 = pj.geometry.join(
          v1, v2, '/tmp/test.sqlite', oformat='SQLite',
          co=['OVERWRITE=YES'], key=['id'], method='OUTER_FULL')
    """
    kwargs.update({'output': output})
    if isinstance(jvec1, _pj.JimVect) and isinstance(jvec2, _pj.JimVect):
        avect = jvec1._jipjimvect.join(jvec2._jipjimvect, kwargs)
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only join two JimVect objects')


def magnify(jim_object, n):
    """Magnify the image.

    :param jim_object: a Jim object
    :param n: a positive integer for magnifying size by pixel replication
    :return: a Jim object containing the magnified image
    """
    if jim_object.properties.nrOfBand() > 1:
        returnJim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if returnJim:
                jimband = _pj.geometry.cropBand(jim_object,
                                                band=band)
                jimband = _pj.Jim(jimband._jipjim.imageMagnify(n))

                returnJim.geometry.stackBand(jimband)
            else:
                returnJim = _pj.Jim(_pj.geometry.cropBand(
                    jim_object, band=band)._jipjim.imageMagnify(n))

        return returnJim
    else:
        return _pj.Jim(jim_object._jipjim.imageMagnify(n))


def plane2band(jim):
    """Convert 3-dimensional single-band object to a 2-dimensional
    multi-band object.

    The result will be a multi-band single plane object

    Example: convert a single band object with 12 planes to a 2-dimensional
    multi-band object with 1 plane::

        jim3d=pj.Jim('/path/to/multi/band/image.tif',band2plane=True)
        jim3d.properties.nrOfBand()
        1
        jim3d.properties.nrOfPlane()
        12
        jim2d=pj.geometry.plane2band(jim3d)
        jim2d.properties.nrOfPlane()
        1
        jim2d.properties.nrOfBand()
        12

    """
    result = None
    for iplane in range(0, jim.properties.nrOfPlane()):
        jim_plane = _pj.geometry.cropPlane(jim, iplane)
        if result is None:
            result = jim_plane
        else:
            result.geometry.stackBand(jim_plane)
    return result


def plotLine(jim_object, x1, y1, x2, y2, val):
    """Draw a line from [x1, y1] to [x2, y2] by setting pixels of Jim to val.

    :param jim_object: a Jim object
    :param x1: an integer for x-coordinate of 1st point
    :param y1: an integer for y-coordinate of 1st point
    :param x2: an integer for x-coordinate of 2nd point
    :param y2: an integer for y-coordinate of 2nd point
    :return: a Jim object
    """
    return _pj.Jim(jim_object._jipjim.plotLine(x1, y1, x2, y2, val))


def polygonize(jim_object, output, **kwargs):
    """Polygonize Jim object based on GDALPolygonize.

    :param jim_object: a Jim object
    :param output: output filename of JimVect object that is returned.
        Use /vsimem for in memory vectors
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
    | mask             | mask with identical geometry as input raster object  |
    |                  | (zero is invalid, non-zero is valid)                 |
    +------------------+------------------------------------------------------+
    """
    kwargs.update({'output': output})

    if isinstance(jim_object, _pj.Jim):
        mask = kwargs.pop('mask', None)
        if mask is not None:
            if isinstance(mask, _pj.Jim):
                avect = jim_object._jipjim.polygonize(kwargs,
                                                      mask._jipjim)
            else:
                raise TypeError('Error: mask should be of Jim type')
        else:
            avect = jim_object._jipjim.polygonize(kwargs)
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only polygonize Jim object')


# def rasterize(jim_object, jim_vect, burnValue=1,eo=['ALL_TOUCHED'],ln=None):
#     """Rasterize Jim object based on GDALRasterizeLayersBuf

#     :param jim_object: a template Jim object
#     :param jim_vect: JimVect object that needs to be polygonized
#     :param burnValue: burn value
#     :param eo: option (default is ALL_TOUCHED)
#     :param ln: layer names (optional)
#     :return: rasterized Jim object

#     .. note::
#        Possible values for the key 'eo' are:

#        ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

#        For instance you can use 'eo':'ATTRIBUTE=fieldname'
#     """
#     if not isinstance(jim_object, _pj.Jim):
#         raise TypeError('Error: template must be a Jim object')
#     if not isinstance(jim_vect, _pj.JimVect):
#         raise TypeError('Error: can only rasterize a JimVect')

#     ajim=_pj.Jim(jim_object)
#     kwargs={}
#     kwargs.update({'burn':float(burnValue)})
#     kwargs.update({'eo':eo})
#     kwargs.update({'ln':ln})
#     ajim._jipjim.d_rasterizeBuf(item._jipjimvect,kwargs)
#     return ajim


def reducePlane(jim, rule, ref_band=None, nodata=None):
    """Reduce planes of Jim object.

    :param jim: jim object on which to reduce planes
    :param rule: rule to reduce (mean, median, min or max) or callback function
    :param ref_band: band on which to apply rule
        (default is to check all bands,
         not supported when rule is callback function)
    :param nodata: value to ignore when applying rule
        (not supported when rule is callback function)
    :return: reduced single plane jim object
    """
    if jim.properties.nrOfPlane() < 2:
        print("Warning: single plane, no reduction is performed")
        return None

    jimreduced = _pj.geometry.cropPlane(jim, 0)

    if isinstance(rule, str):
        if rule in ('mean', 'avg', 'median'):
            theType = jim.properties.getDataType()
            if nodata is not None and theType not in ('GDT_Float32',
                                                    'GDT_Float64'):
                jim = _pj.pixops.convert(jim, otype='GDT_Float32')
            if ref_band is not None:
                mask = _pj.geometry.cropBand(jim, band=ref_band)

                planes = [
                    mask.np()[i] for i in range(
                        jim.properties.nrOfPlane())]
                stacked_planes = numpy.vstack(planes)
                nr_of_row = jim.properties.nrOfRow()
                nr_of_col = jim.properties.nrOfCol()
                same_values = numpy.reshape(
                    numpy.diff(numpy.vstack(planes).reshape(len(planes), -1),
                            axis=0) == 0,
                    (nr_of_row, nr_of_col))

                nodata_mask = same_values & (mask.np()[0] == nodata)

            if rule == 'mean' or rule == 'avg':
                for iband in range(0, jim.properties.nrOfBand()):
                    if nodata is not None:
                        if ref_band is None:
                            planes = [
                                jim.np(iband)[i] for i in range(
                                    jim.properties.nrOfPlane())]
                            stacked_planes = numpy.vstack(planes)
                            nr_of_row = jim.properties.nrOfRow()
                            nr_of_col = jim.properties.nrOfCol()
                            same_values = numpy.reshape(
                                numpy.diff(
                                    numpy.vstack(planes).reshape(len(planes),
                                                                -1),
                                    axis=0) == 0,
                                (nr_of_row, nr_of_col))

                            nodata_mask = same_values & (
                                    jim.np(iband)[0] == nodata)

                            jim.np(iband)[jim.np(iband) == nodata] = numpy.nan
                            jimreduced.np(iband)[:] = numpy.nanmean(jim.np(iband),
                                                                    axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                        else:
                            jim.np(iband)[mask.np() == nodata] = numpy.nan
                            jimreduced.np(iband)[:] = numpy.nanmean(jim.np(iband),
                                                                    axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                    else:
                        jimreduced.np(iband)[:] = numpy.mean(jim.np(iband), axis=0)
            if rule == 'median':
                for iband in range(0, jim.properties.nrOfBand()):
                    if nodata is not None:
                        if ref_band is None:
                            planes = [
                                jim.np(iband)[i] for i in range(
                                    jim.properties.nrOfPlane())]
                            stacked_planes = numpy.vstack(planes)
                            nr_of_row = jim.properties.nrOfRow()
                            nr_of_col = jim.properties.nrOfCol()
                            same_values = numpy.reshape(
                                numpy.diff(
                                    numpy.vstack(planes).reshape(len(planes),
                                                                -1),
                                    axis=0) == 0,
                                (nr_of_row, nr_of_col))

                            nodata_mask = same_values & (
                                    jim.np(iband)[0] == nodata)

                            jim.np(iband)[jim.np(iband) == nodata] = numpy.nan
                            jimreduced.np(iband)[:] = numpy.nanmedian(
                                jim.np(iband), axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                        else:
                            jim.np(iband)[mask.np() == nodata] = numpy.nan
                            jimreduced.np(iband)[:] = numpy.nanmedian(
                                jim.np(iband), axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                    else:
                        jimreduced.np(iband)[:] = numpy.median(jim.np(iband),
                                                            axis=0)
            if nodata is not None:
                if theType not in ('GDT_Float32', 'GDT_Float64'):
                    jimreduced.pixops.convert(otype=theType)
                    jim.pixops.convert(otype=theType)
        else:
            if rule == 'max':
                def rule(reduced,plane):
                    return reduced<plane
            elif rule == 'min':
                def rule(reduced,plane):
                    return reduced>plane
            else:
                raise AttributeError(
                    'Error: rule not supported')

            if ref_band is not None:
                maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)
            for iplane in range(1, jim.properties.nrOfPlane()):
                jimplane = _pj.geometry.cropPlane(jim, iplane)

                if nodata is not None and ref_band is None:
                    raise AttributeError(
                        'Error: use ref_band option for nodata')

                if ref_band is not None:
                    maskplane = _pj.geometry.cropBand(jimplane, ref_band)
                    themask=rule(maskreduced,maskplane)
                    if nodata is not None:
                        themask |= maskreduced == nodata
                        themask &= maskplane != nodata
                    maskreduced[themask] = maskplane
                else:
                    themask=rule(jimreduced,jimplane)

                jimreduced[themask] = jimplane


                if nodata is not None:
                    nodata_mask = (maskreduced == nodata) & \
                                    (maskplane == nodata)
                    jimreduced[nodata_mask] = nodata
    else:
        if nodata is not None or ref_band is not None:
            raise AttributeError(
                'Error: nodata and ref_band are not supported for this rule')
        jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)
            jimreduced=rule(jimreduced,jimplane)
    return jimreduced


def stackBand(jim_object, jim_other=None, band=None):
    """Stack bands from raster datasets into new multiband Jim object.

    :param jim_object: a Jim or JimList object used for stacking the bands
    :param jim_other: a Jim object or jimlist from which to copy bands
        (optional)
    :param band: List of band indices to stack (index is 0 based).
        Default is to stack all bands.
    :return: Jim object with stacked bands

    Append all the bands of jim1 to jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        jim_stacked=pj.geometry.stackBand(jim0,jim1)

    Stack all bands of the JimList, returning a multi-band Jim object::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        jimlist=pj.JimList([jim0,jim1])
        jim_stacked=pj.geometry.stackBand(jimlist)

    Append the first three bands of raster dataset jim1 to the image jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        jim_stacked=pj.geometry.stackBand(jim0,jim1,band=[0,1,2])

    """
    if isinstance(jim_object, _pj.JimList):
        if band:
            retJim = _pj.Jim(jim_object._jipjimlist.stackBand({'band': band}))
        else:
            retJim = _pj.Jim(jim_object._jipjimlist.stackBand())

        if isinstance(jim_other, _pj.Jim):
            if band:
                return retJim.geometry.stackBand(jim_other, band=band)
            else:
                retJim.geometry.stackBand(jim_other)
        elif isinstance(jim_other, list):
            if band:
                jim_to_stack = _pj.Jim(
                    jim_other._jipjimlist.stackBand({'band': band}))
            else:
                jim_to_stack = _pj.Jim(jim_other._jipjimlist.stackBand())

            retJim = _pj.Jim(retJim._jipjim.stackBand(
                jim_to_stack._jipjim))

        return retJim
    elif isinstance(jim_object, _pj.Jim):
        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            if band:
                jim_object = _pj.Jim(jim_object._jipjim.stackBand(
                    jim._jipjim, {'band': band}))
            else:
                jim_object = _pj.Jim(jim_object._jipjim.stackBand(
                    jim._jipjim))

        return jim_object
    else:
        raise TypeError('Error: expected a Jim object')


def stackPlane(jim_object, jim_other=None, plane=None):
    """Stack planes from raster datasets into new multiplane Jim object.

    :param jim_object: a Jim or JimList object used for stacking the planes
    :param jim_other: a Jim object or jimlist from which to copy planes
        (optional)
    :return: Jim object with stacked planes

    Append all the planes of jim1 to jim0::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        jim_stacked=pj.geometry.stackPlane(jim0,jim1)

    Stack all planes of the JimList, returning a multi-plane Jim object::

        jim0=pj.Jim('/path/to/raster0.tif')
        jim1=pj.Jim('/path/to/raster1.tif')
        jimlist=pj.JimList([jim0,jim1])
        jim_stacked=pj.geometry.stackPlane(jimlist)
    """
    if isinstance(jim_object, _pj.JimList):
        retJim = _pj.Jim(jim_object._jipjimlist.stackPlane())
        if isinstance(jim_other, _pj.Jim):
            retJim.geometry.stackPlane(jim_other)
            return retJim
        else:
            return retJim
    elif isinstance(jim_object, _pj.Jim):
        retJim = _pj.Jim(jim_object)
        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            retJim.geometry.stackPlane(jim_other)

        return retJim
    else:
        raise TypeError('Error: expected a Jim object')


def warp(jim_object, t_srs, **kwargs):
    """Warp a raster dataset to a target spatial reference system.

    :param jim_object: a Jim object
    :param t_srs: Target spatial reference system
    :param kwargs: See table below
    :return: Cropped subimage as Jim instance

    +----------+---------------------------------------------------------------------------------------------+
    | key      | value                                                                                       |
    +==========+=============================================================================================+
    | s_srs    | Source spatial reference system (default is to read                                         |
    |          | from input)                                                                                 |
    +----------+---------------------------------------------------------------------------------------------+
    | resample | Resample algorithm used for reading pixel data in                                           |
    |          | case of interpolation                                                                       |
    |          | (default: GRIORA_NearestNeighbour). Check                                                   |
    |          | https://gdal.org/api/raster_c_api.html?highlight=griora_nearestn#_CPPv418GDALRIOResampleAlg |
    |          | for available options.                                                                      |
    +----------+---------------------------------------------------------------------------------------------+
    | nodata   | Nodata value to put in image if out of bounds                                               |
    +----------+---------------------------------------------------------------------------------------------+

    Example:

    Read a raster dataset from disk and warp to the target spatial reference
    system::

        jim = pj.Jim('/path/to/file.tif')
        jim_warped=pj.geometry.warp(jim, 'epsg:3035')

    Read a raster dataset from disk that is in lat lon (epsg:4326), select
    a bounding box in a different spatial reference system (epsg:3035).
    Notice the raster dataset read is still in the original
    projection (epsg:4326). Then warp the raster dataset to the target
    spatial reference system (epsg:3035)::

        jim = pj.Jim('/path/to/file.tif', t_srs='epsg:3035', ulx=1000000,
                     uly=4000000, lrx=1500000, lry=3500000)
        jim_warped=pj.geometry.warp(jim, 'epsg:3035', s_srs='epsg:4326')

    """
    kwargs.update({'t_srs': t_srs})

    jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.warp(kwargs))
    for iplane in range(1, jim_object.properties.nrOfPlane()):
        jimplane = _pj.Jim(_pj.geometry.cropPlane(jim_object,
                                                  iplane)._jipjim.warp(kwargs))
        jim.geometry.stackPlane(jimplane)
    return jim
    # return _pj.Jim(jim_object._jipjim.warp(kwargs))


class _Geometry():
    """Define all Geometry methods."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_object = caller

    def aggregate_vector(self, jvec, rule, output, **kwargs):
        """Extract pixel values from raster image based on a vector dataset.

        :param jvec: reference JimVect instance
        :param rule: Rule how to calculate zonal statistics per feature
            (see list of :ref:`supported rules <extract_rules>`)
        :param output: Name of the output vector dataset in which the zonal
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
        | planename        | List of plane names corresponding to list of     |
        |                  | planes to extract                                |
        +------------------+--------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)           |
        +------------------+--------------------------------------------------+
        | endband          | End band sequence number (0 indexed)             |
        +------------------+--------------------------------------------------+
        | plane            | List of planes to extract (0 indexed). Default is|
        |                  | to use extract all planes                        |
        +------------------+--------------------------------------------------+
        | planename        | List of plane name corresponding to list of      |
        |                  | planes to extract                                |
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

        :Supported key values to mask pixels that must be ignored in \
        the extraction process:

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | srcnodata        | List of nodata values not to extract             |
        +------------------+--------------------------------------------------+
        | buffer           | Buffer (in geometric units of raster dataset).   |
        |                  | Use neg. value to exclude pixels within buffer   |
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

        Extract the mean value of the pixels within the polygon of the provided
        reference vector. Exclude the pixels within a buffer of 10m of
        the polygon boundary. Use a temporary vector in memory for
        the calculation. Then write the result to the final destination
        on disk::

            reference = pj.JimVect('/path/to/reference.sqlite')
            jim0 = pj.Jim('/path/to/raster.tif')
            v = jim0.geometry.aggregateVector(
                reference, buffer=-10, rule=['mean'],
                output='/vsimem/temp.sqlite', oformat='SQLite')
            v.write('/path/to/output.sqlite)

        """
        # make list of rules
        rules = []
        if rule:
            for irule in rule:
                rules.append(irule)

        kwargs.update({'output': output})
        kwargs.update({'rule': rules})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']

        bandnames = kwargs.pop('bandname', None)
        if bandnames is None:
            bandnames = ['b'+str(iband) for iband in range(
                0, self._jim_object.properties.nrOfBand())]

        planenames = kwargs.pop('planename', None)
        if planenames is None:
            planenames = ['t'+str(iplane) for iplane in range(
                0, self._jim_object.properties.nrOfPlane())]

        if jvec is None:
            raise AttributeError('Error: missing jvec option')

        if self._jim_object.properties.nrOfPlane() > 1:
            plane = kwargs.pop('plane', None)
            if plane is None:
                plane = range(0, self._jim_object.properties.nrOfPlane())
            else:
                for iplane in plane:
                    if iplane >= self._jim_object.properties.nrOfPlane():
                        raise AttributeError(
                            'Error: illegal plane {}'.format(iplane))
            if len(planenames) != self._jim_object.properties.nrOfPlane():
                raise AttributeError(
                    'Error: number of planes does not correspond to planename')

            if jvec.properties.getLayerCount() > 1:
                raise AttributeError(
                    'Error: multiple layers not supported when aggregating '
                    'vectors over multi-plane raster datasets, please use '
                    'single layer vector object')
        else:
            plane = [0]
            planenames = None

        if len(plane) == 1:
            if '/' not in kwargs['output']:
                kwargs['output'] = os.path.join('/vsimem/', kwargs['output'])
            if self._jim_object.properties.nrOfPlane() == 1:
                avect = self._jim_object.geometry.extractOgr(jvec, **kwargs)
            else:
                avect = _pj.geometry.cropPlane(
                    self, plane[0])._jim_object.geometry.extractOgr(jvec,
                                                                    **kwargs)
            avect.io.write()
            avect.io.close()

            pjvect = _pj.JimVect(kwargs['output'])
            pjvect._set(avect)
            return pjvect

        intersectfn = '/vsimem/sampleintersect.sqlite'
        sampleintersect = _pj.geometry.intersect(
            jvec, self._jim_object, output=intersectfn, oformat='SQLite',
            co=['OVERWRITE=YES'])
        if sampleintersect.properties.isEmpty():
            sampleintersect.io.close()
            raise ValueError('intersect is empty')
        else:
            sampleintersect.io.write()

        firstExtract = True
        joinfn = None

        vsioutput = os.path.join('/vsimem', 'aggregate_polygon.sqlite')
        for iplane in range(0, len(plane)):
            ibandnames = []
            try:
                for band in bandnames:
                    if len(rules) > 1:
                        # rules are automatically pre-pended in extractogr
                        ibandnames.append(
                            '_' + planenames[iplane] + '_' + band)
                    else:
                        ibandnames.append(
                            rules[0] + '_' + planenames[iplane] + '_' + band)

                fieldnames = sampleintersect.properties.getFieldNames()

                try:
                    if 'buffer' in kwargs:
                        jim = _pj.geometry.cropPlane(self._jim_object,
                                                     plane[iplane])
                        v = jim.geometry.extractOgr(
                            sampleintersect, rule=rules, output=vsioutput,
                            oformat='SQLite', co=['OVERWRITE=YES'],
                            bandname=ibandnames, copy=fieldnames, fid='fid',
                            buffer=kwargs['buffer'])
                    else:
                        jim = _pj.geometry.cropPlane(self._jim_object,
                                                     plane[iplane])
                        v = jim.geometry.extractOgr(
                            sampleintersect, rule=rules, output=vsioutput,
                            oformat='SQLite', co=['OVERWRITE=YES'],
                            bandname=ibandnames, copy=fieldnames, fid='fid')
                    if v.properties.isEmpty():
                        v.io.close()
                except:
                    print("no coverage for plane {}, continue with next "
                          "product".format(planenames[iplane]))
                    if v:
                        v.io.close()
                    continue
                if not v.properties.isEmpty():
                    v.io.write()
                    # join vectors
                    if '/' not in output:
                        joinfn = os.path.join('/vsimem/', output)
                    else:
                        joinfn = output
                    # joinfn='/vsimem/vjoin.sqlite'

                    if firstExtract:
                        vjoin = _pj.JimVect(v, output=joinfn,
                                            co='OVERWRITE=YES')
                        vjoin.io.write()
                        vjoin.io.close()
                        firstExtract = False
                    else:
                        vprev = _pj.JimVect(joinfn)
                        vjoin = _pj.geometry.join(
                            vprev, v, joinfn, oformat='SQLite',
                            co=['OVERWRITE=YES'], key=['fid'],
                            method='OUTER_FULL')
                        vjoin.io.write()
                        vjoin.io.close()
                    v.io.close()
                else:
                    v.io.close()
            except:
                print("raised exception dataset in for plane {}".format(
                    planenames[iplane]))
                continue

        sampleintersect.io.close()
        if joinfn:
            v = _pj.JimVect(joinfn)
            return v
        else:
            print("Error: joinfn is None, no valid features found")
            raise ValueError('Error: joinfn is None, no valid features found')

    def band2plane(self):
        """Convert 2-dimensional multi-band object to a 3-dimensional.

        The result will be a single band multi-plane object

        Example: convert a multi-band object with 12 bands to a 3-dimensional
        single band object with 12 planes::

           jim=pj.Jim('/path/to/multi/band/image.tif')
           jim.properties.nrOfBand()
           12
           jim.properties.nrOfPlane()
           1
           jim.geometry.band2plane()
           jim.properties.nrOfPlane()
           12
           jim.properties.nrOfBand()
           1

        Notice that a multi-band image can also be read directly as
        a multi-plane object::

           jim=pj.Jim('/path/to/multi/band/image.tif',band2plane=True)
           jim.properties.nrOfBand()
           1
           jim.properties.nrOfPlane()
           12

        """
        self._jim_object._jipjim.d_band2plane()

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
        :param dx: spatial resolution in x to crop (stride if nogeo is True)
        :param dy: spatial resolution in y to crop (stride if nogeo is True)
        :param nogeo: use image coordinates if True, default is spatial
            reference system coordinates
        :param nodata: set no data in case the specified bounding box is not
            within the object boundaries (default is 0)

        Example:

        Crop bounding box in georeferenced coordinates::

            jim=pj.Jim('/path/to/raster.tif')
            jim.crop(ulx=1000000,uly=5000000,lrx=2000000,lry=4000000,
                     dx=1000,dy=1000)

        Crop bounding box in image coordinates (starting from upper left pixel
        coordinate 0, 0). For instance, get first 10 columns in first 10 rows::

            jim=pj.Jim('/path/to/raster.tif')
            jim.crop(ulx=0,uly=0,lrx=10,lry=10, nogeo=True)

        Notice that for this case, a more pythonic way to is available
        via :ref:`indexing`::

            jim0[0:10,0:10]

        However, crop can also be used to enlarge a Jim object. For instance,
        to add a border of one pixel use::

            jim=pj.Jim('/path/to/raster.tif')
            jim.geometry.crop(ulx=-1,uly=-1,lrx=jim.properties.nrOfCol()+1,
                              lry=jim.properties.nrOfRow()+1,nogeo=True)

        """
        if ulz is not None or lrz is not None:
            assert len(kwargs) == 0, \
                'It is not supported to use both z coords and special ' \
                'cropping parameters'
            gt = self._jim_object.properties.getGeoTransform()
            nr_of_cols = self._jim_object.properties.nrOfCol()
            nr_of_rows = self._jim_object.properties.nrOfRow()
            if nogeo:
                if ulx is None:
                    ulx = 0
                if uly is None:
                    uly = 0
                if lrx is None:
                    lrx = self._jim_object.properties.nrOfCol()
                if lry is None:
                    lry = self._jim_object.properties.nrOfRow()
                if dx is None:
                    dx = 1
                if dy is None:
                    dy = 1
                uli = ulx
                ulj = uly
                lri = lrx
                lrj = lry
                upperLeft = self._jim_object.geometry.image2geo(float(ulx),
                                                                float(uly))
                lowerRight = self._jim_object.geometry.image2geo(float(lrx),
                                                                 float(lry))
                ulx = upperLeft[0]
                uly = upperLeft[1]
            else:
                upperLeftImage = self._jim_object.geometry.geo2image(ulx, uly)
                uli = upperLeftImage[0]
                ulj = upperLeftImage[1]
                lowerRightImage = self._jim_object.geometry.geo2image(lrx, lry)
                lri = lowerRightImage[0]
                lrj = lowerRightImage[1]
            for iband in range(0, self._jim_object.properties.nrOfBand()):
                self._jim_object._jipjim.d_imageFrameSubtract([
                    uli, nr_of_cols - lri,
                    ulj, nr_of_rows - lrj,
                    ulz, self._jim_object.properties.nrOfPlane() - lrz])
            gt[0] = ulx
            gt[3] = uly
            self._jim_object.properties.setGeoTransform(gt)
        else:
            if nogeo:
                if ulx is None:
                    ulx = 0
                if uly is None:
                    uly = 0
                if lrx is None:
                    lrx = self._jim_object.properties.nrOfCol()
                if lry is None:
                    lry = self._jim_object.properties.nrOfRow()
                if dx is None:
                    dx = 1
                if dy is None:
                    dy = 1
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

            jim = _pj.Jim(_pj.geometry.cropPlane(self._jim_object,
                                                 0)._jipjim.crop(kwargs))
            for iplane in range(1, self._jim_object.properties.nrOfPlane()):
                jimplane = _pj.Jim(_pj.geometry.cropPlane(
                    self._jim_object, iplane)._jipjim.crop(kwargs))
                jim.geometry.stackPlane(jimplane)
            self._jim_object._set(jim._jipjim)
            # self._jim_object._set(self._jim_object._jipjim.crop(kwargs))
        # else:
        #     if nogeo:
        #         uli = ulx
        #         ulj = uly
        #         lri = lrx
        #         lrj = lry
        #         upperLeft = self._jim_object.geometry.image2geo(ulx, uly)
        #         lowerRight = self._jim_object.geometry.image2geo(lrx, lry)
        #         ulx = upperLeft[0]
        #         uly = upperLeft[1]
        #         lrx = lowerRight[0]+self._jim_object.properties.getDeltaX()/2.0
        #         lry = lowerRight[1]-self._jim_object.properties.getDeltaY()/2.0
        #         if dx is None:
        #             dx = 1
        #         if dy is None:
        #             dy = 1
        #     else:
        #         if ulx is None:
        #             ulx = self._jim_object.properties.getUlx()
        #         if uly is None:
        #             uly = self._jim_object.properties.getUly()
        #         if lrx is None:
        #             lrx = self._jim_object.properties.getLrx()
        #         if lry is None:
        #             lry = self._jim_object.properties.getLry()
        #         if dx is None:
        #             dx = self._jim_object.properties.getDeltaX()
        #         if dy is None:
        #             dy = self._jim_object.properties.getDeltaY()
        #     kwargs.update({'ulx': ulx})
        #     kwargs.update({'uly': uly})
        #     kwargs.update({'lrx': lrx})
        #     kwargs.update({'lry': lry})
        #     kwargs.update({'dx': dx})
        #     kwargs.update({'dy': dy})
        #     kwargs.update({'nogeo': nogeo})
        #     self._jim_object._set(self._jim_object._jipjim.crop(kwargs))
        #     # return _pj.Jim(self._jim_object.crop(kwargs))

    def cropBand(self, band):
        """Subset raster dataset.

        Subset raster dataset in spectral domain.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim0.cropBand(band=[0,1,2])

        """
        self._jim_object._jipjim.d_cropBand({'band': band})

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
        jim = _pj.Jim(_pj.geometry.cropPlane(
            self._jim_object, 0)._jipjim.cropOgr(extent._jipjimvect, kwargs))
        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.Jim(_pj.geometry.cropPlane(
                self._jim_object, iplane)._jipjim.cropOgr(extent._jipjimvect,
                                                          kwargs))
            jim.geometry.stackPlane(jimplane)
        self._jim_object._set(jim._jipjim)
        # self._jim_object._set(
        #     self._jim_object._jipjim.cropOgr(extent._jipjimvect, kwargs))

    def cropPlane(self, plane):
        """Subset raster dataset.

        Subset raster dataset in temporal domain.

        Modifies the instance on which the method was called.

        :param plane: List of plane indices to crop (index is 0 based)

        Example:

        Crop the first three planes from raster dataset jim0::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim0.cropPlane(plane=[0,1,2])

        """
        self._jim_object._jipjim.d_cropPlane({'plane': plane})

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

        :Supported key values to mask pixels that must be ignored in \
        the extraction process:

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

            jvec=pj.Jim('/path/to/landcovermap.tif')

            classes=[2,12,25,41,50]
            thresholds=['20%','25%','25%','10%','5%']

            jvec=pj.Jim('/path/to/multiband.tif',
                        dx=jim.getDeltaX(),dy=jim.getDeltaY(),
                        ulx=jim.getUlx(),uly=jim.getUly(),
                        lrx=jim.getLrx(),lry=jim.getLry())

            outputfn='/path/to/output.sqlite'
            sample=jim.extractImg(jvec,srcnodata=[0],output=outputfn,
                                  class=classes,threshold=thresholds,
                                  bandname=['B02','B03','B04','B08'],
                                  band=[0,1,2,3])
        """
        kwargs.update({'output': output})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']
        return self._jim_object._jipjim.extractImg(reference._jipjim, kwargs)

    # deprecated: use aggregate_vector instead
    def extractOgr(self, jvec, rule, output, **kwargs):
        """Extract pixel values from raster image based on a vector dataset.

        :param jvec: reference JimVect instance
        :param rule: Rule how to calculate zonal statistics per feature
            (see list of :ref:`supported rules <extract_rules>`)
        :param output: Name of the output vector dataset in which the zonal
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
        | bandname         | List of band names corresponding to list of      |
        |                  | bands to extract                                 |
        +------------------+--------------------------------------------------+
        | planename        | List of plane names corresponding to list of     |
        |                  | planes to extract                                |
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

        :Supported key values to mask pixels that must be ignored in \
        the extraction process:

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | srcnodata        | List of nodata values not to extract             |
        +------------------+--------------------------------------------------+
        | buffer           | Buffer (in geometric units of raster dataset).   |
        |                  | Use neg. value to exclude pixels within buffer   |
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

        Extract the mean value of the pixels within the polygon of the provided
        reference vector. Exclude the pixels within a buffer of 10m of
        the polygon boundary. Use a temporary vector in memory for
        the calculation. Then write the result to the final destination
        on disk::

            reference = pj.JimVect('/path/to/reference.sqlite')
            jim0 = pj.Jim('/path/to/raster.tif')
            v = jim0.geometry.extractOgr(
                reference, buffer=-10, rule=['mean'],
                output='/vsimem/temp.sqlite', oformat='SQLite')
            v.io.write('/path/to/output.sqlite)

        """
        kwargs.update({'output': output})
        kwargs.update({'rule': rule})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']

        avect = self._jim_object._jipjim.extractOgr(jvec._jipjimvect, kwargs)
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect

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
        |                  | features (in geometric units of raster dataset)  |
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
        on a 3x3 window (buffer value of 100 m) in a vector
        dataset in memory::

            v01 = jim0.extractSample(random=100, buffer=100, rule=['mean'],
                                     output='mem01', oformat='Memory')

        Extract a sample of 100 points using a regular grid sampling scheme.
        For each grid point, calculate the median value based on a 3x3 window
        (buffer value of 100 m neighborhood). Write the result in a SQLite
        vector dataset on disk::

            outputfn = '/path/to/output.sqlite'
            npoint = 100
            gridsize = int(jim.nrOfCol()*jim.getDeltaX()/math.sqrt(npoint))
            v = jim.extractSample(grid=gridsize, buffer=100, rule=['median'],
                                  output=outputfn, oformat='SQLite')
            v.write()

        """
        kwargs.update({'output': output})
        if 'threshold' in kwargs:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']
        return self._jim_object._jipjim.extractSample(kwargs)

    def geo2image(self, x, y):
        """Convert georeferenced coordinates to image coordinates (col and row).

        :param x: georeferenced coordinate in x according to the object spatial
            reference system
        :param y: georeferenced coordinate in y according to the object spatial
            reference system
        :return: image coordinates (row and column, starting from 0)

        Get column and row index (0 based) of some georeferenced coordinates
        x and y (in this case first pixel: 0, 0)::

          jim=pj.Jim('/path/to/raster.tif')
          x=jim.properties.getUlx()
          y=jim.properties.getUly()
          jim.geometry.geo2image(x,y)
        """
        coord = self._jim_object._jipjim.geo2image(x, y)
        return [int(coord[0]), int(coord[1])]

    def image2geo(self, i, j):
        """Convert image coordinates (col and row) to georeferenced coordinates.

        :param i: image column number (starting from 0)
        :param j: image row number (starting from 0)
        :return: georeferenced coordinates according to the object spatial
            reference system

        Get upper left corner in georeferenced coordinates
        (in SRS of the Jim object)::

          jim=pj.Jim('/path/to/raster.tif')
          jim.geometry.image2geo(0,0)

        """
        return self._jim_object._jipjim.image2geo(i, j)

    def imageFrameAdd(self, l=0, r=0, t=0, b=0, u=0, d=0, val=0):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        """
        if self._jim_object.properties.nrOfBand() > 1:
            returnJim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if returnJim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
                    returnJim.geometry.stackBand(jimband)
                else:
                    returnJim = _pj.geometry.cropBand(self._jim_object,
                                                      band=band)
                    returnJim._jipjim.d_imageFrameAdd(
                        [l, r, t, b, u, d], val)
            self._jim_object._set(returnJim._jipjim)
        else:
            self._jim_object._jipjim.d_imageFrameAdd(
                [l, r, t, b, u, d], val)

    def imageFrameSet(self, l=0, r=0, t=0, b=0, u=0, d=0, val=0, band=None):
        """Set the values of the image frame to value val.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        :param val: value of frame
        :param band: value of band
        """
        if band is None:
            self._jim_object._jipjim.d_imageFrameSet(
                [l, r, t, b, u, d], val, -1)
        else:
            self._jim_object._jipjim.d_imageFrameSet(
                [l, r, t, b, u, d], val, band)

    def imageFrameSubtract(self, l=0, r=0, t=0, b=0, u=0, d=0):
        """Subtract an image frame.

        Modifies the instance on which the method was called.

        :param l: width of left frame
        :param r: width of right frame
        :param t: width of top frame
        :param b: width of bottom frame
        :param u: width of upper frame
        :param d: width of lower frame
        """
        if self._jim_object.properties.nrOfBand() > 1:
            returnJim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if returnJim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
                    returnJim.geometry.stackBand(jimband)
                else:
                    returnJim = _pj.geometry.cropBand(self._jim_object,
                                                      band=band)
                    returnJim._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
            self._jim_object._set(returnJim._jipjim)
        else:
            self._jim_object._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])

    def imageInsert(self, sec_jim_object, x, y, z, band=None):
        """Merge Jim instance with values of sec_jim_object in given coords.

        Modifies the instance on which the method was called.

        :param sec_jim_object: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param band: List of band indices to insert (index is 0 based)
        """
        bands = []
        if band is None:
            bands = range(0, self._jim_object.properties.nrOfBand())
        else:
            bands.extend(band)

        for band in bands:
            self._jim_object._jipjim.d_imageInsert(sec_jim_object._jipjim,
                                                   x, y, z, band)

    def imageInsertCompose(self, imlbl, im2, x, y, z, val, band=None):
        """Merge Jim instance with values of im2 if val of imlbl == val.

        Modifies the instance on which the method was called.

        :param imRaster_imlbl: a Jim object
        :param imRaster_im2: a Jim object
        :param x: x coordinate of 1st pixel
        :param y: y coordinate of 1st pixel
        :param z: z coordinate of 1st pixel
        :param val: integer for label value
        :param band: List of band indices to compose (index is 0 based)
        """
        bands = []
        if band is None:
            bands = range(0, self._jim_object.properties.nrOfBand())
        else:
            try:
                bands.extend(band)
            except TypeError:
                bands.append(band)
        for band in bands:
            self._jim_object._jipjim.d_imageInsertCompose(imlbl._jipjim,
                                                          im2._jipjim,
                                                          x, y, z, val, band)

    def magnify(self, n):
        """Magnify the image.

        Modifies the instance on which the method was called.

        :param n: a positive integer for magnifying size by pixel replication
        """
        if self._jim_object.properties.nrOfBand() > 1:
            returnJim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if returnJim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband = _pj.Jim(jimband._jipjim.imageMagnify(n))

                    returnJim.geometry.stackBand(jimband)
                else:
                    returnJim = _pj.Jim(_pj.geometry.cropBand(
                        self._jim_object, band=band)._jipjim.imageMagnify(n))
            self._jim_object._set(returnJim._jipjim)
        else:
            self._jim_object._set(self._jim_object._jipjim.imageMagnify(n))

    def plane2band(self):
        """Convert 3-dimensional single-band object to a 2-dimensional
        multi-band object.

        The result will be a multi-band single plane object

        Example: convert a single band object with 12 planes to a 2-dimensional
        multi-band object with 1 plane::

           jim=pj.Jim('/path/to/multi/band/image.tif',band2plane=True)
           jim.properties.nrOfBand()
           1
           jim.properties.nrOfPlane()
           12
           jim.geometry.plane2band()
           jim.properties.nrOfPlane()
           1
           jim.properties.nrOfBand()
           12

        """
        result = None
        for iplane in range(0, self._jim_object.properties.nrOfPlane()):
            jim_plane = _pj.geometry.cropPlane(self._jim_object, iplane)
            if result is None:
                result = jim_plane
            else:
                result.geometry.stackBand(jim_plane)
        self._jim_object._set(result._jipjim)

    def plotLine(self, x1, y1, x2, y2, val):
        """Draw a line from [x1, y1] to [x2, y2] by setting pixels to a value.

        Modifies the instance on which the method was called.

        :param x1: an integer for x-coordinate of 1st point
        :param y1: an integer for y-coordinate of 1st point
        :param x2: an integer for x-coordinate of 2nd point
        :param y2: an integer for y-coordinate of 2nd point
        """
        self._jim_object._jipjim.d_plotLine(x1, y1, x2, y2, val)

    def polygonize(self, output, **kwargs):
        """Polygonize Jim object based on GDALPolygonize.

        Returns a new JimVect object.

        :param output: output filename of JimVect object that is returned.
            Use /vsimem for in memory vectors
        :param kwargs: See table below
        :return: JimVect object with polygons

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | ln               | Output layer name                                |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+
        | name             | Field name of the output layer (default is DN)   |
        +------------------+--------------------------------------------------+
        | nodata           | Discard this nodata value when creating polygons |
        +------------------+--------------------------------------------------+
        | mask             | mask with identical geometry as input raster     |
        |                  | object (zero is invalid, non-zero is valid)      |
        +------------------+--------------------------------------------------+

        Example: create a polygon vector file from a Sentinel-2 classification
        raster dataset, where clouds are represented by the pixel value 9::

          sclfn = '/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L2A/2018/07/01/065/S2A_MSIL2A_20180701T102021_N0208_R065_T33UUT_20180701T141038.SAFE/GRANULE/L2A_T33UUT_A015792_20180701T102404/IMG_DATA/R20m/T33UUT_20180701T102021_SCL_20m.jp2'
          sclJim = pj.Jim(sclfn)
          sclJim[sclJim != 9] = 0
          sclJim[sclJim == 9] = 1
          vcloud = sclJim.geometry.polygonize('/vsimem/cloud.sqlite',
                                              name='cloud', nodata=0)
          vcloud.io.write('/path/to/cloud.sqlite')
        """
        kwargs.update({'output': output})
        mask = kwargs.pop('mask', None)
        if mask is not None:
            if isinstance(mask, _pj.Jim):
                avect = self._jim_object._jipjim.polygonize(kwargs,
                                                            mask._jipjim)
            else:
                raise TypeError('Error: mask should be of Jim type')
        else:
            avect = self._jim_object._jipjim.polygonize(kwargs)
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect

    #todo: to be tested (we can also use jim[jimVect] instead...)
    # def rasterize(self, jim_vect, burnValue=1,eo=['ALL_TOUCHED'],ln=None):
    #     """Rasterize Jim object based on GDALRasterizeLayersBuf

    # CPLErr Jim::rasterizeBuf(VectorOgr& ogrReader, double burnValue,
    #                          const std::vector<std::string>& eoption,
    #                          const std::vector<std::string>& layernames ){

    #     :param jim_vect: JimVect object that needs to be polygonized
    #     :param burnValue: burn value
    #     :param eo: option (default is ALL_TOUCHED)
    #     :param ln: layer names (optional)

    #     .. note::
    #       Possible values for the key 'eo' are:

    #       ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

    #       For instance you can use 'eo':'ATTRIBUTE=fieldname'
    #     """

    #     print(jim_vect)
    #     print("type of jim_vect: {}".format(type(jim_vect)))
    #     if not isinstance(jim_vect, _pj.JimVect):
    #         raise TypeError('Error: can only rasterize JimVect')

    #     kwargs={}
    #     kwargs.update({'burn':float(burnValue)})
    #     kwargs.update({'eo':eo})
    #     kwargs.update({'ln':ln})
    #     self._jim_object._jipjim.d_rasterizeBuf(jim_vect._jipjimvect,kwargs)

    def reducePlane(self, rule, ref_band=None, nodata=None):
        """Reduce planes of Jim object

        :param rule: rule to reduce (mean, median, min or max) or callback function
        :param ref_band: band on which to apply rule
            (default is to check all bands,
             not supported when rule is callback function)
        :param nodata: value to ignore when applying rule
            (not supported when rule is callback function)

        Stack planes of two single plane jim objects, then reduce by taking
        the means::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim_stacked = pj.geometry.stackPlane(jim0, jim1)
            jim_stacked.geometry.reducePlane('mean')

        Stack planes of two single plane jim objects, then reduce by taking
        the maximum using callback function::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim_stacked = pj.geometry.stackPlane(jim0, jim1)
            def getMax(reduced, plane):
                return pj.pixops.supremum(reduced, plane)
            jim_stacked.geometry._reducePlaneSimple(getMax)

        """
        if self._jim_object.properties.nrOfPlane() < 2:
            print("Warning: single plane, no reduction is performed")
            return None

        jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
        if isinstance(rule, str):
            if rule in ('mean', 'avg', 'median'):
                theType = self._jim_object.properties.getDataType()
                if nodata is not None and theType not in ('GDT_Float32',
                                                        'GDT_Float64'):
                    self._jim_object.pixops.convert(otype='GDT_Float32')
                if ref_band is not None:
                    mask = _pj.geometry.cropBand(self._jim_object, band=ref_band)

                    planes = [
                        mask.np()[i] for i in range(
                            self._jim_object.properties.nrOfPlane())]
                    stacked_planes = numpy.vstack(planes)
                    nr_of_row = self._jim_object.properties.nrOfRow()
                    nr_of_col = self._jim_object.properties.nrOfCol()
                    same_values = numpy.reshape(
                        numpy.diff(stacked_planes.reshape(len(planes), -1),
                                axis=0) == 0,
                        (nr_of_row, nr_of_col))

                    nodata_mask = same_values & (mask.np()[0] == nodata)


                if rule == 'mean' or rule == 'avg':
                    for iband in range(0, self._jim_object.properties.nrOfBand()):
                        if nodata is not None:
                            if ref_band is None:
                                planes = [
                                    self._jim_object.np(iband)[i] for i in range(
                                        self._jim_object.properties.nrOfPlane())]
                                stacked_planes = numpy.vstack(planes)
                                nr_of_row = self._jim_object.properties.nrOfRow()
                                nr_of_col = self._jim_object.properties.nrOfCol()
                                same_values = numpy.reshape(
                                    numpy.diff(
                                        numpy.vstack(planes).reshape(len(planes),
                                                                    -1),
                                        axis=0) == 0,
                                    (nr_of_row, nr_of_col))

                                nodata_mask = same_values & (
                                        self._jim_object.np(iband)[0] == nodata)

                                self._jim_object.np(iband)[self._jim_object.np(
                                    iband) == nodata] = numpy.nan
                                jimreduced.np(iband)[:] = numpy.nanmean(
                                    self._jim_object.np(iband), axis=0)
                                jimreduced.np(iband)[nodata_mask] = nodata
                            else:
                                self._jim_object.np(iband)[mask.np() == nodata] = \
                                    numpy.nan
                                jimreduced.np(iband)[:] = numpy.nanmean(
                                    self._jim_object.np(iband), axis=0)
                                jimreduced.np(iband)[nodata_mask] = nodata
                        else:
                            jimreduced.np(iband)[:] = numpy.mean(
                                self._jim_object.np(iband), axis=0)
                elif rule == 'median':
                    for iband in range(0, self._jim_object.properties.nrOfBand()):
                        if nodata is not None:
                            if ref_band is None:
                                planes = [
                                    self._jim_object.np(iband)[i] for i in range(
                                        self._jim_object.properties.nrOfPlane())]
                                stacked_planes = numpy.vstack(planes)
                                nr_of_row = self._jim_object.properties.nrOfRow()
                                nr_of_col = self._jim_object.properties.nrOfCol()
                                same_values = numpy.reshape(
                                    numpy.diff(
                                        numpy.vstack(planes).reshape(len(planes),
                                                                    -1),
                                        axis=0) == 0,
                                    (nr_of_row, nr_of_col))

                                nodata_mask = same_values & (
                                        self._jim_object.np(iband)[0] == nodata)

                                self._jim_object.np(iband)[self._jim_object.np(
                                        iband) == nodata] = numpy.nan
                                jimreduced.np(iband)[:] = numpy.nanmedian(
                                    self._jim_object.np(iband), axis=0)
                                jimreduced.np(iband)[nodata_mask] = nodata
                            else:
                                self._jim_object.np(iband)[mask.np() == nodata] = \
                                    numpy.nan
                                jimreduced.np(iband)[:] = numpy.nanmedian(
                                    self._jim_object.np(iband), axis=0)
                                jimreduced.np(iband)[nodata_mask] = nodata
                        else:
                            jimreduced.np(iband)[:] = numpy.median(
                                self._jim_object.np(iband), axis=0)
                if nodata is not None:
                    if theType not in ('GDT_Float32', 'GDT_Float64'):
                        jimreduced.pixops.convert(otype=theType)
                        self._jim_object.pixops.convert(otype=theType)
                self._jim_object._set(jimreduced._jipjim)
            else:
                if rule == 'max':
                    def rule(reduced,plane):
                        return reduced<plane
                elif rule == 'min':
                    def rule(reduced,plane):
                        return reduced>plane
                else:
                    raise AttributeError('Error: rule not supported')
                maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)
                for iplane in range(1, self._jim_object.properties.nrOfPlane()):
                    jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)

                    if nodata is not None and ref_band is None:
                        raise AttributeError(
                            'Error: use ref_band option for nodata')

                    if ref_band is not None:
                        maskplane = _pj.geometry.cropBand(jimplane, ref_band)
                        themask=rule(maskreduced,maskplane)
                        if nodata is not None:
                            themask |= maskreduced == nodata
                            themask &= maskplane != nodata
                        maskreduced[themask] = maskplane
                    else:
                        themask=rule(jimreduced,jimplane)

                    jimreduced[themask] = jimplane


                    if nodata is not None:
                        nodata_mask = (maskreduced == nodata) & \
                                    (maskplane == nodata)
                        jimreduced[nodata_mask] = nodata

                self._jim_object._set(jimreduced._jipjim)
        else:
            if nodata is not None or ref_band is not None:
                raise AttributeError(
                    'Error: nodata and ref_band are not supported for this rule')
            jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
            for iplane in range(1, self._jim_object.properties.nrOfPlane()):
                jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)
                jimreduced=rule(jimreduced,jimplane)
            self._jim_object._set(jimreduced._jipjim)


    def _reducePlaneSimple(self, rule):
        """Reduce planes of Jim object using callback function without nodata (for performance reasons).

        :param rule: rule to reduce (mean, median) or callback function
            (e.g., for max composite, use callBack(jimreduced, jimplane): return jimreduced<jimplane)

        Stack planes of two single plane jim objects, then reduce by taking
        the maximum::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim_stacked = pj.geometry.stackPlane(jim0, jim1)
            def getMax(reduced, plane):
                return pj.pixops.supremum(reduced, plane)
            jim_stacked.geometry._reducePlaneSimple(getMax)


        """
        if self._jim_object.properties.nrOfPlane() < 2:
            print("Warning: single plane, no reduction is performed")
            return None

        if rule == 'max':
            def rule(reduced,plane):
                return _pj.pixops.supremum(reduced,plane)
        elif rule == 'min':
            def rule(reduced,plane):
                return _pj.pixops.infimum(reduced,plane)
        jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)
            jimreduced=rule(jimreduced,jimplane)
        self._jim_object._set(jimreduced._jipjim)

    def stackBand(self, jim_other, band=None):
        """Stack the bands of another Jim object to the current Jim object.

        Modifies the instance on which the method was called.

        :param jim_other: a Jim object or jimlist from which to copy bands
        :param band: List of band indices to stack (index is 0 based)

        Example:

        Append all the bands of a multiband Jim object jim1 to the current
        single band Jim object jim0::

            jim0 = pj.Jim('/path/to/singleband.tif')
            jim1 = pj.Jim('/path/to/multiband.tif')
            jim0.geometry.stackBand(jim1)

        Append the first three bands of raster dataset jim1 to the current
        Jim object jim0::

            jim0 = pj.Jim('/path/to/singleband.tif')
            jim1 = pj.Jim('/path/to/multiband.tif')
            jim0.geometry.stackBand(jim1, band=[0, 1, 2])
        """
        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            if band:
                self._jim_object._jipjim.d_stackBand(jim._jipjim,
                                                     {'band': band})
            else:
                self._jim_object._jipjim.d_stackBand(jim._jipjim)

    def stackPlane(self, jim_other):
        """Stack the planes of another Jim object to the current Jim object.

        Modifies the instance on which the method was called.

        :param jim_other: a Jim object or jimlist from which to copy planes

        Example:

        Append all the planes of a multiplane Jim object jim1 to the current
        single plane Jim object jim0::

            jim0 = pj.Jim('/path/to/singleplane.tif')
            jim1 = pj.Jim('/path/to/multiplane.tif')
            jim0.geometry.stackPlane(jim1)
        """
        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            self._jim_object._jipjim.d_stackPlane(jim._jipjim)

    def warp(self, t_srs, **kwargs):
        """Warp a raster dataset to a target spatial reference system.

        :param t_srs: Target spatial reference system
        :param kwargs: See table below

        +----------+-------------------------------------------------------------------+
        | key      | value                                                             |
        +==========+===================================================================+
        | s_srs    | Source spatial reference system                                   |
        |          | (default is to read from input)                                   |
        +----------+-------------------------------------------------------------------+
        | resample | Resample algorithm used for reading pixel data in                 |
        |          | case of interpolation                                             |
        |          | (default: GRIORA_NearestNeighbour). Check                         |
        |          | http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a |
        |          | or available options.                                             |
        +----------+-------------------------------------------------------------------+
        | nodata   | Nodata value to put in image if out of bounds                     |
        +----------+-------------------------------------------------------------------+

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

            jim=pj.Jim('/path/to/file.tif',t_srs='epsg:3035',
                       ulx=1000000,uly=4000000,lrx=1500000,lry=3500000)
            jim.warp('epsg:3035',s_srs='epsg:4326')

        """
        kwargs.update({'t_srs': t_srs})
        jim = _pj.Jim(_pj.geometry.cropPlane(
            self._jim_object, 0)._jipjim.warp(kwargs))
        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.Jim(_pj.geometry.cropPlane(
                self._jim_object, iplane)._jipjim.warp(kwargs))
            jim.geometry.stackPlane(jimplane)
        self._jim_object._set(jim._jipjim)
        # self._jim_object._set(self._jim_object._jipjim.warp(kwargs))


class _GeometryList():
    """Define all Geometry methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_list = caller

    def stackBand(self, jim_other=None, band=None):
        """Stack bands from raster datasets into new multiband Jim object.

        :param band: List of band indices to stack (index is 0 based)
        :return: multiband Jim object

        Create a multiband Jim object from a list of two Jim objects::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim1=pj.Jim('/path/to/raster1.tif')
            jim_stacked=pj.JimList([jim0,jim1]).geometry.stackBand()

        Create a multiband Jim object from a JimList object, selecting
        the first and third band::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim1=pj.Jim('/path/to/raster1.tif')
            jim2=pj.Jim('/path/to/raster2.tif')
            jimlist=pj.JimList([jim0,jim1,jim2])
            jim_stacked=jimlist.geometry.stackBand([0,2])
        """
        if band:
            retJim = _pj.Jim(
                self._jim_list._jipjimlist.stackBand({'band': band}))
        else:
            retJim = _pj.Jim(self._jim_list._jipjimlist.stackBand())

        if isinstance(jim_other, _pj.Jim):
            if band:
                return retJim.geometry.stackBand(jim_other, band=band)
            else:
                retJim.geometry.stackBand(jim_other)
        elif isinstance(jim_other, list):
            if band:
                jim_to_stack = _pj.Jim(
                    jim_other._jipjimlist.stackBand({'band': band}))
            else:
                jim_to_stack = _pj.Jim(jim_other._jipjimlist.stackBand())

            retJim = _pj.Jim(retJim._jipjim.stackBand(
                jim_to_stack._jipjim))

        return retJim

    def stackPlane(self):
        """Stack planes from raster datasets into new multiplane Jim object.

        :return: multiplane Jim object

        Stack planes from another raster dataset to current raster dataset.

        :return: multiplane Jim object

        Create a multiplane Jim object from a list of two Jim objects::

            jim0=pj.Jim('/path/to/raster0.tif')
            jim1=pj.Jim('/path/to/raster1.tif')
            jim_stacked=pj.JimList([jim0,jim1]).geometry.stackPlane()
        """
        return _pj.Jim(self._jim_list._jipjimlist.stackPlane())

    def extractOgr(self, sample, rule, output, **kwargs):
        """Extract pixel values based on a JimVect vector dataset.

        :param sample: reference JimVect instance
        :param rule: Rule how to calculate zonal statistics per feature
            (see list of :ref:`supported rules <extract_rules>`)
        :param output: Name of the output vector dataset in which the zonal
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
        | planename        | List of plane names corresponding to list of     |
        |                  | planes to extract                                |
        +------------------+--------------------------------------------------+
        | startband        | Start band sequence number (0 indexed)           |
        +------------------+--------------------------------------------------+
        | endband          | End band sequence number (0 indexed)             |
        +------------------+--------------------------------------------------+
        | plane            | List of planes to extract (0 indexed). Default is|
        |                  | to use extract all planes                        |
        +------------------+--------------------------------------------------+
        | planename        | List of plane name corresponding to list of      |
        |                  | planes to extract                                |
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

        :Supported key values to mask pixels that must be ignored in \
            the extraction process:

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | srcnodata        | List of nodata values not to extract             |
        +------------------+--------------------------------------------------+
        | buffer           | Buffer (in geometric units of raster dataset).   |
        |                  | Use neg. value to exclude pixels within buffer   |
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

        Extract the mean value of the pixels within the polygon of the provided
        reference vector. Exclude the pixels within a buffer of 10m of
        the polygon boundary. Use a temporary vector in memory for
        the calculation. Then write the result to the final destination
        on disk::

            jiml = pj.JimList([jim0,jim1,jim2])
            v = jiml.geometry.extractOgr(
                reference, bandname, buffer=-10, rule=['mean'],
                output='/vsimem/temp.sqlite', oformat='SQLite')
            v.io.write('/path/to/output.sqlite)

        """
        # make list of rules
        # if not isinstance(rule, Iterable) or isinstance(rule, basestring):
        #     rules = [rule]
        # else:
        #     rules=rule[:]
        # # rules = []
        # # if rule:
        # #     for irule in rule:
        # #         rules.append(irule)

        # kwargs.update({'output': output})
        # kwargs.update({'rule': rules})
        # if 'threshold' in kwargs:
        #     if '%' in kwargs['threshold']:
        #         kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
        #     else:
        #         kwargs['threshold'] = -kwargs['threshold']

        # if sample is None:
        #     raise Exception('Error: missing sample option')

        # firstExtract = True
        # joinfn = None

        # vsioutput = os.path.join('/vsimem', 'aggregate_polygon.sqlite')
        # filenames=[]
        # ijim=0
        # for jim in self._jim_list:
        #     print("ijim: {}".format(ijim))
        #     filenames.append('t'+str(ijim))
        #     print("filenames: {}".format(filenames))
        #     bandnames = kwargs.pop('bandname', None)
        #     if bandnames is None:
        #         bandnames = ['b'+str(iband) for iband in range(
        #             0, jim.properties.nrOfBand())]
        #     if jim.properties.nrOfPlane() > 1:
        #         if sample.properties.getLayerCount() > 1:
        #             raise Exception(
        #                 'Error: multiple layers not supported when '
        #                 'aggregating vectors over multi-plane raster '
        #                 'datasets, please use single layer vector object')

        #     intersectfn = '/vsimem/sampleintersect.sqlite'
        #     sampleintersect = _pj.geometry.intersect(
        #         sample, jim, output=intersectfn, oformat='SQLite',
        #         co=['OVERWRITE=YES'])
        #     if sampleintersect.properties.isEmpty():
        #         sampleintersect.io.close()
        #         raise Exception('intersect is empty')
        #     else:
        #         sampleintersect.io.write()

        #     print("debug 0 ijim: {}".format(ijim))
        #     ibandnames = []
        #     # try:
        #     if True:
        #         print("debug 1 ijim: {}".format(ijim))
        #         print("rules: {}".format(rules))
        #         for band in bandnames:
        #             if len(rules) > 1:
        #                 # rules are automatically pre-pended in extractogr
        #                 ibandnames.append(
        #                     '_' + filenames[ijim] + '_' + band)
        #             else:
        #                 ibandnames.append(
        #                     rules[0] + '_' + filenames[ijim] + '_' + band)

        #         fieldnames = sampleintersect.properties.getFieldNames()

        #         # try:
        #         if True:
        #             if 'buffer' in kwargs:
        #                 v = jim.geometry.extractOgr(
        #                     sampleintersect, rule=rules, output=vsioutput,
        #                     oformat='SQLite', co=['OVERWRITE=YES'],
        #                     bandname=ibandnames, copy=fieldnames, fid='fid',
        #                     buffer=kwargs['buffer'])
        #             else:
        #                 v = jim.geometry.extractOgr(
        #                     sampleintersect, rule=rules, output=vsioutput,
        #                     oformat='SQLite', co=['OVERWRITE=YES'],
        #                     bandname=ibandnames, copy=fieldnames, fid='fid')
        #             if v.properties.isEmpty():
        #                 v.io.close()
        #             sampleintersect.io.close()
        #         print("debug 2 ijim: {}".format(ijim))
        #         try:
        #             print("debug0")
        #         except:
        #             print("no coverage for jim {}, continue with next "
        #                   "product".format(filenames[ijim]))
        #             sampleintersect.io.close()
        #             if v:
        #                 v.io.close()
        #             continue
        #         print("debug 3 ijim: {}".format(ijim))
        #         if not v.properties.isEmpty():
        #             print("debug 4 ijim: {}".format(ijim))
        #             v.io.write()
        #             # join vectors
        #             if '/' not in output:
        #                 joinfn = os.path.join('/vsimem/', output)
        #             else:
        #                 joinfn = output
        #             # joinfn='/vsimem/vjoin.sqlite'

        #             if firstExtract:
        #                 vjoin = _pj.JimVect(v, output=joinfn,
        #                                     co='OVERWRITE=YES')
        #                 vjoin.io.write()
        #                 vjoin.io.close()
        #                 firstExtract = False
        #             else:
        #                 vprev = _pj.JimVect(joinfn)
        #                 vjoin = _pj.geometry.join(
        #                     vprev, v, joinfn, oformat='SQLite',
        #                     co=['OVERWRITE=YES'], key=['fid'],
        #                     method='OUTER_FULL')
        #                 vjoin.io.write()
        #                 vjoin.io.close()
        #             v.io.close()
        #         else:
        #             v.io.close()
        #     try:
        #         print("debug1")
        #     except:
        #         print("raised exception dataset")
        #         continue
        #     print("debug 5 ijim: {}".format(ijim))
        #     ijim=ijim+1
        #     print("debug 6 ijim: {}".format(ijim))

        # if joinfn:
        #     v = _pj.JimVect(joinfn)
        #     return v
        # else:
        #     print("Error: joinfn is None, no valid features found")
        #     raise Exception('Error: joinfn is None, no valid features found')
        #############
        if isinstance(sample, _pj.JimVect):
            kwargs.update({'rule': rule})
            kwargs.update({'output': output})
            if 'threshold' in kwargs:
                if '%' in kwargs['threshold']:
                    kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
                else:
                    kwargs['threshold'] = -kwargs['threshold']

            bandname = kwargs.pop('bandname', None)

            #todo: support multi-band images in JimList...
            if bandname is None:
                bandname = [
                    't'+str(ifile) for ifile in range(0, len(self._jim_list))]
            # elif isinstance(bandname,list):
            #     if len(bandname) != len(self._jim_list):
            #         raise ValueError('Error: len of bandname should be '
            #                          'identical to len of JimList')
            kwargs.update({'bandname': bandname})

            avect = self._jim_list._jipjimlist.extractOgr(sample._jipjimvect,
                                                          kwargs)
            pjvect = _pj.JimVect()
            pjvect._set(avect)
            return pjvect
        else:
            raise TypeError('Error: extractOgr must operate on vector sample'
                            'of type JimVect')


class _GeometryVect():
    """Define all Geometry methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self, caller):
        self._jim_vect = caller

    # def append(self, jvec):
    #     """Append JimVect object with another JimVect object.

    #     :param jvec: JimVect object to append
    #     """
    #     if isinstance(jvec, _pj.JimVect):
    #         self._jim_vect._jipjimvect.append(jvec._jipjimvect)
    #         # return pjvect
    #         # return _pj.JimVect(
    #               # self._jim_vect._jipjimvect.join(jvec._jipjimvect,kwargs))
    #     else:
    #         raise TypeError('Error: can only join with JimVect object')

    def convexHull(self, **kwargs):
        """Create the convex hull on a JimVect object.

        Modifies the instance on which the method was called.

        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+
        """
        non_existing_path = _pj._get_random_path()

        kwargs.update({'output': non_existing_path})
        avect = self._jim_vect._jipjimvect.convexHull(kwargs)
        self._jim_vect._set(avect)

    def intersect(self, jim, **kwargs):
        """Intersect JimVect object with Jim object.

        Keeps only those features with an intersect.

        Modifies the instance on which the method was called.

        :param jim: Jim object with which to intersect
        :param kwargs: See table below
        :return: intersected JimVect object

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+

        Example: intersect a sample with a Jim object::

          jim = pj.Jim('/path/to/raster.tif')
          v = pj.JimVect('/path/to/vector.sqlite')
          sampleintersect = pj.geometry.intersect(
              v, jim, output='/vsimem/intersect', oformat='SQLite',
              co=['OVERWRITE=YES'])
          sampleintersect.io.write('/path/to/output.sqlite')

        """
        non_existing_path = _pj._get_random_path()

        kwargs.update({'output': non_existing_path})
        if isinstance(jim, _pj.Jim):
            avect = self._jim_vect._jipjimvect.intersect(jim._jipjim, kwargs)
            self._jim_vect._set(avect)
        else:
            raise TypeError('Error: can only intersect with Jim object')

    def join(self, jvec2, **kwargs):
        """Join JimVect object with another JimVect object.

        A key field is used to find corresponding features in both objects.

        :param jve2c: second JimVect object to join
        :param kwargs: See table below
        :return: joined JimVect object

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

        .. |inner| image:: figures/join_inner.png
             :width: 20 %
        .. |outer_left| image:: figures/join_outer_left.png
             :width: 20 %
        .. |outer_right| image:: figures/join_outer_right.png
             :width: 20 %
        .. |outer_full| image:: figures/join_outer_full.png
             :width: 20 %

        The join methods currently supported are:

            :INNER |inner|: join two JimVect objects, keeping only those \
                            features for which identical keys in both objects \
                            are found
            :OUTER_LEFT |outer_left|: join two JimVect objects, keeping all \
                                      features from first object
            :OUTER_RIGHT |outer_right|: join two JimVect objects, keeping all \
                                        features from second object
            :OUTER_FULL |outer_full|: join two JimVect objects, keeping all \
                                      features from both objects

        Example: join two vectors, based on the key 'id', which is a common
        field shared between v1 and v2. Use OUTER_FULL as the join method::

          v1 = pj.JimVect('/path/to/vector1.sqlite')
          v2 = pj.JimVect('/path/to/vector2.sqlite')
          v3 = pj.geometry.join(
              v1, v2, '/tmp/test.sqlite', oformat='SQLite',
              co=['OVERWRITE=YES'], key=['id'], method='OUTER_FULL')
        """
        non_existing_path = _pj._get_random_path()

        kwargs.update({'output': non_existing_path})
        if isinstance(jvec2, _pj.JimVect):
            avect = self._jim_vect._jipjimvect.join(jvec2._jipjimvect, kwargs)
            self._jim_vect._set(avect)
        else:
            raise TypeError('Error: can only join two JimVect objects')
