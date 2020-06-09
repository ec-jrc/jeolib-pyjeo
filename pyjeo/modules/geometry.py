"""Module for operations working with the geometry of the Jim objects."""

import os
import numpy as _np
import warnings as _warnings

import pyjeo as _pj


def append(jvec1,
           jvec2,
           output: str,
           **kwargs):
    """Append JimVect object with another JimVect object.

    :param jvec1: first JimVect object to append
    :param jvec2: second JimVect object to append
    :param output: output filename of JimVect object that is returned.
        Use /vsimem for in memory vectors

    Example: append two vectors

      v1 = pj.JimVect('/path/to/vector1.sqlite')
      v2 = pj.JimVect('/path/to/vector2.sqlite')
      v3 = pj.geometry.append(
          v1, v2, '/tmp/test.sqlite', oformat='SQLite',
          co=['OVERWRITE=YES'])
    """
    kwargs.update({'output': output})
    if isinstance(jvec1, _pj.JimVect) and isinstance(jvec2, _pj.JimVect):
        avect = jvec1._jipjimvect.merge(jvec2._jipjimvect, kwargs)
        avect.write()
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only append two JimVect objects')

def band2plane(jim):
    """Convert 2-dimensional multi-band object to a 3-dimensional.

    The result will be a single band multi-plane object

    Example: convert a multi-band object with 12 bands to a 3-dimensional
    single band object with 12 planes::

        jim2d = pj.Jim('/path/to/multi/band/image.tif')
        jim2d.properties.nrOfBand()
        12
        jim2d.properties.nrOfPlane()
        1
        jim3d = pj.geometry.band2plane(jim2d)
        jim3d.properties.nrOfPlane()
        12
        jim3d.properties.nrOfBand()
        1

    Notice that a multi-band image can also be read directly as
    a multi-plane object::

        jim3d = pj.Jim('/path/to/multi/band/image.tif', band2plane=True)
        jim3d.properties.nrOfBand()
        1
        jim3d.properties.nrOfPlane()
        12

    """
    result = _pj.Jim(jim)
    result.geometry.band2plane()
    return result


def convexHull(jim_vect,
               output: str,
               **kwargs):
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
    avect.write()

    pjvect = _pj.JimVect()
    pjvect._set(avect)

    return pjvect


def crop(jim_object,
         ulx: float = None,
         uly: float = None,
         ulz: float = None,
         lrx: float = None,
         lry: float = None,
         lrz: float = None,
         dx: float = None,
         dy: float = None,
         nogeo: bool = False, **kwargs):
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


def cropBand(jim_object,
             band: int):
    """Subset raster dataset.

    Subset raster dataset in spectral domain.

    :param jim_object: a Jim object
    :param band: List of band indices to crop (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three bands from raster dataset jim::

        jim = pj.Jim('/path/to/raster.tif')
        jim3 = pj.geometry.cropBand(jim, band=[0, 1, 2])

    """
    if isinstance(band, list):
        bands = band
    else:
        bands = [band]
    band = [jim_object.properties.nrOfBand() + b if b < 0 else b
            for b in bands]
    return _pj.Jim(jim_object._jipjim.cropBand({'band': band}))


def cropOgr(jim_object,
            extent,
            **kwargs):
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

       For instance you can use 'eo':'ATTRIBUTE=fieldname' to burn
           the (numeric) fieldname in the pixel value'
    """
    jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.cropOgr(
        extent._jipjimvect, kwargs))
    for iplane in range(1, jim_object.properties.nrOfPlane()):
        jimplane = _pj.Jim(_pj.geometry.cropPlane(
            jim_object, iplane)._jipjim.cropOgr(extent._jipjimvect, kwargs))
        jim.geometry.stackPlane(jimplane)
    return jim


def cropPlane(jim_object,
              plane: int):
    """Subset raster dataset.

    Subset raster dataset in temporal domain.

    :param jim_object: a Jim object
    :param plane: List of plane indices to crop (index is 0 based)
    :return: Cropped subimage as Jim instance

    Example:

    Crop the first three planes from raster dataset jim::

        jim = pj.Jim('/path/to/raster.tif')
        jim3 = pj.geometry.cropPlane(jim, plane=[0, 1, 2])

    """
    if isinstance(plane, list):
        planes = plane
    else:
        planes = [plane]
    plane = [jim_object.properties.nrOfPlane() + b if b < 0 else b
             for b in planes]
    return _pj.Jim(jim_object._jipjim.cropPlane({'plane': plane}))


def extract(jvec,
            jim,
            output: str,
            rule: str = None,
            **kwargs):
    """Extract pixel values from raster image based on a vector dataset.

    :param jvec: overlay JimVec object or Jim thematic raster
    :param jim: Jim object (list) on which vector is overlaid to extract from
    :param rule: Rule how to calculate zonal statistics per feature
        (see list of :ref:`supported rules <extract_rules>`)
    :param output: Name of the output vector dataset in which the zonal
        statistics will be saved
    :param kwargs: See table below
    :return: A JimVect with the same geometry as the sample vector
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
    | classes          | Only when overlaying Jim thematic raster dataset |
    |                  | dataset.                                         |
    +------------------+--------------------------------------------------+
    | fid              | Create extra field named 'fid' with this field   |
    |                  | identifier (sequence of features)                |
    +------------------+--------------------------------------------------+
    | bandname         | List of band names corresponding to list of      |
    |                  | bands to extract                                 |
    +------------------+--------------------------------------------------+
    | planename        | List of plane names corresponding to list of     |
    |                  | planes to extract                                |
    +------------------+--------------------------------------------------+
    | oformat          | Output vector dataset format                     |
    +------------------+--------------------------------------------------+
    | co               | Creation option for output vector dataset        |
    +------------------+--------------------------------------------------+

    .. note::
        To ignore some pixels from the extraction process, see list
        of :ref:`mask key values <extract_mask>`

    Example:

    Extract the mean value of the pixels within the polygon of the provided
    reference vector. Exclude the pixels within a buffer of 10m of
    the polygon boundary. Use a temporary vector in memory for
    the calculation. Then write the result to the final destination
    on disk::

        reference = pj.JimVect('/path/to/reference.sqlite')
        jim0 = pj.Jim('/path/to/raster.tif')
        v = reference.geometry.extract(jim0,
            buffer=-10, rule=['mean'],
            output='/vsimem/temp.sqlite', oformat='SQLite')
        v.io.write('/path/to/output.sqlite)

    Open a raster sample dataset based on land cover map (e.g., Corine) and
    use it to extract a stratified sample of 100 points from an input
    raster dataset with four spectral bands ('B02', 'B03', 'B04', 'B08').
    Only sample classes 2 (urban), 12 (agriculture), 25 (forest),
    41 (water) and an aggregated (rest) class 50::

        reference=pj.Jim('/path/to/landcovermap.tif')

        classes=[2,12,25,41,50]
        thresholds=['20%','25%','25%','10%','5%']

        jim=pj.Jim('/path/to/s2_multiband.tif')

        outputfn='/path/to/output.sqlite'
        sample=pj.geometry.extract(jim, reference, srcnodata=[0],
                                    output=outputfn,
                                    classes=classes,
                                    threshold=thresholds,
                                    bandname=['B02','B03','B04','B08'])
    """
    kwargs.update({'output': output})
    if 'threshold' in kwargs:
        if kwargs['threshold'] is not None:
            if not isinstance(kwargs['threshold'],list):
                kwargs['threshold']=[kwargs['threshold']]
            kwargs['threshold']=[float(t.strip('%')) if isinstance(t,str) else -t for t in kwargs['threshold']]

    if 'classes' in kwargs:
        classes = kwargs.pop('classes')
        kwargs['class'] = classes

    pjvect = _pj.JimVect()
    if isinstance(jvec, _pj.JimVect):
        kwargs.update({'rule': rule})
        if isinstance(jim, _pj.Jim):
            avect = jim._jipjim.extractOgr(jvec._jipjimvect, kwargs)
        elif isinstance(jim, _pj.JimList):
            bandname = kwargs.pop('bandname', None)
            #todo: support multi-band images in JimList...
            if bandname is None:
                bandname = [
                    't'+str(ifile) for ifile in range(0, len(jim._jim_list))]
            kwargs.update({'bandname': bandname})

            avect = jim._jipjimlist.extractOgr(jvec._jipjimvect, kwargs)
        else:
            raise TypeError('Error: extract must operate on Jim or JimList')
    elif isinstance(jvec, _pj.Jim):
        avect = jim._jipjim.extractImg(jvec._jipjim, kwargs)
    else:
        raise TypeError('Error: sample to extract must be either JimVect or Jim')
    avect.write()
    pjvect._set(avect)
    return pjvect


def geo2image(jim_object,
              x: float,
              y: float):
    """Convert georeferenced coordinates to image coordinates (col, row).

    :param jim_object: a Jim object
    :param x: georeferenced coordinate in x according to the object spatial
        reference system
    :param y: georeferenced coordinate in y according to the object spatial
        reference system
    :return: image coordinates (row and column, starting from 0)

    Get column and row index (0 based) of some georeferenced coordinates
    x and y (in this case first pixel: 0, 0)::

      jim = pj.Jim('/path/to/raster.tif')
      x = jim.properties.getUlx()
      y = jim.properties.getUly()
      pj.geometry.geo2image(jim, x, y)
    """
    coords = jim_object._jipjim.geo2image(x, y)
    return [int(coords[0]), int(coords[1])]


def image2geo(jim_object,
              i: int,
              j: int):
    """Convert image coordinates (col, row) to georeferenced coordinates.

    :param jim_object: a Jim object
    :param i: image column number (starting from 0)
    :param j: image row number (starting from 0)
    :return: georeferenced coordinates according to the object spatial
        reference system

    Get upper left corner in georeferenced coordinates
    (in SRS of the Jim object)::

      jim = pj.Jim('/path/to/raster.tif')
      pj.geometry.image2geo(jim, 0, 0)

    """
    return jim_object._jipjim.image2geo(i, j)


def imageFrameAdd(jim_object,
                  l: int = 0,
                  r: int = 0,
                  t: int = 0,
                  b: int = 0,
                  u: int = 0,
                  d: int = 0,
                  val: int = 0):
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
        ret_jim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if ret_jim:
                jimband = _pj.geometry.cropBand(jim_object, band=band)
                jimband._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
                ret_jim.geometry.stackBand(jimband)
            else:
                ret_jim = _pj.geometry.cropBand(jim_object, band=band)
                ret_jim._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
        return _pj.Jim(ret_jim)
    else:
        return _pj.Jim(jim_object._jipjim.imageFrameAdd(
            [l, r, t, b, u, d], val))


def imageFrameSet(jim_object,
                  l: int = 0,
                  r: int = 0,
                  t: int = 0,
                  b: int = 0,
                  u: int = 0,
                  d: int = 0,
                  val: int = 0,
                  band: int = None):
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


def imageFrameSubtract(jim_object,
                       l: int = 0,
                       r: int = 0,
                       t: int = 0,
                       b: int = 0,
                       u: int = 0,
                       d: int = 0):
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
        ret_jim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if ret_jim:
                jimband = _pj.geometry.cropBand(jim_object,
                                                band=band)
                jimband._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
                ret_jim.geometry.stackBand(jimband)
            else:
                ret_jim = _pj.geometry.cropBand(jim_object,
                                                band=band)
                ret_jim._jipjim.d_imageFrameSubtract(
                    [l, r, t, b, u, d])
        return _pj.Jim(ret_jim)
    else:
        return _pj.Jim(jim_object._jipjim.imageFrameSubtract(
            [l, r, t, b, u, d]))


def imageInsert(jim_object,
                sec_jim_object,
                x: float, y: float,
                z: float,
                band: list = None):
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

    ret_jim = None

    for band in bands:
        if ret_jim:
            jimband = jim_object._jipjim.imageInsert(
                sec_jim_object._jipjim, x, y, z, band)
            ret_jim.geometry.stackBand(jimband)
        else:
            ret_jim = jim_object._jipjim.imageInsert(
                sec_jim_object._jipjim, x, y, z, band)

    return _pj.Jim(ret_jim)


def imageInsertCompose(jim_object,
                       imlbl,
                       im2,
                       x: float,
                       y: float,
                       z: float,
                       val: int,
                       band: list = None):
    """Merge Jim instance with values of im2 if val of imlbl == val.

    :param jim_object: a Jim object
    :param imlbl: a Jim object
    :param im2: a Jim object
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

    ret_jim = None

    for band in bands:
        if ret_jim:
            jimband = jim_object._jipjim.imageInsertCompose(
                imlbl._jipjim, im2._jipjim, x, y, z, val, band)
            ret_jim.geometry.stackBand(jimband)
        else:
            ret_jim = jim_object._jipjim.imageInsertCompose(
                imlbl._jipjim, im2._jipjim, x, y, z, val, band)

    return _pj.Jim(ret_jim)


def intersect(jvec,
              jim,
              output: str,
              **kwargs):
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
        avect.write()
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only intersect with Jim object')


def join(jvec1,
         jvec2,
         output: str,
         **kwargs):
    """Join JimVect object with another JimVect object.

    A key field is used to find corresponding features in both objects.

    :param jvec1: first JimVect object to join
    :param jvec2: second JimVect object to join
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
        :OUTER_LEFT |outer_left|: join two JimVect objects, keeping \
                                  all features from first object
        :OUTER_RIGHT |outer_right|: join two JimVect objects, keeping all \
                                    features from second object
        :OUTER_FULL |outer_full|: join two JimVect objects, keeping \
                                  all features from both objects

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
        avect.write()
        pjvect = _pj.JimVect()
        pjvect._set(avect)
        return pjvect
    else:
        raise TypeError('Error: can only join two JimVect objects')


def magnify(jim_object,
            n: int):
    """Magnify the image.

    :param jim_object: a Jim object
    :param n: a positive integer for magnifying size by pixel replication
    :return: a Jim object containing the magnified image
    """
    if jim_object.properties.nrOfBand() > 1:
        ret_jim = None
        for band in range(0, jim_object.properties.nrOfBand()):
            if ret_jim:
                jimband = _pj.geometry.cropBand(jim_object,
                                                band=band)
                jimband = _pj.Jim(jimband._jipjim.imageMagnify(n))

                ret_jim.geometry.stackBand(jimband)
            else:
                ret_jim = _pj.Jim(_pj.geometry.cropBand(
                    jim_object, band=band)._jipjim.imageMagnify(n))

        return ret_jim
    else:
        return _pj.Jim(jim_object._jipjim.imageMagnify(n))


def plane2band(jim):
    """Convert 3-dimensional 1-band Jim to a 2-dimensional multi-band one.

    The result will be a multi-band single plane object

    Example: convert a single band object with 12 planes to a 2-dimensional
    multi-band object with 1 plane::

        jim3d = pj.Jim('/path/to/multi/band/image.tif', band2plane=True)
        jim3d.properties.nrOfBand()
        1
        jim3d.properties.nrOfPlane()
        12
        jim2d = pj.geometry.plane2band(jim3d)
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


def plotLine(jim_object,
             x1: float,
             y1: float,
             x2: float,
             y2: float,
             val: float):
    """Draw a line from [x1, y1] to [x2, y2] by setting pixels of Jim to val.

    Works only for 1-plane Jims.

    :param jim_object: a Jim object
    :param x1: an integer for x-coordinate of 1st point
    :param y1: an integer for y-coordinate of 1st point
    :param x2: an integer for x-coordinate of 2nd point
    :param y2: an integer for y-coordinate of 2nd point
    :param val: value to be used for line pixels
    :return: a Jim object
    """
    if jim_object.properties.nrOfPlane() > 1:
        raise TypeError(
            'Error: plotLine does not support multi-plane Jim objects')
    return _pj.Jim(jim_object._jipjim.plotLine(x1, y1, x2, y2, val))


def polygonize(jim_object,
               output: str,
               **kwargs):
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


def rasterize(jim_object,
              jim_vect,
              burn_value: int = 1,
              eo: list = None,
              ln: str = None):
    """Rasterize Jim object based on GDALRasterizeLayersBuf.

    :param jim_object: a template Jim object
    :param jim_vect: JimVect object that needs to be rasterized
    :param burn_value: burn value
    :param eo: option (default is ALL_TOUCHED)
    :param ln: layer names (optional)
    :return: rasterized Jim object

    .. note::
       Possible values for the key 'eo' are:

       ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

       For instance you can use 'eo':'ATTRIBUTE=fieldname to burn
       the (numeric) fieldname in the pixel value'
    """
    if not isinstance(jim_object, _pj.Jim):
        raise TypeError('Error: template must be a Jim object')
    if not isinstance(jim_vect, _pj.JimVect):
        raise TypeError('Error: can only rasterize a JimVect')

    kwargs = {}
    if burn_value is not None:
        kwargs.update({'burn': float(burn_value)})
    if ln is not None:
        kwargs.update({'ln': ln})

    if eo is not None:
        kwargs.update({'eo': eo})
    else:
        kwargs.update({'eo': ['ALL_TOUCHED']})

    ajim = _pj.Jim(jim_object, copy_data=False)
    ajim._jipjim.d_rasterizeBuf(jim_vect._jipjimvect, kwargs)
    return ajim


def reducePlane(jim,
                rule: str = 'overwrite',
                ref_band: int = None,
                nodata: float = None):
    """Reduce planes of Jim object.

    :param jim: jim object on which to reduce planes
    :param rule: rule to reduce (mean, median, min or max)
        or callback function
    :param ref_band: band on which to apply rule
        (default is to check all bands,
        not supported when rule is callback function)
    :param nodata: value to ignore when applying rule
        (not supported when rule is callback function)
    :return: reduced single plane jim object
    """
    nr_of_planes = jim.properties.nrOfPlane()
    if jim.properties.nrOfPlane() < 2:
        _warnings.warn(
            'Single-plane Jim: No plane reduction performed', SyntaxWarning
        )
        return jim

    jimreduced = _pj.geometry.cropPlane(jim, 0)

    if isinstance(rule, str):
        nr_of_row = jim.properties.nrOfRow()
        nr_of_col = jim.properties.nrOfCol()

        if rule in ('mean', 'avg', 'median'):
            if rule == 'median':
                nan_func = _np.nanmedian
                func = _np.median
            else:
                nan_func = _np.nanmean
                func = _np.mean

            d_type = jim.properties.getDataType()
            nr_of_bands = jim.properties.nrOfBand()
            if nodata is not None and d_type not in ('GDT_Float32',
                                                     'GDT_Float64'):
                jim = _pj.pixops.convert(jim, otype='GDT_Float32')
            if ref_band is not None:
                # create plane-wise mask for the ref_band
                mask = _pj.geometry.cropBand(jim,
                                             band=ref_band)

                # create one-plane boolean mask containing True where
                # nodata value is in all planes
                planes = [mask.np()[i] for i in range(nr_of_planes)]
                stacked_planes = _np.vstack(planes)

                diffs = _np.abs(
                    _np.diff(stacked_planes.reshape(len(planes), -1),
                             axis=0)) == 0

                same_values = diffs[0]
                for i in range(1, len(diffs)):
                    same_values = _np.logical_and(same_values, diffs[i])

                same_values = _np.reshape(same_values,
                                          (nr_of_row, nr_of_col))

                nodata_mask = same_values & (mask.np()[0] == nodata)

            for iband in range(0, nr_of_bands):
                if nodata is not None:
                    if ref_band is None:
                        # create one-plane boolean mask containing True
                        # where nodata value is in all planes
                        planes = [jim.np(iband)[i] for i
                                  in range(nr_of_planes)]
                        stacked_planes = _np.vstack(planes)

                        diffs = _np.abs(
                            _np.diff(
                                stacked_planes.reshape(len(planes), -1),
                                axis=0)) == 0

                        same_values = diffs[0]
                        for i in range(1, len(diffs)):
                            same_values = _np.logical_and(same_values,
                                                          diffs[i])

                        same_values = _np.reshape(same_values,
                                                  (nr_of_row, nr_of_col))

                        nodata_mask = same_values & (
                            jim.np(iband)[0] == nodata)

                        # compute the reduction function
                        jim.np(iband)[jim.np(iband) == nodata] = _np.nan
                        jimreduced.np(iband)[:] = nan_func(
                            jim.np(iband), axis=0)
                        jimreduced.np(iband)[nodata_mask] = nodata
                    else:
                        jim.np(iband)[
                            mask.np() == nodata] = _np.nan
                        jimreduced.np(iband)[:] = nan_func(
                            jim.np(iband), axis=0)
                        jimreduced.np(iband)[nodata_mask] = nodata
                else:
                    jimreduced.np(iband)[:] = func(
                        jim.np(iband), axis=0)
            if nodata is not None:
                if d_type not in ('GDT_Float32', 'GDT_Float64'):
                    jimreduced.pixops.convert(otype=d_type)
                    jim.pixops.convert(otype=d_type)
        else:
            if rule == 'max':
                def rule(reduced, plane):
                    return reduced < plane
            elif rule == 'min':
                def rule(reduced, plane):
                    return reduced > plane
            elif rule == 'overwrite':
                def rule(reduced, plane):
                    return _pj.pixops.setData(plane, 1)
            else:
                raise AttributeError(
                    'Error: rule not supported')

            if ref_band is not None:
                maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)

            for iplane in range(1, nr_of_planes):
                jimplane = _pj.geometry.cropPlane(jim, iplane)

                if nodata is not None and ref_band is None:
                    raise AttributeError(
                        'Error: use ref_band option for nodata')

                if ref_band is not None:
                    maskplane = _pj.geometry.cropBand(jimplane, ref_band)
                    themask = rule(maskreduced, maskplane)
                    if nodata is not None:
                        themask |= maskreduced == nodata
                        themask &= maskplane != nodata
                    maskreduced[themask] = maskplane
                else:
                    themask = rule(jimreduced, jimplane)

                jimreduced[themask] = jimplane

                if nodata is not None:
                    nodata_mask = (maskreduced == nodata) & \
                                  (maskplane == nodata)
                    jimreduced[nodata_mask] = nodata
    else:
        if nodata is not None or ref_band is not None:
            # TODO: Allow nodata and ref_band parameters (see #48)
            raise AttributeError('Error: nodata and ref_band are not '
                                 'supported for callback rules')
        jimreduced = _pj.geometry.cropPlane(jim, 0)
        for iplane in range(1, nr_of_planes):
            jimplane = _pj.geometry.cropPlane(jim, iplane)
            jimreduced = rule(jimreduced, jimplane)
        # else:
        #     if ref_band is not None:
        #         maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)

        #     for iplane in range(1, nr_of_planes):
        #         jimplane = _pj.geometry.cropPlane(jim, iplane)

        #         if nodata is not None and ref_band is None:
        #             raise AttributeError(
        #                 'Error: use ref_band option for nodata')

        #         if ref_band is not None:
        #             maskplane = _pj.geometry.cropBand(jimplane, ref_band)
        #             themask = rule(maskreduced, maskplane)
        #             if nodata is not None:
        #                 themask |= maskreduced == nodata
        #                 themask &= maskplane != nodata
        #             maskreduced[themask] = maskplane
        #         else:
        #             themask = rule(jimreduced, jimplane)

        #         jimreduced[themask] = jimplane

        #         if nodata is not None:
        #             nodata_mask = (maskreduced == nodata) & \
        #                           (maskplane == nodata)
        #             jimreduced[nodata_mask] = nodata

    return jimreduced


def sample(jim,
           output: str,
           **kwargs):
    """Extract a random or grid sample from a raster dataset.

    :jim: Jim object to sample (multi-band supported, but multi-plane not yet)
    :output: Name of the output vector dataset in which the zonal
        statistics will be saved
    :param output: path to the output
    :param kwargs: See table below
    :return: A JimVect with sample

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

        v01 = jim0.sample(random=100, buffer=100, rule=['mean'],
                          output='mem01', oformat='Memory')

    Extract a sample of 100 points using a regular grid sampling scheme.
    For each grid point, calculate the median value based on a 3x3 window
    (buffer value of 100 m neighborhood). Write the result in a SQLite
    vector dataset on disk::

        outputfn = '/path/to/output.sqlite'
        npoint = 100
        gridsize = int(jim.nrOfCol()*jim.getDeltaX()/math.sqrt(npoint))
        v = jim.sample(grid=gridsize, buffer=100, rule=['median'],
                       output=outputfn, oformat='SQLite')
        v.write()

    """
    kwargs.update({'output': output})
    if 'threshold' in kwargs:
        if kwargs['threshold'] is not None:
            if '%' in kwargs['threshold']:
                kwargs['threshold'] = float(kwargs['threshold'].strip('%'))
            else:
                kwargs['threshold'] = -kwargs['threshold']

    pjvect = _pj.JimVect()
    pjvect._set(jim._jipjim.extractSample(kwargs))
    return pjvect

def stackBand(jim_object,
              jim_other=None,
              band: int = None):
    """Stack bands from raster datasets into new multiband Jim object.

    :param jim_object: a Jim or JimList object used for stacking the bands
    :param jim_other: a Jim object or jimlist from which to copy bands
        (optional)
    :param band: List of band indices to stack (index is 0 based).
        Default is to stack all bands.
    :return: Jim object with stacked bands

    Append all the bands of jim1 to jim0::

        jim0 = pj.Jim('/path/to/raster0.tif')
        jim1 = pj.Jim('/path/to/raster1.tif')
        jim_stacked = pj.geometry.stackBand(jim0, jim1)

    Stack all bands of the JimList, returning a multi-band Jim object::

        jim0 = pj.Jim('/path/to/raster0.tif')
        jim1 = pj.Jim('/path/to/raster1.tif')
        jimlist = pj.JimList([jim0, jim1])
        jim_stacked = pj.geometry.stackBand(jimlist)

    Append the first three bands of raster dataset jim1 to the image jim0::

        jim0 = pj.Jim('/path/to/raster0.tif')
        jim1 = pj.Jim('/path/to/raster1.tif')
        jim_stacked = pj.geometry.stackBand(jim0, jim1, band=[0, 1, 2])

    """
    def _check_number_of_bands(queried_band: int, nr_of_band: int):
        """Raise an error if the user wants to use a band out of bounds."""
        if queried_band is not None and queried_band > nr_of_band:
            raise AttributeError('Band out of bounds')

    if isinstance(jim_object, _pj.JimList):
        if band:
            nr_of_band = jim_object[0].properties.nrOfBand()
            _check_number_of_bands(band, nr_of_band)

            ret_jim = _pj.Jim(
                jim_object._jipjimlist.stackBand({'band': band}))
        else:
            ret_jim = _pj.Jim(jim_object._jipjimlist.stackBand())

        if isinstance(jim_other, _pj.Jim):
            if band:
                ret_jim.geometry.stackBand(jim_other, band=band)
            else:
                ret_jim.geometry.stackBand(jim_other)
        elif isinstance(jim_other, list):
            if band:
                jim_to_stack = _pj.Jim(
                    jim_other._jipjimlist.stackBand({'band': band}))
            else:
                jim_to_stack = _pj.Jim(jim_other._jipjimlist.stackBand())

            ret_jim = _pj.Jim(ret_jim._jipjim.stackBand(
                jim_to_stack._jipjim))
    elif isinstance(jim_object, _pj.Jim):
        _check_number_of_bands(band, jim_object.properties.nrOfBand())

        ret_jim = jim_object
        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            if band:
                ret_jim = _pj.Jim(ret_jim._jipjim.stackBand(
                    jim._jipjim, {'band': band}))
            else:
                ret_jim = _pj.Jim(ret_jim._jipjim.stackBand(
                    jim._jipjim))
    else:
        raise TypeError('Error: expected a Jim or JimList object as '
                        'the first argument')

    return ret_jim


def stackPlane(jim_object,
               jim_other=None,
               *args):
    """Stack planes from raster datasets into new multiplane Jim object.

    :param jim_object: a Jim or JimList object used for stacking the planes
    :param jim_other: a Jim object or jimlist from which to copy planes
        (optional)
    :return: Jim object with stacked planes

    Append all the planes of jim1 to jim0::

        jim0 = pj.Jim('/path/to/raster0.tif')
        jim1 = pj.Jim('/path/to/raster1.tif')
        jim_stacked = pj.geometry.stackPlane(jim0, jim1)

    Stack all planes of the JimList, returning a multi-plane Jim object::

        jim0 = pj.Jim('/path/to/raster0.tif')
        jim1 = pj.Jim('/path/to/raster1.tif')
        jimlist = pj.JimList([jim0, jim1])
        jim_stacked = pj.geometry.stackPlane(jimlist)
    """
    if not isinstance(jim_other, list):
        args_list = [jim_other]
    else:
        args_list = jim_other

    if args:
        args_list.extend(args)

    if isinstance(jim_object, _pj.JimList):
        ret_jim = _pj.Jim(jim_object._jipjimlist.stackPlane())

        if jim_other is not None:
            for jim in args_list:
                ret_jim._jipjim.d_stackPlane(jim._jipjim)
    elif isinstance(jim_object, _pj.Jim):
        if jim_other is None:
            return _pj.Jim(jim_object)

        ret_jim = _pj.Jim(jim_object)

        for jim in args_list:
            ret_jim._jipjim.d_stackPlane(jim._jipjim)
    else:
        raise TypeError('Error: expected a Jim or JimList object as '
                        'the first argument')

    return ret_jim


def warp(jim_object,
         t_srs,
         **kwargs):
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
        jim_warped = pj.geometry.warp(jim, 'epsg:3035')

    Read a raster dataset from disk that is in lat lon (epsg:4326), select
    a bounding box in a different spatial reference system (epsg:3035).
    Notice the raster dataset read is still in the original
    projection (epsg:4326). Then warp the raster dataset to the target
    spatial reference system (epsg:3035)::

        jim = pj.Jim('/path/to/file.tif', t_srs='epsg:3035', ulx=1000000,
                     uly=4000000, lrx=1500000, lry=3500000)
        jim_warped = pj.geometry.warp(jim, 'epsg:3035', s_srs='epsg:4326')

    """
    kwargs.update({'t_srs': t_srs})

    jim = _pj.Jim(_pj.geometry.cropPlane(jim_object, 0)._jipjim.warp(kwargs))
    for iplane in range(1, jim_object.properties.nrOfPlane()):
        jimplane = _pj.Jim(_pj.geometry.cropPlane(jim_object,
                                                  iplane)._jipjim.warp(kwargs))
        jim.geometry.stackPlane(jimplane)
    return jim
    # return _pj.Jim(jim_object._jipjim.warp(kwargs))


class _Geometry(_pj.modules.JimModuleBase):
    """Define all Geometry methods."""

    def band2plane(self):
        """Convert 2-dimensional multi-band object to a 3-dimensional.

        The result will be a single band multi-plane object

        Example: convert a multi-band object with 12 bands to a 3-dimensional
        single band object with 12 planes::

           jim = pj.Jim('/path/to/multi/band/image.tif')
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

           jim = pj.Jim('/path/to/multi/band/image.tif', band2plane=True)
           jim.properties.nrOfBand()
           1
           jim.properties.nrOfPlane()
           12

        """
        self._jim_object._jipjim.d_band2plane()

    def crop(self,
             ulx: float = None,
             uly: float = None,
             ulz: float = None,
             lrx: float = None,
             lry: float = None,
             lrz: float = None,
             dx: float = None,
             dy: float = None,
             nogeo: bool = False,
             **kwargs):
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

            jim = pj.Jim('/path/to/raster.tif')
            jim.crop(ulx=1000000, uly=5000000, lrx=2000000, lry=4000000,
                     dx=1000, dy=1000)

        Crop bounding box in image coordinates (starting from upper left pixel
        coordinate 0, 0). For instance, get first 10 columns in first 10 rows::

            jim = pj.Jim('/path/to/raster.tif')
            jim.crop(ulx=0, uly=0, lrx=10, lry=10, nogeo=True)

        Notice that for this case, a more pythonic way to is available
        via :ref:`indexing`::

            jim0[0:10, 0:10]

        However, crop can also be used to enlarge a Jim object. For instance,
        to add a border of one pixel use::

            jim = pj.Jim('/path/to/raster.tif')
            jim.geometry.crop(ulx=-1, uly=-1, lrx=jim.properties.nrOfCol()+1,
                              lry=jim.properties.nrOfRow()+1, nogeo=True)

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

    def cropBand(self,
                 band: int):
        """Subset raster dataset.

        Subset raster dataset in spectral domain.

        Modifies the instance on which the method was called.

        :param band: List of band indices to crop (index is 0 based)

        Example:

        Crop the first three bands from raster dataset jim0::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim0.cropBand(band=[0, 1, 2])

        """
        if isinstance(band, list):
            bands = band
        else:
            bands = [band]
        band = [self._jim_object.properties.nrOfBand() + b if b < 0 else b
                for b in bands]
        self._jim_object._jipjim.d_cropBand({'band': band})

    def cropOgr(self,
                extent,
                **kwargs):
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

           For instance you can use 'eo':'ATTRIBUTE=fieldname' to burn
           the (numeric) fieldname in the pixel value'
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

    def cropPlane(self,
                  plane: int):
        """Subset raster dataset.

        Subset raster dataset in temporal domain.

        Modifies the instance on which the method was called.

        :param plane: List of plane indices to crop (index is 0 based)

        Example:

        Crop the first three planes from raster dataset jim0::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim0.cropPlane(plane=[0, 1, 2])

        """
        if isinstance(plane, list):
            planes = plane
        else:
            planes = [plane]
        plane = [self._jim_object.properties.nrOfPlane() + b if b < 0 else b
                 for b in planes]
        self._jim_object._jipjim.d_cropPlane({'plane': plane})

    def geo2image(self,
                  x: float,
                  y: float):
        """Convert georeferenced coordinates to image coordinates (col, row).

        :param x: georeferenced coordinate in x according to the object spatial
            reference system
        :param y: georeferenced coordinate in y according to the object spatial
            reference system
        :return: image coordinates (row and column, starting from 0)

        Get column and row index (0 based) of some georeferenced coordinates
        x and y (in this case first pixel: 0, 0)::

          jim = pj.Jim('/path/to/raster.tif')
          x = jim.properties.getUlx()
          y = jim.properties.getUly()
          jim.geometry.geo2image(x, y)
        """
        coord = self._jim_object._jipjim.geo2image(x, y)
        return [int(coord[0]), int(coord[1])]

    def image2geo(self,
                  i: int,
                  j: int):
        """Convert image coordinates (col, row) to georeferenced coordinates.

        :param i: image column number (starting from 0)
        :param j: image row number (starting from 0)
        :return: georeferenced coordinates according to the object spatial
            reference system

        Get upper left corner in georeferenced coordinates
        (in SRS of the Jim object)::

          jim = pj.Jim('/path/to/raster.tif')
          jim.geometry.image2geo(0, 0)

        """
        return self._jim_object._jipjim.image2geo(i, j)

    def imageFrameAdd(self,
                      l: int = 0,
                      r: int = 0,
                      t: int = 0,
                      b: int = 0,
                      u: int = 0,
                      d: int = 0,
                      val: float = 0):
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
            ret_jim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if ret_jim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband._jipjim.d_imageFrameAdd([l, r, t, b, u, d], val)
                    ret_jim.geometry.stackBand(jimband)
                else:
                    ret_jim = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    ret_jim._jipjim.d_imageFrameAdd(
                        [l, r, t, b, u, d], val)
            self._jim_object._set(ret_jim._jipjim)
        else:
            self._jim_object._jipjim.d_imageFrameAdd(
                [l, r, t, b, u, d], val)

    def imageFrameSet(self,
                      l: int = 0,
                      r: int = 0,
                      t: int = 0,
                      b: int = 0,
                      u: int = 0,
                      d: int = 0,
                      val: float = 0,
                      band: int = None):
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

    def imageFrameSubtract(self,
                           l: int = 0,
                           r: int = 0,
                           t: int = 0,
                           b: int = 0,
                           u: int = 0,
                           d: int = 0):
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
            ret_jim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if ret_jim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
                    ret_jim.geometry.stackBand(jimband)
                else:
                    ret_jim = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    ret_jim._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])
            self._jim_object._set(ret_jim._jipjim)
        else:
            self._jim_object._jipjim.d_imageFrameSubtract([l, r, t, b, u, d])

    def imageInsert(self,
                    sec_jim_object,
                    x: float,
                    y: float,
                    z: float,
                    band: list = None):
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

    def imageInsertCompose(self,
                           imlbl,
                           im2,
                           x: float,
                           y: float,
                           z: float,
                           val: float,
                           band: list = None):
        """Merge Jim instance with values of im2 if val of imlbl == val.

        Modifies the instance on which the method was called.

        :param imlbl: a Jim object
        :param im2: a Jim object
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

    def magnify(self,
                n: int):
        """Magnify the image.

        Modifies the instance on which the method was called.

        :param n: a positive integer for magnifying size by pixel replication
        """
        if self._jim_object.properties.nrOfBand() > 1:
            ret_jim = None
            for band in range(0, self._jim_object.properties.nrOfBand()):
                if ret_jim:
                    jimband = _pj.geometry.cropBand(self._jim_object,
                                                    band=band)
                    jimband = _pj.Jim(jimband._jipjim.imageMagnify(n))

                    ret_jim.geometry.stackBand(jimband)
                else:
                    ret_jim = _pj.Jim(_pj.geometry.cropBand(
                        self._jim_object, band=band)._jipjim.imageMagnify(n))
            self._jim_object._set(ret_jim._jipjim)
        else:
            self._jim_object._set(self._jim_object._jipjim.imageMagnify(n))

    def plane2band(self):
        """Convert 3-dimensional 1-band Jim to a 2-dimensional multi-band one.

        The result will be a multi-band single plane object

        Example: convert a single band object with 12 planes to a 2-dimensional
        multi-band object with 1 plane::

           jim = pj.Jim('/path/to/multi/band/image.tif', band2plane=True)
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

    def plotLine(self,
                 x1: int,
                 y1: int,
                 x2: int,
                 y2: int,
                 val: float):
        """Draw a line from [x1, y1] to [x2, y2] by setting pixels to a value.

        Works only for 1-plane Jims.

        Modifies the instance on which the method was called.

        :param x1: an integer for x-coordinate of 1st point
        :param y1: an integer for y-coordinate of 1st point
        :param x2: an integer for x-coordinate of 2nd point
        :param y2: an integer for y-coordinate of 2nd point
        :param val: value to be used for line pixels
        """
        if self._jim_object.properties.nrOfPlane() > 1:
            raise TypeError(
                'Error: plotLine does not support multi-plane Jim objects')
        self._jim_object._jipjim.d_plotLine(x1, y1, x2, y2, val)

    def polygonize(self,
                   output: str,
                   **kwargs):
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

    # todo: to be tested (we can also use jim[jimVect] instead...)
    # def rasterize(self, jim_vect, burn_value=1,eo=['ALL_TOUCHED'],ln=None):
    #     """Rasterize Jim object based on GDALRasterizeLayersBuf

    def rasterize(self,
                  jim_vect,
                  burn_value: int = 1,
                  eo: list = None,
                  ln: str = None):
        """Rasterize Jim object based on GDALRasterizeLayersBuf.

        :param jim_vect: JimVect object that needs to be rasterized
        :param burn_value: burn value
        :param eo: option (default is ALL_TOUCHED)
        :param ln: layer names (optional)

        .. note::
          Possible values for the key 'eo' are:

          ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG.

          For instance you can use 'eo':'ATTRIBUTE=fieldname to burn
          the (numeric) fieldname in the pixel value'

        Example: rasterize vector using the label attribute by creating first
        a mask from an existing raster::

          jim0 = pj.Jim(rasterfn, band=0)
          sample = pj.JimVect(vectorfn)
          mask = pj.Jim(jim0, copy_data=False)
          mask.pixops.convert('GDT_Byte')
          mask.geometry.rasterize(sample, eo=['ATTRIBUTE=label'])
        """
        if not isinstance(jim_vect, _pj.JimVect):
            raise TypeError('Error: can only rasterize a JimVect')

        kwargs = {}
        if burn_value is not None:
            kwargs.update({'burn': float(burn_value)})
        if ln is not None:
            kwargs.update({'ln': ln})

        if eo is not None:
            kwargs.update({'eo': eo})
        else:
            kwargs.update({'eo': ['ALL_TOUCHED']})

        self._jim_object.pixops.setData(0)
        self._jim_object._jipjim.d_rasterizeBuf(jim_vect._jipjimvect, kwargs)

    def reducePlane(self,
                    rule: str = 'overwrite',
                    ref_band: int = None,
                    nodata: float = None):
        """Reduce planes of Jim object.

        :param rule: rule to reduce (mean, median, min or max)
            or callback function
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
            jim_stacked.geometry.reducePlane(getMax)

        """
        nr_of_planes = self._jim_object.properties.nrOfPlane()
        if nr_of_planes < 2:
            _warnings.warn(
                'Single-plane Jim: No plane reduction performed', SyntaxWarning
            )
            return None

        jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
        if isinstance(rule, str):
            nr_of_row = self._jim_object.properties.nrOfRow()
            nr_of_col = self._jim_object.properties.nrOfCol()

            if rule in ('mean', 'avg', 'median'):
                if rule == 'median':
                    nan_func = _np.nanmedian
                    func = _np.median
                else:
                    nan_func = _np.nanmean
                    func = _np.mean

                d_type = self._jim_object.properties.getDataType()
                nr_of_bands = self._jim_object.properties.nrOfBand()
                if nodata is not None and d_type not in ('GDT_Float32',
                                                         'GDT_Float64'):
                    self._jim_object.pixops.convert(otype='GDT_Float32')
                if ref_band is not None:
                    # create plane-wise mask for the ref_band
                    mask = _pj.geometry.cropBand(self._jim_object,
                                                 band=ref_band)

                    # create one-plane boolean mask containing True where
                    # nodata value is in all planes
                    planes = [mask.np()[i] for i in range(nr_of_planes)]
                    stacked_planes = _np.vstack(planes)

                    diffs = _np.abs(
                        _np.diff(stacked_planes.reshape(len(planes), -1),
                                 axis=0)) == 0

                    same_values = diffs[0]
                    for i in range(1, len(diffs)):
                        same_values = _np.logical_and(same_values, diffs[i])

                    same_values = _np.reshape(same_values,
                                              (nr_of_row, nr_of_col))

                    nodata_mask = same_values & (mask.np()[0] == nodata)

                for iband in range(0, nr_of_bands):
                    if nodata is not None:
                        if ref_band is None:
                            # create one-plane boolean mask containing True
                            # where nodata value is in all planes
                            planes = [self._jim_object.np(iband)[i] for i
                                      in range(nr_of_planes)]
                            stacked_planes = _np.vstack(planes)

                            diffs = _np.abs(
                                _np.diff(
                                    stacked_planes.reshape(len(planes), -1),
                                    axis=0)) == 0

                            same_values = diffs[0]
                            for i in range(1, len(diffs)):
                                same_values = _np.logical_and(
                                    same_values, diffs[i])

                            same_values = _np.reshape(same_values,
                                                      (nr_of_row, nr_of_col))

                            nodata_mask = same_values & (
                                self._jim_object.np(iband)[0] == nodata)

                            # compute the reduction function
                            self._jim_object.np(iband)[self._jim_object.np(
                                iband) == nodata] = _np.nan
                            jimreduced.np(iband)[:] = nan_func(
                                self._jim_object.np(iband), axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                        else:
                            self._jim_object.np(iband)[
                                mask.np() == nodata] = _np.nan
                            jimreduced.np(iband)[:] = nan_func(
                                self._jim_object.np(iband), axis=0)
                            jimreduced.np(iband)[nodata_mask] = nodata
                    else:
                        jimreduced.np(iband)[:] = func(
                            self._jim_object.np(iband), axis=0)
                if nodata is not None:
                    if d_type not in ('GDT_Float32', 'GDT_Float64'):
                        jimreduced.pixops.convert(otype=d_type)
                        self._jim_object.pixops.convert(otype=d_type)
            else:
                if rule == 'max':
                    def rule(reduced, plane):
                        return reduced < plane
                elif rule == 'min':
                    def rule(reduced, plane):
                        return reduced > plane
                elif rule == 'overwrite':
                    def rule(reduced, plane):
                        return _pj.pixops.setData(plane, 1)
                else:
                    raise AttributeError('Error: rule not supported')

                if ref_band is not None:
                    maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)

                for iplane in range(1, nr_of_planes):
                    jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)

                    if nodata is not None and ref_band is None:
                        raise AttributeError(
                            'Error: use ref_band option for nodata')

                    if ref_band is not None:
                        maskplane = _pj.geometry.cropBand(jimplane, ref_band)
                        themask = rule(maskreduced, maskplane)
                        if nodata is not None:
                            themask.pixops.convert('GDT_Byte')
                            themask |= maskreduced == nodata
                            themask &= maskplane != nodata
                        maskreduced[themask] = maskplane
                    else:
                        themask = rule(jimreduced, jimplane)

                    jimreduced[themask] = jimplane

                    if nodata is not None:
                        nodata_mask = (maskreduced == nodata) & \
                                      (maskplane == nodata)
                        jimreduced[nodata_mask] = nodata
        else:
            if nodata is not None or ref_band is not None:
                # TODO: Allow nodata and ref_band parameters (see #48)
                raise AttributeError('Error: nodata and ref_band are not '
                                     'supported for callback rules')
            jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
            for iplane in range(1, nr_of_planes):
                jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)
                jimreduced = rule(jimreduced, jimplane)
            # else:
            #     if ref_band is not None:
            #         maskreduced = _pj.geometry.cropBand(jimreduced, ref_band)

            #     for iplane in range(1, nr_of_planes):
            #         jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)

            #         if nodata is not None and ref_band is None:
            #             raise AttributeError(
            #                 'Error: use ref_band option for nodata')

            #         if ref_band is not None:
            #             maskplane = _pj.geometry.cropBand(jimplane, ref_band)
            #             themask = rule(maskreduced, maskplane)
            #             if nodata is not None:
            #                 themask |= maskreduced == nodata
            #                 themask &= maskplane != nodata
            #             maskreduced[themask] = maskplane
            #         else:
            #             themask = rule(jimreduced, jimplane)

            #         jimreduced[themask] = jimplane

            #         if nodata is not None:
            #             nodata_mask = (maskreduced == nodata) & \
            #                         (maskplane == nodata)
            #             jimreduced[nodata_mask] = nodata

        self._jim_object._set(jimreduced._jipjim)

    def _reducePlaneSimple(self,
                           rule: str):
        """Reduce planes of Jim object using callback function without nodata.

        (for performance reasons)

        :param rule: rule to reduce (mean, median) or callback function
            (e.g., for max composite,
            use callBack(jimreduced, jimplane): return jimreduced < jimplane)

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
            def rule(reduced, plane):
                return _pj.pixops.supremum(reduced, plane)
        elif rule == 'min':
            def rule(reduced, plane):
                return _pj.pixops.infimum(reduced, plane)
        jimreduced = _pj.geometry.cropPlane(self._jim_object, 0)
        for iplane in range(1, self._jim_object.properties.nrOfPlane()):
            jimplane = _pj.geometry.cropPlane(self._jim_object, iplane)
            jimreduced = rule(jimreduced, jimplane)
        self._jim_object._set(jimreduced._jipjim)

    def stackBand(self,
                  jim_other,
                  band: int = None):
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
        if band is not None and band > self._jim_object.properties.nrOfBand():
            raise AttributeError('Band out of bounds')

        if not isinstance(jim_other, list):
            jim_other = [jim_other]

        for jim in jim_other:
            if band:
                self._jim_object._jipjim.d_stackBand(jim._jipjim,
                                                     {'band': band})
            else:
                self._jim_object._jipjim.d_stackBand(jim._jipjim)

    def stackPlane(self,
                   jim_other=None,
                   *args):
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
        if jim_other is None:
            return None

        if not isinstance(jim_other, list):
            args_list = [jim_other, *args]

        if args:
            args_list.extend(args)

        for jim in args_list:
            self._jim_object._jipjim.d_stackPlane(jim._jipjim)

    def warp(self,
             t_srs,
             **kwargs):
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

            jim = pj.Jim('/path/to/file.tif')
            jim.warp('epsg:3035')

        Read a raster dataset from disk that is in lat lon (epsg:4326), select
        a bounding box in a different spatial reference system (epsg:3035).
        Notice the raster dataset read is still in the original projection
        (epsg:4326). Then warp the raster dataset to the target spatial
        reference system (epsg:3035)::

            jim = pj.Jim('/path/to/file.tif', t_srs='epsg:3035',
                         ulx=1000000, uly=4000000, lrx=1500000, lry=3500000)
            jim.warp('epsg:3035', s_srs='epsg:4326')

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


class _GeometryList(_pj.modules.JimListModuleBase):
    """Define all Geometry methods for JimLists."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self,
                    caller):
        self._jim_list = caller

    def stackBand(self,
                  jim_other=None,
                  band: int = None):
        """Stack bands from raster datasets into new multiband Jim object.

        :param jim_other: a Jim object or jimlist from which to copy bands
            (optional)
        :param band: List of band indices to stack (index is 0 based)
        :return: multiband Jim object

        Create a multiband Jim object from a list of two Jim objects::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim_stacked = pj.JimList([jim0, jim1]).geometry.stackBand()

        Create a multiband Jim object from a JimList object, selecting
        the first and third band::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim2 = pj.Jim('/path/to/raster2.tif')
            jimlist = pj.JimList([jim0, jim1, jim2])
            jim_stacked = jimlist.geometry.stackBand([0, 2])
        """
        if band:
            nr_of_band = self._jim_list[0].properties.nrOfBand()
            if band is not None and band > nr_of_band:
                raise AttributeError('Band out of bounds')

            ret_jim = _pj.Jim(
                self._jim_list._jipjimlist.stackBand({'band': band}))
        else:
            ret_jim = _pj.Jim(self._jim_list._jipjimlist.stackBand())

        if isinstance(jim_other, _pj.Jim):
            if band:
                ret_jim.geometry.stackBand(jim_other, band=band)
            else:
                ret_jim.geometry.stackBand(jim_other)
        elif isinstance(jim_other, list):
            if band:
                jim_to_stack = _pj.Jim(
                    jim_other._jipjimlist.stackBand({'band': band}))
            else:
                jim_to_stack = _pj.Jim(jim_other._jipjimlist.stackBand())

            ret_jim = _pj.Jim(ret_jim._jipjim.stackBand(
                jim_to_stack._jipjim))

        return ret_jim

    def stackPlane(self, jim_other=None, *args):
        """Stack planes from raster datasets into new multiplane Jim object.

        :param jim_other: Either Jim or JimList to be stacked on top of
            the JimList on which the method was called
        :return: multiplane Jim object

        Stack planes from another raster dataset to current raster dataset.

        Create a multiplane Jim object from a list of two Jim objects::

            jim0 = pj.Jim('/path/to/raster0.tif')
            jim1 = pj.Jim('/path/to/raster1.tif')
            jim_stacked = pj.JimList([jim0, jim1]).geometry.stackPlane()
        """
        if not isinstance(jim_other, list):
            args_list = [jim_other]
        else:
            args_list = jim_other

        if args:
            args_list.extend(args)

        ret_jim = _pj.Jim(self._jim_list._jipjimlist.stackPlane())

        if jim_other is not None:
            for jim in args_list:
                ret_jim._jipjim.d_stackPlane(jim._jipjim)
        return ret_jim


class _GeometryVect(_pj.modules.JimVectModuleBase):
    """Define all Geometry methods for JimVects."""

    def __init__(self):
        """Initialize the module."""
        pass

    def _set_caller(self,
                    caller):
        self._jim_vect = caller

    def append(self,
              jvec,
              **kwargs):
        """Append JimVect object with another JimVect object.

        :param jvec: JimVect object to append
        :param kwargs: See table below
        :return: joined JimVect object

        Modifies the instance on which the method was called.

        Example: append two vectors::

          v1 = pj.JimVect('/path/to/vector1.sqlite')
          appendedv1.geometry.append(v2, '/path/to/appended.sqlite')
        """
        non_existing_path = _pj._get_random_path()
        non_existing_path = os.path.join('/vsimem',
                                         os.path.basename(non_existing_path))
        kwargs.update({'output': non_existing_path})
        if isinstance(jvec, _pj.JimVect):
            avect = self._jim_vect._jipjimvect.merge(jvec._jipjimvect, kwargs)
            avect.write()
            self._jim_vect._set(avect)
        else:
            raise TypeError('Error: can only append two JimVect objects')

    def convexHull(self,
                   **kwargs):
        """Create the convex hull on a JimVect object.

        Modifies the instance on which the method was called.

        :param kwargs: See table below

        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+
        """
        non_existing_path = _pj._get_random_path()
        non_existing_path = os.path.join('/vsimem',
                                         os.path.basename(non_existing_path))
        kwargs.update({'output': non_existing_path})
        avect = self._jim_vect._jipjimvect.convexHull(kwargs)
        avect.write()
        self._jim_vect._set(avect)

    def extract(self,
                jim,
                output: str,
                rule: str = None,
                **kwargs):
        """Extract pixel values from raster image based on a vector dataset.

        :param jim: Jim object (list) on which vector is overlaid to extract from
        :param rule: Rule how to calculate zonal statistics per feature
            (see list of :ref:`supported rules <extract_rules>`)
        :param output: Name of the output vector dataset in which the zonal
            statistics will be saved
        :param kwargs: See table below
        :return: A JimVect with the same geometry as the sample vector
            dataset and an extra field for each of the calculated raster value
            (zonal) statistics. The same layer name(s) of the sample will be
            used for the output vector dataset


        +------------------+--------------------------------------------------+
        | key              | value                                            |
        +==================+==================================================+
        | copy             | Copy these fields from the sample vector dataset |
        |                  | (default is to copy all fields)                  |
        |                  | set copy to single field if rule is 'allpoints'  |
        +------------------+--------------------------------------------------+
        | label            | Create extra field named 'label' with this value |
        +------------------+--------------------------------------------------+
        | fid              | Create extra field named 'fid' with this field   |
        |                  | identifier (sequence of features)                |
        +------------------+--------------------------------------------------+
        | bandname         | List of band names corresponding to list of      |
        |                  | bands to extract                                 |
        +------------------+--------------------------------------------------+
        | planename        | List of plane names corresponding to list of     |
        |                  | planes to extract                                |
        +------------------+--------------------------------------------------+
        | oformat          | Output vector dataset format                     |
        +------------------+--------------------------------------------------+
        | co               | Creation option for output vector dataset        |
        +------------------+--------------------------------------------------+

        .. _extract_rules:

        :Supported rules to extract <extract_rules>:

        +------------------+--------------------------------------------------+
        | rule             | description                                      |
        +==================+==================================================+
        | point            | extract a single pixel within the polygon or on  |
        |                  | each point feature                               |
        +------------------+--------------------------------------------------+
        | allpoints        | Extract all pixel values covered by the polygon  |
        |                  | (copy field must be set to one field)            |
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
        |                  | (classes must be set with the option classes)    |
        +------------------+--------------------------------------------------+
        | proportion       | Extract proportion of class(es) within           |
        |                  | the polygon                                      |
        |                  | (classes must be set with the option classes)    |
        +------------------+--------------------------------------------------+
        | count            | Extract count of class(es) within the polygon    |
        |                  | (classes must be set with the option classes)    |
        +------------------+--------------------------------------------------+
        | percentile       | Extract percentile as defined by option perc     |
        |                  | (e.g, 95th percentile of values covered by       |
        |                  | polygon)                                         |
        +------------------+--------------------------------------------------+


        .. note::
            To ignore some pixels from the extraction process, see list
            of :ref:`mask key values <extract_mask>`

        .. _extract_mask:

        :Supported key values to mask pixels that must be ignored in \
        the extraction process <extract_mask>:

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
        | mask             | Use the specified file as a validity mask        |
        +------------------+--------------------------------------------------+
        | mskband          | Use the specified band of the mask file          |
        |                  | defined                                          |
        +------------------+--------------------------------------------------+
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
            v = reference.geometry.extract(
                jim0, buffer=-10, rule=['mean'],
                output='/vsimem/temp.sqlite', oformat='SQLite')
            v.io.write('/path/to/output.sqlite)

        """
        kwargs.update({'output': output})
        kwargs.update({'rule': rule})
        if 'threshold' in kwargs:
            if kwargs['threshold'] is not None:
                if not isinstance(kwargs['threshold'],list):
                    kwargs['threshold']=[kwargs['threshold']]
                kwargs['threshold']=[float(t.strip('%')) if isinstance(t,str) else -t for t in kwargs['threshold']]

        if 'classes' in kwargs:
            classes = kwargs.pop('classes')
            kwargs['class'] = classes

        if isinstance(jim, _pj.Jim):
            avect = jim._jipjim.extractOgr(self._jim_vect._jipjimvect, kwargs)
        elif isinstance(jim, _pj.JimList):
            bandname = kwargs.pop('bandname', None)
            #todo: support multi-band images in JimList...
            if bandname is None:
                bandname = [
                    't'+str(ifile) for ifile in range(0, len(self._jim_list))]
            kwargs.update({'bandname': bandname})

            avect = jim._jipjimlist.extractOgr(self._jim_vect._jipjimvect, kwargs)
        else:
            raise TypeError('Error: extract must operate on Jim or JimList')
        avect.write()
        self._jim_vect._set(avect)

    def intersect(self,
                  jim,
                  **kwargs):
        """Intersect JimVect object with Jim object.

        Keeps only those features with an intersect.

        Modifies the instance on which the method was called.

        :param jim: Jim object with which to intersect

        Example: intersect a sample with a Jim object::

          jim = pj.Jim('/path/to/raster.tif')
          v = pj.JimVect('/path/to/vector.sqlite')
          sampleintersect = pj.geometry.intersect(
              v, jim, output='/vsimem/intersect', oformat='SQLite',
              co=['OVERWRITE=YES'])
          sampleintersect.io.write('/path/to/output.sqlite')

        """
        non_existing_path = _pj._get_random_path()
        non_existing_path = os.path.join('/vsimem',
                                         os.path.basename(non_existing_path))
        kwargs.update({'output': non_existing_path})
        if isinstance(jim, _pj.Jim):
            avect = self._jim_vect._jipjimvect.intersect(jim._jipjim, kwargs)
            avect.write()
            self._jim_vect._set(avect)
        else:
            raise TypeError('Error: can only intersect with Jim object')

    def join(self,
             jvec2,
             **kwargs):
        """Join JimVect object with another JimVect object.

        A key field is used to find corresponding features in both objects.

        Modifies the instance on which the method was called.

        :param jvec2: second JimVect object to join
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

        The join methods currently supported are:

            :INNER |inner|: join two JimVect objects, keeping only those \
                            features for which identical keys in both objects \
                            are found
            :OUTER_LEFT |outer_left|: join two JimVect objects, keeping \
                                      all features from first object
            :OUTER_RIGHT |outer_right|: join two JimVect objects, keeping all \
                                        features from second object
            :OUTER_FULL |outer_full|: join two JimVect objects, keeping \
                                      all features from both objects

        Example: join two vectors, based on the key 'id', which is a common
        field shared between v1 and v2. Use OUTER_FULL as the join method::

          v1 = pj.JimVect('/path/to/vector1.sqlite')
          v2 = pj.JimVect('/path/to/vector2.sqlite')
          v3 = pj.geometry.join(
              v1, v2, '/tmp/test.sqlite', oformat='SQLite',
              co=['OVERWRITE=YES'], key=['id'], method='OUTER_FULL')
        """
        non_existing_path = _pj._get_random_path()
        non_existing_path = os.path.join('/vsimem',
                                         os.path.basename(non_existing_path))
        kwargs.update({'output': non_existing_path})
        if isinstance(jvec2, _pj.JimVect):
            avect = self._jim_vect._jipjimvect.join(jvec2._jipjimvect, kwargs)
            avect.write()
            self._jim_vect._set(avect)
        else:
            raise TypeError('Error: can only join two JimVect objects')

