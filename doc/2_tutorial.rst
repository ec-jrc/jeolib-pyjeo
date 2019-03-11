########
Tutorial
########

.. toctree::
   :maxdepth: 4

*******************************************************
How to create basic data structures (Jim) and (JimVect)
*******************************************************

=================================
Creating a Jim raster data object
=================================

To create a new Jim raster data object, use the constructor :py:meth:`Jim`. A Jim object can be created by opening an existing raster data :ref:`from file <create_Jim_from_file>`. Alternatively, a :ref:`new Jim <create_Jim_new>` object can be created by defining image attributes.

   Examples:

   Create Jim image object by opening an existing file (file content will
   automatically be read in memory)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=pj.Jim(ifn)
        #do stuff with jim ...

   Create Jim image object by opening an existing file for specific region of
   interest and spatial resolution using cubic convolution resampling::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim0=pj.Jim(ifn,'noread'=True)
        ULX=jim0.getUlx()
        ULY=jim0.getUly()
        LRX=jim0.getUlx()+100*jim0.getDeltaX()
        LRY=jim0.getUly()-100*jim0.getDeltaY()
        jim=pj.Jim(ifn,ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,dx=5,dy=5,resample='GRIORA_Cubic')
        #do stuff with jim ...


   Create a new georeferenced Jim image object by defining the projection epsg
   code, bounding box, and pixel size::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=pj.Jim(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',dx=100,dy=100)
        #do stuff with jim ...

   Create a new georeferenced Jim image object for writing by defining
   the projection epsg code, bounding box and number of rows and columns::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=pj.Jim(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',ncol=1098,nrow=1098)
        #do stuff with jim ...

=================================
Creating a JimVect data object
=================================

To create a new JimVect vector data object, use the constructor :py:meth:`JimVect`. A JimVect object is typically created by opening an existing vector dataset from file or as the result from a function or method, e.g., the method :py:meth:`~geometry._Geometry.extractOgr` in module :py:mod:`geometry` on a Jim object.

   Examples:

   Open a vector and read all layers::

     v0=pj.JimVect('/path/to/vector.sqlite')

   Open a vector and read layer named lodi::

     v0=pj.JimVect('/path/to/nuts.sqlite', ln='lodi')

   Extract the mean value of the pixels within the polygon of the provided reference vector. Exclude the pixels within a buffer of 10m of the polygon boundary. Use a temporary vector in memory for the calculation. Write the result to the final destination on disk::

     reference = pj.JimVect('/path/to/reference.sqlite')
     jim0 = pj.Jim('/path/to/raster.tif')
     v = jim0.extractOgr(reference, buffer=-10, rule=['mean'], output='/vsimem/temp.sqlite', oformat='SQLite')
     v.write('/path/to/output.sqlite)
