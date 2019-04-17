.. _Tutorial:

########
Tutorial
########

.. toctree::
   :maxdepth: 4

*******************************************************
Sources of documentation and help
*******************************************************

This documentation is written using the `sphinx <http://www.sphinx-doc.org>`_ tool.

====================
Online documentation
====================

The primary source of documentation is available online via this `website <https://cidportal.jrc.ec.europa.eu/services/processing/pyjeohelp>`_.
The documentation is divided in three main parts:

* :ref:`Introduction`: Introduction explaining the organization of the pyjeo package, its modules, design, installation and usage.

* :ref:`Tutorial`: A tutorial to guide you through the first steps of using pyjeo, introducing the main data structures :py:class:`Jim`, :py:class:`JimList`, and :py:class:`JimVect`.

* :ref:`Reference`: A manual describing all functions and methods available in pyjeo


====================
Inline documentation
====================

In addition to the online help pages, there is the inline documetation. 

To get help on a specific module, e.g., :py:mod:`geometry`::

  help(pj.geometry)

To get help on a class, e.g., :py:class:`Jim`::

  help(pj.Jim)

To get help on a function, e.g., :py:func:`geometry.warp`::

  help(pj.geometry.warp)

To get help on a class method, e.g., :py:meth:`~geometry._Geometry.warp`::

  jim=pj.Jim()
  help(jim.geometry.warp)

*******************************************************
How to create basic data structures (Jim) and (JimVect)
*******************************************************

=================================
Creating a Jim raster data object
=================================

To create a new Jim object, use the constructor :py:meth:`Jim`. A Jim object can be created by opening an existing raster data :ref:`from file <create_Jim_from_file>`. Alternatively, a :ref:`new Jim <create_Jim_new>` object can be created by defining image attributes.

Creating a Jim raster data object from file
===========================================

   Examples:

   Create Jim object by opening an existing file (file content will
   automatically be read in memory)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=pj.Jim(ifn)
        #do stuff with jim ...

   Create Jim object by opening an existing file for specific region of
   interest and spatial resolution using cubic convolution resampling::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim0=pj.Jim(ifn,noread=True)
        ULX=jim0.properties.getUlx()
        ULY=jim0.properties.getUly()
        LRX=jim0.properties.getUlx()+100*jim0.properties.getDeltaX()
        LRY=jim0.properties.getUly()-100*jim0.properties.getDeltaY()
        jim=pj.Jim(ifn,ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,dx=5,dy=5,resample='GRIORA_Cubic')
        #do stuff with jim ...

   Create Jim object by opening an existing file, reading 10 columns and row, starting from the sixth (index=5) row and column. Here we use  pixel coordinates instead of georeferenced coordinates to define the bounding box (nogeo=True)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=pj.Jim(ifn,ulx=5,uly=5,lrx=14,lry=14,nogeo=True)
        #do stuff with jim ...


Creating a new Jim raster data object (not from file)
=====================================================

   Create a new georeferenced Jim image object by defining the projection epsg
   code, bounding box, and pixel size::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=pj.Jim(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='uint16',dx=100,dy=100)
        #do stuff with jim ...

   Create a new georeferenced Jim image object for writing by defining
   the projection epsg code, bounding box and number of rows and columns::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=pj.Jim(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='uint16',ncol=1098,nrow=1098)
        #do stuff with jim ...

        
   To initialize the pixel values to some value, either use the methode :py:meth:`Jim:setData`::

        jim.pixops.setData(100)

   Alternatively, you can initialize the pixel values when creating the Jim object, selecting values for mean, stdev or uniform. For instance, create a new Jim object with data type float, 256 rows, 256 columns and set all values to 10::

        im=pj.Jim(otype='float32',ncol=256,nrow=256,mean=10)

   Random values using a Gaussian distribution function can be generated by setting a value for stdev (e.g., 2)::

        im=pj.Jim(otype='float32',ncol=256,nrow=256,mean=10, stdev=2)

   To generate random values with a uniform distribution with values between 100 and 200::

        im=pj.Jim(otype='float32',ncol=256,nrow=256,uniform=[100,200])

Creating a new copy of a Jim raster data object
===============================================

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
