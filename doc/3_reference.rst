.. _Reference:

################
Reference Manual
################

.. toctree::
   :maxdepth: 4


****************************************************
Data structures: rasters (Jim) and vectors (JimVect)
****************************************************

.. method:: Jim(filename=None, **kwargs)

  Creates a new Jim object, either :ref:`from file <create_Jim_from_file>`
  or :ref:`create new <create_Jim_new>`

:param filename: Path to a raster dataset or another Jim object
:param: see supported keys in table below
:return: a Jim object


.. _create_Jim_from_file:

===========================
Create Jim object from file
===========================

   Supported keys as arguments:

    ======== ===================================================
    band     Bands to open, index starts from 0
    ulx      Upper left x value bounding box
    uly      Upper left y value bounding box
    lrx      Lower right x value bounding box
    lry      Lower right y value bounding box
    dx       Resolution in x
    dy       Resolution in y
    resample Resample algorithm used for reading pixel data
    extent   filename of a vector dataset from which to get boundary constraint
    nodata   Nodata value to put in image
    nogeo    Use image coordinates (index starts from 0,0 for upper left pixel)
    noread   Set this flag to True to not read data when opening
    ======== ===================================================

   .. note::
        You can specify a different spatial reference system to define the
        region of interest to read set with keys ulx, uly, lrx, and lry with
        the extra key 't_srs'. Notice this will not re-project the resulting
        image. You can use the function :py:func:`geometry.warp` or the corresponding
        method on a Jim object :py:meth:`~geometry._Geometry.warp` in the :py:mod:`geometry` module.
   ..

   .. note::
      for resample values, please check `gdal site <http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a>`_.

   Example:

   Create Jim image object by opening an existing file (file content will
   automatically be read in memory)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=pj.Jim(ifn)
        #do stuff with jim ...

   Create Jim image object by opening an existing file for specific region of
   interest and spatial resolution using cubic convolution resampling::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim0=pj.Jim(ifn,'noread'=True)
        ULX=jim0.properties.getUlx()
        ULY=jim0.properties.getUly()
        LRX=jim0.properties.getUlx()+100*jim0.properties.getDeltaX()
        LRY=jim0.properties.getUly()-100*jim0.properties.getDeltaY()
        jim=pj.Jim(ifn,ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,dx=5,dy=5,resample='GRIORA_Cubic')
        #do stuff with jim ...

   Create Jim image object by opening an existing file, reading 10 columns and row, starting from the sixth (index=5) row and column::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        jim=pj.Jim(ifn,ulx=5,uly=5,lrx=14,lry=14)
        #do stuff with jim ...

.. _create_Jim_new:

===============================================================================
Create a new Jim image object by defining image attributes (not read from file)
===============================================================================

Supported keys as arguments:

===== =================
ncol  Number of columns
nrow  Number of rows
nband (default: 1) Number of bands
otype (default: Byte) Data type
a_srs Assign the spatial reference for the output file, e.g., psg:3035 to
      use European projection and force to European grid
===== =================

Supported keys used to initialize random pixel values in new Jim image
object:

======= ============================================
seed    (default: 0) seed value for random generator
mean    (default: 0) Mean value for random generator
stdev   (default: 0) Standard deviation for Gaussian random generator
uniform (default: 0) Start and end values for random value with uniform
        distribution
======= ============================================

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


=====================================
Converting Jim object to numpy arrays
=====================================

.. method:: Jim.np()

   Return numpy array from Jim object
   :return: numpy array representation

=================================
Creating a JimVect data object
=================================

.. method:: JimVect(filename, **kwargs)

  Create a new JimVect object from file

:param filename:        Path to a vector dataset
:param ln:              Layer name to read (default is to read all layers)
:param attributeFilter: Set an attribute filter in restricted SQL WHERE format
:ulx:                   Upper left x value bounding box
:uly:                   Upper left y value bounding box
:lrx:                   Lower right x value bounding box
:lry:                   Lower right y value bounding box
:param noread:          Set this flag to True to not read data when opening
:return: a JimVect object

Example:

Open a vector and read all layers::

  v=pj.JimVect('/path/to/vector.sqlite')

Open a vector and read layer named lodi::

  v=pj.JimVect('/path/to/nuts.sqlite', ln='lodi')

Open a vector and use an attribute filter (the field intern_id must be between 10000 and 10500)::

  v=pj.JimVect('/path/to/vector.sqlite', attributeFilter='(intern_id>10000) AND (intern_id<10500)')


.. _indexing:

********************
Indexing Jim objects
********************

=============
Get Jim items
=============

   .. method:: Jim[item]

      Get subset of the raster dataset. Item can be of type:

      :tuple: get all pixels defined by tuple (e.g., [0:10,0:10] for first 10 rows and columns in a single band image)

      :Jim: get all pixels where the specified raster dataset object is > 0

      :JimVect: get spatial subset of all pixels covered by the specified vector dataset object

      :returns: subset of the raster dataset

      Example:

      Get first 10 columns in first 10 rows::

        jim0[0:10,0:10]

      Get a binary mask of all values not within [0,250] (notice the parentheses to take care of the precedence of the operators!) ::

        jim0[(jim0<0) | (jim0>250)]

      Get a binary mask identifying  pixels where jim0<jim1::

        jim0[(jim0<jim1)]

      Crop a raster dataset according to the extent of a vector dataset (crop_to_cutline), set all pixels not covered to 0 (or value defined as no data)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        cfn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/QI_DATA/MSK_CLOUDS_B00.gml'
        jim=pj.Jim(ifn)
        v=pj.JimVect(cfn)
        jimcloud=jim[v]

=============
Set Jim items
=============

   .. method:: Jim[item]=

        Set items of the raster dataset. Item can be of type:

        :tuple: set all pixels defined by tuple (e.g., [0:10,0:10] for first 10 rows and columns in a single band image)
        :Jim: set all pixels where the specified raster dataset object is > 0
        :JimVect: set all pixels covered by the specified vector dataset object

        Modifies the instance on which the method was called.

        Example:

        Set first 10 columns in first 10 rows to 255::

          jim0[0:10,0:10]=255

        Set all edge pixels to 0::

          jim[0,:]=0
          jim[:,0]=0
          jim[-1,:]=0
          jim[:,-1]=0

        Mask all values not within [0,250] and set to 255 (no data)::

          jim0[(jim0<0) | (jim0>250)]=255

        Select the maximum of two Jim images::

          jim0[(jim0<jim1)]=jim1

        Set a gml cloud mask to a Jim image (setting all cloudy pixels to 255)::

          ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
          cfn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/QI_DATA/MSK_CLOUDS_B00.gml'
          jim=pj.Jim(ifn)
          v=pj.JimVect(cfn)
          jim[v]=255
   
*********
Operators
*********

.. _comparison_operators:

====================
Comparison operators
====================

      Jim objects support comparison operations performed at pixel level. The result is a new binary Jim object (of type Byte) with value 1 if the comparison result is True and value 0 if the comparison result is False.

========   ============
Operator   Meaning
========   ============
``==``     Pixel wise check for equality
``!=``     Pixel wise check for inequality
``<``      Pixel wise check for less than
``>``      Pixel wise check for greater than
``<=``     Pixel wise check for less than or equal to
``>=``     Pixel wise check for greater than or equal to
========   ============

      Examples:

      Pixel wise check for equality. The result is a binary Jim object of type Byte: 1 if pixel values are equal and 0 if objects differ::

        result = jim1==jim2

      Set all pixel to 0 where Jim objects are equal::

        jim1[jim1==jim2] = 0

      Pixel wise check for inequality. The result is a binary Jim object of type Byte: 0 if pixel values are equal and 1 if objects differ::

        result = jim1!=jim2

      Set all pixel to 0 where Jim objects differ::

        jim1[jim1!=jim2] = 0

      Pixel wise check if an image is less than another::

        result = jim1<jim2

      Set all pixel values less than 0 to 0::

        jim0[jim0<0] = 0

      Pixel wise check if an image is less than or equal to another::

        result = jim1<=jim2

      Set all pixel values less than or equal to 0 to 0::

        jim0[jim0<=0] = 0

      Pixel wise check if an image is greater than another::

        result = jim1>jim2

      Set all pixel values greater than 0 to 0::

        jim0[jim0>0] = 0

      Pixel wise check if an image is greater than or equal to another::

        result = jim1>=jim2

      Set all pixel values greater than or equal to 0 to 0::

        jim0[jim0>=0] = 0


.. _boolean_operators:

=================
Boolean operators
=================

Jim objects support boolean operations performed at pixel level. Both input and result are assumed to be binary Jim objects of type Byte (0 is False, 1 is True).

========   ============
Operator   Meaning
========   ============
``|``      bitwise or
``^``      bitwise exclusive or
``&``      bitwise and
========   ============


      Examples:

      Calculate the bitwise or value of to Jim object::

        result=jim1 | jim2

      Get a binary mask of all values not within [0,250] (notice the parentheses to take care of the precedence of the operators!) ::

        result = jim0[(jim0<0) | (jim0>250)]

      Calculate the bitwise exclusive or (xor) value of to Jim object::

        result=jim1 ^ jim2

      Calculate the bitwise and value of to Jim object::

        result=jim1 & jim2

      Get a binary mask of all values within [100,200] (notice the parentheses to take care of the precedence of the operatands!) ::

        result = jim0[(jim0>=100) & (jim0<=200)]


.. _unary_operators:

==========================
Arithmetic unary operators
==========================

      Jim objects support unary arithmetic operations performed at pixel level.

========   ============
Operator   Meaning
========   ============
``abs``    absolute value
``-``      negation
========   ============

      Examples:

      Calculate the absolute value of Jim object::

        jimabs=abs(jim)

      Calculate the negation of a jim object. Notice that the output data type can be changed if the input was not signed. A warning message is given::

        jimneg=-jim

.. _binary_operators:

===========================
Arithmetic binary operators
===========================

      Jim objects support binary arithmetic operations performed at pixel level.

========   ============
Operator   Meaning
========   ============
``+``      addition
``+=``     addition and assignment
``-``      subtraction
``-=``     subtraction and assignment
``*``      multiplication
``*=``     multiplication and assignment
``/``      division
``/=``     division and assignment
``%``      modulo
``%=``     modulo and assignment
``<<``     bitwise left shift
``<<=``    bitwise left shift and assignment
``>>``     bitwise right shift
``>>=``    bitwise right shift and assignment
========   ============

      Examples:

      Add two Jim objects and return the result as a new Jim object::

        jim=jim1+jim2

      Replace the content of jim1 with the sum of jim1 and jim2::

        jim1+=jim2

      Subtract two Jim objects and return the result as a new Jim object::

        jim=jim1-jim2

      Replace the content of jim1 with the difference of jim1 and jim2::

        jim1-=jim2

      Multiply two Jim objects and return the result as a new Jim object::

        jim=jim1*jim2

      Replace the content of jim1 with the multiplication of jim1 and jim2::

        jim1*=jim2

      Divide two Jim objects and return the result as a new Jim object::

        jim=jim1/jim2

      Replace the content of jim1 with the division of jim1 and jim2::

        jim1/=jim2

      Calculate the modulo of 2 (remainder after division of 2) and return the result as a new Jim object::

        jim=jim0%2

      Replace the content of jim1 with the modulo 2 (remainder after division of 2)::

        jim%=2

      Multiply the pixel values by 2 and return as a new Jim object (by calculating the bitwise left shift)::

        jim2 = jim1 << 2

      Multiply the pixel values by 2 and replace the current Jim object (by calculating the bitwise left shift)::

        jim<<=2

      Divide the pixel values by 2 and return as a new Jim object (by calculating the bitwise right shift)::

        jim2 = jim1 >> 2

      Divide the pixel values by 2 and replace the current Jim object (by calculating the bitwise right shift)::

        jim>>=2

********************
Accessing properties
********************

.. automodule:: properties
   :members:

==============
Jim properties
==============

.. autoclass:: _Properties
   :members:

==================
JimVect properties
==================

.. autoclass:: _PropertiesVect
   :members:

==================
JimList properties
==================

.. autoclass:: _PropertiesList
   :members:

********************
Input/Output methods
********************

.. automodule:: pjio
   :members:

========================
Jim Input/Output methods
========================

.. autoclass:: _IO
   :members:

============================
JimVect Input/Output methods
============================

.. autoclass:: _IOVect
   :members:

============================
JimList Input/Output methods
============================

.. autoclass:: _IOList
   :members:

****************
Pixel operations
****************

=========================
Pixel operation functions
=========================

.. automodule:: pixops
   :members:
   :exclude-members: simpleArithOp, simpleBitwiseOp, simpleThreshold, setLevel

==============================
Pixel operation methods on Jim
==============================

.. autoclass:: _PixOps
   :members:
   :exclude-members: simpleArithOp, simpleBitwiseOp, simpleThreshold, setLevel

==================================
Pixel operation methods on JimList
==================================

.. autoclass:: _PixOpsList
   :members:

***********************
Neighborhood operations
***********************

===========================================
Neighborhood operation from scipy (ndimage)
===========================================

The neighborhood operations from scipy ndimage can be applied to a :py:class:`Jim` object by using its numpy representation (:py:meth:`Jim.np`)

Perform a Gaussian filter using a standard deviation (sigma) of 2::

  jim=pj.Jim('/path/to/image.tif')
  npfiltered=ndimage.gaussian_filter(jim.np(),2)
  jim_filtered=pj.np2jim(npfiltered)

The jim_filtered object does not contain projection and geotransform information. These can be easily obtained from the original jim object::

  jim_filtered.properties.setProjection(jim.properties.getProjection())
  jim_filtered.properties.setGeoTransform(jim.properties.getGeoTransform())


=======================================
Native neighborhood operation functions
=======================================

.. automodule:: ngbops
   :members:

============================================
Native neighborhood operation methods on Jim
============================================

.. autoclass:: _NgbOps
   :members:
     

*******************
Geometry operations
*******************

============================
Geometry operation functions
============================

.. automodule:: geometry
   :members:
   :exclude-members: imageInsert, imageInsertCompose, imageFrameSet, imageFrameAdd, magnify 

=================================
Geometry operation methods on Jim
=================================

.. autoclass:: _Geometry
   :members:
   :exclude-members: imageInsert, imageInsertCompose, imageFrameSet, imageFrameAdd, magnify 

=====================================
Geometry operation methods on JimVect
=====================================

.. autoclass:: _GeometryVect
   :members:


******************************
Connected component operations
******************************

=======================================
Connected component operation functions
=======================================

.. automodule:: ccops
   :members:

============================================
Connected component operation methods on Jim
============================================

.. autoclass:: _CCOps
   :members:

**************
Classification
**************

========================
Classification functions
========================

.. automodule:: classify
   :members:

=============================
Classification methods on Jim
=============================

.. autoclass:: _Classify
   :members:

=================================
Classification methods on JimVect
=================================

.. autoclass:: _ClassifyVect
   :members:

*****************
Digital elevation
*****************

===========================
Digital elevation functions
===========================

.. automodule:: demops
   :members:

================================
Digital elevation methods on Jim
================================

.. autoclass:: _DEMOps
   :members:
