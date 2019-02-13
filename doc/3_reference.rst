################
Reference Manual
################

.. toctree::
   :maxdepth: 4

****************************************************
Data structures: rasters (Jim) and vectors (JimVect)
****************************************************

=================================
Creating a Jim raster data object
=================================

.. method:: Jim(filename=None, **kwargs)

  Create a new Jim object, either :ref:`from file <create_Jim_from_file>`
  or :ref:`create new <create_Jim_new>`

:param filename: Path to a raster dataset or another Jim object
:param: see supported keys in table below
:return: a Jim object


.. _create_Jim_from_file:

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
    resample Resample algorithm used for reading pixel data in case of interpolation
    extent   filename of a vector dataset from which to get boundary constraint
    nodata   Nodata value to put in image
    noread   Set this flag to True to not read data when opening
    ======== ===================================================

   .. note::
        You can specify a different spatial reference system to define the
        region of interest to read set with keys ulx, uly, lrx, and lry with
        the extra key 't_srs'. Notice this will not re-project the resulting
        image. You can use the function :py:func:`geometry.warp` or the corresponding
        method on a Jim object :py:meth:`geometry._Geometry.warp` for this.
   ..

   .. note::
        resample values: please check
        http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a

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
        ULX=jim0.getUlx()
        ULY=jim0.getUly()
        LRX=jim0.getUlx()+100*jim0.getDeltaX()
        LRY=jim0.getUly()-100*jim0.getDeltaY()
        jim=pj.Jim(ifn,ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,dx=5,dy=5,resample='GRIORA_Cubic')
        #do stuff with jim ...

.. _create_Jim_new:

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
        jim=pj.Jim.createImg(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',dx=100,dy=100)
        #do stuff with jim ...
        jim.close()

Create a new georeferenced Jim image object for writing by defining
the projection epsg code, bounding box and number of rows and columns::

        projection='epsg:32612'
        ULX=600000.0
        ULY=4000020.0
        LRX=709800.0
        LRY=3890220.0
        jim=pj.Jim.createImg(ulx=ULX,uly=ULY,lrx=LRX,lry=LRY,a_srs=projection,otype='GDT_UInt16',ncol=1098,nrow=1098)
        #do stuff with jim ...
        jim.close()

*********
Operators
*********

======================
Pixel access operators
======================

   .. method:: Jim[item]

      Get subset of the raster dataset.

      :param items can be of type:

      :tuple: get all pixels defined by tuple (e.g., [0:10,0:10] for first 10 rows and columns in a single band image)

      :Jim object: get all pixels where the specified raster dataset object is > 0

      :VectorOgr: get spatial subset of all pixels covered by the specified vector dataset object

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

   .. method:: Jim[item]=

        Set items of the raster dataset.

        :param items can be of type:

        :tuple: set all pixels defined by tuple (e.g., [0:10,0:10] for first 10 rows and columns in a single band image)
        :Jim object: set all pixels where the specified raster dataset object is > 0
        :VectorOgr: set all pixels covered by the specified vector dataset object

        Modifies the instance on which the method was called.

        Example:

        Set first 10 columns in first 10 rows to 255::

          jim0[0:10,0:10]=255

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
   

====================
Comparison operators
====================

      Jim objects support comparison operations performed at pixel level. The result is a new binary Jim object (of type Byte) with value 1 if the comparison result is True and value 0 if the comparison result is False.

      :equality: ``==``

      Example:

      Pixel wise check for equality. The result is a binary Jim object of type Byte: 1 if pixel values are equal and 0 if objects differ::

        result = jim1==jim2

      Set all pixel to 0 where Jim objects are equal::

        jim1[jim1==jim2] = 0

      :inequality: ``!=``

      Example:

      Pixel wise check for inequality. The result is a binary Jim object of type Byte: 0 if pixel values are equal and 1 if objects differ::

        result = jim1!=jim2

      Set all pixel to 0 where Jim objects differ::

        jim1[jim1!=jim2] = 0

      :less than: ``<``

      Example:

      Pixel wise check if an image is less than another::

        result = jim1<jim2

      Set all pixel values less than 0 to 0::

        jim0[jim0<0] = 0

      :less or equal than: ``<=``

      Example:

      Pixel wise check if an image is less than or equal to another::

        result = jim1<=jim2

      Set all pixel values less than or equal to 0 to 0::

        jim0[jim0<=0] = 0

      :greater than: ``>``

      Example:

      Pixel wise check if an image is greater than another::

        result = jim1>jim2

      Set all pixel values greater than 0 to 0::

        jim0[jim0>0] = 0

      :greater or equal than: ``>=``

      Example:

      Pixel wise check if an image is greater than or equal to another::

        result = jim1>=jim2

      Set all pixel values greater than or equal to 0 to 0::

        jim0[jim0>=0] = 0

=================
Boolean operators
=================

      Jim objects support boolean operations performed at pixel level. Both input and result are assumed to be binary Jim objects of type Byte (0 is False, 1 is True).

      :or: ``|`` (bitwise or)

      Example:

      Calculate the bitwise or value of to Jim object::

        result=jim1 | jim2

      Get a binary mask of all values not within [0,250] (notice the parentheses to take care of the precedence of the operators!) ::

        result = jim0[(jim0<0) | (jim0>250)]

      :xor: ``^`` (bitwise exclusive or)

      Example:

      Calculate the bitwise exclusive or (xor) value of to Jim object::

        result=jim1 ^ jim2

      :and: ``&`` (bitwise and)

      Example:

      Calculate the bitwise and value of to Jim object::

        result=jim1 & jim2

      Get a binary mask of all values within [100,200] (notice the parentheses to take care of the precedence of the operatands!) ::

        result = jim0[(jim0>=100) & (jim0<=200)]


==========================
Arithmetic unary operators
==========================

      Jim objects support unary arithmetic operations performed at pixel level.

      :abs: ``abs`` (absolute value)

      Example:

      Calculate the absolute value of Jim object::

        jimabs=abs(jim)

      :neg: ``-`` (negation). Notice that the output data type can be changed if the input was not signed. A warning message is given.

      Example:

      Calculate the negation of a jim object::

        jimneg=-jim

===========================
Arithmetic binary operators
===========================

      Jim objects support binary arithmetic operations performed at pixel level.

      :addition: ``+``

      Example:

      Add two Jim objects and return the result as a new Jim object::

        jim=jim1+jim2

      :addition and assignment operator: ``+=``

      Example:

      Replace the content of jim1 with the sum of jim1 and jim2::

        jim1+=jim2

      :subtraction: ``-``

      Example:

      Subtract two Jim objects and return the result as a new Jim object::

        jim=jim1-jim2

      :subtraction and assignment operator: ``-=``

      Example:

      Replace the content of jim1 with the difference of jim1 and jim2::

        jim1-=jim2

      :multiplication: ``*``

      Example:

      Multiply two Jim objects and return the result as a new Jim object::

        jim=jim1*jim2

      :multiplication and assignment operator: ``*=``

      Example:

      Replace the content of jim1 with the multiplication of jim1 and jim2::

        jim1*=jim2

      :division: ``/``

      Example:

      Divide two Jim objects and return the result as a new Jim object::

        jim=jim1/jim2

      :division and assignment operator: ``/=``

      Example:

      Replace the content of jim1 with the division of jim1 and jim2::

        jim1/=jim2

      :modulo: ``%`` (return the remainder of a division)

      Example:

      Calculate the modulo of 2 (remainder after division of 2) and return the result as a new Jim object::

        jim=jim0%2

      :modulo and assignment operator: ``%=``

      Example:

      Replace the content of jim1 with the modulo 2 (remainder after division of 2)::

        jim%=2

      :bitwise left shift: ``<<`` (Calculate the bitwise left shift (shifts the bits of the pixel values by the specified number)

      Example:

      Multiply the pixel values by 2 and return as a new Jim object::

        jim2 = jim1 << 2

      :bitwise left shift and assignment operator: ``<<=`` Replace the current object by the bitwise left shifted object (shifts the bits of the pixel values by the specified number)

      Example:

      Multiply the pixel values by 2 and replace the current Jim object::

        jim<<=2

      :bitwise right shift: ``<<`` (Calculate the bitwise right shift (shifts the bits of the pixel values by the specified number)

      Example:

      Divide the pixel values by 2 and return as a new Jim object::

        jim2 = jim1 >> 2

      :bitwise right shift and assignment operator: ``<<=`` Replace the current object by the bitwise right shifted object (shifts the bits of the pixel values by the specified number)

      Example:

      Divide the pixel values by 2 and replace the current Jim object::

        jim>>=2

********************
Accessing properties
********************

.. automodule:: properties
   :members:
                
.. autoclass:: _Properties
   :members:

.. autoclass:: _PropertiesList
   :members:

.. autoclass:: _PropertiesVect
   :members:

************
Input/Output
************

.. automodule:: pjio
   :members:

.. autoclass:: _IO
   :members:

.. autoclass:: _IOList
   :members:

.. autoclass:: _IOVect
   :members:

****************
Pixel operations
****************

.. automodule:: pixops
   :members:

.. autoclass:: _PixOps
   :members:

.. autoclass:: _PixOpsList
   :members:

.. autoclass:: _PixOpsVect
   :members:

***********************
Neighborhood operations
***********************

.. automodule:: ngbops
   :members:

.. autoclass:: _NgbOps
   :members:
     
.. autoclass:: _NgbOpsList
   :members:
     
.. autoclass:: _NgbOpsVect
   :members:
     

*******************
Geometry operations
*******************

.. automodule:: geometry
   :members:

.. autoclass:: _Geometry
   :members:

.. autoclass:: _GeometryList
   :members:

.. autoclass:: _GeometryVect
   :members:


******************************
Connected component operations
******************************

.. automodule:: ccops
   :members:

.. autoclass:: _CCOps
   :members:

.. autoclass:: _CCOpsList
   :members:

.. autoclass:: _CCOpsVect
   :members:

**************
Classification
**************

.. automodule:: classify
   :members:

.. autoclass:: _Classify
   :members:

.. autoclass:: _ClassifyList
   :members:

.. autoclass:: _ClassifyVect
   :members:

*****************
Digital elevation
*****************

.. automodule:: demops
   :members:

.. autoclass:: _DEMOps
   :members:
.. autoclass:: _DEMOpsList
   :members:
.. autoclass:: _DEMOpsVect
   :members:
