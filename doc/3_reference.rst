Reference Manual
================

.. toctree::
   :maxdepth: 4

Global functions
################

.. automodule:: pyjeo
   :members:

.. autoclass:: Jim
   :members:

   .. method:: boil(time=10)

      Boil the noodle *time* minutes.

   .. method:: Jim[ item ]

      Get subset of the raster dataset.

      :param items can be of type:

      :tuple: get all pixels defined by tuple (e.g., [0:10,0:10] for first 10 rows and columns in a single band image)

      :Jim object: get all pixels where the specified raster dataset object is > 0

      :VectorOgr: get spatial subset of all pixels covered by the specified vector dataset object

      :returns: subset of the raster dataset

      Example:

      Get first 10 columns in first 10 rows::

        jim0[0:10,0:10]

      Get a binary mask of all values not within [0,250]::

        jim0[(jim0<0) | (jim0>250)]

      Get a binary mask identifying  pixels where jim0<jim1::

        jim0[(jim0<jim1)]

      Crop a raster dataset according to the extent of a vector dataset (crop_to_cutline), set all pixels not covered to 0 (or value defined as no data)::

        ifn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/IMG_DATA/T32TNR_20170805T102031_B08.jp2'
        cfn='/eos/jeodpp/data/SRS/Copernicus/S2/scenes/source/L1C/2017/08/05/065/S2A_MSIL1C_20170805T102031_N0205_R065_T32TNR_20170805T102535.SAFE/GRANULE/L1C_T32TNR_A011073_20170805T102535/QI_DATA/MSK_CLOUDS_B00.gml'
        jim=pj.io.createJim(ifn)
        v=pj.io.createVector(cfn)
        jimcloud=jim[v]

   .. method:: Jim[ item ]=

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
          jim=pj.io.createJim(ifn)
          v=pj.io.createVector(cfn)
          jim[v]=255
   
Accessing properties
####################

.. automodule:: properties
   :members:
                
.. autoclass:: _Properties
   :members:

Input/Output
############

.. automodule:: pjio
   :members:

Pixel operations
################

.. automodule:: pixops
   :members:

Neighborhood operations
#######################

.. automodule:: ngbops
   :members:

Geometry operations
###################

.. automodule:: geometry
   :members:

.. autoclass:: _Geometry
   :members:


CC operations
#############

.. automodule:: ccops
   :members:

Classification
##############

.. automodule:: clssfy
   :members:

Digital elevation
#################

.. automodule:: demops
   :members:
