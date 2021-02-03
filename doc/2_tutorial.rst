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

The primary source of documentation is available online via this `website <https://jeodpp.jrc.ec.europa.eu/services/processing/pyjeohelp>`_.
The documentation is divided in three main parts:

* :ref:`Introduction`: Introduction explaining the organization of the pyjeo package, its modules, design, installation and usage.

* :ref:`Tutorial`: A tutorial to guide you through the first steps of using pyjeo, introducing the main data structures :py:class:`Jim`, :py:class:`JimList`, and :py:class:`JimVect`.

* :ref:`Reference_manual`: A manual describing all functions and methods available in pyjeo


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

  jim = pj.Jim()
  help(jim.geometry.warp)

.. _Tutorial_jim:

=======================
Tutorial on Jim objects
=======================

Jim is the main class to represent raster data objects. You can create a new Jim either :ref:`from file <create_Jim_from_file>`
  or :ref:`create new <create_Jim_new>` Jim by specifying all the attributes (e.g., data type, columns, rows, etc.).

When creating an geospatial image from file (e.g., in GeoTIFF format), the attributes for the geotransform and projection is set automatically and stored in the Jim object. These attributes can be retrieved with the methods in the :py:mod:`properties` module::


  import pyjeo as pj

  jim = pj.Jim('/path/to/raster.tif'))
  jim.properties.getBBox()
  
  [399960.0, 5100000.0, 405080.0, 5094880.0]


Where the list returned represents the upper left x, upper left y, lower right x, and lower right y coordinates respectively. 

.. _data_model:

Data model
==========

It is important to understand the data model that is used for a Jim object. For a more detailed description, please refer to :cite:`kempeneers2019`.

The data model used a multi-band three dimensional (3D) model. Each band represents a 3D contiguous array in memory, where data are organized as *[plane][row][column]* (see also :py:meth:`Jim.np`). Two dimensions refer to the spatial domain (x and y) and pyjeo refers to the third dimension as *plane* (see :numref:`cube`). This plane is typically used for either the temporal or spectral dimension, but can also be used to address volumetic data.  The data cube has a single georeference.

.. _cube:

.. figure:: figures/cube.png
   :width: 100 %

   Data model used for Jim objects.

When reading a multi-band raster dataset into a Jim raster data object, the user can choose if the bands are to be considered as planes, resulting in a single band 3D object, or as bands (resulting in a multi-band 2D raster object). Planes and bands can be stacked :py:meth:`~geometry._Geometry.stackPlane`, :py:meth:`~geometry._Geometry.stackBand` and subset (:py:meth:`~geometry._Geometry.cropPlane`, :py:meth:`~geometry._Geometry.cropBand`) as shown in :numref:`cube`. Dimensions must correspond across all bands within the same *Jim* object. To collect objects with different dimensions, a *JimList* can be used (which inherits from a plain Python list). Multi-dimensional data with dimensions above three are not supported in pyjeo.

Functions on the individual bands can easily be processed in parallel as different data pointers are used for each band. Warping multi-band images take advantage of the single georeference object they have in common (only a single function call to *GDALReprojectImage* (`gdalwarper.h <https://gdal.org/doxygen/gdalwarper_8h.html>`_ is needed).

The simple data model with contiguous arrays allows to combine pyjeo with other software such as GDAL, GSL, and Python packages that are compatible with Numpy arrays (e.g., `xarray <http://xarray.pydata.org>`_ , `SciPy <https://www.scipy.org/>`_). A direct bridge to Numpy arrays and xarray (if installed) is included in pyjeo, without the generation of an extra copy in memory.

In particular, when dealing with geospatial data that have a large memory footprint, careful memory handling is important. To this end, users can choose if pyjeo functions modify objects in-place or return a new object (see also :ref:`functions_methods`).  

Jim conversions and bridge to third party packages
==================================================

Jim objects can be easily converted to Numpy array and xarray objects either with or without duplicating the memory (see also :ref:`jim_conversions`). A Numpy array object derived from a Jim object without a memory copy references to the same data in memory as the original Jim object. This reduces the memory footprint, but can lead to memory errors. The Jim object should remain the owner of the data and the referenced Numpy array object should not be altered in shape nor destroyed.

However, if handled with care, this can be a powerful technique.
As shown in :ref:`ndimage`, third party libraries operating on Numpy arrays can directly written into Jim objects. For instance, to Gaussian filter a Jim object using `SciPy <https://www.scipy.org/>`_, simply use::

  jim.np()[:] = ndimage.gaussian_filter(jim.np(), 2)[:]

.. _Tutorial_classification:

==========================
Tutorial on classification
==========================

.. _random_forest_classifier:

---------------------------------
Random Forest ensemble classifier
---------------------------------

Import relevant modules::

  from sklearn.ensemble import RandomForestClassifier
  from sklearn.model_selection import train_test_split
  from sklearn.metrics import confusion_matrix
  from sklearn.metrics import accuracy_score

We will use a vector file that contains reference data (in a numerical field 'label')

.. image:: figures/labels.png
   :width: 100 %

Load the vector file in a JimVect object::

  reference = pj.JimVect('training.sqlite')

Create a 3D Jim object that loads the raster data containing all features (read bands as planes)::

  jim = pj.Jim('/path/to/raster.tif', band2plane=True)

Extract features from the Jim object::

  featurevect = jim.geometry.extract(reference, rule=['allpoints'],
                                     output='/vsimem/features.sqlite',
                                     oformat='SQLite',
                                     co=['OVERWRITE=YES'],
                                     classes=[1, 2, 3, 4],
                                     copy='label',
                                     fid='fid')

Use the Numpy representation of the vector feature data to create arrays for the features (x) and label data (y)::

  x = featurevect.np()[:, 1:]
  y = featurevect.np()[:, 0:1]

Split the data in a training and test set::

  x_train, x_test, y_train, y_test = train_test_split(x, y,
                                                      test_size=0.33,
                                                      random_state=42)

Create a Random Forest classifier from sklearn and fit the model using the training data ::

  rfModel = RandomForestClassifier(n_estimators=100,
                                   max_depth=9,
                                   min_samples_leaf=5,
                                   min_samples_split=3,
                                   criterion='gini')
  rfModel.fit(x_train, y_train.ravel())

Use the unseen test set to perform an accuracy assessment::

  y_predict = rfModel.predict(x_test)
  print(confusion_matrix(y_test, y_predict))
  print('accuracy score: {}'.format(accuracy_score(y_test, y_predict)))

  [[29  8  2  6  0]
  [ 2 14  2  5  2]
  [ 0  2 13  3  9]
  [ 6  2  1 24  0]
  [ 0  0 15  0 20]]

  accuracy score: 0.606060606061


Classify the image using the Numpy representation of the Jim object::

  x = jim.np()
  x = x.reshape(jim.properties.nrOfPlane(), jim.properties.nrOfRow() * \
                jim.properties.nrOfCol()).T

  jim_class = pj.Jim(ncol=jim.properties.nrOfCol(),
                     nrow=jim.properties.nrOfRow(),
                     otype='Byte')
  jim_class.properties.copyGeoReference(jim)
  jim_class.np()[:] = rfModel.predict(x).astype(np.dtype(np.uint8)).\
      reshape(jim.properties.nrOfRow(), jim.properties.nrOfCol())

Show the classified map with matplotlib::

  import matplotlib.pyplot as plt

  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.imshow(jim_class.np())
  plt.show()


.. image:: figures/rf_class.png
   :width: 100 %

.. _svm_classifier:

---------------------------------
Support Vector Machine classifier
---------------------------------

Import relevant modules::

  from sklearn.svm import SVC
  from sklearn import preprocessing
  from sklearn.model_selection import train_test_split
  from sklearn.metrics import confusion_matrix
  from sklearn.metrics import accuracy_score

Refer to random_forest_classifier_ to create a training and test data set

Create a support vector machine (SVM) classifier from sklearn and fit the model using the training data ::

  svmModel = SVC(gamma='auto')
  svmModel.fit(preprocessing.MinMaxScaler().fit_transform(x_train,
                                                          y_train.ravel())

Classify the image using the Numpy representation of the Jim object::

  x = jim.np()
  x = x.reshape(jim.properties.nrOfPlane(), jim.properties.nrOfRow() * \
                jim.properties.nrOfCol()).T


  jim_class = pj.Jim(ncol=jim.properties.nrOfCol(),
                     nrow=jim.properties.nrOfRow(), otype='Byte')
  jim_class.properties.copyGeoReference(jim)
  jim_class.np()[:] = svmModel.predict(preprocessing.MinMaxScaler().\
      fit_transform(x)).astype(np.dtype(np.uint8)).reshape(
          jim.properties.nrOfRow(), jim.properties.nrOfCol())

.. _Symbolic machine learning:

-------------------------
Symbolic machine learning
-------------------------

For the reference for this method please refer to :py:meth:`~classify._Classify.sml` in the :py:mod:`classify` module.

Import Numpy, Path and pyjeo::

  import numpy as np
  from pathlib import Path
  import pyjeo as pj

Locate the the input and reference data (Corine Land Cover)::

  datadir = Path.home() / 'pyjeo/tests/data'
  clc = datadir / 'clc_32632.tif'
  testFile = datadir / 'modis_ndvi_2010.tif'

The input image is multi-band image based on the NDVI values based on MODIS acquisitions in the year 2010. Each band corresponds to a monthly composite (12 bands in total). We open a spatial subset of the input image, making sure the reference image entirely covers the input image. The bands are loaded in a 3D image as planes::

  bbox = [4246000, 2547000, 4349500, 2441000]
  jim = pj.Jim(testFile, band2plane=True,
               ulx=bbox[0], uly=bbox[1], lrx=bbox[2], lry=bbox[3])

Show the first three months of the NDVI image image. The dimension of the 3D image is organized as [plane][row][column]. We need to roll the first axis due to obtain a [row][col][colour] image (see also :py:meth:`Jim.np`) ::

  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.imshow(np.rollaxis(jim[0:3,:,:].np(),0,3))

.. image:: figures/modis_ndvi_2010.png
   :width: 100 %

Create the reference dataset that contains only the target classes (2, 12, 25, 41, 50)::

  class_dict = {'urban': 2, 'agriculture': 12, 'forest': 25,
                'water': 41, 'rest': 50}


  class_from = range(0, 50)
  class_to = [50] * 50
  for i in range(0, 50):
      if 1 <= i < 10:
          class_to[i] = class_dict['urban']
      elif 11 <= i < 22:
          class_to[i] = class_dict['agriculture']
      elif 23 <= i < 25:
          class_to[i] = class_dict['forest']
      elif 40 <= i < 45:
          class_to[i] = class_dict['water']
      else:
          class_to[i] = class_dict['rest']

The spatial resolution of the input image is 500 m. The SML algorithm typically expects a reference dataset with a coarser spatial resolution, for example 1000 m::

  jim_ref = pj.Jim(clc, dx=1000, dy=1000)
  jim_ref.classify.reclass(classes=list(class_from), reclasses=class_to)

The reference image should be in the same projection as the input image::

  jim_ref.geometry.warp(jim.properties.getProjection())

.. image:: figures/reference_coarse.png
   :width: 100 %

The SML algorithm supports multiple reference images in a JimList object. Here, only a single reference is provided (list of a single image)::

  reflist = pj.JimList([jim_ref])

Train the SML model and save it::

  model = datadir / 'model_sml.dat'
  jim.classify.trainSML(reflist, output=model,
                        classes=sorted(class_dict.values()))

Classify the input image using the trained model::

  sml = pj.classify.sml(jim, model=model)

The result is a multi-band Jim object where the number of bands equals the number of classes and each band represents the probability for the respective class. To obtain a discrete classification result, based on the maximum probability, create a 3D multi-plane image from the multi-band image::

  sml.geometry.band2plane()

Use Numpy to select the maximum probability class. Put the result in the first plane (overwriting the probability plane of the first class)::

  sml.np(0)[:] = np.argmax(sml.np(), axis=0)

Keep only the first plane::

  sml.geometry.cropPlane(0)

Value 0 is a valid class (do not consider 0 as nodata)::

  sml.properties.clearNoData()

Map the classes [0-4] to the original class values::
  
  sml.classify.reclass(classes=[0, 1, 2, 3, 4],
                       reclasses=[2, 12, 25, 41, 50])

.. image:: figures/sml.png
   :width: 100 %
