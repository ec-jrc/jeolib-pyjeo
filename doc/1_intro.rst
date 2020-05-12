.. _Introduction:

============
Introduction
============

Organization of pyJEO
---------------------

Modules
^^^^^^^

The pyJEO package is grouped in modules. A module combines a number of operations that belong together. The models defined pyJEO are:

 * properties: get properties of a class object
 * io: input/output operations
 * geometry: operations that relate to geometry
 * pixops: pixel operations
 * ngbops: neighborhood operations
 * clssfy: classification operations
 * ccops: segmentation operations
 * demops: operations on digital elevation models


Functions and methods
^^^^^^^^^^^^^^^^^^^^^

Operations on pyJEO objects are distinguished in functions and methods. Methods directly operate on objects, i.e., instances of a class. A method implicitly has access to the object attributes on which it was called. For instance, methodX that operates on myobject that belongs to moduleA and takes a single argument (arg1) is typically called like this::

  myobject.moduleA.methodX(arg1)
     
Functions that operate on objects must have the objects passed as arguments. An example of a function (functionY) that belongs to the same module (moduleA) and operates on myobject using an extra boolean argument (flag1) is::

  newobject = pyjeo.moduleA.functionY(myobject, flag1=True)

In this case a new object (newobject) is returned. The functions that return a new object are non-destructive, i.e., they do not alter the object that was passed as an argument. This follows the principle of command-query-separation. 

Methods can be either destructive or non-destructive. Destructive methods do not return any object. For instance, in the following example::

  object.moduleB.destructiveMethod()
  
the object is altered and no object is returned. Methods that are non-destructive can return a new object, for instance, the method getMethodZ in module properties returns the attributeZ from myobject::

   myattributeZ = myobject.properties.getAttributeZ()


Design of pyJEO
---------------

Design in C++
^^^^^^^^^^^^^

The pyJEO package largely depends of the JIPlib library, which is implemented in C++ and contains three main classes: Jim, JimList and VectorOgr (see :numref:`Fig:jiplib_classes`). Each of these classes is represented by a Python proxy class in pyJEO: Jim, JimList and JimVec. In addition to the methods from the C++ class JimList, the Python class JimList inherits all methods from a Python list (e.g., :code:`len`, :code:`append`, :code:`extend`, :code:`insert`, :code:`remove`, :code:`count`)


.. _Fig:jiplib_classes:
.. figure:: figures/jiplib_classes.png
    :align: center

    Overview of the classes in the JIPlib library

Building a Python interface via SWIG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Simplified Wrapper and Interface Generator (SWIG) produces a wrapper file that is compiled and linked into a dynamic library (see :numref:`Fig:jiplib_swig`). 

.. _Fig:jiplib_swig:
.. figure:: figures/jiplib_swig.png
    :align: center

    Using SWIG to build a Python interface from C++

   
Installation
------------

From the directory of the repository, run::

  sudo python setup.py install

To test the installation, run::

  python -W ignore -m unittest -v tests

To build the documentation, go to directory doc and run make::

  cd doc
  make html

Usage
-----

Usage in your local environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In your local environment, import the pyjeo module::

  import pyjeo as pj

Usage in the JEOdesk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The pyjeo module is can be imported as::

  import pyjeo as pj

An environment variable has been automatically set for all users::

  PYTHONPATH=/eos/jeodpp/shared/prod/lib/python/jeodesk-16

In case you should change this environment variable (e.g., /foo/bar), please make sure to append it as follows::


  PYTHONPATH=${PYTHONPATH}:/foo/bar

Usage in the execute function in the interactive processing JEOlab
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The pyjeo module is automatically imported as pj

Usage on the cluster
^^^^^^^^^^^^^^^^^^^^
In your condor submit file, use this docker file::
  
  docker_image    =  jeoreg.cidsn.jrc.it:5000/jeodpp-proc/jeodpp_jupyter_inter_py2_deb9


In your execution script launched by the condor submit file, define the following environment variables::

   export PYTHONPATH="${PYTHONPATH}:/eos/jeodpp/shared/dev/lib/python:/eos/jeodpp/shared/dev/lib/python/jeodpp:/eos/jeodpp/shared/dev/lib/python/pyjeo"
   export LD_LIBRARY_PATH=/eos/jeodpp/shared/dev/lib/cpp/jeodpp
