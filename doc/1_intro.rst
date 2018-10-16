============
Introduction
============

.. toctree::
   :maxdepth: 3

Design in C++
-------------

The JIPlib library is implemented in C++ and contains three main classes: Jim, JimList and VectorOgr (see :numref:`Fig:jiplib_classes`)


.. _Fig:jiplib_classes:
.. figure:: figures/jiplib_classes.png
    :align: center

    Overview of the classes in the JIPlib library

Building a Python interface via SWIG
------------------------------------

The Simplified Wrapper and Interface Generator (SWIG) produces a wrapper file that is compiled and linked into a dynamic library (see :numref:`Fig:jiplib_swig`). 

.. _Fig:jiplib_swig:
.. figure:: figures/jiplib_swig.png
    :align: center

    Using SWIG to build a Python interface from C++

   
Installation
------------


