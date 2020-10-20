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

  jim = pj.Jim()
  help(jim.geometry.warp)
