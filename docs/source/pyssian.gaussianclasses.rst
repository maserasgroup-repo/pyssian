gaussianclasses
===============

.. automodule:: pyssian.gaussianclasses


To do:
   * Create a Class to represent the basis and core electronic pontentials


GaussianInFile
--------------

*   Accepts a context manager usage similar to 'with open(file) as F:...'
*   Does not parse ONIOM calculation inputs
*   Does not parse MM calculation inputs

.. autoclass:: pyssian.GaussianInFile
   :members:

GaussianOutFile
---------------

*   Accepts a context manager usage similar to 'with open(file) as F:...'

.. autoclass:: pyssian.GaussianOutFile
   :members: FileStructure, GetLinks, update, read, clean


InternalJob
-----------

.. autoclass:: pyssian.InternalJob
   :members:
