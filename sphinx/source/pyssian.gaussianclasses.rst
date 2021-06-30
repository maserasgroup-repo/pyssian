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
*   Not tested for z-matrix inputs

.. autoclass:: pyssian.GaussianInFile
   :members:

GaussianOutFile
---------------

*   Accepts a context manager usage similar to 'with open(file) as F:...'
*   Will only parse a Gaussian Output File that has the "#p". 

.. autoclass:: pyssian.GaussianOutFile
   :members: print_file_structure, get_links, update, read, clean


InternalJob
-----------

.. autoclass:: pyssian.InternalJob
   :members:
