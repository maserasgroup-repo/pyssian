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

MultiGaussianInFile
-------------------

*  Allows the concatenation of multiple GaussianInFile instances to create 
   a linked (--Link1--) gaussian input.
*  Provides some convenience methods to enforce same resources across each 
   concatenated job
*  It was conceived for writing gaussian input files rather than reading 
   gaussian input files with linked calculations, so no much testing was 
   put into this task.

.. autoclass:: pyssian.MultiGaussianInFile
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
