Updating the version
--------------------

Places where the version has to be updated: 

- ./setup.py
- ./pyssian/__init__.py

Notes: 

Bugfixes -> version += 0.0.1
Non-breaking updates -> version += 0.1.0 (resets minor)
Major updates/breaking API updates -> version += 1.0.0 (resets minor and middle) 

Building and updating the Docs
------------------------------

To build the docs for publication go to the docs folder where you have the 
source code then run:

.. code:: shell-session

   $ cd sphinx/
   $ python -m pip install -r requirements.txt
   $ make github

To build locally previous to release: 

.. code:: shell-session

   $ cd sphinx/
   $ python -m pip install -r requirements.txt
   $ make html

Now if you go to the _docs/html folder you can open the index.html file in your 
folder and navigate through the documentation easily. 

Building the distribution
-------------------------

.. code:: shell-session

   $ # Go to where the pyproject.toml file is. 
   $ python -m build 
   $ # The distribution should be at the dist/ folder. 

Uploading to PyPI
-----------------

TODO