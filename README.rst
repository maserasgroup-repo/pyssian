==============
pyssian
==============

----------------------------------------
A python interface to Gaussian i/o files
----------------------------------------

.. project-description-start

This project focuses on writing an object oriented library for parsing Gaussian
io files and provide usefull command line tools for computational chemists.


.. setup-instructions

Getting Started
---------------

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Prerequisites
.............

- python >= 3.6
- python library: setuptools

Installing the dependencies
...........................

python3.7 installation in Ubuntu 18.04 LTS

.. code:: shell-session

   $ sudo apt-get update
   $ sudo apt-get install python3.7 python3.7-dev


If for any reason it is not reachable:

.. code:: shell-session

   $ sudo add-apt-repository ppa:deadsnakes/ppa
   $ sudo apt-get update
   $ sudo apt-get install python3.7 python3.7-dev

Now you can skip this if you don't want to set up a virtual environment
(Remember to change NewFolder for the actual path of the directory where you
want the virtual environment)

.. code:: shell-session

   $ sudo apt-get install python3.7-venv
   $ python3.7 -m venv NewFolder
   $ source NewFolder/bin/activate

Now we install the python default installer pip

.. code:: shell-session

   $ python -m pip install pip
   $ python -m pip install --upgrade pip
   $ python -m pip install setuptools

If it proceeded without any errors (pip and setuptools should already be installed)

Now to install the optional dependencies (Skip if you don't desire them):

.. code:: shell-session

   $ python -m pip install numpy
   $ python -m pip install matplotlib

Installing pyssian
..................


Get the source code either git or download and unpack it into "pyssian"
(Currently you need the developer permission to access the code)

.. code:: shell-session

   $ git clone https://github.com/maserasgroup-repo/pyssian.git pyssian

Now install pyssian

.. code:: shell-session

   $ python -m pip install pyssian


Installing with the -e option before NewDir will make that
all the changes in the source file will have a visible effect when using pyssian


Running the tests
-----------------

To run the tests run the following command:

.. code:: shell-session

   $ python -m unittest pyssian.tests -v


Developed with
--------------

- python 3.7.3
- Ubuntu 16.04 LTS

.. examples-msg

Examples
--------

Please open the Examples.rst in github to visualize the basic usage examples
or read the documentation.

.. project-author-license

Authors
-------

* **Raúl Pérez-Soto** - [rperezsoto](https://github.com/rperezsoto)


License
-------

(None currently)
