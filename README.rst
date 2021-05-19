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
- python library: numpy

Installing pyssian
..................

(Work In Progress) Installing from PyPI: 
++++++++++++++++++++++++++++++++++++++++

.. code:: shell-session

   $ python -m pip install pyssian

Installing from source: 
+++++++++++++++++++++++

Check your python version: 

.. code:: shell-session

   $ python --version

If the version is not one of the required substitute in the following commands 
"python" by the path to the python with the appropiate version. Some systems 
do not have the python symlink but come with the python3 symlink.

If working within a virtual environment with the appropiate python version, the 
installation should proceed as normal using "python". 

Get the source code either git or download and unpack it into "pyssian"
(Currently you need the developer permission to access the code)

.. code:: shell-session

   $ git clone https://github.com/maserasgroup-repo/pyssian.git pyssian

Now install pyssian

.. code:: shell-session

   $ python -m pip install pyssian

*Installing with the -e option before pyssian will make that
all the changes in the source file will be effective without havin to reinstall*

Or if you prefer you can download the source code in a compressed folder and 
install directly: 

.. code:: shell-session
 
   $ python -m pip install pyssian.tar.gz


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
   $ python -m pip install numpy

If it proceeded without any errors (pip and setuptools should already be installed)

Running the tests
-----------------

To run the tests run the following command:

.. code:: shell-session

   $ python -m unittest -v pyssian.tests 


Building the docs
-----------------

To build the docs go to the docs folder where you have the source code then run:

.. code:: shell-session

   $ cd docs/
   $ python -m pip install -r requirements.txt 
   $ make html

Now if you go to the build/html folder you can open the index.html file in your 
folder and navigate through the documentation easily. 

Developed with
--------------

- python 3.7
- Ubuntu 16.04 LTS, 18.04 LTS and 20.04 LTS

.. examples-msg

Examples
--------

Please open the Examples.rst in github to visualize the basic usage examples
or read the documentation.

.. project-author-license

Authors
-------

* **Raúl Pérez-Soto** - [rperezsoto](https://github.com/rperezsoto)
* **Maria Besora** - [MaBeBo](https://github.com/MaBeBo)
* **Feliu Maseras** - [maserasgroup](https://github.com/maserasgroup)

License
-------

pyssian is freely available under an [MIT](https://opensource.org/licenses/MIT)
