import setuptools

__version__ = '0.0.0'

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
  name = 'pyssian',
  packages = setuptools.find_packages(), #['pyssian']
  version = __version__,
  description = 'Parser Library for Gaussian Files',
  author = 'Raúl Pérez-Soto',
  author_email = 'rperezsoto.research@gmail.com',
  long_description=long_description,
  long_description_content_type="text/x-rst",
  url = 'https://github.com/rperezsoto/pyssian',
  keywords = ['compchem', 'gaussian','parser'],
  classifiers = ["Programming Language :: Python :: 3",],
  install_requires=['setuptools','pathlib','numpy'],
  python_requires='>=3.7',
  include_package_data=True
)
