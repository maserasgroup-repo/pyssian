import setuptools

__version__ = '1.0.0'

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
  name = 'pyssian',
  version = __version__,
  description = 'Parser Library for Gaussian Files',
  long_description=long_description,
  long_description_content_type="text/x-rst",
  url = 'https://github.com/maserasgroup-repo/pyssian',
  author = 'Raúl Pérez-Soto',
  author_email = 'rperezsoto.research@gmail.com',
  classifiers = ['License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.6'
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 ],
  keywords = ['compchem, gaussian, parser'],
  packages = setuptools.find_packages(),
  python_requires='>=3.6, <4',
  install_requires=['setuptools','pathlib','numpy'],
  include_package_data=True,
  package_data = {'test_files': ['pyssian/tests/test_files/*.txt'],
                  'tests' : ['pyssian/tests/*.py']},
  project_urls={'Bug Reports': 'https://github.com/maserasgroup-repo/pyssian/issues',
                'Source': 'https://github.com/maserasgroup-repo/pyssian',
                'Docs' : 'https://pyssian.readthedocs.io/en/latest/'
               },
)
