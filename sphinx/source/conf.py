# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import configparser
sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'pyssian'
copyright = '2021, Pérez-Soto, R.; Besora, M.; Maseras, F.'
author = 'Pérez-Soto, R.; Besora, M.; Maseras, F.'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.napoleon',
    'sphinx_copybutton'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = {'.rst': 'restructuredtext'}

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# The master toctree document.
master_doc = 'index'

# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'python'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'source/modules.rst']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme' #'nature'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'one-dark' # 'paraiso-dark' #'default', 'emacs', 'friendly', 'colorful', 'autumn', 'murphy', 'manni', 'material', 'monokai', 'perldoc', 'pastie', 'borland', 'trac', 'native', 'fruity', 'bw', 'vim', 'vs', 'tango', 'rrt', 'xcode', 'igor', 'paraiso-light', 'paraiso-dark', 'lovelace', 'algol', 'algol_nu', 'arduino', 'rainbow_dash', 'abap', 'solarized-dark', 'solarized-light', 'sas', 'stata', 'stata-light', 'stata-dark', 'inkpot', 'zenburn', 'gruvbox-dark', 'gruvbox-light'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'pyssian-doc'

# -- Extension configuration -------------------------------------------------
copybutton_prompt_text = "$ "

# ---Multi version support ---------------------------------------------------
# get the environment variable build_all_docs and pages_root
build_all_docs = os.environ.get("build_all_docs")
pages_root = os.environ.get("pages_root", "")

# if not there, we dont call this
if build_all_docs is not None:
    # we get the current language and version
    current_version = os.environ.get("current_version")

    # we set the html_context with current version 
    html_context = {
        'current_version' : current_version,
        'versions' : [],
    }

    # and we append all versions accordingly 
    # we treat the master branch as latest 

    html_context['versions'].append(['latest', pages_root])

    # and loop over all other versions from our versions.ini file
    # to set versions

    versions = configparser.ConfigParser(allow_no_value=True)
    versions.read(['versions.ini'])

    for version in versions.options('versions'):
        html_context['versions'].append([version, f'{pages_root}/{version}'])
