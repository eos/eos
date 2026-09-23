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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

import eos
import os

project = 'EOS'
copyright = '2019-2026, The EOS Authors'
author = 'The EOS Authors'
master_doc = 'index'

# The full version, including alpha/beta/rc tags
release = eos.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinxarg.ext',
    'sphinxcontrib.contentui',
    'nbsphinx'
]

# Without these, the sections that sphinxarg generates per subcommand all claim
# the same label, and equally named sections collide across documents.
autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 2

# Resolve references to the types that EOS' API exposes from these projects. The
# inventories are pinned in '_inventories/', so that the build needs no network
# access and does not change when one of these projects publishes new documentation.
# Refresh an inventory by downloading '<url>/objects.inv' over the pinned file.
intersphinx_mapping = {
    'python':     ('https://docs.python.org/3',             '_inventories/python.inv'),
    'numpy':      ('https://numpy.org/doc/stable',           '_inventories/numpy.inv'),
    'scipy':      ('https://docs.scipy.org/doc/scipy',       '_inventories/scipy.inv'),
    'matplotlib': ('https://matplotlib.org/stable',          '_inventories/matplotlib.inv'),
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

todo_include_todos = True

# Suppress 'WARNING: Citation [...] is not referenced.' messages.
# Implemented references need to be listed as documentation itself and should
# also be available for citation in other parts of the documentation.
# Many references, however, remain intentionally uncited.
suppress_warnings = ['ref.citation']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
html_js_files = [
    'js/version-banner.js',
]
html_baseurl = 'https://eos.github.io/doc/'

# Release builds are archived below 'releases/' in the documentation repository,
# two levels below the root of the documentation site.
eos_release = os.environ.get('EOS_RELEASE', '')
html_context = {
    'eos_release':  eos_release,
    'eos_doc_root': '../../' if eos_release else '',
}

html_show_sourcelink = False

html_logo = '_static/github-eos-logo.png'
logo_only = True


# The classes that the bindings expose under a private name, mapped to the names
# under which they are documented.
documented_names = {
    '_Parameters': 'eos.Parameters',
}

# Boost.Python ends the signature it generates with a colon, which autodoc then
# reads as part of the return annotation.
def strip_signature_colon(app, what, name, obj, options, signature, return_annotation):
    if return_annotation and return_annotation.endswith(' :'):
        return_annotation = return_annotation[:-2]

    return (signature, documented_names.get(return_annotation, return_annotation))


def setup(app):
    app.connect('autodoc-process-signature', strip_signature_colon)
