# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'atooms'
copyright = '2022, Daniele Coslovich'
author = 'Daniele Coslovich'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# import sphinx_rtd_theme
# extensions = [
#    'sphinx_rtd_theme',
# ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The root document.
root_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
# html_theme_options = {
#     'display_version': True,
#     'vcs_pageview_mode': '',
#     # Toc options
#     'collapse_navigation': False,
#     'sticky_navigation': True,
#     'navigation_depth': 4,
#     'includehidden': True,
#     'titles_only': False
# }

pygments_style = 'friendly'
html_theme = 'alabaster'
html_static_path = ['_static/custom.css']
html_theme_options = {
    'description': 'A Python framework for simulations of interacting particles',
    'fixed_sidebar': True,
    'sidebar_collapse': True,
    'extra_nav_links': {'Run this tutorial on Binder': 'https://atooms.frama.io/atooms/',
                        'Org-mode and jupyter notebooks on Framagit': 'https://atooms.frama.io/atooms/'},
    'gray_2': '#F4F4F4ED',
    'sidebar_width': '380px',
    'body_max_width': 'auto',
    'page_width': '1200px',
#    'code_highlight_bg': '#111',
}

html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'searchbox.html',
        'relations.html',
        'donate.html',
    ]
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
