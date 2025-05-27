import os
import sys
# Add the path to your project's root directory
sys.path.insert(0, os.path.abspath('../../')) # Adjust this path if your 'my_package' is elsewhere

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MDDB_Workflow'
copyright = '2025, IRB Barcelona'
author = 'IRB Barcelona'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = [
    'sphinx.ext.autodoc',       # To automatically pull in docstrings
    'sphinx.ext.napoleon',      # To parse Google/NumPy style docstrings
    'sphinx_rtd_theme',         # ReadTheDocs theme
    'myst_parser',              # For Markdown support
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Source suffix (this is for input files, still RST if that's what you're writing)
source_suffix = '.rst'

def setup(app):
    app.add_css_file('custom.css')
# Add the ReadTheDocs theme
html_theme = 'sphinx_rtd_theme'
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-html_favicon
html_favicon = '_static/MDDB_favicon.png'
html_logo = '_static/MDDB_Logo_colour.png'
html_theme_options = {'style_nav_header_background': "#FFFFFF",}
html_static_path = ['_static']
