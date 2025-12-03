import os
import sys
import inspect

# Add the path to mddb_workflow location
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../scripts'))


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MDDB_Workflow'
copyright = '2025, IRB Barcelona'
author = 'IRB Barcelona'
github_user = 'mmb-irb'
github_repo = 'MDDB-workflow'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

extensions = [
    'sphinx.ext.autodoc',       # To automatically pull in docstrings
    'sphinx.ext.napoleon',      # To parse Google/NumPy style docstrings
    'sphinx_rtd_theme',         # ReadTheDocs theme
    'myst_parser',              # For Markdown support
    'sphinx.ext.linkcode',      # For linking to source code
    'nbsphinx',                 # For Jupyter Notebook support
    'format_substitutions',
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
html_theme_options = {'style_nav_header_background': "#FFFFFF", }
html_static_path = ['_static']
html_context = {
    'display_github': True,
    'github_user': github_user,
    'github_repo': github_repo,
    'github_version': 'master/docs/source/',
}


def linkcode_resolve(domain, info):
    """Get links to source code in the documentation."""
    if domain != 'py':
        return None
    if not info['module']:
        return None

    # Get the module and object
    module = sys.modules.get(info['module'])
    if module is None:
        return None

    obj = module
    for part in info['fullname'].split('.'):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None

    # Try to get the file path and line number
    try:
        # Get the absolute path to the source file
        source_file = inspect.getsourcefile(obj)
        if source_file is None:
            return None

        # Get the path relative to the project root
        # This works both locally and on ReadTheDocs
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
        filepath = os.path.relpath(source_file, start=project_root)

        lineno = inspect.getsourcelines(obj)[1]
    except (TypeError, OSError):
        return None

    return f"https://github.com/{github_user}/{github_repo}/blob/master/{filepath}#L{lineno}"
