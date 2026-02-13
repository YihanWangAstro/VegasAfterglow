# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
import os
import sys
from unittest.mock import MagicMock

# Mock C++ extension modules for documentation building
class BetterMock(MagicMock):
    __all__ = []

    @classmethod
    def __getattr__(cls, name):
        mock = MagicMock()
        mock.__doc__ = f"Mocked {name} class for documentation."
        if name == '__init__':
            mock.__signature__ = None
        return mock

# Create fake module structure
systems_mock = type('VegasAfterglowC', (), {'__all__': []})

sys.modules['VegasAfterglow.VegasAfterglowC'] = systems_mock
sys.path.insert(0, os.path.abspath('../../'))

project = 'VegasAfterglow'
copyright = '2024, VegasAfterglow Team'
author = 'VegasAfterglow Team'
release = '0.1.0'  # Update this with your actual version

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.graphviz',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.ifconfig',
    'sphinx.ext.extlinks',  # For easily linking to external sites
    'breathe',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

# Enable todos
todo_include_todos = True

# Default role for inline markup
default_role = 'any'

# Define common links
extlinks = {
    'doxygen': ('doxygen/%s', '%s'),
    'source': ('doxygen/files.html#%s', '%s')
}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '../../assets/logo.svg'
html_favicon = '../../assets/logo.svg'
html_theme_options = {
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    'globaltoc_collapse': False,
    'globaltoc_maxdepth': 4,
}
html_css_files = [
    'css/custom.css',
]

# Add javascript for source code toggling
html_js_files = [
    'js/custom.js',
]

# Add syntax highlighting style
pygments_style = 'material'
highlight_language = 'python'  # Ensuring Python is the default language
# GitHub Pages settings
html_baseurl = 'https://yihanwangastro.github.io/VegasAfterglow/docs/'

# Create a custom css file for basic styling
css_dir = os.path.join(os.path.dirname(__file__), '_static', 'css')
os.makedirs(css_dir, exist_ok=True)
with open(os.path.join(css_dir, 'custom.css'), 'w') as f:
    f.write("""
/* Basic styling for documentation */
dl.cpp.function {
    margin-bottom: 15px;
    padding: 10px;
    border-radius: 5px;
    background-color: #f7f7f7;
}

dl.cpp.class {
    padding: 10px;
    margin: 10px 0;
    border: 1px solid #eee;
    border-radius: 5px;
}
""")

# -- Breathe configuration ---------------------------------------------------
breathe_projects = {
    "VegasAfterglow": "../doxygen/xml"
}
# -- Breathe project and member defaults ------------------------------------
breathe_default_project               = "VegasAfterglow"                      # valid config key  [oai_citation:0‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_default_members               = ('members', 'undoc-members')           # valid config key  [oai_citation:1‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)

# -- What to show in the documentation --------------------------------------
breathe_show_include                  = True                                  # valid config key  [oai_citation:2‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_show_enumvalue_initializer    = True                                  # valid config key  [oai_citation:3‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_show_define_initializer       = True                                  # valid config key  [oai_citation:4‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)

# -- Parameter list placement (replaces breathe_separate_parameterlist) ------
breathe_order_parameters_first        = False                                  # valid config key  [oai_citation:5‡GitHub](https://github.com/breathe-doc/breathe/blob/main/breathe/renderer/sphinxrenderer.py?utm_source=chatgpt.com)

# -- Project-wide linking and file handling ---------------------------------
breathe_use_project_refids            = True                                  # valid config key  [oai_citation:6‡Breathe](https://breathe.readthedocs.io/en/latest/directives.html?utm_source=chatgpt.com)
breathe_implementation_filename_extensions = ['.c', '.cc', '.cpp']

# Additional Breathe customization for better detailed documentation
breathe_domain_by_extension = {
    "h": "cpp",
    "hpp": "cpp",
    "c": "c",
    "cpp": "cpp",
    "cc": "cpp",
}
breathe_domain_by_file_pattern = {
    '*/include/*': 'cpp',
    '*/src/*': 'cpp',
}
breathe_ordered_classes = True
breathe_show_include_files = True
breathe_doxygen_mapping = {
    # Map implementation file elements to their header declarations
    'function': 'function',
    'define': 'define',
    'property': 'variable',
    'variable': 'variable',
    'enum': 'enum',
    'enumvalue': 'enumvalue',
    'method': 'method',
    'typedef': 'typedef',
    'class': 'class',
    'struct': 'struct',
}

# Improved debug options for troubleshooting
breathe_debug_trace_directives = True
breathe_debug_trace_doxygen_ids = True
breathe_debug_trace_qualification = True

# Enhanced options for template and inline function documentation
breathe_template_relations = True
breathe_inline_details = True
breathe_show_define_initializer = True
breathe_show_template_parameters = True
breathe_show_templateparams = True

# -- intersphinx configuration -----------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# -- autodoc configuration ---------------------------------------------------
autodoc_member_order = 'groupwise'  # Changed to groupwise for better organization
autodoc_typehints = 'both'  # Changed to 'both' to show in signature and description
autodoc_default_options = {
    'members': True,
    'member-order': 'groupwise',
    'undoc-members': True,
    'private-members': True,  # Show private members
    'special-members': True,  # Show special members
    'show-inheritance': True,
    'inherited-members': True
}

# -- napoleon configuration --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_attr_annotations = True

# -- Source code link configuration ------------------------------------------
viewcode_follow_imported_members = True
viewcode_enable_epub = True

# Simple setup function
def setup(app):
    # Add custom CSS file for styling
    app.add_css_file('css/custom.css')

# Create a custom css file for styling code blocks
css_dir = os.path.join(os.path.dirname(__file__), '_static', 'css')
os.makedirs(css_dir, exist_ok=True)
with open(os.path.join(css_dir, 'custom.css'), 'w') as f:
    f.write("""
/* Basic styling for C++ documentation */
dl.cpp.function {
    margin-bottom: 15px;
    padding: 10px;
    border-radius: 5px;
    background-color: #f7f7f7;
}

dl.cpp.class {
    padding: 10px;
    margin: 10px 0;
    border: 1px solid #eee;
    border-radius: 5px;
}
""")
