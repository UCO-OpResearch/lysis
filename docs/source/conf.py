# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src', 'python')))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Lysis'
copyright = '2022, Brittany Bannish <bbannish@uco.edu> & Brad Paynter <bpaynter@uco.edu>'
author = 'Brittany Bannish <bbannish@uco.edu> & Brad Paynter <bpaynter@uco.edu>'
release = '0.1'
version = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',  # Add source code links
    'sphinx.ext.githubpages',  # Create .nojekyll file for GitHub Pages
]

templates_path = ['_templates']
exclude_patterns = ['*.ipynb', '.ipynb_checkpoints']

# -- Autodoc configuration ---------------------------------------------------
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'undoc-members': True,
    'show-inheritance': True,
}

# Mock imports for optional modules that aren't needed for documentation
# These are specialized dependencies only used in specific contexts:
# - cupy: GPU array library (only for CUDA acceleration)
# - nvtx: NVIDIA profiling tools (only for performance profiling)
# - GooseSLURM: SLURM cluster management (only for HPC job submission)
autodoc_mock_imports = [
    'cupy',
    'cp',
    'nvtx',
    'GooseSLURM',
    'gs',
]

# -- Napoleon settings -------------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
