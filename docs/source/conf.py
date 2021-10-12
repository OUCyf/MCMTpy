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
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'MCMTpy-test'
copyright = '2021, Fu Yin'
author = 'Fu Yin'

# The short X.Y version
version = '0.1.0a1'
# The full version, including alpha/beta/rc tags
release = '0.1.0a1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
]

latex_engine = 'xelatex'



latex_use_xindy = True
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    'papersize': 'a4paper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '12pt',
    'pointsize': '10pt',
    # 'classoptions': ',english',
    'inputenc': '',
    'utf8extra': '',
    'extraclassoptions': 'openany',

    # Additional stuff for the LaTeX preamble.
    'preamble': r'''
\usepackage{xeCJK}
\usepackage{indentfirst}
\setlength{\parindent}{2em}
\setCJKmainfont[BoldFont=SimHei, ItalicFont=STKaiti]{SimSun}
\setCJKmonofont[Scale=0.9]{SimSun}
\setCJKfamilyfont{song}[BoldFont=SimSun]{SimSun}
\setCJKfamilyfont{sf}[BoldFont=SimSun]{SimSun}
''',
    # 'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    # 'printindex': r'\footnotesize\raggedright\printindex',
}

#
# https://juejin.im/post/5c7253c2e51d4512543327b4
#
latex_elements['preamble'] = r"""
\usepackage{xeCJK}
\usepackage{indentfirst}
\setlength{\parindent}{2em}
\setCJKmainfont[BoldFont=SimHei, ItalicFont=STKaiti]{SimSun}
\setCJKmonofont[Scale=0.9]{SimSun}
\setCJKfamilyfont{song}[BoldFont=SimSun]{SimSun}
\setCJKfamilyfont{sf}[BoldFont=SimSun]{SimSun}
\XeTeXlinebreaklocale "zh"
\XeTeXlinebreakskip = 0pt plus 1pt
\parindent 2em
\definecolor{VerbatimColor}{rgb}{0.95,0.95,0.95}
\setcounter{tocdepth}{3}
\renewcommand\familydefault{\ttdefault}
\renewcommand\CJKfamilydefault{\CJKrmdefault}
"""
latex_logo = '../figures/logo/logo-small-cut.png'
latex_show_urls = 'footnote'





# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'
# pygments_style = "monokai"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_logo = '../figures/logo/logo-small-cut.png'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']