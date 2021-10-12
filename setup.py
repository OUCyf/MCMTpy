#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCMTpy is Created on Sun Apr 18 17:24:43 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    Nowadays, it’s usually do the research based on a python’s workflow, with 
    the help of packages such as obspy, numpy, pyasdf and mpi4py. MCMTpy is 
    Python Package for Simultaneous Inversion of Source Location, Focal 
    Mechanism, and Rupture Directivity.


Modify history:
    1) Apr 18 17:24:43 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""

import os
import time
from setuptools import find_packages, setup


############################################
project_root = os.path.dirname(os.path.realpath(__file__))                      # Gets the absolute path to the currently executing script
DOCSTRING = __doc__.strip().split("\n")                                         # The above comment document



def make_info_module(project_root,version):
    '''Put version and revision information into file MCMTpy/info.py.'''
    datestr = time.strftime('%Y-%m-%d_%H:%M:%S')
    s = '''# This module is automatically created from setup.py
project_root = %s
version = %s
installed_date = %s
''' % tuple([repr(x) for x in (
        project_root, version,  datestr)])

    try:
        f = open(os.path.join('MCMTpy', 'info.py'), 'w')
        f.write(s)
        f.close()
    except Exception:
        pass



def readme(project_root):
    """
    read 'README.rst'
    """
    README_file = os.path.join(project_root, 'README.rst')
    with open(README_file) as f:
        return f.read()



def get_requires(project_root):
    """
    read requirements.txt
    """
    requirements_file = os.path.join(project_root, 'requirements.txt')
    with open(requirements_file) as f:
        install_reqs = f.read().splitlines()
        return install_reqs



def get_version(project_root):
    """
    get __version__.py version number
    """
    version_file = os.path.join(project_root,'MCMTpy', '__version__.py')
    with open(version_file,'r') as f:
        for line in f:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")



def get_package_data(project_root):
    """
    Read package_data info
    Returns a list of all files needed for the installation relative to the
    'MCMTpy' subfolder.
    """
    filenames = []
    # Recursively include all files in these folders:
    folders = [
        os.path.join(project_root, "data"),
    ]
    for folder in folders:
        for directory, _, files in os.walk(folder):
            for filename in files:
                # Exclude hidden files.
                if filename.startswith("."):
                    continue
                filenames.append(
                    os.path.relpath(os.path.join(
                        directory, filename), project_root)
                )
    return filenames



setup_config = dict(
    name="MCMTpy",
    version=get_version(project_root),
    description=DOCSTRING[0],
    long_description=readme(project_root),
    author="Fu Yin",
    author_email="yinfu@mail.ustc.edu.cn",
    url="https://www.google.com.hk/",                                           # change later
    license="MIT",
    packages=find_packages(),
    package_data={"MCMTpy": get_package_data(project_root)},
    entry_points={
    'console_scripts':
        ['MCMTpy = MCMTpy.apps.MCMTpy:main']
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=get_requires(project_root),
)





############################################
if __name__ == "__main__":
    version=get_version(project_root)
    make_info_module(project_root,version)
    
    setup(
        **setup_config
    )



