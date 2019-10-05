#!/usr/bin/env python

import setuptools
from numpy.distutils.core import setup, Extension
from os import listdir

scripts = ['Scripts/' + i for i in listdir('Scripts/')]

setup(
    name                = 'splitpy',
    version             = '0.1.0',
    description         = 'Software for teleseismic shear-wave splitting analysis',
    author              = 'Pascal Audet, Andrew Schaeffer',
    maintainer          = 'Pascal Audet, Andrew Schaeffer',
    maintainer_email    = 'pascal.audet@uottawa.ca, andrew.schaeffer@canada.ca',
    classifiers         = [
                            'Development Status :: 3 - Alpha',
                            'License :: OSI Approved :: MIT License',
                            'Programming Language :: Python :: 3.5',
                            'Programming Language :: Python :: 3.6',
                            'Programming Language :: Python :: 3.7'
                            ],
    install_requires    = ['obspy', 'dill', 'PyQt5'],
    python_requires     = '>=3.5',
    packages            = ['splitpy'],
    package_dir         = {'SplitPy': 'SplitPy'},
    scripts             = scripts,
    url                 = 'https://github.com/paudetseis/SplitPy'
    )
