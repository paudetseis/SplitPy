#!/usr/bin/env python

import setuptools
from numpy.distutils.core import setup, Extension
from os import listdir

version = open('version.txt').read().split()[0]
scripts = ['Scripts/' + i for i in listdir('Scripts/')]

setup(
    name =              'splitpy',
    version =           version,
    description =       'Python Module for computing SKS splitting orientations',
    author =            'Pascal Audet, Andrew Schaeffer',
    maintainer =        'Andrew Schaeffer',
    maintainer_email =  'pascal.audet@uottawa.ca',
    packages =          ['splitpy'],
    package_dir =       {'SplitPy': 'SplitPy'},
    scripts =           scripts,
    url =               'http://gitlab.com/uottawa-geophysics')
