import os.path
from os import listdir
import re
from setuptools import setup
# from numpy.distutils.core import setup
from pathlib import Path


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


scripts = [str(x) for x in Path('Scripts').iterdir() if x.is_file()]

setup(
    name='splitpy',
    version=find_version('splitpy', '__init__.py'),
    description='Software for teleseismic shear-wave splitting analysis',
    author='Pascal Audet, Andrew Schaeffer',
    author_email='pascal.audet@uottawa.ca',
    maintainer='Pascal Audet, Andrew Schaeffer',
    maintainer_email='pascal.audet@uottawa.ca, andrew.schaeffer@canada.ca',
    url='https://github.com/paudetseis/SplitPy',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'],
    install_requires=['stdb', 'obspy', 'PyQt5'],
    python_requires='>=3.6',
    packages=['splitpy'],
    scripts=scripts)
