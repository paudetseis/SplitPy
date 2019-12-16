import os.path
from os import listdir
import re
from numpy.distutils.core import setup


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


scripts = ['Scripts/' + i for i in listdir('Scripts/')]

setup(
    name='splitpy',
    version=find_version('splitpy', '__init__.py'),
    description='Software for teleseismic shear-wave splitting analysis',
    author='Pascal Audet, Andrew Schaeffer',
    maintainer='Pascal Audet, Andrew Schaeffer',
    maintainer_email='pascal.audet@uottawa.ca, andrew.schaeffer@canada.ca',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    install_requires=['stdb', 'obspy', 'dill', 'PyQt5'],
    python_requires='>=3.6',
    packages=['splitpy'],
    scripts=scripts,
    url='https://github.com/paudetseis/SplitPy')
