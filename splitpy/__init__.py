'''
Module: SplitPy

Python Module for computing SKS splitting parameters using both the RC and SC methods.

Contains SubModules:
  classes.py -> defines main splitting classes
  gui.py    -> gui preparation and interface operation
  utils.py  -> shared utilities
'''

#__version__=0.1
__author__ = "Pascal Audet"

from . import conf, io, calc, utils
from .classes import Split, SeisPlot, DiagPlot
from .gui import Pick, Keep, Save, Repeat