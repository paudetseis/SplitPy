import numpy as np
from splitpy import arguments
from pkg_resources import resource_filename
from pathlib import Path
from . import get_meta 


dbfile = resource_filename('splitpy',
                           'examples/data/MMPY.pkl')

def test_get_args():
    args = arguments.get_arguments([dbfile])
    return args

def test_get_args_prep():
    args = arguments.get_arguments_prep([dbfile])
    return args

def test_get_args_offline():
    args = arguments.get_arguments_offline([dbfile])
    return args

def test_get_args_plot():
    args = arguments.get_arguments_plot([dbfile])
    return args

