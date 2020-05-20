import numpy as np
from splitpy import arguments
from pkg_resources import resource_filename
from pathlib import Path
from . import get_meta 


dbfile = resource_filename('splitpy',
                           'examples/data/MMPY.pkl')

def test_get_args_calc_auto():
    args = arguments.get_arguments_calc_auto([dbfile])
    return args

def test_get_args_calc_manual():
    args = arguments.get_arguments_calc_manual([dbfile])
    return args

def test_get_args_average():
    args = arguments.get_arguments_average([dbfile])
    return args

