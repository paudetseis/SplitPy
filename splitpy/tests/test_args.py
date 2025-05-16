import numpy as np
from pkg_resources import resource_filename
from splitpy.scripts import split_calc_auto as sca
from splitpy.scripts import split_calc_manual as scm
from splitpy.scripts import split_average as sa
from pathlib import Path
from . import get_meta 


dbfile = resource_filename('splitpy',
                           'examples/data/MMPY.pkl')

def test_get_args_calc_auto():
    args = sca.get_arguments_calc_auto([dbfile])
    return args

def test_get_args_calc_manual():
    args = scm.get_arguments_calc_manual([dbfile])
    return args

def test_get_args_average():
    args = sa.get_arguments_average([dbfile])
    return args

