def test_stdb_import():
    import stdb

def test_obspy_import():
    import obspy

def test_splitpy_modules():
    import splitpy
    from splitpy import io, calc, classes, arguments, gui
    from splitpy.classes import Meta, Data, Result, Split
    from splitpy import Pick, Keep, Save, Repeat
    from splitpy import PickPlot, DiagPlot
    import matplotlib
    matplotlib.use('Qt5Agg')
