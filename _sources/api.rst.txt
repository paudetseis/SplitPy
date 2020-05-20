.. figure:: ../splitpy/examples/figures/SplitPy_logo.png
   :align: center

Classes
=======

:mod:`~splitpy` defines the following base classes:

- :class:`~splitpy.classes.Split`
- :class:`~splitpy.classes.PickPlot`
- :class:`~splitpy.classes.DiagPlot`

The class :class:`~splitpy.classes.Split` contains attributes
and methods for the analysis of teleseismic shear-wave splitting 
from three-component seismograms. 

The class :class:`~splitpy.classes.PickPlot` defines figure handles 
for a picking window showing the seismograms and the predicted teleseismic
shear-wave phase arrivals. This figure is interactive and new picks can
be generated to refine the analysis.

The class :class:`~splitpy.classes.DiagPlot` defines figure handles
for a diagnostic figure showing a summary of the splitting results. It can
be called after each application of the `split.analyze` method to show 
the summary of the analysis as a figure. This figure can also be saved as
a .png file.

Split
-----

.. autoclass:: splitpy.classes.Split
   :members:

Meta
----

.. autoclass:: splitpy.classes.Meta
   :members:

Result
------

.. autoclass:: splitpy.classes.Result
   :members:

PickPlot
--------

.. autoclass:: splitpy.classes.PickPlot
   :members:

DiagPlot
--------

.. autoclass:: splitpy.classes.DiagPlot
   :members:

GUI classes
===========

.. automodule:: splitpy.gui
   :members:

Modules
=======

calc
----

.. automodule:: splitpy.calc
   :members:

utils
-----

.. automodule:: splitpy.utils
   :members:
