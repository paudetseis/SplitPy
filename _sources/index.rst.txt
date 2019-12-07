.. figure:: ../splitpy/examples/figures/SplitPy_logo.png
   :align: center

SplitPy Project documentation
=============================

SplitPy is a teleseismic shear-wave (SKS) Splitting Toolbox based on the 
Matlab Tool SplitLab, developed by 
`Wustefeld et al (2008) <https://doi.org/doi:10.1016/j.cageo.2007.08.002>`_.

The software is bundled with executable ``Python`` scripts that enable 
SKS splitting analysis in two different environments:

1) On-the-fly processing where waveforms are downloaded and
   analyzed for SKS splitting (:ref:`split`).
2) Offline processing where waveforms are first downloaded
   and saved to disk (:ref:`prep`), and later processed offline (:ref:`offline`).

Once the splitting analysis has been done, the results can be 
aggregated to produce averaged splitting parameters and plotted 
(:ref:`plot`). The tutorials below show an example using data
from `one station <http://ds.iris.edu/mda/TA/EPYK/>`_ of the 
USArray Transportable Array network.

.. note::

    The software ``SplitPy`` uses a `StDb <https://github.com/schaefferaj/StDb>`_
    database for processing. Check out `StDb <https://github.com/schaefferaj/StDb>`_ 
    for more details. An example is shown in the program :ref:`split`.

.. image:: https://zenodo.org/badge/211722700.svg
   :target: https://zenodo.org/badge/latestdoi/211722700

Quick links
"""""""""""

* `StDb Git repository <https://github.com/schaefferaj/StDb>`_
* `SplitPy Git repository <https://github.com/paudetseis/SplitPy>`_
* `SplitLab <http://splitting.gm.univ-montp2.fr>`_

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   api

.. toctree::
   :maxdepth: 1
   :caption: API Documentation

   classes
   calc
   utils
   gui

.. toctree::
   :maxdepth: 1
   :caption: Scripts & Tutorials

   split
   prep
   offline
   plot

