.. figure:: ../splitpy/examples/figures/SplitPy_logo.png
   :align: center

Documentation
=============

SplitPy is a teleseismic shear-wave (SKS) Splitting Toolbox based on the 
Matlab Tool SplitLab, developed by 
`Wustefeld et al (2008) <https://doi.org/doi:10.1016/j.cageo.2007.08.002>`_.

The software is bundled with executable ``Python`` scripts that enable 
teleseismic shear-wave splitting analysis in two different modes:

1) Automated processing where waveforms are downloaded and
   analyzed using default parameters (:ref:`splitauto`).
2) Manual processing where picking window can be refined (:ref:`splitmanual`).

Once the splitting analysis has been done, the results can be 
aggregated to produce averaged splitting parameters and plotted 
(:ref:`splitaverage`). The tutorials below show an example using data
from `one station <http://ds.iris.edu/mda/NY/TGTN/>`_ of the 
Yukon-Northwest Seismograph Network.

.. note::

    The software ``SplitPy`` uses a `StDb <https://github.com/schaefferaj/StDb>`_
    database for processing. Check out `StDb <https://github.com/schaefferaj/StDb>`_ 
    for more details. An example is shown in the tutorials.

.. image:: https://zenodo.org/badge/211722700.svg
    :target: https://zenodo.org/badge/latestdoi/211722700
.. image:: https://travis-ci.org/paudetseis/SplitPy.svg?branch=master
    :target: https://travis-ci.org/paudetseis/

.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   splitpy

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   api

.. toctree::
   :maxdepth: 2
   :caption: Scripts & Tutorials

   scripts
   tutorials

