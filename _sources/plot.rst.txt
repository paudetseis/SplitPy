.. _plot:

Plotting Average results
========================

Program ``sks_plot_results.py``
-------------------------------

Plots the results for a given station based on ``.pkl`` files present in the
running directory.

Usage
-----

.. code-block::

   $ sks_plot_results.py -h
    Usage: sks_plot_results.py [options] <Folder> [Folder 2]

    Script to plot the average splitting results for a given station. Loads the
    available pkl files in the specified Station Directory.

    Options:
      -h, --help            show this help message and exit
      --no-figure           Specify to prevent plots from opening during
                            processing; they are still saved to disk. [Default
                            plots open and save]

      Null Selection Settings:
        Settings associated with selecting which Null or Non-Null data is
        included

        --nulls, --Nulls    Specify this flag to include Null Values in the
                            average. [Default Non-Nulls only]
        --no-nons, --No-Nons
                            Specify this flag to exclude Non-Nulls from the
                            average [Default False]

      Quality Selection Settings:
        Settings associated with selecting the qualities to include in the
        selection.

        --No-Good, --no-good
                            Specify to exclude 'Good' measurements from the
                            average. [Default Good + Fair]
        --No-Fair, --no-fair
                            Specify to exclude 'Fair' measurements from the
                            average [Default Good + Fair]
        --Poor, --poor      Specify to include 'Poor' measurements in the average
                            [Default No Poors]

      Split Type Settings:
        Settings to Select which Split types are included in the selection.

        --RC-Only, --rc-only, --RC-only
                            Specify to only include RC splits in the average.
                            [Default RC + SC]
        --SC-Only, --sc-only, --SC-only
                            Specify to only include SC splits in the average.
                            [Default RC + SC]

Example
-------

In this example we rapidly analyzed SKS data for station EPYK for approximately
1 year of data. To view the aggregate results and obtain vector averages
of the fast axis direction, simply execute the program:

.. code-block::

   $ sks_plot_results.py ./RESULTS/TA.EPYK

   ---------------------------
    Selection Criteria 
     Null Value: 
        Non Nulls: True
        Nulls:     False
     Quality Value: 
        Goods:  True
        Fairs:  True
        Poors:  False
    ---------------------------
    Working on Station: ./RESULTS/TA.EPYK
      Processing 10 Events...
      1) ./RESULTS/TA.EPYK/Split.EPYK.2012.294.230032.pkl
          Good Null -> Skipped
      2) ./RESULTS/TA.EPYK/Split.EPYK.2012.345.165309.pkl
          Poor Null -> Skipped
      3) ./RESULTS/TA.EPYK/Split.EPYK.2013.096.044236.pkl
          Good Non-Null -> Retained
      4) ./RESULTS/TA.EPYK/Split.EPYK.2013.106.225527.pkl
          Good Non-Null -> Retained
      5) ./RESULTS/TA.EPYK/Split.EPYK.2013.116.065328.pkl
          Good Non-Null -> Retained
      6) ./RESULTS/TA.EPYK/Split.EPYK.2013.185.171557.pkl
          Good Null -> Skipped
      7) ./RESULTS/TA.EPYK/Split.EPYK.2013.187.050506.pkl
          Good Null -> Skipped
      8) ./RESULTS/TA.EPYK/Split.EPYK.2013.213.200142.pkl
          Good Non-Null -> Retained
      9) ./RESULTS/TA.EPYK/Split.EPYK.2013.284.212459.pkl
          Good Non-Null -> Retained
      10) ./RESULTS/TA.EPYK/Split.EPYK.2013.327.074832.pkl
          Good Non-Null -> Retained

    *** Station Average from 6 measurements ***
       ./RESULTS/TA.EPYK
       Loc: -136.7191, 66.3701
       PHI:  60.655 d +- 28.987
       DT:    0.907 s +- 0.664
       Saved to: ./RESULTS/TA.EPYK_RC-SC_Nons_G-F_results

This produces a summary to prompt, and generates a Figure as in the example below

.. figure:: ../splitpy/examples/figures/Figure_3.png
   :align: center

This Figure shows the fast axis directions and delay times 
between the fast and slow components as a function 
of the event back-azimuths (panels **A** and **B**). 
The program calculates a vector average of the spliting 
results (panel **C**) and estimates uncertainties, shown in 
the prompt (above).

The results and Figure are saved to disk in the files with prefix
``./RESULTS/TA.EPYK_RC-SC_Nons_G-F_results`` and suffix ``.dat`` 
or ``.png``.

