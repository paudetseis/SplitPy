.. _split:

SKS splitting on the fly
========================

Program ``sks_split.py``
------------------------

Calculates single station SKS splitting results by downloading 
waveforms on the fly using
an FDSN client. Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
`StDb` dictionary.

.. _splitusage:

Usage
-----

.. code-block::

    $ sks_split.py -h
    Usage: sks_split.py [options] <station database>

    Script wrapping together the python-based implementation of SplitLab by
    Wustefeld and others. This version requests data on the fly for a given date
    range. Data is requested from the internet using the client services framework
    or from data provided on a local disk. The stations are processed one by one
    with the SKS Splitting parameters measured individually using both the
    Rotation-Correlation (RC) and Silver & Chan (SC) methods.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing Split results.
                            Default behaviour prompts for those that already
                            exist. Selecting overwrite and skip (ie, both flags)
                            negate each other, and both are set to false (every
                            repeat is prompted). [Default False]
      -K, --skip-existing   Skip any event for which existing splitting results
                            are saved to disk. Default behaviour prompts for each
                            event. Selecting skip and overwrite (ie, both flags)
                            negate each other, and both are set to False (every
                            repeat is prompted). [Default False]

      Server Settings:
        Settings associated with which datacenter to log into.

        -S SERVER, --Server=SERVER
                            Specify the server to connect to. Options include:
                            BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU,
                            NCEDC, NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS,
                            USP. [Default IRIS]
        -U USERAUTH, --User-Auth=USERAUTH
                            Enter your IRIS Authentification Username and Password
                            (--User-Auth='username:authpassword') to access and
                            download restricted data. [Default no user and
                            password]

      Local Data Settings:
        Settings associated with defining and using a local data base of pre-
        downloaded day-long SAC files.

        --local-data=LOCALDATA
                            Specify a comma separated list of paths containing
                            day-long sac files of data already downloaded. If data
                            exists for a seismogram is already present on disk, it
                            is selected preferentially over downloading the data
                            using the Client interface
        --no-data-zero      Specify to force missing data to be set as zero,
                            rather than default behaviour which sets to nan.
        --no-local-net      Specify to prevent using the Network code in the
                            search for local data (sometimes for CN stations the
                            dictionary name for a station may disagree with that
                            in the filename. [Default Network used]

      Event Settings:
        Settings associated with refining the events to include in matching
        station pairs

        --start-time=STARTT
                            Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default more recent
                            start date for each station pair]
        --end-time=ENDT     Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station end times [Default older end date
                            for each the pair of stations]
        -R, --reverse-order
                            Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
        --min-mag=MINMAG    Specify the minimum magnitude of event for which to
                            search. [Default 6.0]
        --max-mag=MAXMAG    Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]

      Geometry Settings:
        Settings associatd with the event-station geometries

        --min-dist=MINDIST  Specify the minimum great circle distance (degrees)
                            between the station and event. [Default 85]
        --max-dist=MAXDIST  Specify the maximum great circle distance (degrees)
                            between the station and event. [Default 120]

      Parameter Settings:
        Miscellaneous default values and settings

        --Vp=VP             Specify default P velocity value. [Default 6.0 km/s]
        --SNR=MSNR          Specify the SNR threshold used to determine whether
                            events are processedc. [Default 7.5]
        --window=DTS        Specify time window length before and after the SKS
                            arrival. The total window length is 2*dst. [Default
                            120 s]
        --max-delay=MAXDT   Specify the maximum delay time. [Default 4 s]
        --time-increment=DDT
                            Specify the time increment. [Default 0.1 s]
        --angle-increment=DPHI
                            Specify the angle increment. [Default 1 d]
        --transverse-SNR=SNRTLIM
                            Specify the minimum SNR Threshold for the Transverse
                            component to be considered Non-Null. [Default 1.]

Example
-------

Let's first create a ``StDb`` database and send the prompt to a logfile

.. code-block::

   $ query_fdsn_stdb.py -N TA -C BH? -S EPYK epyk.pkl > logfile

To see the station info for EPYK, use the program ``ls_stdb.py``:

.. code-block::

   $ ls_stdb.py epyk.pkl
    Listing Station Pickle: epyk.pkl
    TA.EPYK
    --------------------------------------------------------------------------
    1) TA.EPYK
         Station: TA EPYK 
          Alternate Networks: None
          Channel: BH ;  Location: --
          Lon, Lat, Elev:  66.37010, -136.71910,   0.717
          StartTime: 2012-10-10 00:00:00
          EndTime:   2599-12-31 23:59:59
          Status:    open
          Polarity: 1
          Azimuth Correction: 0.000000

Then simply run ``sks_split.py`` with ``epyk.pkl`` and you will be able
to perform SKS analysis, refine the analysis window using a GUI and
save the results to disk.

.. code-block::

   $ sks_split.py --keys=TA.EPYK --local-data=/mnt/datadisk/DaySac/ epyk.pkl

This uses all default settings for window lengths, magnitude criteria, etc. 
In this example, data will be used from both IRIS as well as any local data 
on disk (defined with the --local-data flag). If no data exists on disk, then 
the program will search on the specific data sever (through ``obspy`` clients).

Based on the criteria specified (see :ref:`splitusage`), seismograms will be 
processed. The analysis will proceed for an event where the minimum SNR 
threshold is exceeded. Two Figure windows will pop up. ``Figure 1`` shows the three 
rotated component waveforms (LQT), along with lines representing the SKS, SKKS, 
S and ScS arrivals from model ``iasp91``. Red vertical lines denote the analysis 
window. This figure is interactive and the picking red can be refined by clicking
at the two x-positions of the new analysis window.

.. figure:: ../splitpy/examples/figures/Figure_1.png
   :align: center

``Figure 2`` summarizes the results of the splitting calculation. The top left "Q,T" 
frame shows the un-corrected radial (Q) and tangential (T) components within the 
time window. The second row of panels correspond to the 'Rotation-Correlation' 
results, and the third row of panels is for the 'Silver-Chan' results. In each 
case, the first column shows the corrected Q and T fast and slow components, the 
second column the corrected Q and T components, the third column the before and after 
particle motion, and the fourth column the map of the error surfaces. A text box 
prints out the summary of the results, including whether or not the estimate is a 
Null, and the quality of the estimate ('good', 'fair', 'poor').

.. figure:: ../splitpy/examples/figures/Figure_2.png
   :align: center

A message box will pop up asking whether to Re-pick the window in ``Figure 1``. 
This can be done to refine the signal window in which the measurements are made 
in order to improve the measurements.

The terminal will show a summary of the processing, including an examination 
for the Null/Non-Null classification as well as the quality of the estimate.

Once ``No`` is selected for the picking/re-picking of the window, a second box 
will pop up asking whether to keep the estimates. Click ``Yes`` to save the results, 
or ``No`` to discard the measurement.

The results of processing are saved into a ``./RESULTS`` folder in the current 
working directory, in a subfolder named after the station key. In this example, 
``./RESULTS/TA.EPYK``.

Each measurement is stored in a separate pickle (.pkl) file named 
``Split.STN.YYYY.JJJ.HHMMSS.pkl``.

Plotting and subsequent processing of splitting results is carried out using 
:ref:`plot`, where options are present to control selection of nulls 
and quality settings, as well as which methods are used. The path to a station 
``./RESULT`` folder is provided and all Split.*.pkl files are loaded. The final 
average splits are then saved in a text file for future use.