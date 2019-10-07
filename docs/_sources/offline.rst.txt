.. _offline:

SKS splitting offline
=====================

Program ``sks_offline.py``
--------------------------

Calculates single station SKS splitting results from data previously
downloaded using the program :ref:`prep`. Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ sks_offline.py -h
    Usage: sks_offline.py [options] <station database>

    Script to process and calculate the spliting parmaters for a dataset that has
    already been downloaded by sks_prep.py.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]

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
                            search. [Default None, ie no limit]

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

Assume we have followed the steps in :ref:`prep` and saved some data to disk
in a folder called ``DATA_D85-120_M6.0+_S7.5+``. Now we can simply run 
``sks_offline.py`` by specifying the data folder containing the prepared 
waveforms, which execute the SKS splitting analysis and generate
Figures and picking windows similar to :ref:`split` (see the program
for details).

.. code-block::

    $ sks_offline.py ./DATA_D85-120_M6.0+_S7.5+
     
     
    |===============================================|
    |===============================================|
    |                    TA.EPYK                    |
    |===============================================|
    |===============================================|
    |  Working on    13 saved events                |
    |===============================================|
     
    ****************************************************
    * #1 (2/13):  20130409_115250
    *   Origin Time: 2013-04-09 11:52:50
    *   Lat:  28.45; Lon:   51.63
    *   Dep:  12.80; Mag: 6.3
    *     EPYK  -> Ev: 9467.98 km;   85.15 deg; 352.65;   3.35
    * SNR Passed: 20.13 >= 7.5
    * --> Calculating Rotation-Correlation (RC) Splitting
    * --> Calculating Silver-Chan (SC) Splitting

    ...

This uses all default settings for window lengths, magnitude criteria, etc. 
See the program :ref:`split` for details on the interactions with the code
to produce splitting results.