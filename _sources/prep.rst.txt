.. _prep:

SKS data preparation
====================

Program ``sks_prep.py``
-----------------------

Downloads seismic waveforms for SKS splitting analysis, performs
component rotation and determines whether or not to save the
waveforms to disk for subsequent analysis offline (see program
:ref:`offline`). Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
`StDb` dictionary.

Usage
-----

.. code-block::

    $ sks_prep.py -h
    Usage: sks_prep.py [options] <station database>

    Script to download and prepare datasets for SKS splitting processing. This
    script downloads and prepares event and station data, so that splitting can
    then be calculated offline.

    Options:
      -h, --help            show this help message and exit
      --keys=STKEYS         Specify a comma separated list of station keys for
                            which to perform analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.
      --local-data=LOCALDATA
                            Specify a comma separated list of paths containing
                            day-long sac files of data already downloaded. If data
                            exists for a seismogram is already present on disk, it
                            is selected preferentially over downloading the data
                            using the Client interface
      --no-data-zero        Specify to force missing data to be set as zero,
                            rather than default behaviour which sets to nan.
      --no-local-net        Specify to prevent using the Network code in the
                            search for local data (sometimes for CN stations the
                            dictionary name for a station may disagree with that
                            in the filename. [Default Network used]
      -D DATADIR, --data-directory=DATADIR
                            Specify the directory prefix in which the prepared
                            data is stored. [Default 'DATA']. The start and end
                            time and date as well as min and max magnitudes are
                            included in the final folder name.

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

  Event Settings:
    Settings associated with refining the events to include in matching
    station pairs

    --start-time=STARTT
                        Specify a UTCDateTime compatible string representing
                        the start time for the event search. This will
                        override any station start times. [Default more recent
                        start date for each station pair]
    -R, --reverse-order
                        Reverse order of events. Default behaviour starts at
                        oldest event and works towards most recent. Specify
                        reverse order and instead the program will start with
                        the most recent events and work towards older
    --end-time=ENDT     Specify a UTCDateTime compatible string representing
                        the start time for the event search. This will
                        override any station end times [Default older end date
                        for each the pair of stations]
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

Then simply run ``sks_prep.py`` with ``epyk.pkl`` to obtain the waveforms
from an FDSN client and save them to disk. 

.. code-block::

   $ sks_prep.py --keys=TA.EPYK --local-data=/mnt/datadisk/DaySac/ epyk.pkl

This uses all default settings for window lengths, magnitude criteria, etc. 
In this example, data will be used from both IRIS as well as any local data 
on disk (defined with the --local-data flag). If no data exists on disk, then 
the program will search on the specific data sever (through ``obspy`` clients).

Based on the criteria specified (see :ref:`splitusage`), seismograms will be 
processed and saved to disk in a default local folder labelled ``DATA_`` where
the label also contains search criteria: e.g., ``DATA_D85-120_M6.0+_S7.5+``, where
``D85-120`` shows the epicentral distance searched, ``M6.0`` is the minimum 
earthquake magnitude, and ``S7.5`` is the minimum SNR threshold for the radial
component to save the data to disk. Within this folder, the data are stored
by station-key, such as:

.. code-block::

   $ ls -l DATA*
   total 0
   drwxr-xr-x  15 username  staff  480  3 Oct 21:31 TA.EPYK

.. code-block::

   $ ls -l DATA*/TA.EPYK
   total 0
   drwxr-xr-x  12 username  staff  384  3 Oct 21:30 20121020_230032
   drwxr-xr-x  12 username  staff  384  3 Oct 21:30 20121110_145750
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20121210_165309
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20121217_091631
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20121221_222808
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130130_230347
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130208_111213
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130208_152638
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130209_210222
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130216_043737
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130310_225151
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130406_044236
   drwxr-xr-x  12 username  staff  384  3 Oct 21:31 20130409_115250

Once you have downloaded all the waveforms for this station, you can process
them for SKS splitting using the code :ref:`offline`.
