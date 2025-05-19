.. figure:: ../splitpy/examples/figures/SplitPy_logo.png
   :align: center

Scripts
=======

There are several Python scripts that accompany :mod:`~splitpy`, which can be used
in bash scripts to automate data processing. These include scripts to download 
three-component seismogram data and automatically calculate splitting parameters, 
manually refine the analysis window for more precise estimates, and final averaging
of splitting parameters at single stations. All of them use 
a station database provided as a :class:`~stdb.StDb` dictionary. 

.. _splitauto:

``split_calc_auto``
+++++++++++++++++++

Description
-----------

This script can be used in one of two ways: 1) to download/collect available
teleseismic shear-wave data for later processing; and 2) to further perform
an automated processing for the shear-wave splitting estimates using default
parameters. Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
:class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ split_calc_auto -h

    ###################################################################
    #            _ _ _                _                     _         #
    #  ___ _ __ | (_) |_     ___ __ _| | ___     __ _ _   _| |_ ___   #
    # / __| '_ \| | | __|   / __/ _` | |/ __|   / _` | | | | __/ _ \  #
    # \__ \ |_) | | | |_   | (_| (_| | | (__   | (_| | |_| | || (_) | #
    # |___/ .__/|_|_|\__|___\___\__,_|_|\___|___\__,_|\__,_|\__\___/  #
    #     |_|          |_____|             |_____|                    #
    #                                                                 #
    ###################################################################

    usage: split_calc_auto [arguments] <station database>

    Script wrapping together the python-based implementation of SplitLab by
    Wustefeld and others. This version requests data on the fly for a given date
    range. Data is requested from the internet using the client services framework
    or from data provided on a local disk. The stations are processed one by one
    with the SKS Splitting parameters measured individually using both the
    Rotation-Correlation (RC) and Silver & Chan (SC) methods.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for
                            which to perform the analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -V, --verbose         Specify to increase verbosity.
      -O, --overwrite       Force the overwriting of pre-existing Split results.
                            Default behaviour prompts for those that already
                            exist. Selecting overwrite and skip (ie, both flags)
                            negate each other, and both are set to false (every
                            repeat is prompted). [Default False]
      --zcomp ZCOMP         Specify the Vertical Component Channel Identifier. [Default Z].
      --calc                Analyze data for shear-wave splitting. [Default saves
                            data to folders for subsequent analysis]
      --plot-diagnostic
                            Plot diagnostic window at end of process. [Default
                            False]
      --recalc              Re-calculate estimates and overwrite existing splitting results
                            without re-downloading data. [Default False]

    Server Settings:
      Settings associated with which datacenter to log into.

      --server SERVER       Base URL of FDSN web service compatible server (e.g.
                            “http://service.iris.edu”) or key string for recognized server (one
                            of 'AUSPASS', 'BGR', 'EARTHSCOPE', 'EIDA', 'EMSC', 'ETH', 'GEOFON',
                            'GEONET', 'GFZ', 'ICGC', 'IESDMC', 'INGV', 'IPGP', 'IRIS', 'IRISPH5',
                            'ISC', 'KNMI', 'KOERI', 'LMU', 'NCEDC', 'NIEP', 'NOA', 'NRCAN',
                            'ODC', 'ORFEUS', 'RASPISHAKE', 'RESIF', 'RESIFPH5', 'SCEDC',
                            'TEXNET', 'UIB-NORSAR', 'USGS', 'USP'). [Default 'IRIS']
      --user-auth USERAUTH  Authentification Username and Password for the waveform server
                            (--user-auth='username:authpassword') to access and download
                            restricted data. [Default no user and password]
      --eida-token TOKENFILE
                            Token for EIDA authentication mechanism, see http://geofon.gfz-
                            potsdam.de/waveform/archive/auth/index.php. If a token is provided,
                            argument --user-auth will be ignored. This mechanism is only
                            available on select EIDA nodes. The token can be provided in form of
                            the PGP message as a string, or the filename of a local file with the
                            PGP message in it. [Default None]

    Local Data Settings:
      Settings associated with defining and using a local data base of pre-downloaded day-long
      SAC or MSEED files.

      --local-data LOCALDATA
                            Specify absolute path to a SeisComP Data Structure (SDS) archive
                            containing day-long SAC or MSEED files(e.g., --local-
                            data=/Home/username/Data/SDS). See
                            https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html
                            for details on the SDS format. If this option is used, it takes
                            precedence over the --server settings.
      --dtype DTYPE         Specify the data archive file type, either SAC or MSEED. Note the
                            default behaviour is to search for SAC files. Local archive files
                            must have extensions of '.SAC' or '.MSEED'. These are case dependent,
                            so specify the correct case here.

    Parameter Settings:
      Miscellaneous default values and settings

      --sampling-rate NEW_SAMPLING_RATE
                            Specify new sampling rate in Hz. [Default 10.]
      --min-snr MSNR        Minimum SNR value calculated on the radial (Q)
                            component to proceed with analysis (dB). [Default 5.]
      --window DTS          Specify time window length before and after the SKS
                            arrival. The total window length is 2*dst (sec).
                            [Default 120]
      --max-delay MAXDT     Specify the maximum delay time in search (sec).
                            [Default 4]
      --dt-delay DDT        Specify the time delay increment in search (sec).
                            [Default 0.1]
      --dphi DPHI           Specify the fast angle increment in search (degree).
                            [Default 1.]
      --snrT SNRTLIM        Specify the minimum SNR Threshold for the Transverse
                            component to be considered Non-Null. [Default 1.]
      --fmin FMIN           Specify the minimum frequency corner for SNR filter
                            (Hz). [Default 0.02]
      --fmax FMAX           Specify the maximum frequency corner for SNR filter
                            (Hz). [Default 0.5]

    Event Settings:
      Settings associated with refining the events to include in matching
      station pairs

      --start STARTT        Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default start date
                            of each station]
      --end ENDT            Specify a UTCDateTime compatible string representing
                            the end time for the event search. This will override
                            any station end times [Default end date of each
                            station]
      --reverse             Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older
      --min-mag MINMAG      Specify the minimum magnitude of event for which to
                            search. [Default 6.0]
      --max-mag MAXMAG      Specify the maximum magnitude of event for which to
                            search. [Default None, i.e. no limit]

    Geometry Settings:
      Settings associatd with the event-station geometries

      --min-dist MINDIST    Specify the minimum great circle distance (degrees)
                            between the station and event. [Default 85]
      --max-dist MAXDIST    Specify the maximum great circle distance (degrees)
                            between the station and event. [Default 120]
      --phase PHASE         Specify the phase name to use. Be careful with the
                            distance. setting. Options are 'SKS' or 'SKKS'.
                            [Default 'SKS']

.. _splitmanual:

``split_calc_manual``
+++++++++++++++++++++

Description
-----------

This script is used if the user desires manual re-picking of the analysis window
for refined estimates. Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
:class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ split_calc_manual -h

    ###################################################################################
    #            _ _ _                _                                            _  #
    #  ___ _ __ | (_) |_     ___ __ _| | ___     _ __ ___   __ _ _ __  _   _  __ _| | #
    # / __| '_ \| | | __|   / __/ _` | |/ __|   | '_ ` _ \ / _` | '_ \| | | |/ _` | | #
    # \__ \ |_) | | | |_   | (_| (_| | | (__    | | | | | | (_| | | | | |_| | (_| | | #
    # |___/ .__/|_|_|\__|___\___\__,_|_|\___|___|_| |_| |_|\__,_|_| |_|\__,_|\__,_|_| #
    #     |_|          |_____|             |_____|                                    #
    #                                                                                 #
    ###################################################################################

    usage: split_calc_manual [arguments] <station database>

    Script to process and calculate the spliting parameters for a dataset that has
    already been downloaded by split_calc_auto.

    positional arguments:
      indb                  Station Database to process from.

    optional arguments:
      -h, --help            show this help message and exit
      --keys STKEYS         Specify a comma separated list of station keys for
                            which to perform analysis. These must be contained
                            within the station database. Partial keys will be used
                            to match against those in the dictionary. For
                            instance, providing IU will match with all stations in
                            the IU network [Default processes all stations in the
                            database]
      -v, -V, --verbose     Specify to increase verbosity.

    Parameter Settings:
      Miscellaneous default values and settings

      --window DTS          Specify time window length before and after the SKS
                            arrival. The total window length is 2*dst (sec).
                            [Default 120]
      --max-delay MAXDT     Specify the maximum delay time. [Default 4 s]
      --time-increment DDT  Specify the time increment. [Default 0.1 s]
      --angle-increment DPHI
                            Specify the angle increment. [Default 1 d]
      --transverse-SNR SNRTLIM
                            Specify the minimum SNR Threshold for the Transverse
                            component to be considered Non-Null. [Default 1.]

    Event Settings:
      Settings associated with refining the events to include in matching
      station pairs

      --start STARTT        Specify a UTCDateTime compatible string representing
                            the start time for the event search. This will
                            override any station start times. [Default more recent
                            start date for each station pair]
      --end ENDT            Specify a UTCDateTime compatible string representing
                            the end time for the event search. This will override
                            any station end times [Default older end date for each
                            the pair of stations]
      --reverse-order       Reverse order of events. Default behaviour starts at
                            oldest event and works towards most recent. Specify
                            reverse order and instead the program will start with
                            the most recent events and work towards older

.. _splitaverage:

``split_average``
+++++++++++++++++

Description
-----------

This script is used for producing station average shear-wave splitting estimates obtained 
from either the automated or manual mode. 
Station selection is specified by a network and 
station code. The data base is provided in a pickled file as a 
:class:`~stdb.StDb` dictionary.

Usage
-----

.. code-block::

    $ split_average -h

    ###############################################################
    #            _ _ _                                            #
    #  ___ _ __ | (_) |_     __ ___   _____ _ __ __ _  __ _  ___  #
    # / __| '_ \| | | __|   / _` \ \ / / _ \ '__/ _` |/ _` |/ _ \ #
    # \__ \ |_) | | | |_   | (_| |\ V /  __/ | | (_| | (_| |  __/ #
    # |___/ .__/|_|_|\__|___\__,_| \_/ \___|_|  \__,_|\__, |\___| #
    #     |_|          |_____|                        |___/       #
    #                                                             #
    ###############################################################

    usage: split_average [arguments] <station database>

    Script to plot the average splitting results for a given station. Loads the available .pkl
    files in the specified Station Directory.

    positional arguments:
      indb           Station Database to process from.

    options:
      -h, --help     show this help message and exit
      --keys STKEYS  Specify a comma separated list of station keys for which to perform
                     analysis. These must be contained within the station database. Partial keys
                     will be used to match against those in the dictionary. For instance,
                     providing IU will match with all stations in the IU network [Default
                     processes all stations in the database]
      -V, --verbose  Specify to increase verbosity.
      --show-fig     Specify show plots during processing - they are still saved to disk.
                     [Default only saves]
      --auto         Specify to use automatically processed split results. [Default uses refined
                     ('manual') split results]

    Null Selection Settings:
      Settings associated with selecting which Null or Non-Null data is included

      --nulls        Specify this flag to include Null Values in the average. [Default Non-Nulls
                     only]
      --no-nons      Specify this flag to exclude Non-Nulls from the average [Default False]

    Quality Selection Settings:
      Settings associated with selecting the qualities to include in the selection.

      --no-good      Specify to exclude 'Good' measurements from the average. [Default Good +
                     Fair]
      --no-fair      Specify to exclude 'Fair' measurements from the average [Default Good +
                     Fair]
      --poor         Specify to include 'Poor' measurements in the average [Default No Poors]
