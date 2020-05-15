# Copyright 2019 Pascal Audet & Andrew Schaeffer
#
# This file is part of SplitPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""

Module containing the main utility functions used in the `SplitPy` scripts
that accompany this package.

"""

# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from os.path import exists as exist
from obspy import UTCDateTime
from numpy import nan


def get_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for data processing on-the-fly (requires web connection)

    """

    parser = ArgumentParser(
        usage="Usage: %prog [arguments] <station database>",
        description="Script wrapping "
        "together the python-based implementation of SplitLab by " +
        "Wustefeld and others. This version " +
        "requests data on the fly for a given date range. Data is " +
        "requested from the internet using " +
        "the client services framework or from data provided on a " +
        "local disk. The stations are processed " +
        "one by one with the SKS Splitting parameters measured " +
        "individually using both the " +
        "Rotation-Correlation (RC) and Silver & Chan (SC) methods.")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys for " +
        "which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match with " +
        "all stations in the IU network [Default processes " +
        "all stations in the database]")
    parser.add_argument(
        "-v", "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing Split results. " +
        "Default behaviour prompts for those that " +
        "already exist. Selecting overwrite and skip (ie, both flags) " +
        "negate each other, and both are set to " +
        "false (every repeat is prompted). [Default False]")
    parser.add_argument(
        "-K", "--skip-existing",
        action="store_true",
        dest="skip",
        default=False,
        help="Skip any event for which existing splitting results are " +
        "saved to disk. Default behaviour prompts for " +
        "each event. Selecting skip and overwrite (ie, both flags) " +
        "negate each other, and both are set to " +
        "False (every repeat is prompted). [Default False]")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which " +
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")

    # Database Settings
    DataGroup = parser.add_argument_group(
        title="Local Data Settings",
        description="Settings associated with defining and using a " +
        "local data base of pre-downloaded day-long SAC files.")
    DataGroup.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify a comma separated list of paths containing " +
        "day-long sac files of data already downloaded. " +
        "If data exists for a seismogram is already present on " +
        "disk, it is selected preferentially over downloading " +
        "the data using the Client interface")
    DataGroup.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, rather " +
        "than default behaviour which sets to nan.")
    DataGroup.add_argument(
        "--no-local-net",
        action="store_false",
        dest="useNet",
        default=True,
        help="Specify to prevent using the Network code in the " +
        "search for local data (sometimes for CN stations " +
        "the dictionary name for a station may disagree with that " +
        "in the filename. [Default Network used]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--Vp",
        action="store",
        type=float,
        dest="vp",
        default=6.,
        help="Specify default P velocity value. [Default 6.0 km/s]")
    ConstGroup.add_argument(
        "--SNR",
        action="store",
        type=float,
        dest="msnr",
        default=7.5,
        help="Specify the SNR threshold used to determine whether " +
        "events are processedc. [Default 7.5]")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="dts",
        default=120.,
        help="Specify time window length before and after the SKS "
        "arrival. The total window length is 2*dst. [Default 120 s]")
    ConstGroup.add_argument(
        "--max-delay",
        action="store",
        type=float,
        dest="maxdt",
        default=4.,
        help="Specify the maximum delay time. [Default 4 s]")
    ConstGroup.add_argument(
        "--time-increment",
        action="store",
        type=float,
        dest="ddt",
        default=0.1,
        help="Specify the time increment. [Default 0.1 s]")
    ConstGroup.add_argument(
        "--angle-increment",
        action="store",
        type=float,
        dest="dphi",
        default=1.,
        help="Specify the angle increment. [Default 1 d]")
    ConstGroup.add_argument(
        "--transverse-SNR",
        action="store",
        type=float,
        dest="snrTlim",
        default=1.,
        help="Specify the minimum SNR Threshold for the Transverse " +
        "component to be considered Non-Null. [Default 1.]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with refining "
        "the events to include in matching station pairs")
    EventGroup.add_argument(
        "--start-time",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station start times. [Default start date of each station]")
    EventGroup.add_argument(
        "--end-time",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the end time for the event search. This will override any " +
        "station end times [Default end date of each station]")
    EventGroup.add_argument(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour starts at " +
        "oldest event and works towards most recent. " +
        "Specify reverse order and instead the program will start " +
        "with the most recent events and work towards older")
    EventGroup.add_argument(
        "--min-mag",
        action="store",
        type=float,
        dest="minmag",
        default=6.0,
        help="Specify the minimum magnitude of event for which to " +
        "search. [Default 6.0]")
    EventGroup.add_argument(
        "--max-mag",
        action="store",
        type=float,
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event for which to " +
        "search. [Default None, i.e. no limit]")

    # Geometry Settings
    GeomGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with the "
        "event-station geometries")
    GeomGroup.add_argument(
        "--min-dist",
        action="store",
        type=float,
        dest="mindist",
        default=85.,
        help="Specify the minimum great circle distance (degrees) " +
        "between the station and event. [Default 85]")
    GeomGroup.add_argument(
        "--max-dist",
        action="store",
        type=float,
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance (degrees) " +
        "between the station and event. [Default 120]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for " +
                "User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # Check existing file behaviour
    if args.skip and args.ovr:
        args.skip = False
        args.ovr = False

    # Parse Local Data directories
    if args.localdata is not None:
        args.localdata = args.localdata.split(',')
    else:
        args.localdata = []

    # Check NoData Value
    if args.ndval:
        args.ndval = 0.0
    else:
        args.ndval = nan

    return args


def get_arguments_prep(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for preparation of SKS data for offline processing

    """

    parser = ArgumentParser(
        usage="Usage: %prog [arguments] <station database>",
        description="Script to " +
        "download and prepare datasets for SKS splitting processing. " +
        "This script downloads and prepares event and station data, " +
        "so that splitting can then be calculated offline.")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys for which " +
        "to perform analysis. These must be contained within the " +
        "station database. Partial keys will be used to match against " +
        "those in the dictionary. For instance, providing IU will match " +
        "with all stations in the IU network [Default " +
        "processes all stations in the database]")
    parser.add_argument(
        "-v", "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify a comma separated list of paths containing " +
        "day-long sac files of data already downloaded. If data exists " +
        "for a seismogram is already present on disk, it is selected " +
        "preferentially over downloading the data using the Client interface")
    parser.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, rather " +
        "than default behaviour which sets to nan.")
    parser.add_argument(
        "--no-local-net",
        action="store_false",
        dest="useNet",
        default=True,
        help="Specify to prevent using the Network code in the search " +
        "for local data (sometimes for CN stations the dictionary name " +
        "for a station may disagree with that in the filename. " +
        "[Default Network used]")
    parser.add_argument(
        "-D", "--data-directory",
        action="store",
        type=str,
        dest="datadir",
        default="DATA",
        help="Specify the directory prefix in which the prepared data " +
        "is stored. [Default 'DATA']. The start and end time and date " +
        "as well as min and max magnitudes are included in the final " +
        "folder name.")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--Vp",
        action="store",
        type=float,
        dest="vp",
        default=6.,
        help="Specify default P velocity value. [Default 6.0 km/s]")
    ConstGroup.add_argument(
        "--SNR",
        action="store",
        type=float,
        dest="msnr",
        default=7.5,
        help="Specify the SNR threshold used to determine whether " +
        "events are processedc. [Default 7.5]")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="dts",
        default=120.,
        help="Specify time window length before and after the " +
        "SKS arrival. The total window length is "
        "2*dst. [Default 120 s]")
    ConstGroup.add_argument(
        "--max-delay",
        action="store",
        type=float,
        dest="maxdt",
        default=4.,
        help="Specify the maximum delay time. [Default 4 s]")
    ConstGroup.add_argument(
        "--time-increment",
        action="store",
        type=float,
        dest="ddt",
        default=0.1,
        help="Specify the time increment. [Default 0.1 s]")
    ConstGroup.add_argument(
        "--angle-increment",
        action="store",
        type=float,
        dest="dphi",
        default=1.,
        help="Specify the angle increment. [Default 1 d]")
    ConstGroup.add_argument(
        "--transverse-SNR",
        action="store",
        type=float,
        dest="snrTlim",
        default=1.,
        help="Specify the minimum SNR Threshold for the Transverse " +
        "component to be considered Non-Null. [Default 1.]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with "
        "refining the events to include in matching station pairs")
    EventGroup.add_argument(
        "--start-time",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing the " +
        "start time for the event search. This will override any station " +
        "start times. [Default more recent start date for each station pair]")
    EventGroup.add_argument(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour starts at " +
        "oldest event and works towards most " +
        "recent. Specify reverse order and instead the program will " +
        "start with the most recent events and "
        "work towards older")
    EventGroup.add_argument(
        "--end-time",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the end time for the event search. " +
        "This will override any station end times [Default older end " +
        "date for each the pair of stations]")
    EventGroup.add_argument(
        "--min-mag",
        action="store",
        type=float,
        dest="minmag",
        default=6.0,
        help="Specify the minimum magnitude of event for which to " +
        "search. [Default 6.0]")
    EventGroup.add_argument(
        "--max-mag",
        action="store",
        type=float,
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event for which to " +
        "search. [Default None, ie no limit]")

    # Geometry Settings
    GeomGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with the "
        "event-station geometries")
    GeomGroup.add_argument(
        "--min-dist",
        action="store",
        type=float,
        dest="mindist",
        default=85.,
        help="Specify the minimum great circle distance (degrees) " +
        "between the station and event. [Default 85]")
    GeomGroup.add_argument(
        "--max-dist",
        action="store",
        type=float,
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance (degrees) " +
        "between the station and event. [Default 120]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings " +
                "for User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # Parse Local Data directories
    if args.localdata is not None:
        args.localdata = args.localdata.split(',')
    else:
        args.localdata = []

    # Check NoData Value
    if args.ndval:
        args.ndval = 0.0
    else:
        args.ndval = nan

    return args


def get_arguments_offline(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    This function is used for processing SKS data offline 

    """

    parser = ArgumentParser(
        usage="Usage: %prog [arguments] <station database>",
        description="Script to process "
        "and calculate the spliting parameters for a dataset " +
        "that has already been downloaded by sks_prep.py. ")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network [Default " +
        "processes all stations in the database]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--Vp",
        action="store",
        type=float,
        dest="vp",
        default=6.,
        help="Specify default P velocity value. [Default 6.0 km/s]")
    ConstGroup.add_argument(
        "--SNR",
        action="store",
        type=float,
        dest="msnr",
        default=7.5,
        help="Specify the SNR threshold used to determine whether " +
        "events are processedc. [Default 7.5]")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="dts",
        default=120.,
        help="Specify time window length before and after the " +
        "SKS arrival. The total window length is 2*dst. " +
        "[Default 120 s]")
    ConstGroup.add_argument(
        "--max-delay",
        action="store",
        type=float,
        dest="maxdt",
        default=4.,
        help="Specify the maximum delay time. [Default 4 s]")
    ConstGroup.add_argument(
        "--time-increment",
        action="store",
        type=float,
        dest="ddt",
        default=0.1,
        help="Specify the time increment. [Default 0.1 s]")
    ConstGroup.add_argument(
        "--angle-increment",
        action="store",
        type=float,
        dest="dphi",
        default=1.,
        help="Specify the angle increment. [Default 1 d]")
    ConstGroup.add_argument(
        "--transverse-SNR",
        action="store",
        type=float,
        dest="snrTlim",
        default=1.,
        help="Specify the minimum SNR Threshold for the Transverse " +
        "component to be considered Non-Null. [Default 1.]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with " +
        "refining the events to include in matching station pairs")
    EventGroup.add_argument(
        "--start-time",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing the " +
        "start time for the event search. This will override any station " +
        "start times. [Default more recent start date for each station pair]")
    EventGroup.add_argument(
        "--end-time",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing the " +
        "end time for the event search. This will override any station " +
        "end times [Default older end date for each the pair of stations]")
    EventGroup.add_argument(
        "--reverse-order", "-R",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour starts at oldest " +
        "event and works towards most recent. Specify reverse order and " +
        "instead the program will start with the most recent events and " +
        "work towards older")
    EventGroup.add_argument(
        "--min-mag",
        action="store",
        type=float,
        dest="minmag",
        default=6.0,
        help="Specify the minimum magnitude of event for which to search. " +
        "[Default 6.0]")
    EventGroup.add_argument(
        "--max-mag",
        action="store",
        type=float,
        dest="maxmag",
        default=None,
        help="Specify the maximum magnitude of event for which to search. " +
        "[Default None, ie no limit]")

    # Geometry Settings
    GeomGroup = parser.add_argument_group(
        title="Geometry Settings",
        description="Settings associatd with "
        "the event-station geometries")
    GeomGroup.add_argument(
        "--min-dist",
        action="store",
        type=float,
        dest="mindist",
        default=85.,
        help="Specify the minimum great circle distance (degrees) " +
        "between the station and event. [Default 85]")
    GeomGroup.add_argument(
        "--max-dist",
        action="store",
        type=float,
        dest="maxdist",
        default=120.,
        help="Specify the maximum great circle distance (degrees) " +
        "between the station and event. [Default 120]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from start time: " + args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except:
            parser.error(
                "Cannot construct UTCDateTime from end time: " + args.endT)
    else:
        args.endT = None

    return args


def get_arguments_plot(argv=None):

    parser = ArgumentParser(
        usage="Usage: %prog [arguments] <station database>",
        description="Script to plot the average splitting results for a " +
        "given station. Loads the available pkl files in the specified " +
        "Station Directory.")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network [Default " +
        "processes all stations in the database]")
    parser.add_argument(
        "--no-figure",
        action="store_false",
        dest="showfig",
        default=True,
        help="Specify to prevent plots from opening during processing; " +
        "they are still saved to disk. [Default plots open and save]")

    # Null Settings
    NullGroup = parser.add_argument_group(
        title="Null Selection Settings",
        description="Settings "
        "associated with selecting which Null or Non-Null data is included")
    NullGroup.add_argument(
        "--nulls", "--Nulls",
        action="store_true",
        dest="nulls",
        default=False,
        help="Specify this flag to include Null Values in the average. " +
        "[Default Non-Nulls only]")
    NullGroup.add_argument(
        "--no-nons", "--No-Nons",
        action="store_false",
        dest="nons",
        default=True,
        help="Specify this flag to exclude Non-Nulls from the average " +
        "[Default False]")

    # Quality Settings
    QualGroup = parser.add_argument_group(
        title="Quality Selection Settings",
        description="Settings associated with selecting the qualities " +
        "to include in the selection.")
    QualGroup.add_argument(
        "--No-Good", "--no-good",
        action="store_false",
        dest="goods",
        default=True,
        help="Specify to exclude 'Good' measurements from the average. " +
        "[Default Good + Fair]")
    QualGroup.add_argument(
        "--No-Fair", "--no-fair",
        action="store_false",
        dest="fairs",
        default=True,
        help="Specify to exclude 'Fair' measurements from the average " +
        "[Default Good + Fair]")
    QualGroup.add_argument(
        "--Poor", "--poor",
        action="store_true",
        dest="poors",
        default=False,
        help="Specify to include 'Poor' measurements in the average " +
        "[Default No Poors]")

    # Split Type Settings
    SpTypGroup = parser.add_argument_group(
        title="Split Type Settings",
        description="Settings to Select "
        "which Split types are included in the selection.")
    SpTypGroup.add_argument(
        "--RC-Only", "--rc-only", "--RC-only",
        action="store_false",
        dest="SCinc",
        default=True,
        help="Specify to only include RC splits in the average. " +
        "[Default RC + SC]")
    SpTypGroup.add_argument(
        "--SC-Only", "--sc-only", "--SC-only",
        action="store_false",
        dest="RCinc",
        default=True,
        help="Specify to only include SC splits in the average. " +
        "[Default RC + SC]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # Check Nulls
    if not args.nons and not args.nulls:
        parser.error("One of Non-Nulls or Nulls must be included.")

    # Check Quality
    if not args.goods and not args.fairs and not args.poors:
        parser.error("At least one Quality must be included.")

    # Check Types
    if not args.RCinc and not args.SCinc:
        parser.error("At leat one Splitting Tyhpe must be included.")

    # Construct Null FileName Components
    NullName = ""
    if args.nons:
        NullName = "_Nons"
        if args.nulls:
            NullName = NullName + "-Nulls"
    else:
        if args.nulls:
            NullName = "_Nulls"
    args.NullName = NullName

    # Construct Quality FileName Components
    QualName = ""
    if args.goods:
        QualName = "_G"
        if args.fairs:
            QualName = QualName + "-F"
        if args.poors:
            QualName = QualName + "-P"
    else:
        if args.fairs:
            QualName = "_F"
            if args.poors:
                QualName = QualName + "-P"
        else:
            if args.poors:
                QualName = "_P"
    args.QualName = QualName

    # Construct Type FileName Components
    TypeName = ""
    if args.RCinc and args.SCinc:
        TypeName = "_RC-SC"
    elif args.RCinc and not args.SCinc:
        TypeName = "_RC"
    elif not args.RCinc and args.SCinc:
        TypeName = "_SC"
    args.TypeName = TypeName

    return args

