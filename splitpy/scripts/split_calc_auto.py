#!/usr/bin/env python

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

# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pickle
import stdb

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

from splitpy import utils
from splitpy import Split, DiagPlot

from argparse import ArgumentParser
from os.path import exists as exist
from pathlib import Path

matplotlib.use('Qt5Agg')


def get_arguments_calc_auto(argv=None):

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
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
        "-V", "--verbose",
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
        "--skip-existing",
        action="store_true",
        dest="skip",
        default=False,
        help="Skip any event for which existing splitting results are " +
        "saved to disk. Default behaviour prompts for " +
        "each event. Selecting skip and overwrite (ie, both flags) " +
        "negate each other, and both are set to " +
        "False (every repeat is prompted). [Default False]")
    parser.add_argument(
        "--calc",
        action="store_true",
        dest="calc",
        default=False,
        help="Analyze data for shear-wave splitting. [Default saves data "+
        "to folders for subsequent analysis]")
    parser.add_argument(
        "--plot-diagnostic",
        action="store_true",
        dest="diagplot",
        default=False,
        help="Plot diagnostic window at end of process. [Default False]")
    parser.add_argument(
        "--recalc",
        action="store_true",
        dest="recalc",
        default=False,
        help="Re-calculate estimates and overwrite existing splitting "+
        "results without re-downloading data. [Default False]")

    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which " +
        "datacenter to log into.")
    ServerGroup.add_argument(
        "--server",
        action="store",
        type=str,
        dest="server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. [Default IRIS]")
    ServerGroup.add_argument(
        "--user-auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--user-auth='username:authpassword') to access and download " +
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
        "--dtype",
        action="store",
        type=str,
        dest="dtype",
        default='SAC',
        help="Specify the data archive file type, either SAC " +
        " or MSEED. Note the default behaviour is to search for " +
        "SAC files. Local archive files must have extensions of '.SAC' "+
        " or '.MSEED. These are case dependent, so specify the correct case"+
        "here.")
    # DataGroup.add_argument(
    #     "--no-data-zero",
    #     action="store_true",
    #     dest="ndval",
    #     default=False,
    #     help="Specify to force missing data to be set as zero, rather " +
    #     "than default behaviour which sets to nan.")
    # DataGroup.add_argument(
    #     "--no-local-net",
    #     action="store_false",
    #     dest="useNet",
    #     default=True,
    #     help="Specify to prevent using the Network code in the " +
    #     "search for local data (sometimes for CN stations " +
    #     "the dictionary name for a station may disagree with that " +
    #     "in the filename. [Default Network used]")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
    ConstGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=10.,
        help="Specify new sampling rate in Hz. [Default 10.]")
    ConstGroup.add_argument(
        "--min-snr",
        action="store",
        type=float,
        dest="msnr",
        default=5.,
        help="Minimum SNR value calculated on the radial (Q) component "+
        "to proceed with analysis (dB). [Default 5.]")
    ConstGroup.add_argument(
        "--window",
        action="store",
        type=float,
        dest="dts",
        default=120.,
        help="Specify time window length before and after the SKS "
        "arrival. The total window length is 2*dst (sec). [Default 120]")
    ConstGroup.add_argument(
        "--max-delay",
        action="store",
        type=float,
        dest="maxdt",
        default=4.,
        help="Specify the maximum delay time in search (sec). "+
        "[Default 4]")
    ConstGroup.add_argument(
        "--dt-delay",
        action="store",
        type=float,
        dest="ddt",
        default=0.1,
        help="Specify the time delay increment in search (sec). "+
        "[Default 0.1]")
    ConstGroup.add_argument(
        "--dphi",
        action="store",
        type=float,
        dest="dphi",
        default=1.,
        help="Specify the fast angle increment in search (degree). "+
        "[Default 1.]")
    ConstGroup.add_argument(
        "--snrT",
        action="store",
        type=float,
        dest="snrTlim",
        default=1.,
        help="Specify the minimum SNR Threshold for the Transverse " +
        "component to be considered Non-Null. [Default 1.]")
    ConstGroup.add_argument(
        "--fmin",
        action="store",
        type=float,
        dest="fmin",
        default=0.02,
        help="Specify the minimum frequency corner for bandpass " +
        "filter (Hz). [Default 0.02]")
    ConstGroup.add_argument(
        "--fmax",
        action="store",
        type=float,
        dest="fmax",
        default=0.5,
        help="Specify the maximum frequency corner for bandpass " +
        "filter (Hz). [Default 0.5]")

    # Event Selection Criteria
    EventGroup = parser.add_argument_group(
        title="Event Settings",
        description="Settings associated with refining "
        "the events to include in matching station pairs")
    EventGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station start times. [Default start date of each station]")
    EventGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the end time for the event search. This will override any " +
        "station end times [Default end date of each station]")
    EventGroup.add_argument(
        "--reverse",
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
    GeomGroup.add_argument(
        "--phase",
        action="store",
        type=str,
        dest="phase",
        default='SKS',
        help="Specify the phase name to use. Be careful with the distance. " +
        "setting. Options are 'SKS' or 'SKKS'. [Default 'SKS']")

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

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    # Check selected phase
    if args.phase not in ['SKS', 'SKKS', 'PKS']:
        parser.error(
            "Error: choose between 'SKS', 'SKKS and 'PKS'.")

    # Check distances for all phases
    if not args.mindist:
        if args.phase == 'SKS':
            args.mindist = 85.
        elif args.phase == 'SKKS':
            args.mindist = 90.
        elif args.phase == 'PKS':
            args.mindist = 130.
    if not args.maxdist:
        if args.phase == 'SKS':
            args.maxdist = 120.
        elif args.phase == 'SKKS':
            args.maxdist = 130.
        elif args.phase == 'PKS':
            args.maxdist = 150.
    if args.mindist < 85. or args.maxdist > 180.:
        parser.error(
            "Distances should be between 85 and 180 deg. for " +
            "teleseismic 'SKS', 'SKKS' and 'PKS' waves.")

    return args

def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_arguments_calc_auto()

    # Load Database
    # stdb=0.1.4
    try:
        db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # stdb=0.1.3
    except:
        db = stdb.io.load_db(fname=args.indb)

        # Construct station key loop
        allkeys = db.keys()
        sorted(allkeys)

        # Extract key subset
        if len(args.stkeys) > 0:
            stkeys = []
            for skey in args.stkeys:
                stkeys.extend([s for s in allkeys if skey in s])
        else:
            stkeys = db.keys()
            sorted(stkeys)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Output directory
        datapath = Path('DATA') / stkey
        if not datapath.is_dir():
            datapath.mkdir(parents=True)

        # Establish client
        if len(args.UserAuth) == 0:
            data_client = Client(args.server)
        else:
            data_client = Client(
                args.server,
                user=args.UserAuth[0],
                password=args.UserAuth[1])

        # Establish client for events
        event_client = Client()

        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT

        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT
        if tstart > sta.enddate or tend < sta.startdate:
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|"+"="*50+"|")
        print("|                   {0:>8s}                       |".format(
            sta.station))
        print("|"+"="*50+"|")
        print("|  Station: {0:>2s}.{1:5s}                               |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}     |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                   |".format(
            sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}             |".format(
            sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}             |".format(
            sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|"+"-"*50+"|")
        print("| Searching Possible events:                       |")
        print("|   Start: {0:19s}                     |".format(
            tstart.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                     |".format(
            tend.strftime("%Y-%m-%d %H:%M:%S")))
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}".format(args.minmag) +
                  "                                     |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(args.minmag) + \
                " - {0:3.1f}".format(args.maxmag) + \
                "                           |"
            print(msg)

        print("| ...                                              |")

        # Get catalogue using deployment start and end
        cat = event_client.get_events(
            starttime=tstart,
            endtime=tend,
            minmagnitude=args.minmag,
            maxmagnitude=args.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d}".format(nevtT) +
            " possible events                     |")
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(args.localdata) > 0:
            print("|"+"-"*50+"|")
            print("| Cataloging Local Data...                         |")
            if args.useNet:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata,
                    sta=sta.station,
                    net=sta.network,
                    dtype=args.dtype,
                    altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files              ".format(
                          sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata,
                    dtype=args.dtype,
                    sta=sta.station)
                print("|   {0:5s}: {1:6d} files                      ".format(
                          sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|"+"="*50+"|")

        # Select order of processing
        if args.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            # Initialize Split object with station info
            split = Split(sta)

            # Add event to split object
            accept = split.add_event(
                ev,
                gacmin=args.mindist,
                gacmax=args.maxdist,
                phase=args.phase,
                returned=True)

            # Define time stamp
            yr = str(split.meta.time.year).zfill(4)
            jd = str(split.meta.time.julday).zfill(3)
            hr = str(split.meta.time.hour).zfill(2)

            # If event is accepted (data exists)
            if accept:

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("|"+"*"*50+"|")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s} {4}".format(
                    nevK, inum, nevtT, split.meta.time.strftime(
                        "%Y%m%d_%H%M%S"), stkey))
                if args.verb:
                    print("*   Phase: {}".format(args.phase))
                    print("*   Origin Time: " +
                          split.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                    print(
                        "*   Lat: {0:6.2f};        Lon: {1:7.2f}".format(
                            split.meta.lat, split.meta.lon))
                    print(
                        "*   Dep: {0:6.2f} km;     Mag: {1:3.1f}".format(
                            split.meta.dep, split.meta.mag))
                    print("*   Dist: {0:7.2f} km;".format(split.meta.epi_dist) +
                          "   Epi dist: {0:6.2f} deg\n".format(split.meta.gac) +
                          "*   Baz:  {0:6.2f} deg;".format(split.meta.baz) +
                          "   Az: {0:6.2f} deg".format(split.meta.az))

                # Event Folder
                timekey = split.meta.time.strftime("%Y%m%d_%H%M%S")
                datadir = datapath / timekey
                ZNEfile = datadir / 'ZNE_data.pkl'
                LQTfile = datadir / 'LQT_data.pkl'
                metafile = datadir / 'Meta_data.pkl'
                stafile = datadir / 'Station_data.pkl'
                splitfile = datadir / 'Split_results_auto.pkl'

                # Check if RF data already exist and overwrite has been set
                if datadir.exists():
                    if splitfile.exists():
                        if not args.ovr:
                            continue

                if args.recalc:
                    if np.sum([file.exists() for file in
                               [ZNEfile, metafile, stafile]]) < 3:
                        continue
                    sta = pickle.load(open(stafile, "rb"))
                    split = Split(sta)
                    meta = pickle.load(open(metafile, "rb"))
                    split.meta = meta
                    dataZNE = pickle.load(open(ZNEfile, "rb"))
                    split.dataZNE = dataZNE

                    # Rotate from ZNE to 'LQT'
                    split.rotate(align='LQT')

                    # Filter rotated traces
                    split.dataLQT.filter(
                        'bandpass',
                        freqmin=args.fmin,
                        freqmax=args.fmax)

                    # Calculate snr over dt_snr seconds
                    split.calc_snr()

                    # Save LQT Traces
                    pickle.dump(split.dataLQT, open(LQTfile, "wb"))

                else:

                    # Get data
                    has_data = split.download_data(
                        client=data_client,
                        dts=args.dts,
                        stdata=stalcllist,
                        dtype=args.dtype,
                        ndval=args.ndval,
                        new_sr=args.new_sampling_rate,
                        returned=True,
                        verbose=args.verb)

                    if not has_data:
                        continue

                    # Rotate from ZNE to 'LQT'
                    split.rotate(align='LQT')

                    # Filter rotated traces
                    split.dataLQT.filter(
                        'bandpass',
                        freqmin=args.fmin,
                        freqmax=args.fmax)

                    # Calculate snr over dt_snr seconds
                    split.calc_snr()

                    # If SNR lower than user-specified threshold, continue
                    if split.meta.snrq < args.msnr:
                        if args.verb:
                            print(
                                "* SNRQ < {0:.1f}, continuing".format(args.msnr))
                            print("*"*50)
                        continue

                    # Make sure no processing happens for NaNs
                    if np.isnan(split.meta.snrq):
                        if args.verb:
                            print("* SNR NaN, continuing")
                            print("*"*50)
                        continue

                    # Create Folder if it doesn't exist
                    if not datadir.exists():
                        datadir.mkdir(parents=True)

                    # Save ZNE Traces
                    pickle.dump(split.dataZNE, open(ZNEfile, "wb"))

                    # Save LQT Traces
                    pickle.dump(split.dataLQT, open(LQTfile, "wb"))

                if args.verb:
                    print("* SNRQ: {}".format(split.meta.snrq))
                    print("* SNRT: {}".format(split.meta.snrt))

                if args.calc or args.recalc:

                    # Analyze
                    split.analyze(verbose=args.verb)

                    # Continue if problem with analysis
                    if split.RC_res.edtt is None or split.SC_res.edtt is None:
                        if args.verb:
                            print("* !!! DOF Error. --> Skipping...")
                            print("*"*50)
                        continue

                    # Determine if Null and Quality of estimate
                    split.is_null(args.snrTlim, verbose=args.verb)
                    split.get_quality(verbose=args.verb)

                # Display results
                if args.verb:
                    split.display_meta()
                    if args.calc or args.recalc:
                        split.display_results()
                        split.display_null_quality()

                # Save event meta data
                pickle.dump(split.meta, open(metafile, "wb"))

                # Save Station Data
                pickle.dump(split.sta, open(stafile, "wb"))

                if args.calc or args.recalc:
                    # Save Split Data
                    file = open(splitfile, "wb")
                    pickle.dump(split.SC_res, file)
                    pickle.dump(split.RC_res, file)
                    pickle.dump(split.null, file)
                    pickle.dump(split.quality, file)
                    file.close()

                    # Initialize diagnostic figure and plot it
                    if args.diagplot:
                        dplot = DiagPlot(split)
                        dplot.plot_diagnostic()
                        plt.figure(dplot.axes[0].number)
                        plt.show()


if __name__ == "__main__":

    # Run main program
    main()
