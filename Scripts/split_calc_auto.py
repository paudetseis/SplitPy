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
from pathlib import Path
from splitpy import arguments, utils
from splitpy import Split, DiagPlot
import matplotlib.pyplot as plt
import numpy as np
import pickle
import stdb
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib
matplotlib.use('Qt5Agg')


def main():

    # Run Input Parser
    args = arguments.get_arguments_calc_auto()

    # Load Database
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
            data_client = Client(args.Server)
        else:
            data_client = Client(
                args.Server, user=args.UserAuth[0], password=args.UserAuth[1])

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
            print("|   Mag:   >{0:3.1f}", format(args.minmag) +
                  "                           |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(args.minmag) + \
                " - {0:3.1f}".format(args.maxmag) + \
                "                           |"
            print(msg)

        print("| ...                                              |")

        # Get catalogue using deployment start and end
        cat = event_client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=args.minmag, maxmagnitude=args.maxmag)

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
                    lcldrs=args.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files              " +
                      "        |".format(
                          sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = utils.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                      " +
                      "   |".format(
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
                ev, gacmin=args.mindist, gacmax=args.maxdist,
                phase=args.phase, returned=True)

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
                    split.dataLQT.filter('bandpass', freqmin=args.fmin,
                                         freqmax=args.fmax)
                    
                    # Calculate snr over dt_snr seconds
                    split.calc_snr()

                    # Save LQT Traces
                    pickle.dump(split.dataLQT, open(LQTfile, "wb"))

                else:

                    # Get data
                    has_data = split.download_data(
                        client=data_client, dts=args.dts, stdata=stalcllist,
                        ndval=args.ndval, new_sr=args.new_sampling_rate,
                        returned=True, verbose=args.verb)

                    if not has_data:
                        continue

                    # Rotate from ZNE to 'LQT'
                    split.rotate(align='LQT')

                    # Filter rotated traces
                    split.dataLQT.filter('bandpass', freqmin=args.fmin,
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
