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

from splitpy import Pick, Keep, Save, Repeat
from splitpy import PickPlot, DiagPlot
from splitpy import utils
from splitpy import Split

from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication

from argparse import ArgumentParser
from os.path import exists as exist
from pathlib import Path

matplotlib.use('Qt5Agg')


def get_arguments_calc_manual(argv=None):

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script to process "
        "and calculate the spliting parameters for a dataset " +
        "that has already been downloaded by split_calc_auto. ")

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
        "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")

    # Constants Settings
    ConstGroup = parser.add_argument_group(
        title='Parameter Settings',
        description="Miscellaneous default values and settings")
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
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing the " +
        "start time for the event search. This will override any station " +
        "start times. [Default more recent start date for each station pair]")
    EventGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing the " +
        "end time for the event search. This will override any station " +
        "end times [Default older end date for each the pair of stations]")
    EventGroup.add_argument(
        "--reverse-order",
        action="store_true",
        dest="reverse",
        default=False,
        help="Reverse order of events. Default behaviour starts at oldest " +
        "event and works towards most recent. Specify reverse order and " +
        "instead the program will start with the most recent events and " +
        "work towards older")

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

def main(args=None):

    if args is None:
        # Run Input Parser
        args = get_arguments_calc_manual()

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

        # Data directory
        datapath = Path("DATA") / stkey

        # Get only directories for which the key is available
        evs = [str(x) for x in datapath.iterdir() if x.is_dir()]

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

        # Get List of events to process
        evs.sort()
        nevs = len(evs)

        # Update Display
        if args.verb:
            print(" ")
            print(" ")
            print("|"+"="*50+"|")
            print(
                "|                   {0:>8s}                    |".format(stkey))
            print("|"+"="*50+"|")
            print(
                "|  Working on {0:5d} saved events                |".format(nevs))
            print("|"+"="*50+"|")

        # select order of processing
        if args.reverse:
            ievs = range(0, nevs)
        else:
            ievs = range(nevs - 1, -1, -1)
        nevK = 0

        # Read through catalogue
        for iev in ievs:

            datekey = evs[iev].split('/')[-1].split('_')[0]
            timekey = evs[iev].split('/')[-1].split('_')[-1]
            evdate = UTCDateTime(
                datekey[0:4]+'-'+datekey[4:6]+'-'+datekey[6:8]+
                'T'+timekey[0:2]+':'+timekey[2:4]+':'+timekey[4:6])
            if not (evdate > tstart and evdate < tend):
                continue

            # Event Name
            evSTR = evs[iev]

            # Load Relevant Files
            # Station data
            stafile = Path(evSTR) / "Station_data.pkl"
            sta = pickle.load(open(stafile, "rb"))
            split = Split(sta)

            # Event data
            metafile = Path(evSTR) / "Meta_data.pkl"
            meta = pickle.load(open(metafile, "rb"))
            split.meta = meta

            # ZNE data
            ZNEfile = Path(evSTR) / "ZNE_data.pkl"
            dataZNE = pickle.load(open(ZNEfile, "rb"))
            split.dataZNE = dataZNE

            # LQT data
            LQTfile = Path(evSTR) / "LQT_data.pkl"
            dataLQT = pickle.load(open(LQTfile, "rb"))
            split.dataLQT = dataLQT

            # Split results
            splitfile = Path(evSTR) / "Split_results_auto.pkl"
            if not splitfile.exists():

                print("* Split results not available... calculating")
                split.analyze(verbose=args.verb)
                split.is_null(args.snrTlim, verbose=args.verb)
                split.get_quality(verbose=args.verb)

            else:
                file = open(splitfile, "rb")
                split.SC_res = pickle.load(file)
                split.RC_res = pickle.load(file)
                split.null = pickle.load(file)
                split.quality = pickle.load(file)
                file.close()

            if args.verb:
                split.display_results()
                split.display_meta()
                split.display_null_quality()

            # Define time stamp
            yr = str(split.meta.time.year).zfill(4)
            jd = str(split.meta.time.julday).zfill(3)
            hr = str(split.meta.time.hour).zfill(2)

            # Initialize LQT seismogram figure and plot it
            pplot = PickPlot(split)
            pplot.plot_LQT_phases(args.dts)

            # Initialize diagnostic figure and plot it
            dplot = DiagPlot(split)
            dplot.plot_diagnostic()

            # Choose whether to re-pick window times
            iselect = 'a'
            while iselect == 'a':

                # Call interactive window for picking
                pick = Pick()
                ans = pick.reply

                # If user clicks yes:
                if ans:

                    print("* Refining Window in Pick Plot")

                    # Make sure pick plot figure is active
                    plt.figure(pplot.axes[0].number)

                    # Get clicks from LQT figure
                    xc = plt.ginput(2, show_clicks=True)

                    # Extract times from clicks
                    tp1 = [xx for xx, yy in xc][0]
                    tp2 = [xx for xx, yy in xc][1]

                    # Make sure tp2 > tp1
                    if tp2 < tp1:
                        tp11 = tp2
                        tp2 = tp1
                        tp1 = tp11

                    # Re-define time window
                    t1 = split.meta.time + split.meta.ttime + tp1
                    t2 = split.meta.time + split.meta.ttime + tp2

                    # Update LQT figure
                    pplot.update_LQT(tp1, tp2)

                    # Re-analyze splits
                    split.analyze(verbose=args.verb, t1=t1, t2=t2)

                    # Check for fault result
                    if (split.RC_res.edtt is None or
                            split.SC_res.edtt is None):
                        print("* !!! DOF Error. --> Skipping...")
                        print("*"*50)
                        continue

                    # Determine if estimate is Null and quality
                    # of estimate
                    split.calc_snr(t1=t1, dt=(t2-t1))
                    split.is_null(args.snrTlim, verbose=args.verb)
                    split.get_quality(verbose=args.verb)

                    if args.verb:
                        split.display_results()
                        split.display_meta()
                        split.display_null_quality()

                    # Re-initialize diagnostic figure and plot
                    dplot = DiagPlot(split)
                    dplot.plot_diagnostic(t1, t2)

                # If user clicks no:
                else:
                    iselect = 'c'

                    # Call interactive window for decision on
                    # estimate
                    keep = Keep()
                    ans2 = keep.reply

                    # If user keeps estimate:
                    if ans2:

                        # Print estimates to screen
                        print("* Split Accepted")
                        split.display_null_quality()
                        split.display_results()

                        # Save split results
                        splitfile = Path(evSTR) / "Split_results_manual.pkl"
                        file = open(splitfile, "wb")
                        split.SC_res = pickle.dump(split.SC_res, file)
                        split.RC_res = pickle.dump(split.RC_res, file)
                        split.null = pickle.dump(split.null, file)
                        split.quality = pickle.dump(split.quality, file)
                        file.close()

                        # Save pick plot
                        pplotfile = Path(evSTR) / "Plot_window_manual.png"
                        pplot.save(pplotfile)

                        # Save diagnostic plot
                        dplotfile = Path(evSTR) / "Plot_diagnostic_manual.png"
                        dplot.save(dplotfile)

                        print("* Estimate Saved")
                        print("*"*50)
                    else:
                        print("* Estimate Discarded ")
                        print("*"*50)
                        continue


if __name__ == "__main__":

    # Run main program
    main()
