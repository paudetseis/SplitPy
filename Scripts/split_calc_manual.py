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
import numpy as np
import pickle
import stdb
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from splitpy import Split
from splitpy import arguments, utils
from splitpy import PickPlot, DiagPlot
from splitpy import Pick, Keep, Save, Repeat
from pathlib import Path


def main():

    # Run Input Parser
    args = arguments.get_arguments_calc_manual()

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

        # Data directory
        datapath = Path("DATA") / stkey

        # Get only directories for which the key is available
        evs = [str(x) for x in datapath.iterdir() if x.is_dir()]

        # Add condition for data range
        #############################

        # Get List of events to process
        evs.sort()
        nevs = len(evs)

        # Update Display
        if args.verb:
            print(" ")
            print(" ")
            print("|"+"="*50+"|")
            print("|                   {0:>8s}                    |".format(stkey))
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
