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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#!/usr/bin/env python


"""
PROGRAM sks_OL_proc.py
----------------------

Calculates single station SKS splitting results using and offline process
where data exists on disk. Station selection is specified by a network and 
station code. The data base is provided in stations_db.pkl as a 
`StDb` dictionary.

"""

# -*- coding: utf-8 -*-
# Import splitpy, its classes and the conf module
import splitpy
from splitpy import Split, PickPlot, DiagPlot
from splitpy import Pick, Keep, Save, Repeat
from splitpy import conf as cf

# Import miscellaneous
import sys
import stdb
import dill
from os import path, listdir
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

# Import Obspy Modules
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client

# Main function
def main():

    # Run Input Parser
    (opts, indr) = splitpy.utils.get_options_OL_proc()

    # Set Global Search Variables
    cf.maxdt = opts.maxdt #   maxdt is Max delay time
    cf.ddt = opts.ddt     #   ddt is time increment
    cf.dphi = opts.dphi   #   dphi is angle increment

    # Generate list of available Station Folders
    stnkeys = listdir(indr)

    # Sort available stations
    stnkeys.sort()

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = stnkeys

    # Loop over station keys
    for ik in range(len(stkeys)):

        # Station Key
        stkey = stkeys[ik]

        # Get List of events to process
        evs = listdir(path.join(indr, stkey))
        evs.sort()
        nevs = len(evs)

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(stkey))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Working on {0:5d} saved events                |".format(nevs))
        print("|===============================================|")      

        # select order of processing
        if opts.reverse:
            ievs = range(0,nevs)
        else:
            ievs = range(nevs - 1, -1, -1)
        nevK=0

        # Read through catalogue
        for iev in ievs:

            # Event Name
            evSTR = evs[iev]

            # Load Relevant Files
            # Station Data
            sta = dill.load(open(path.join(indr, stkey, evSTR, "Station_Data.pkl"),"rb"))
            split = Split(sta)

            # Event Data
            ll = dill.load(open(path.join(indr, stkey, evSTR, "Event_Data.pkl"),"rb"))
            ev = ll[0]; tt=ll[1]
            split.add_event(ev)

            # Trace Filenames
            trpref = evSTR + "_" + sta.network + "." + sta.station

            # NEZ Trace files
            trs = dill.load(open(path.join(indr, stkey, evSTR, "NEZ_Data.pkl"),"rb"))
            split.add_NEZ(trs)

            # LQT Trace files
            trs = dill.load(open(path.join(indr, stkey, evSTR, "LQT_Data.pkl"),"rb"))
            split.add_LQT(trs)

            # Output directory
            outdir = path.join('RESULTS', sta.network + "." + sta.station)
            if not path.isdir(outdir):
                makedirs(outdir)

            # Define time stamp
            yr = str(split.meta["time"].year).zfill(4)
            jd = str(split.meta["time"].julday).zfill(3)
            hr = str(split.meta["time"].hour).zfill(2)

            # If distance between 85 and 120 deg:
            if (split.meta["gac"] > opts.mindist and split.meta["gac"] < opts.maxdist):

                # Display Event Info
                nevK = nevK + 1
                if opts.reverse:
                    inum = iev + 1
                else:
                    inum = nevs - iev + 1
                print(" ")
                print("****************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(nevK, inum, nevs, split.meta["time"].strftime("%Y%m%d_%H%M%S")))
                print("*   Origin Time: " + split.meta["time"].strftime("%Y-%m-%d %H:%M:%S"))
                print("*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(split.meta["lat"], split.meta["lon"]))
                print("*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(split.meta["dep"]/1000., split.meta["mag"]))
                print("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; {3:6.2f}; {4:6.2f}".format(\
                    split.sta.station, split.meta["epi_dist"], split.meta["gac"], split.meta["baz"], split.meta["az"]))


                # Loop over all times in tt
                for t in tt:

                    # Extract SKS arrival
                    if t.name=='SKS':

                        # Add SKS phase to Split
                        split.add_phase(t, opts.vp)

                        # Break out of loop 
                        break

                # Calculate snr over 30 seconds 
                split.calc_snrq(dt=30.)

                # Make sure no processing happens for NaNs
                if np.isnan(split.snrq): 
                    print("* SNR NaN...Skipping")
                    print("****************************************************")
                    continue 

                # SNR below threshold
                elif split.snrq < opts.msnr:
                    print("* SNR Failed: {0:.2f} < {1:.2f}...Skipping".format(split.snrq, opts.msnr))
                    print("****************************************************")

                # If SNR is higher than threshold
                elif split.snrq >= opts.msnr:
                    print("* SNR Passed: {0:4.2f} >= {1:3.1f}".format(split.snrq, opts.msnr))

                    # Output file
                    outfile = outdir + '/Split' + '.' + split.sta.station + '.' + \
                                    split.meta["time"].strftime("%Y.%j.%H%M%S") + '.pkl'
                    outfig = outdir + '/Split' + '.' + split.sta.station + '.' + \
                                    split.meta["time"].strftime("%Y.%j.%H%M%S") + '.png'

                    # Analyze
                    split.analyze()

                    # Continue if problem with analysis
                    if split.RC_res["edtt"] is None or split.SC_res["edtt"] is None:
                        print("* !!! DOF Error. --> Skipping...")
                        print("****************************************************")
                        continue

                    # Determine is Null and quality of estimate
                    split.calc_snrt()
                    split.is_null(opts.snrTlim, 5)
                    split.get_quality(5)

                    # Initialize LQT seismogram figure and plot it
                    pplot = PickPlot(split)
                    pplot.plot_LQT_phases(tt, opts.dts)

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

                            print("* Refine Window in Figure 1")

                            # Get clicks from LQT figure
                            plt.figure(1)
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
                            t1 = split.meta["time"] + split.ts + tp1
                            t2 = split.meta["time"] + split.ts + tp2

                            # Update LQT figure
                            pplot.update_LQT(tp1, tp2)

                            # Re-analyze splits
                            split.analyze(t1=t1, t2=t2)

                            # Check for fault result
                            if split.RC_res["edtt"] is None or split.SC_res["edtt"] is None:
                                print("* !!! DOF Error. --> Skipping...")
                                print("****************************************************")
                                continue

                            # Determine if estimate is Null and quality of estimate
                            split.calc_snrt(t1, t2 - t1)
                            split.is_null(opts.snrTlim, 5)
                            split.get_quality(5)

                            # Re-initialize diagnostic figure and plot
                            dplot = DiagPlot(split)
                            dplot.plot_diagnostic(t1, t2)

                        # If user clicks no:
                        else:
                            iselect = 'c'

                            # Call interactive window for decision on estimate
                            keep = Keep()
                            ans2 = keep.reply

                            # If user keeps estimate:
                            if ans2:

                                # Print estimates to screen
                                print("* Split Accepted")
                                split.display_null_quality(5)
                                split.display_results(5)

                                # Save split results
                                split.save(outfile)

                                # Save diagnostic plot
                                dplot.save(outfig)

                                print("* Estimate Saved")
                                print("****************************************************")
                            else:
                                print("* Estimate Discarded ")
                                print("****************************************************")
                                continue


if __name__ == "__main__":

    # Run main program
    main()
