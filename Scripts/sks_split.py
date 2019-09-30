#!/usr/bin/env python


'''
PROGRAM sks_split.py

Calculates single station SKS splitting results.
Station selection is specified by a network and station code.
The data base is provided in stations_db.py as a dictionary.

'''

# Import splitpy, its classes and the conf module
import splitpy
from splitpy import Split, PickPlot, DiagPlot
from splitpy import Pick, Keep, Save, Repeat
from splitpy import conf as cf

# Import miscellaneous
import sys
import stdb
import os.path
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

    app = QApplication([])

    # Run Input Parser
    (opts, indb) = splitpy.utils.get_options()

    # Set Global Search Variables
    cf.maxdt = opts.maxdt #   maxdt is Max delay time
    cf.ddt = opts.ddt     #   ddt is time increment
    cf.dphi = opts.dphi   #   dphi is angle increment

    # Load Database
    db = stdb.io.load_db(fname=indb)

    # Construct station key loop
    allkeys = db.keys()
    sorted(allkeys)

    # Extract key subset
    if len(opts.stkeys) > 0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s])
    else:
        stkeys = db.keys()
        sorted(stkeys)
        # stkeys.sort()

    # Initialize Taup Model
    tpmodel = TauPyModel(model='iasp91')

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Initialize Split object
        split = Split(sta)

        # Output directory
        outdir = 'RESULTS/' + stkey
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        # Establish client
        if len(opts.UserAuth) == 0:
            client = Client(opts.Server)
        else:
            client = Client(opts.Server, user=opts.UserAuth[0], password=opts.UserAuth[1])

        # Get catalogue search start time
        if opts.startT is None:
            evst = sta.dstart
        else:
            evst = opts.startT

        # Get catalogue search end time
        if opts.endT is None:
            evet = sta.dend
        else:
            evet = opts.endT
        if evst > sta.dend or evet < sta.dstart:
            continue

        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0: tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0: tlocs[il] = "--"
        sta.location = tlocs

        # Update Display
        print(" ")
        print(" ")
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")
        print("| Searching Possible events:                    |")
        print("|   Start: {0:19s}                  |".format(evst.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(evet.strftime("%Y-%m-%d %H:%M:%S")))
        if opts.maxmag is None:
            print("|   Mag:   >{0:3.1f}                                 |".format(opts.minmag))
        else:
            print("|   Mag:   {0:3.1f} - {1:3.1f}                            |".format(opts.minmag, opts.maxmag))
        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = client.get_events(starttime=evst, endtime=evet,    \
                    minmagnitude=opts.minmag, maxmagnitude=opts.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print("|  Found {0:5d} possible events                  |".format(nevtT))

        # Get Local Data Availabilty
        if len(opts.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if opts.useNet:
                stalcllist = splitpy.io.list_local_data_stn(lcldrs=opts.localdata, sta=sta.station, \
                                 net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files                      |".format( \
                    sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = splitpy.io.list_local_data_stn(lcldrs=opts.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                         |".format(sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|===============================================|")      

        # Select order of processing
        if opts.reverse:
            ievs = range(0, nevtT)
        else:
            ievs = range(nevtT-1, -1, -1)

        # Read through catalogue
        for iev in ievs:

            # Extract event
            ev = cat[iev]

            # Add event to Split object
            split.add_event(ev)

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
                    inum = nevtT - iev + 1
                print(" ")
                print("****************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(nevK, inum, nevtT, split.meta["time"].strftime("%Y%m%d_%H%M%S")))
                print("*   Origin Time: " + split.meta["time"].strftime("%Y-%m-%d %H:%M:%S"))
                print("*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(split.meta["lat"], split.meta["lon"]))
                print("*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(split.meta["dep"]/1000., split.meta["mag"]))
                print("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; {3:6.2f}; {4:6.2f}".format(\
                    split.sta.station, split.meta["epi_dist"], split.meta["gac"], split.meta["baz"], split.meta["az"]))

                # Get Travel times (Careful: here dep is in meters)
                tt = tpmodel.get_travel_times(distance_in_degree=split.meta["gac"], \
                    source_depth_in_km=split.meta["dep"]/1000.)

                # Loop over all times in tt
                for t in tt:

                    # Extract SKS arrival
                    if t.name == 'SKS':

                        # Add SKS phase to Split
                        split.add_phase(t, opts.vp)

                        # Break out of loop 
                        break

                # Get data
                split.get_data_NEZ(client=client, dts=opts.dts, stdata=stalcllist, ndval=opts.ndval)

                if split.err: continue

                # Rotate from ZEN to LQT (Longitudinal, Radial, Transverse)
                split.rotate_ZEN_LQT()

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

                    # Check if result exists already
                    splitExists = False
                    if os.path.exists(outfile):
                        splitExists = True

                    # Should we Re-Pick
                    if splitExists:
                        if (not opts.ovr and not opts.skip) or (opts.ovr and opts.skip):
                            repeat = Repeat()
                            if not repeat.reply:
                            # if not splitpy.gui.repeat():
                                print("* Split Results Exist --> Skipping")
                                print("****************************************************")
                                continue
                            else:
                                print("* Split Results Exist --> Repeating")

                        elif not opts.ovr and opts.skip:
                            print("* Split Results Exist --> Skipping")
                            print("****************************************************")
                            continue
                        elif opts.ovr and not opts.skip:
                            print("* Split Results Exist --> Overwriting")

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

                            print("*   Time Range: {0:.2f}, {1:.2f}".format(tp1, tp2))

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

                                # Check if result exists
                                writeSplit = True
                                if splitExists:
                                    if not opts.ovr:
                                        save_obj = Save()
                                        ans3 = save_obj.reply
                                        # ans = splitpy.gui.save()
                                        if not ans3:
                                            writeSplit = False

                                # Write Output
                                if writeSplit:

                                    # Save split results
                                    split.save(outfile)

                                    # Save diagnostic plot
                                    dplot.save(outfig)
                                    
                                    print("* Estimate Saved to File: " + outfile)
                                else:
                                    print("* Existing Estimate Retained (" + outfile + ")")
                                print("****************************************************")
                            else:
                                print("* Estimate Discarded ")
                                print("****************************************************")
                                continue


if __name__ == "__main__":

    # Run main program
    main()
