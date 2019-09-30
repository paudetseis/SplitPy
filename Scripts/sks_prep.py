#!/usr/bin/env python


'''
PROGRAM sks_prep.py

First of 2 programs for SKS Split to prepare data for later offline 
calculation. This script downloads the event catalogue for a given 
station and downloads the viable data. The second code 
sks_offline_split.py processes the prepared data.

'''

# Import splitpy, its classes and the conf module
import splitpy
from splitpy import Split
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

    # Run Input Parser
    (opts, indb) = splitpy.utils.get_options_prep_offline()

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
    if len(opts.stkeys)>0:
        stkeys = []
        for skey in opts.stkeys:
            stkeys.extend([s for s in allkeys if skey in s] )
    else:
        stkeys = db.keys()
        sorted(stkeys)

    # Initialize Taup Model
    tpmodel = TauPyModel(model='iasp91')

    # Output directory
    if opts.startT is not None and opts.endT is not None:
        Dtrange = "_{0:s}-{1:s}".format(opts.startT.strftime("%Y%m%d"), opts.endT.strftime("%Y%m%d"))
    elif opts.startT is None and opts.endT is not None:
        Dtrange = "_-{0:s}".format(opts.endT.strftime("%Y%m%d"))
    elif opts.startT is not None and opts.endT is None:
        Dtrange = "_{0:s}-".format(opts.startT.strftime("%Y%m%d"))
    else:
        Dtrange = ""

    if opts.maxmag is None:
        mxmag = "+"
    else:
        mxmag = "-{0:.1f}".format(opts.maxmag)

    dtdir = "{0:s}{1:s}_D{2:.0f}-{3:.0f}_M{4:.1f}{5:s}_S{6:.1f}+".format(opts.datadir, Dtrange, \
                    opts.mindist, opts.maxdist, opts.minmag, mxmag, opts.msnr)
    if not os.path.isdir(dtdir): os.makedirs(dtdir)

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Initialize Split object
        split = Split(sta)

        # Create Station Data Folder
        outdir = dtdir + "/" + stkey

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

        #- Temporary print locations
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
        print ("|  Found {0:5d} possible events                  |".format(nevtT))

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
            yr = str(time.year).zfill(4)
            jd = str(time.julday).zfill(3)
            hr = str(time.hour).zfill(2)

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
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(nevK, inum,nevtT, split.meta["time"].strftime("%Y%m%d_%H%M%S")))
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
                    if t.name=='SKS':

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

                    # Create Event Folder
                    evtdir = outdir + "/" + time.strftime("%Y%m%d_%H%M%S")

                    # Create Folder
                    if not os.path.isdir(evtdir): os.makedirs(evtdir)

                    # Event Data
                    dill.dump([ev, tstart, tend, tt],open(evtdir + "/Event_Data.pkl","wb"))

                    # Station Data
                    dill.dump(sta, open(evtdir + "/Station_Data.pkl","wb"))

                    # Trace Filenames
                    trpref=time.strftime("%Y%m%d_%H%M%S") + "_" + sta.network + "." + sta.station

                    # Raw Trace files
                    dill.dump([split.data["trN"], split.data["trE"], split.data["trZ"]], \
                        open(evtdir + "/NEZ_Data.pkl","wb"))
                    split.data["trN"].write(os.path.join(evtdir, trpref + ".N.mseed"), format='MSEED')
                    split.data["trE"].write(os.path.join(evtdir, trpref + ".E.mseed"), format='MSEED')
                    split.data["trZ"].write(os.path.join(evtdir, trpref + ".Z.mseed"), format='MSEED')

                    # LQT Traces
                    dill.dump([split.data["trL"], split.data["trQ"], split.data["trT"]], \
                        open(evtdir + "/LQT_Data.pkl","wb"))
                    split.data["trL"].write(os.path.join(evtdir, trpref + ".L.mseed"), format='MSEED')
                    split.data["trQ"].write(os.path.join(evtdir, trpref + ".Q.mseed"), format='MSEED')
                    split.data["trT"].write(os.path.join(evtdir, trpref + ".T.mseed"), format='MSEED')

                    # Update
                    print ("* Wrote Output Files to: ")
                    print ("*     "+evtdir)
                    print ("****************************************************")




###############################
# Choose one station to process

if __name__ == "__main__":

    # Run main program
    main()

###################
