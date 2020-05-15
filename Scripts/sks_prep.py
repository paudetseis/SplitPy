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

from pathlib import Path
from splitpy import Split
from splitpy import arguments, io
"""
Program sks_prep.py
-------------------

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

"""

# -*- coding: utf-8 -*-
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication
import matplotlib.pyplot as plt
import stdb
import dill
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')


def main():

    # Run Input Parser
    args = arguments.get_arguments_prep()

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

    # Initialize Taup Model
    tpmodel = TauPyModel(model='iasp91')

    # Output directory
    if args.startT is not None and args.endT is not None:
        Dtrange = "_{0:s}-{1:s}".format(args.startT.strftime("%Y%m%d"),
                                        args.endT.strftime("%Y%m%d"))
    elif args.startT is None and args.endT is not None:
        Dtrange = "_-{0:s}".format(args.endT.strftime("%Y%m%d"))
    elif args.startT is not None and args.endT is None:
        Dtrange = "_{0:s}-".format(args.startT.strftime("%Y%m%d"))
    else:
        Dtrange = ""

    if args.maxmag is None:
        mxmag = "+"
    else:
        mxmag = "-{0:.1f}".format(args.maxmag)

    dtdir = Path("{0:s}{1:s}_D{2:.0f}-{3:.0f}_M{4:.1f}{5:s}_S{6:.1f}+".format(
        args.datadir, Dtrange, args.mindist, args.maxdist,
        args.minmag, mxmag, args.msnr))
    if not dtdir.is_dir():
        dtdir.mkdir()

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        # Initialize Split object
        split = Split(sta, args.maxdt, args.ddt, args.dphi)

        # Create Station Data Folder
        outdir = dtdir / stkey

        # Establish client
        if len(args.UserAuth) == 0:
            client = Client(args.Server)
        else:
            client = Client(
                args.Server, user=args.UserAuth[0], password=args.UserAuth[1])

        # Get catalogue search start time
        if args.startT is None:
            evst = sta.dstart
        else:
            evst = args.startT

        # Get catalogue search end time
        if args.endT is None:
            evet = sta.dend
        else:
            evet = args.endT
        if evst > sta.dend or evet < sta.dstart:
            continue

        # - Temporary print locations
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
        print("|===============================================|")
        print("|===============================================|")
        print("|                   {0:>8s}                    |".format(
            sta.station))
        print("|===============================================|")
        print("|===============================================|")
        print("|  Station: {0:>2s}.{1:5s}                            |".format(
            sta.network, sta.station))
        print("|      Channel: {0:2s}; Locations: {1:15s}  |".format(
            sta.channel, ",".join(tlocs)))
        print("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(
            sta.longitude, sta.latitude))
        print("|      Start time: {0:19s}          |".format(
            sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|      End time:   {0:19s}          |".format(
            sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
        print("|-----------------------------------------------|")
        print("| Searching Possible events:                    |")
        print("|   Start: {0:19s}                  |".format(
            evst.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(
            evet.strftime("%Y-%m-%d %H:%M:%S")))
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}                              " +
                  "   |".format(
                      args.minmag))
        else:
            print(
                "|   Mag:   {0:3.1f} - {1:3.1f}                       " +
                "     |".format(args.minmag, args.maxmag))
        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = client.get_events(
            starttime=evst, endtime=evet,
            minmagnitude=args.minmag, maxmagnitude=args.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d} possible events                  |".format(nevtT))

        # Get Local Data Availabilty
        if len(args.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if args.useNet:
                stalcllist = io.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files                " +
                      "      |".format(
                          sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = io.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                         " +
                      "|".format(
                          sta.station, len(stalcllist)))
        else:
            stalcllist = []
        print("|===============================================|")

        # Select order of processing
        if args.reverse:
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
            yr = str(split.meta.time.year).zfill(4)
            jd = str(split.meta.time.julday).zfill(3)
            hr = str(split.meta.time.hour).zfill(2)

            # If distance between 85 and 120 deg:
            if (split.meta.gac > args.mindist and
                    split.meta.gac < args.maxdist):

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("****************************************************")
                print("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(
                    nevK, inum, nevtT, split.meta.time.strftime(
                        "%Y%m%d_%H%M%S")))
                print("*   Origin Time: " +
                      split.meta.time.strftime("%Y-%m-%d %H:%M:%S"))
                print(
                    "*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(
                        split.meta.lat, split.meta.lon))
                print(
                    "*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(
                        split.meta.dep/1000., split.meta.mag))
                print("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; " +
                      "{3:6.2f}; {4:6.2f}".format(
                          split.sta.station, split.meta.epi_dist,
                          split.meta.gac, split.meta.baz, split.meta.az))

                # Get Travel times (Careful: here dep is in meters)
                tt = tpmodel.get_travel_times(
                    distance_in_degree=split.meta.gac,
                    source_depth_in_km=split.meta.dep/1000.)

                # Loop over all times in tt
                for t in tt:

                    # Extract SKS arrival
                    if t.name == 'SKS':

                        # Add SKS phase to Split
                        split.add_phase(t, args.vp)

                        # Break out of loop
                        break

                # Get data
                split.get_data_NEZ(
                    client=client, dts=args.dts,
                    stdata=stalcllist, ndval=args.ndval)

                if split.err:
                    continue

                # Rotate from ZEN to LQT (Longitudinal, Radial, Transverse)
                split.rotate_ZEN_LQT()

                # Calculate snr over 30 seconds
                split.calc_snrq(dt=30.)

                # Make sure no processing happens for NaNs
                if np.isnan(split.snrq):
                    print("* SNR NaN...Skipping")
                    print("*******************************************" +
                          "*********")
                    continue

                # SNR below threshold
                elif split.snrq < args.msnr:
                    print(
                        "* SNR Failed: {0:.2f} < {1:.2f}...Skipping".format(
                            split.snrq, args.msnr))
                    print("********************************************" +
                          "********")

                # If SNR is higher than threshold
                elif split.snrq >= args.msnr:
                    print(
                        "* SNR Passed: {0:4.2f} >= {1:3.1f}".format(
                            split.snrq, args.msnr))

                    # Create Event Folder
                    evtdir = outdir / \
                        split.meta.time.strftime("%Y%m%d_%H%M%S")

                    # Create Folder
                    if not evtdir.is_dir():
                        evtdir.mkdir()

                    # Event Data
                    filename = evtdir / "Event_Data.pkl"
                    dill.dump([ev, tt], open(filename, "wb"))

                    # Station Data
                    filename = evtdir / "Station_Data.pkl"
                    dill.dump(sta, open(filename, "wb"))

                    # Trace Filenames
                    trpref = split.meta.time.strftime(
                        "%Y%m%d_%H%M%S") + "_" + sta.network + "." + \
                        sta.station

                    # Raw Trace files
                    filename = evtdir / "NEZ_Data.pkl"
                    dill.dump([split.data.trN, split.data.trE, split.data.trZ],
                              open(filename, "wb"))

                    # split.data.trN.write(os.path.join(
                    #     evtdir, trpref + ".N.mseed"), format='MSEED')
                    # split.data.trE.write(os.path.join(
                    #     evtdir, trpref + ".E.mseed"), format='MSEED')
                    # split.data.trZ.write(os.path.join(
                    #     evtdir, trpref + ".Z.mseed"), format='MSEED')

                    # LQT Traces
                    filename = evtdir / "LQT_Data.pkl"
                    dill.dump([split.data.trL, split.data.trQ, split.data.trT],
                              open(filename, "wb"))
                    # split.data.trL.write(os.path.join(
                    #     evtdir, trpref + ".L.mseed"), format='MSEED')
                    # split.data.trQ.write(os.path.join(
                    #     evtdir, trpref + ".Q.mseed"), format='MSEED')
                    # split.data.trT.write(os.path.join(
                    #     evtdir, trpref + ".T.mseed"), format='MSEED')

                    # Update
                    print("* Wrote Output Files to: ")
                    print("*     "+str(evtdir))
                    print("******************************************" +
                          "**********")


if __name__ == "__main__":

    # Run main program
    main()
