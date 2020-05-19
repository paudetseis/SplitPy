from obspy import UTCDateTime
import numpy as np
import pickle
import stdb
from obspy.clients.fdsn import Client
from splitpy import Split
from . import test_args, get_meta
from pathlib import Path

def test_split(tmp_path):

    args = test_args.test_get_args()
    args.startT = UTCDateTime('2016-08-24')
    args.endT = UTCDateTime('2016-08-25')
    args.verb = True
    db = get_meta.get_stdb()

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
        datapath = tmp_path / stkey
        if not datapath.is_dir():
            datapath.mkdir()

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
            tstart.strftime("%Y-%m-%d %H:%M:%S")))
        print("|   End:   {0:19s}                  |".format(
            tend.strftime("%Y-%m-%d %H:%M:%S")))
        if args.maxmag is None:
            print("|   Mag:   >{0:3.1f}", format(args.minmag) +
                  "                                 |")
        else:
            msg = "|   Mag:   {0:3.1f}".format(args.minmag) + \
                " - {0:3.1f}".format(args.maxmag) + \
                "                            |"
            print(msg)

        print("| ...                                           |")

        # Get catalogue using deployment start and end
        cat = event_client.get_events(
            starttime=tstart, endtime=tend,
            minmagnitude=args.minmag, maxmagnitude=args.maxmag)

        # Total number of events in Catalogue
        nevK = 0
        nevtT = len(cat)
        print(
            "|  Found {0:5d}".format(nevtT) +
            " possible events                  |")
        ievs = range(0, nevtT)

        # Get Local Data Availabilty
        if len(args.localdata) > 0:
            print("|-----------------------------------------------|")
            print("| Cataloging Local Data...                      |")
            if args.useNet:
                stalcllist = io.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station,
                    net=sta.network, altnet=sta.altnet)
                print("|   {0:>2s}.{1:5s}: {2:6d} files              " +
                      "        |".format(
                          sta.network, sta.station, len(stalcllist)))
            else:
                stalcllist = io.list_local_data_stn(
                    lcldrs=args.localdata, sta=sta.station)
                print("|   {0:5s}: {1:6d} files                      " +
                      "   |".format(
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

                assert accept

                # Display Event Info
                nevK = nevK + 1
                if args.reverse:
                    inum = iev + 1
                else:
                    inum = nevtT - iev + 1
                print(" ")
                print("**************************************************")
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
                evtdir = datapath / timekey
                splitfile = evtdir / 'Split_Data.pkl'
                ZNEfile = evtdir / 'ZNE_Data.pkl'
                metafile = evtdir / 'Meta_Data.pkl'
                stafile = evtdir / 'Station_Data.pkl'

                # Check if RF data already exist and overwrite has been set
                if evtdir.exists():
                    if RFfile.exists():
                        if not args.ovr:
                            continue

                # Get data
                has_data = split.download_data(
                    client=data_client, dts=args.dts, stdata=stalcllist,
                    ndval=args.ndval, new_sr=args.new_sampling_rate,
                    returned=True, verbose=args.verb)

                if not has_data:
                    continue

                assert has_data 

                # Create Folder if it doesn't exist
                if evtdir.exists():
                    evtdir.mkdir()

                # Save ZNE Traces
                pickle.dump(split.dataZNE, open(tmp_path / 'ZNEfile.pkl', "wb"))

                # Rotate from ZNE to 'LQT'
                split.rotate(align='LQT')

                # Calculate snr over dt_snr seconds
                split.calc_snr(
                    dt=args.dt_snr, fmin=args.fmin, fmax=args.fmax)
                if args.verb:
                    print("* SNRQ: {}".format(split.meta.snrq))
                    print("* SNRT: {}".format(split.meta.snrt))

                # Make sure no processing happens for NaNs
                if np.isnan(split.meta.snrq):
                    if args.verb:
                        print("* SNR NaN...Skipping")
                    print("**************************************************")
                    continue

                # Analyze
                split.analyze(verbose=args.verb)

                # Continue if problem with analysis
                if split.RC_res.edtt is None or split.SC_res.edtt is None:
                    if args.verb:
                        print("* !!! DOF Error. --> Skipping...")
                        print("*********************************************" +
                              "*******")
                    continue

                # Determine if Null and quality of estimate
                split.is_null(args.snrTlim, 5)
                split.get_quality(5)

                # Display results
                split.display_results()
                split.display_meta()
                split.display_null_quality()


