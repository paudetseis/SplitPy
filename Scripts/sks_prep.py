#!/usr/bin/env python


'''
PROGRAM sks_prep.py

First of 2 programs for SKS Split to prepare data for later offline calculation.
This script downloads the event catalogue for a given station and downloads the viable data. The second code sks_offline_split.py processes the prepared data.

'''

from sys import exit as sexit

#-- Import Splitting Module
import SplitPy
import SplitPy.classes as spc
import StDb
from SplitPy import conf as cf

#-- Import pickle
try:
  import cPickle as pickle
except:
  import pickle
import dill

# Import NumPy and Matplotlib modules and functions
import numpy as np
import numpy.fft as ft
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#-- Os Modules
import os.path

#-- Import Obspy Modules
from obspy import UTCDateTime
from obspy.core import read, Stream, Trace
from obspy.taup import TauPyModel
#-- 1.x
from obspy.clients.fdsn import Client
#-- 1.x
from obspy.geodetics.base import gps2dist_azimuth as epi
#-- 1.x
from obspy.geodetics import kilometer2degrees as k2d



# Main function
def main():

  #-- Run Input Parser
  (opts,indb) = SplitPy.utils.get_options_prep_offline()

  #-- Set Global Search Variables
  cf.maxdt = opts.maxdt #   maxdt is Max delay time
  cf.ddt = opts.ddt     #   ddt is time increment
  cf.dphi = opts.dphi   #   dphi is angle increment


  #-- Load Database
  db=StDb.io.load_db(fname=indb)

  #-- construct station key loop
  allkeys=db.keys()
  allkeys.sort()

  #-- Extract key subset
  if len(opts.stkeys)>0:
    stkeys=[]
    for skey in opts.stkeys:
      stkeys.extend([s for s in allkeys if skey in s] )
  else:
      stkeys=db.keys()
      stkeys.sort()

  #-- Initialize Taup Model
  TTmod=TauPyModel(model='iasp91')

  # Output directory
  if opts.startT is not None and opts.endT is not None:
    Dtrange="_{0:s}-{1:s}".format(opts.startT.strftime("%Y%m%d"),opts.endT.strftime("%Y%m%d"))
  elif opts.startT is None and opts.endT is not None:
    Dtrange="_-{0:s}".format(opts.endT.strftime("%Y%m%d"))
  elif opts.startT is not None and opts.endT is None:
    Dtrange="_{0:s}-".format(opts.startT.strftime("%Y%m%d"))
  else:
    Dtrange=""

  if opts.maxmag is None:
    mxmag="+"
  else:
    mxmag="-{0:.1f}".format(opts.maxmag)

  dtdir = "{0:s}{1:s}_D{2:.0f}-{3:.0f}_M{4:.1f}{5:s}_S{6:.1f}+".format(opts.datadir, Dtrange, \
          opts.mindist, opts.maxdist, opts.minmag, mxmag, opts.msnr)
  if not os.path.isdir(dtdir): os.makedirs(dtdir)

  #-- Loop over station keys
  for ik in range(len(stkeys)):

    # Extract station information from dictionary
    sta = db[stkeys[ik]]

    # Create Station Data Folder
    outdir = dtdir + "/" + stkeys[ik]

    # Establish client
    if len(opts.UserAuth)==0:
      client = Client(opts.Server)
    else:
      client=Client(opts.Server,user=opts.UserAuth[0],password=opts.UserAuth[1])

    # Get catalogue search start time
    if opts.startT is None:
      evst=sta.dstart
    else:
      evst=opts.startT
    # Get catalogue search end time
    if opts.endT is None:
      evet=sta.dend
    else:
      evet=opts.endT
    if evst > sta.dend or evet<sta.dstart:
      continue

    #- Temporary print locations
    tlocs=sta.location
    if len(tlocs)==0: tlocs=['']
    for il in range(0,len(tlocs)):
      if len(tlocs[il])==0: tlocs[il]="--"
    sta.location=tlocs

    #-- Update Display
    print (" ")
    print (" ")
    print ("|===============================================|")
    print ("|===============================================|")
    print ("|                   {0:>8s}                    |".format(sta.station))
    print ("|===============================================|")
    print ("|===============================================|")
    print ("|  Station: {0:>2s}.{1:5s}                            |".format(sta.network,sta.station))
    print ("|      Channel: {0:2s}; Locations: {1:15s}  |".format(sta.channel,",".join(tlocs)))
    print ("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(sta.longitude,sta.latitude))
    print ("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
    print ("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
    print ("|-----------------------------------------------|")
    print ("| Searching Possible events:                    |")
    print ("|   Start: {0:19s}                  |".format(evst.strftime("%Y-%m-%d %H:%M:%S")))
    print ("|   End:   {0:19s}                  |".format(evet.strftime("%Y-%m-%d %H:%M:%S")))
    if opts.maxmag is None:
      print ("|   Mag:   >{0:3.1f}                                 |".format(opts.minmag))
    else:
      print ("|   Mag:   {0:3.1f} - {1:3.1f}                            |".format(opts.minmag,opts.maxmag))
    print ("| ...                                           |")


    # Get catalogue using deployment start and end
    cat = client.get_events(starttime=evst, endtime=evet,    \
          minmagnitude=opts.minmag, maxmagnitude=opts.maxmag)

    #-- Total number of events in Catalogue
    nevK=0
    nevtT=len(cat)
    print ("|  Found {0:5d} possible events                  |".format(nevtT))

    #-- Get Local Data Availabilty
    if len(opts.localdata)>0:
      print ("|-----------------------------------------------|")
      print ("| Cataloging Local Data...                      |")
      if opts.useNet:
        stalcllist=SplitPy.io.list_local_data_stn(lcldrs=opts.localdata, sta=sta.station, \
                                       net=sta.network, altnet=sta.altnet)
        print ("|   {0:>2s}.{1:5s}: {2:6d} files                      |".format(sta.network,sta.station,len(stalcllist)))
      else:
        stalcllist=SplitPy.io.list_local_data_stn(lcldrs=opts.localdata, sta=sta.station)
        print ("|   {0:5s}: {1:6d} files                         |".format(sta.station,len(stalcllist)))
    else:
        stalcllist=[]
    print ("|===============================================|")


    ievs=range(0,nevtT)

    # Read through catalogue
    for iev in ievs:
      ev = cat[iev]

      # Extract time, coordinates and depth of events
      time = ev.origins[0].time
      lat = ev.origins[0].latitude
      lon = ev.origins[0].longitude
      dep = ev.origins[0].depth
      mag = ev.magnitudes[0].mag

      # Define time stamp
      yr = str(time.year).zfill(4)
      jd = str(time.julday).zfill(3)
      hr = str(time.hour).zfill(2)

      # Calculate epicentral distance
      epi_dist, az, baz = epi(lat, lon, sta.stla, sta.stlo)
      epi_dist /= 1000
      gac = k2d(epi_dist)

      # If distance between 85 and 120 deg:
      if (gac>opts.mindist and gac<opts.maxdist):

        #-- Display Event Info
        nevK=nevK+1
        inum=nevtT-iev+1
        print (" ")
        print ("****************************************************")
        print ("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(nevK,inum,nevtT,time.strftime("%Y%m%d_%H%M%S")))
        print ("*   Origin Time: "+time.strftime("%Y-%m-%d %H:%M:%S"))
        print ("*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(lat,lon))
        print ("*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(dep/1000.,mag))
        print ("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; {3:6.2f}; {4:6.2f}".format(sta.station,epi_dist,gac,baz,az))

        # Get Travel times (Careful: here dep is in meters)
        tt = TTmod.get_travel_times(distance_in_degree=gac, source_depth_in_km=dep/1000.)

        # Loop over all times in tt
        for t in tt:
          # Extract time of SKS arrival
          if t.name=='SKS':
            # Extract time, phase name and slowness
            ts = t.time
            ph = t.name
            slow = t.ray_param_sec_degree/111.
            inc = np.arcsin(opts.vp * slow)*180./np.pi
            # Break out of loop
            break

        # Define start and end times for requests
        tstart = time + ts - opts.dts
        tend = time + ts + opts.dts

        # Get waveforms
        print ("* Requesting Waveforms: ")
        print ("*    Startime: "+tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print ("*    Endtime:  "+tend.strftime("%Y-%m-%d %H:%M:%S"))
        err,trN,trE,trZ = SplitPy.utils.get_data_NEZ(client=client, sta=sta, start=tstart, \
                            end=tend, stdata=stalcllist, ndval=opts.ndval)
        if err: continue

        # Rotate from ZEN to LQT (Longitudinal, Radial, Transverse)
        trL, trQ, trT = SplitPy.calc.rotate_ZEN_LQT(trZ, trE, trN, inc/2., baz)

        # Calculate snr over XX seconds
        snrq = SplitPy.calc.snr(trQ, time, ts, ts + 20.)

        # Create instance of event
        sp_meta = SplitPy.classes.meta_class(snrq, sta.station, time, dep, mag, \
                                             lon, lat, gac, baz, inc)

        #-- Make sure no processing happens for NaNs
        if np.isnan(sp_meta.snr):
          print("* SNR NaN...Skipping")
          print ("****************************************************")
          continue
        #-- SNR below threshold
        elif sp_meta.snr < opts.msnr:
          print("* SNR Failed: {0:.2f} < {1:.2f}...Skipping".format( \
                                                                         sp_meta.snr,opts.msnr))
          print ("****************************************************")

        # If SNR is higher than threshold
        elif sp_meta.snr >= opts.msnr:
          print ("* SNR Passed: {0:4.2f} >= {1:3.1f}".format(sp_meta.snr, opts.msnr))

          #-- Create Event Folder
          evtdir=outdir + "/" + time.strftime("%Y%m%d_%H%M%S")

          #-- Create Folder
          if not os.path.isdir(evtdir): os.makedirs(evtdir)

          #-- Event Data
          dill.dump([ev,tstart,tend,tt],open(evtdir+"/Event_Data.pkl","wb"))

          #-- Station Data
          dill.dump(sta,open(evtdir+"/Station_Data.pkl","wb"))

	  #-- Trace Filenames
          trpref=time.strftime("%Y%m%d_%H%M%S")+"_"+sta.network+"."+sta.station

          #-- Raw Trace files
          dill.dump([trN, trE, trZ],open(evtdir+"/NEZ_Data.pkl","wb"))
	  trN.write(os.path.join(evtdir,trpref+".N.mseed"),format='MSEED')
	  trE.write(os.path.join(evtdir,trpref+".E.mseed"),format='MSEED')
          trZ.write(os.path.join(evtdir,trpref+".Z.mseed"),format='MSEED')

          #-- LQT Traces
          dill.dump([trL,trQ,trT],open(evtdir+"/LQT_Data.pkl","wb"))
	  trL.write(os.path.join(evtdir,trpref+".L.mseed"),format='MSEED')
	  trQ.write(os.path.join(evtdir,trpref+".Q.mseed"),format='MSEED')
          trT.write(os.path.join(evtdir,trpref+".T.mseed"),format='MSEED')

          #-- Update
          print ("* Wrote Output Files to: ")
          print ("*     "+evtdir)
          print ("****************************************************")




###############################
# Choose one station to process

if __name__ == "__main__":

  #-- Run main program
  main()

###################
