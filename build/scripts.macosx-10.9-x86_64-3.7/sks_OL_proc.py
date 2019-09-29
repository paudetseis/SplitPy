#!/Users/pascalaudet/anaconda3/envs/split/bin/python


'''
PROGRAM sks_split.py

Calculates single station SKS splitting results.
Station selection is specified by a network and station code.
The data base is provided in stations_db.py as a dictionary.

'''

from sys import exit as sexit

#-- Import Splitting Module
import SplitPy
import SplitPy.classes as spc
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
from os import path,listdir,makedirs

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
  (opts,indr) = SplitPy.utils.get_options_OL_proc()

  #-- Set Global Search Variables
  cf.maxdt = opts.maxdt #   maxdt is Max delay time
  cf.ddt = opts.ddt     #   ddt is time increment
  cf.dphi = opts.dphi   #   dphi is angle increment


  #-- Generate list of available Station Folders
  stnkeys=listdir(indr)

  #-- Sort available stations
  stnkeys.sort()

  #-- Extract key subset
  if len(opts.stkeys)>0:
    stkeys=[]
    for skey in opts.stkeys:
      stkeys.extend([s for s in stnkeys if skey in s] )
  else:
      stkeys=stnkeys


  #-- Loop over station keys
  for ik in range(len(stkeys)):
    #-- Station Key
    stkey=stkeys[ik]

    #-- Get List of events to process
    evs=listdir(path.join(indr,stkey))
    evs.sort()
    nevs=len(evs)


    #-- Update Display
    print (" ")
    print (" ")
    print ("|===============================================|")
    print ("|===============================================|")
    print ("|                   {0:>8s}                    |".format(stkey))
    print ("|===============================================|")
    print ("|===============================================|")
#    print ("|  Station: {0:>2s}.{1:5s}                            |".format(sta.network,sta.station))
#    print ("|      Channel: {0:2s}; Locations: {1:15s}  |".format(sta.channel,",".join(tlocs)))
#    print ("|      Lon: {0:7.2f}; Lat: {1:6.2f}                |".format(sta.longitude,sta.latitude))
#    print ("|      Start time: {0:19s}          |".format(sta.startdate.strftime("%Y-%m-%d %H:%M:%S")))
#    print ("|      End time:   {0:19s}          |".format(sta.enddate.strftime("%Y-%m-%d %H:%M:%S")))
#    print ("|-----------------------------------------------|")
#    print ("| Searching Possible events:                    |")
#    print ("|   Start: {0:19s}                  |".format(evst.strftime("%Y-%m-%d %H:%M:%S")))
#    print ("|   End:   {0:19s}                  |".format(evet.strftime("%Y-%m-%d %H:%M:%S")))
#    if opts.maxmag is None:
#      print ("|   Mag:   >{0:3.1f}                                 |".format(opts.minmag))
#    else:
#      print ("|   Mag:   {0:3.1f} - {1:3.1f}                            |".format(opts.minmag,opts.maxmag))
#    print ("| ...                                           |")


#    #-- Total number of events in Catalogue
#    nevK=0
#    nevtT=len(cat)
    print ("|  Working on {0:5d} saved events                |".format(nevs))

#    #-- Get Local Data Availabilty
#    if len(opts.localdata)>0:
#      print ("|-----------------------------------------------|")
#      print ("| Cataloging Local Data...                      |")
#      if opts.useNet:
#        stalcllist=list_local_data_stn(lcldrs=opts.localdata, sta=sta.station, \
#                                       net=sta.network, altnet=sta.altnet)
#        print ("|   {0:>2s}.{1:5s}: {2:6d} files                      |".format(sta.network,sta.station,len(stalcllist)))
#      else:
#        stalcllist=list_local_data_stn(lcldrs=opts.localdata, sta=sta.station)
#        print ("|   {0:5s}: {1:6d} files                         |".format(sta.station,len(stalcllist)))
#    else:
#        stalcllist=[]
    print ("|===============================================|")      



    # Initialize figure and axis handles
    fp = []
    fd = []

    #-- select order of processing
    if opts.reverse:
      ievs=range(0,nevs)
    else:
      ievs=range(nevs-1,-1,-1)
    nevK=0

    # Read through catalogue
    for iev in ievs:
      #-- Event Name
      evSTR=evs[iev]

      #-- Load Relevant Files
      #-- Event Data
      ll=dill.load(open(path.join(indr,stkey,evSTR,"Event_Data.pkl"),"rb"))
      ev=ll[0]; tstart=ll[1]; tend=ll[2]; tt=ll[3]

      #-- Station Data
      sta=dill.load(open(path.join(indr,stkey,evSTR,"Station_Data.pkl"),"rb"))

      #-- Trace Filenames
      trpref=evSTR+"_"+sta.network+"."+sta.station

      #-- NEZ Trace files
      trs=dill.load(open(path.join(indr,stkey,evSTR,"NEZ_Data.pkl"),"rb"))
      trN=trs[0]; trE=trs[1]; trZ=trs[2]
      #trN=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".N.mseed"))[0]
      #trE=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".E.mseed"))[0]
      #trZ=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".Z.mseed"))[0]

      #-- LQT Trace files
      trs=dill.load(open(path.join(indr,stkey,evSTR,"LQT_Data.pkl"),"rb"))
      trL=trs[0]; trQ=trs[1]; trT=trs[2]
      #trL=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".L.mseed"))[0]
      #trQ=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".Q.mseed"))[0]
      #trT=read(format='MSEED',path.join(indr,stkey,evSTR,trpref+".T.mseed"))[0]

      # Output directory
      outdir = path.join('RESULTS',sta.network+"."+sta.station)
      if not path.isdir(outdir):
        makedirs(outdir)

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
        if opts.reverse:
          inum=iev+1
        else:
          inum=nevs-iev
        print (" ")
        print ("****************************************************")
        print ("* #{0:d} ({1:d}/{2:d}):  {3:13s}".format(nevK,inum,nevs,time.strftime("%Y%m%d_%H%M%S")))
        print ("*   Origin Time: "+time.strftime("%Y-%m-%d %H:%M:%S"))
        print ("*   Lat: {0:6.2f}; Lon: {1:7.2f}".format(lat,lon))
        print ("*   Dep: {0:6.2f}; Mag: {1:3.1f}".format(dep/1000.,mag))
        print ("*     {0:5s} -> Ev: {1:7.2f} km; {2:7.2f} deg; {3:6.2f}; {4:6.2f}".format(sta.station,epi_dist,gac,baz,az))

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

          # Create instance of data
          sp_data = SplitPy.classes.data_class(trE, trN, trZ, trL, trQ, trT)

          # Output file
          outfile = outdir + '/Split' + '.' + sp_meta.sta + '.' + \
                  sp_meta.time.strftime("%Y.%j.%H%M%S") + '.pkl'

          # Default window times
          t1 = sp_meta.time + ts - 5.
          t2 = sp_meta.time + ts + 25.

          # Calculate Rotation-Correlation splitting estimate
          print ("* --> Calculating Rotation-Correlation (RC) Splitting")
          Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                  SplitPy.proc.split_RotCorr(sp_data.trQ, sp_data.trT, \
                  sp_meta.baz, t1, t2)

          # Calculate error
          edtt, ephi, errc = SplitPy.calc.split_errorRC(sp_data.trT, \
                  t1, t2, 0.05, Emat)
          
          # Create instance of split for RC
          sp_RC = SplitPy.classes.split_RC(Emat, trQ_c, trT_c, trFast, trSlow,\
                  phi, dtt, phi_min, edtt, ephi, errc)

          # Calculate Silver and Chan splitting estimate
          print ("* --> Calculating Silver-Chan (SC) Splitting")
          Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                  SplitPy.proc.split_SilverChan(sp_data.trQ, sp_data.trT, \
                  sp_meta.baz, t1, t2)
          
          # Calculate errors
          edtt, ephi, errc = SplitPy.calc.split_errorSC(sp_data.trT, \
                  t1, t2, 0.05, Emat)
          
          # Create instance of split for SC
          sp_SC = SplitPy.classes.split_SC(Emat, trQ_c, trT_c, trFast, trSlow,\
                  phi, dtt, phi_min, edtt, ephi, errc)

          # Determine if Null and quality of estimate
          snrt = SplitPy.calc.snr(sp_data.trT, sp_meta.time, ts - 5., ts + 25.)
          null = SplitPy.calc.isnull(snrt, sp_RC.phi, sp_SC.phi,opts.snrTlim)
          quality = SplitPy.calc.quality(null, sp_RC.dtt, sp_SC.dtt, \
                  sp_RC.phi, sp_SC.phi)

          # Create instance of results
          sp_res = SplitPy.classes.result_class(sp_RC, sp_SC, null, quality)

          # Initialize LQT figure and plot
          fp = SplitPy.gui.init_fig_LQT(fp=fp, st1=sp_meta.sta)
          ll = SplitPy.plot.plot_LQT_phases(fp, sp_data.trL, sp_data.trQ, \
                  sp_data.trT, tt, ts, opts.dts, t1, t2)

          # Initialize diagnostic figure and plot
          fd = SplitPy.gui.init_fig_diagnostic(fd=fd)
          SplitPy.plot.plot_diagnostic(fd, t1, t2, sp_data, sp_meta, sp_res)
          
          # Choose whether to re-pick window times
          iselect = 'a'
          while iselect=='a':

            # Call interactive window for picking
            ans = SplitPy.gui.pick()

            # If user clicks yes:
            if ans:
              print ("* Refine Window in Figure 1")
              # Get clicks from LQT figure
              plt.figure(1)
              xc = fp[0].ginput(2, show_clicks=True)

              # Extract times from clicks
              tp1 = [xx for xx, yy in xc][0]
              tp2 = [xx for xx, yy in xc][1]
              
              if tp2 < tp1:
                tp11=tp2
                tp2=tp1
                tp1=tp11

              # Update LQT figure
              SplitPy.plot.update_LQT(fp, ll, tp1, tp2)

              # Re-define time window
              t1 = sp_meta.time + ts + tp1
              t2 = sp_meta.time + ts + tp2

              # Calculate Rotation-Correlation splitting estimate
              print ("*    --> Re-calculating RC Splitting")
              Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                      SplitPy.proc.split_RotCorr(sp_data.trQ, sp_data.trT, \
                      sp_meta.baz, t1, t2)

              # Calculate error
              edtt, ephi, errc = SplitPy.calc.split_errorRC(sp_data.trT, \
                      t1, t2, 0.05, Emat)

              # Create instance of split for RC
              sp_RC = SplitPy.classes.split_RC(Emat, trQ_c, trT_c, trFast, \
                      trSlow, phi, dtt, phi_min, edtt, ephi, errc)

              # Calculate Silver and Chan splitting estimate
              print ("*    --> Re-calculating SC Splitting")
              Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                      SplitPy.proc.split_SilverChan(sp_data.trQ, sp_data.trT, \
                      sp_meta.baz, t1, t2)

              # Calculate error
              edtt, ephi, errc = SplitPy.calc.split_errorSC(sp_data.trT, \
                      t1, t2, 0.05, Emat)

              # Create instance of split for SC
              sp_SC = SplitPy.classes.split_SC(Emat, trQ_c, trT_c, trFast, \
                      trSlow, phi, dtt, phi_min, edtt, ephi, errc)
              
              # Determine if Null and quality of estimate
              snrt = SplitPy.calc.snr(sp_data.trT, sp_meta.time, \
                      ts + tp1, ts + tp2)
              null = SplitPy.calc.isnull(snrt, sp_RC.phi, sp_SC.phi,opts.snrTlim)
              quality = SplitPy.calc.quality(null, sp_RC.dtt, sp_SC.dtt, \
                      sp_RC.phi, sp_SC.phi)

              # Create instance of results
              sp_res = SplitPy.classes.result_class(sp_RC, sp_SC, \
                      null, quality)

              # Re-initialize diagnostic figure and plot
              fd = SplitPy.gui.init_fig_diagnostic(fd=fd)
              SplitPy.plot.plot_diagnostic(fd, t1, t2, sp_data, sp_meta, sp_res)

            # If user clicks no:
            else:
              iselect = 'c'

              # Call interactive window for decision on estimate
              ans2 = SplitPy.gui.keep()

              # If user keeps estimate:
              if ans2:
                # Print estimates to screen
                sp_meta.display()
                sp_RC.display()
                sp_SC.display()
                sp_res.display()

                output = open(outfile, 'wb')
                pickle.dump(sp_data, output)
                pickle.dump(sp_meta, output)
                pickle.dump(sp_res, output)
                output.close()
                print ("* Estimate Saved")
                print ("****************************************************")
              else:
                print ("* Estimate Discarded ")
                print ("****************************************************")
                continue





def update_stats(rf, stla, stlo, baz, slow, snr, type):
    rf.stats.sac = obspy.core.AttribDict()
    rf.stats.sac.stla = stla
    rf.stats.sac.stlo = stlo
    rf.stats.sac.baz = baz
    rf.stats.sac.user0 = slow
    rf.stats.sac.user1 = snr
    rf.stats.channel = type
    return rf

## Loads station db and builds attribute dict of station stats
#def load_db(fname):
#    db = pickle.load(open(fname, 'rb'))
#    for k, v in db.items():
#        db[k] = meta_data(v)
#    return db

# Attribute dict class
class meta_data(dict):
    def __init__(self, stats):
        self.__dict__ = self
        self.network = stats[0]
        self.station = stats[1]
        self.stla = stats[2]
        self.stlo = stats[3]
        self.cha = stats[4]
        self.dstart = stats[5]
        self.dend = stats[6]

###############################
# Choose one station to process

if __name__ == "__main__":

  #-- Run main program
  main()

###################

