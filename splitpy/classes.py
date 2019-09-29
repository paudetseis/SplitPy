'''
SUBMODULE split_classes.py

Module containing class definitions for shear-wave splitting analysis

'''
import numpy as np
import splitpy
from splitpy import conf as cf
from obspy.core import Trace
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
# plt.ion()


class Split(object):
    """Split class 

    Attributes:
        Emat:   Error minimization matrix
        trQ_c:  Corrected radial (Q) component
        trT_c:  Corrected transverse (T) component
        trFast: Corrected Fast component
        trSlow: Corrected Slow component
        phi:    Azimuth of fast axis (deg)
        dtt:    Delay time between fast and slow axes (sec)
        phi_min:    Azimuth of... (deg)

    """

    def __init__(self, sta):
        self.sta = sta


    def add_event(self, event):

        time = event.origins[0].time
        dep = event.origins[0].depth
        lon = event.origins[0].longitude
        lat = event.origins[0].latitude

        # Problem with mag
        mag = event.magnitudes[0].mag
        if mag is None: mag = -9.

        # Calculate epicentral distance
        epi_dist, az, baz = epi(lat, lon, self.sta.stla, self.sta.stlo)
        epi_dist /= 1000
        gac = k2d(epi_dist)

        # Store as dictionary attribute
        self.meta = {"time": time, "dep": dep, "lon": lon, "lat": lat, \
                    "mag": mag, "gac": gac, "epi_dist": epi_dist, "baz": baz, "az": az}


    def add_phase(self, t, vp):

        # Store as attributes
        self.ts = t.time
        self.ph = t.name
        self.meta["slow"] = t.ray_param_sec_degree/111.
        self.meta["inc"] = np.arcsin(vp*self.meta["slow"])*180./np.pi


    def get_data_NEZ(self, client, dts, stdata, ndval):

        # Define start and end times for requests
        tstart = self.meta["time"] + self.ts - dts
        tend = self.meta["time"] + self.ts + dts

        # Get waveforms
        print ("* Requesting Waveforms: ")
        print ("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print ("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        err, trN, trE, trZ = splitpy.utils.get_data_NEZ(client=client, \
            sta=self.sta, start=tstart, \
            end=tend, stdata=stdata, ndval=ndval)

        # Store as attributes with traces in dictionay
        self.err = err
        self.data = {"trN": trN, "trE": trE, "trZ": trZ}


    def rotate_ZEN_LQT(self):
        """rotate_ZEN_LQT: 
        
        Rotates 3-component seismograms from vertical,
        east and nort to longitudinal, radial and tangential components
        of motion

        Returns:
            trL:    Trace of longitudinal component of motion
            trQ:    Trace of radial component of motion
            trT:    Trace of tangential component of motion

        """

        inc = self.meta["inc"]*np.pi/180.
        baz = self.meta["baz"]*np.pi/180.

        M = np.zeros((3,3))
        M[0,0] = np.cos(inc)
        M[0,1] = -np.sin(inc) * np.sin(baz)
        M[0,2] = -np.sin(inc) * np.cos(baz)
        M[1,0] = np.sin(inc)
        M[1,1] = np.cos(inc) * np.sin(baz)
        M[1,2] = np.cos(inc) * np.cos(baz)
        M[2,0] = 0.
        M[2,1] = -np.cos(baz)
        M[2,2] = np.sin(baz)

        # Perform 3-D rotation
        LQT = np.dot(np.array(M), np.array(\
            [self.data["trZ"].data, self.data["trE"].data, self.data["trN"].data]))

        # Store into traces and add as new items in attribute dictionary
        self.data["trL"] = Trace(data=LQT[0], header=self.data["trZ"].stats)
        self.data["trQ"] = Trace(data=LQT[1], header=self.data["trN"].stats)
        self.data["trT"] = Trace(data=LQT[2], header=self.data["trE"].stats)


    def calc_snrq(self, t1=None, dt=30.):

        if t1 is None:
            t1=self.meta["time"] + self.ts - 5.
        self.snrq = _calc_snr(self.data["trQ"], t1=t1, dt=dt)


    def calc_snrt(self, t1=None, dt=30.):

        if t1 is None:
            t1=self.meta["time"] + self.ts - 5.
        self.snrt = _calc_snr(self.data["trT"], t1=t1, dt=dt)


    def analyze(self, t1=None, t2=None):

        if t1 is None and t2 is None:
            t1 = self.meta["time"] + self.ts - 5.
            t2 = self.meta["time"] + self.ts + 25.

        # Calculate Silver and Chan splitting estimate
        print ("* --> Calculating Rotation-Correlation (RC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                        splitpy.calc.split_RotCorr(self.data["trQ"], self.data["trT"], \
                        self.meta["baz"], t1, t2)

        # Calculate error
        edtt, ephi, errc = splitpy.calc.split_errorRC(self.data["trT"], \
                        t1, t2, 0.05, Emat)

        # Store dictionary as attribute
        self.RC_res = {
            "Emat": Emat, "trQ_c": trQ_c, "trT_c": trT_c, "trFast": trFast,
            "trSlow": trSlow, "phi": phi, "dtt": dtt, "phi_min": phi_min, \
            "edtt": edtt, "ephi": ephi, "errc": errc
            }

        # Calculate Silver and Chan splitting estimate
        print ("* --> Calculating Silver-Chan (SC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
                        splitpy.calc.split_SilverChan(self.data["trQ"], self.data["trT"], \
                        self.meta["baz"], t1, t2)
        
        # Calculate errors
        edtt, ephi, errc = splitpy.calc.split_errorSC(self.data["trT"], \
                        t1, t2, 0.05, Emat)

        self.SC_res = {
            "Emat": Emat, "trQ_c": trQ_c, "trT_c": trT_c, "trFast": trFast,
            "trSlow": trSlow, "phi": phi, "dtt": dtt, "phi_min": phi_min, \
            "edtt": edtt, "ephi": ephi, "errc": errc
            }


    def is_null(self, snrTlim=3., ds=-1):
        """isnull:

        Determines if splitting result is a Null result

        Returns:
            null:   Boolean for Null

        """

        self.null = False

        # Calculate Angular Difference for Null Measurement
        dphi = max(abs(self.RC_res["phi"] - self.SC_res["phi"]), \
            abs(self.SC_res["phi"] - self.RC_res["phi"]))
        if dphi > 90.: dphi = 180. - dphi


        # Summarize Null Measurement
        if ds>=0:
          print("*" + " "*ds + "Null Classification: ")
          if self.snrt < snrTlim:
            print("*" + " "*ds + "  SNR T Fail: {0:.2f} < {1:.2f}".format(self.snrt, snrTlim))
          else:
            print("*" + " "*ds + "  SNR T Pass: {0:.2f} > {1:.2f}".format(self.snrt, snrTlim))
          if 22. < dphi < 68.:
            print("*" + " "*ds + "  dPhi Fail: {0:.2f} within 22. < X < 68.".format(dphi))
          else:
            print("*" + " "*ds + "  dPhi Pass:  {0:.2f} outside 22. < X < 68.".format(dphi))
       
        # Check snr on tangential component
        if self.snrt < snrTlim: 
            self.null = True

        if 22. < dphi < 68.: 
            self.null = True

    def get_quality(self, ds):
        """quality:

        Determines the quality of the estimate (either Null or non-Null)
        based on ratio of delay times and difference in fast axis directions
        between Rotation-Correlation and Silver-Chan methods

        Returns:
            quality:    String representing quality of estimate (Good, Fair, Poor)

        """

        # Ratio of delay times
        rho = self.RC_res["dtt"]/self.SC_res["dtt"]
       
        # Test based on difference in fast axis directions
        dphi = max(abs(self.RC_res["phi"] - self.SC_res["phi"]), \
            abs(self.SC_res["phi"] - self.RC_res["phi"]))
        if dphi > 90.: dphi = 180. - dphi

        # If estimate is Null
        if self.null:
            if rho < 0.2 and (37. < dphi < 53): 
                self.quality = 'Good'
            elif rho < 0.3 and (32 < dphi < 58):
                self.quality = 'Fair'
            else:
                self.quality = 'Poor'
            if ds >= 0:
              print("*" + " "*ds + "Quality Estimate: Null -- {0:s}".format(self.quality))
              print("*" + " "*ds + "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
              print("*" + " "*ds + "      Good: rho < 0.2  &&  37 < dphi < 53")
              print("*" + " "*ds + "      Fair: rho < 0.3  &&  32 < dphi < 58")
              print("*" + " "*ds + "      Poor: rho > 0.3  &&  dphi < 32 | dphi > 58")

        # If estimate is non-Null
        else:
            if (0.8 < rho < 1.1) and dphi < 8.:
                self.quality = 'Good'
            elif (0.7 < rho < 1.2) and dphi < 15.:
                self.quality = 'Fair'
            else:
                self.quality = 'Poor'
            if ds >=0:
              print("*" + " "*ds + "Quality Estimate: Non-Null -- {0:s}".format(self.quality))
              print("*" + " "*ds + "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
              print("*" + " "*ds + "      Good: 0.8 < rho < 1.1  &&  dphi < 8")
              print("*" + " "*ds + "      Fair: 0.7 < rho < 1.2  &&  dphi < 15")
              print("*" + " "*ds + "      Poor: rho < 0.7 | rho > 1.3 &&  dphi > 15")


    def display_results(self, ds=0):
        """Prints out best fitting results"""
        print(" "*ds + ' ======= Best-fit splitting results ========')
        print() 
        print(" "*ds + ' Best fit values: RC method')
        print(" "*ds + ' Phi = ' + \
                str("{:3d}").format(int(self.RC_res["phi"])) + \
                ' degrees +/- ' + str("{:2d}").format(int(self.RC_res["ephi"])))
        print(" "*ds + ' dt = ' + str("{:.1f}").format(self.RC_res["dtt"]) + \
                ' seconds +/- ' + str("{:.1f}").format(self.RC_res["edtt"]))
        print()
        print(" "*ds + ' Best fit values: SC method')
        print(" "*ds + ' Phi = '+\
                str("{:3d}").format(int(self.SC_res["phi"])) + \
                ' degrees +/- ' + str("{:2d}").format(int(self.SC_res["ephi"])))
        print(" "*ds + ' dt = ' + str("{:.1f}").format(self.SC_res["dtt"]) + \
                ' seconds +/- ' + str("{:.1f}").format(self.SC_res["edtt"]))
        print()
        
    def display_meta(self,  ds=0):
        """Prints out content of meta-data"""
        print(" "*ds + ' ======= Meta data ========')
        print()
        print(" "*ds + 'SNR (dB):            ' + str("{:.0f}").format(self.snrq))
        print(" "*ds + 'Station:             ' + self.sta.station)
        print(" "*ds + 'Time:                ' + str(self.meta["time"]))
        print(" "*ds + 'Event depth (km):    ' + str("{:.0f}").format(self.meta["dep"]/1000.))
        print(" "*ds + 'Magnitude (Mw):      ' + str("{:.1f}").format(self.meta["mag"]))
        print(" "*ds + 'Longitude (deg):     ' + str("{:.2f}").format(self.meta["lon"]))
        print(" "*ds + 'Latitude (deg):      ' + str("{:.2f}").format(self.meta["lat"]))
        print(" "*ds + 'GAC (deg):           ' + str("{:.2f}").format(self.meta["gac"]))
        print(" "*ds + 'Backazimuth deg):    ' + str("{:.2f}").format(self.meta["baz"]))
        print(" "*ds + 'Incidence(deg):      ' + str("{:.2f}").format(self.meta["inc"]))
        print()

    def display_null_quality(self, ds=0):
        """Prints out results"""
        print(" "*ds + ' ======= Nulls and quality ========')
        print()
        print(" "*ds + ' Is Null?     ', self.null)
        print(" "*ds + ' Quality:     ', self.quality)
        print()
        
    def save(self, file):

        import pickle
        output = open(file, 'wb')
        pickle.dump(self.meta, output)
        pickle.dump(self.RC_res, output)
        pickle.dump(self.SC_res, output)
        output.close()


class SeisPlot(object):

    def __init__(self, split):
        """
            function [fig, ax1, ax2, ax3] \
                    = init_fig_LQT (fp=list, st1=str)
       
            Function to initalize the plotting window, and return handles to 
            the different plotting axes
        """

        def init_pickw(ax, title, ylab):
            """init_pickw: Initialize picking window"""

            ax.clear()
            ax.set_title(title)
            ax.set_ylabel(ylab)
            ax.grid(which='major',axis='x')
            ax.grid(which='both',axis='y')
            ax.set_ylim((-1.1, 1.1))
            return ax

        if plt.fignum_exists(1):
            plt.figure(1).clf()

        # Store split as attribute
        self.split = split

        gs = gspec.GridSpec(24,1)

        # Figure handle
        fig = plt.figure(num=1,facecolor='w')
        
        # L component seismogram
        ax1 = fig.add_subplot(gs[0:7])   #0:3
        ax1 = init_pickw(ax=ax1, title=self.split.sta.station, ylab='L')

        # Q component seismogram
        ax2 = fig.add_subplot(gs[8:15])   #3:7
        ax2 = init_pickw(ax=ax2, title='', ylab='Q')

        # T component seismogram
        ax3 = fig.add_subplot(gs[16:23])  #8:11
        ax3 = init_pickw(ax=ax3, title='', ylab='T')
        ax3.set_xlabel('Time (sec)')
        fp = [fig, ax1, ax2, ax3]

        # Ensure Figure is open
        fp[0].show()

        self.fp = fp


    def plot_LQT_phases(self, tt, dts, t1=None, t2=None):
        """plot_LQT_phases:

        Plots rotated three-components of motion for picking

        """

        if t1 is None and t2 is None:
            t1 = self.split.meta["time"] + self.split.ts - 5.
            t2 = self.split.meta["time"] + self.split.ts + 25.

        # Define time axis
        taxis = np.arange(self.split.data["trL"].stats.npts)/self.split.data["trL"].stats.sampling_rate - dts
        tstart = self.split.data["trL"].stats.starttime

        # Set uniform vertical scale
        maxL = np.abs(self.split.data["trL"].data).max()
        maxQ = np.abs(self.split.data["trQ"].data).max()
        maxT = np.abs(self.split.data["trT"].data).max()
        max = np.amax([maxL, maxQ, maxT])

        # Plot traces
        self.fp[1].plot(taxis, self.split.data["trL"].data/max)
        self.fp[2].plot(taxis, self.split.data["trQ"].data/max)
        self.fp[3].plot(taxis, self.split.data["trT"].data/max)

        self.fp[1].set_xlim((taxis[0], taxis[-1]))
        self.fp[2].set_xlim((taxis[0], taxis[-1]))
        self.fp[3].set_xlim((taxis[0], taxis[-1]))

        # Add vertical lines for picked times
        ll = list(range(6))
        ll[0] = self.fp[1].axvline(t1 - tstart - dts, color='r')
        ll[1] = self.fp[1].axvline(t2 - tstart - dts, color='r')
        ll[2] = self.fp[2].axvline(t1 - tstart - dts, color='r')
        ll[3] = self.fp[2].axvline(t2 - tstart - dts, color='r')
        ll[4] = self.fp[3].axvline(t1 - tstart - dts, color='r')
        ll[5] = self.fp[3].axvline(t2 - tstart - dts, color='r')

        # Store lines as attribute
        self.ll = ll

        # Add vertical lines and text for various phases
        for t in tt:

            name = t.name
            time = t.time
            if not name[0] == 'S': continue
            if name == 'S':
                self.fp[1].axvline(time - self.split.ts, color='k') 
                self.fp[1].text(time - self.split.ts + 5., -1., 'S', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k') 
                self.fp[2].text(time - self.split.ts + 5., -1., 'S', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k') 
                self.fp[3].text(time - self.split.ts + 5., -1., 'S', rotation=90, ha='center', va='bottom')
            elif name == 'SKS':
                self.fp[1].axvline(time - self.split.ts, color='k') 
                self.fp[1].text(time - self.split.ts + 5., -1., 'SKS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k') 
                self.fp[2].text(time - self.split.ts + 5., -1., 'SKS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k') 
                self.fp[3].text(time - self.split.ts + 5., -1., 'SKS', rotation=90, ha='center', va='bottom')
            elif name == 'SKKS':
                self.fp[1].axvline(time - self.split.ts, color='k') 
                self.fp[1].text(time - self.split.ts + 5., -1., 'SKKS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k') 
                self.fp[2].text(time - self.split.ts + 5., -1., 'SKKS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k') 
                self.fp[3].text(time - self.split.ts + 5., -1., 'SKKS', rotation=90, ha='center', va='bottom')
            elif name == 'ScS':
                self.fp[1].axvline(time - self.split.ts, color='k') 
                self.fp[1].text(time - self.split.ts + 5., -1., 'ScS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k') 
                self.fp[2].text(time - self.split.ts + 5., -1., 'ScS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k') 
                self.fp[3].text(time - self.split.ts + 5., -1., 'ScS', rotation=90, ha='center', va='bottom')

        # Update plot
        self.fp[0].canvas.draw()


    def update_LQT(self, tp1, tp2):
        """update_LQT:

        Updates LQT figure with newly picked times

        """

        # Remove lines from previously picked times
        self.ll[0].remove()
        self.ll[1].remove()
        self.ll[2].remove()
        self.ll[3].remove()
        self.ll[4].remove()
        self.ll[5].remove()

        # Add new vertical lines with new picks
        self.ll[0] = self.fp[1].axvline(tp1, color='r')
        self.ll[1] = self.fp[1].axvline(tp2, color='r')
        self.ll[2] = self.fp[2].axvline(tp1, color='r')
        self.ll[3] = self.fp[2].axvline(tp2, color='r')
        self.ll[4] = self.fp[3].axvline(tp1, color='r')
        self.ll[5] = self.fp[3].axvline(tp2, color='r')

        # Update figure
        self.fp[0].canvas.draw()


class DiagPlot(object):

    def __init__(self, split):

        """
            function [fig, ax0, axt, axRC1, axRC2, axRC3, axRC4, \
                    axSC1, axSC2, axSC3, axSC4] = init_fig (fd=list)
       
            Function to initalize the diagnostic window, and return handles to 
            the different plotting axes
        """

        def init_splitw(ax, title):
            """init_splitw: Initialize window for splitting results"""

            ax.clear()
            ax.set_ylim((-1.2, 1.2))
            #ax.set_xlabel('Time (sec)')
            ax.set_title(title, fontsize=12)
            return ax

        def init_pmotion(ax):
            """init_pmotion: Initialize window for particle motion"""

            ax.clear()
            ax.set_ylim((-1.2, 1.2))
            ax.set_xlim((-1.2, 1.2))
            ax.set_yticks(())
            ax.set_xticks(())
            ax.set_xlabel(r'W $\longleftrightarrow$ E')
            ax.set_ylabel(r'S $\longleftrightarrow$ N')
            ax.set_title('Particle Motion', fontsize=12)
            return ax

        def init_emap(ax, title):
            """init_emap: Initialize window for matrix minimization"""

            ax.clear()
            ax.set_title(title, fontsize=12)
            ax.set_ylabel(r'$\phi$ (deg)', labelpad=-10)
            ax.set_xlabel(r'$\delta t$ (sec)', labelpad=-2)
            ax.set_xticks([0,1,2,3,4])
            return ax

        if plt.fignum_exists(2):
            plt.figure(2).clf()

        self.split = split

        # Figure handle
        fig = plt.figure(num=2, figsize=(10,7), facecolor='w')
    
        # Q, T component seismograms
        ax0 = fig.add_axes([0.05, 0.73, 0.2, 0.2])
        ax0 = init_splitw(ax=ax0, title='Q, T')

        # Text box
        axt = fig.add_axes([0.45, 0.73, 0.3, 0.2])
        axt.axis('off')

        #
        # Corrected Fast, Slow window for Rotation-Correlation
        axRC1 = fig.add_axes([0.05, 0.4, 0.2, 0.25])
        axRC1 = init_splitw(ax=axRC1, title='Corrected Fast, Slow')

        # Corrected Q, T window
        axRC2 = fig.add_axes([0.3, 0.4, 0.2, 0.25])
        axRC2 = init_splitw(ax=axRC2, title='Corrected Q, T')

        # Particle motion
        axRC3 = fig.add_axes([0.5375, 0.4, 0.175, 0.25])
        axRC3 = init_pmotion(ax=axRC3)

        # Energy map
        axRC4 = fig.add_axes([0.775, 0.4, 0.2, 0.25])
        axRC4 = init_emap(ax=axRC4, title='Map of correlation coeff')

        # Corrected Fast, Slow window for Silver-Chan
        axSC1 = fig.add_axes([0.05, 0.07, 0.2, 0.25])
        axSC1 = init_splitw(ax=axSC1, title='Corrected Fast, Slow')

        # Corrected Q, T window
        axSC2 = fig.add_axes([0.3, 0.07, 0.2, 0.25])
        axSC2 = init_splitw(ax=axSC2, title='Corrected Q, T')

        # Particle motion
        axSC3 = fig.add_axes([0.5375, 0.07, 0.175, 0.25])
        axSC3 = init_pmotion(ax=axSC3)

        # Energy map
        axSC4 = fig.add_axes([0.775, 0.07, 0.2, 0.25])
        axSC4 = init_emap(ax=axSC4, title='Energy map of T')

        fd = [fig, ax0, axt, axRC1, axRC2, axRC3, axRC4,\
                axSC1, axSC2, axSC3, axSC4]

        # Make sure figure is open
        fd[0].show()

        self.fd = fd


    def plot_diagnostic(self, t1=None, t2=None):
        """plot_diagnostic:

        Plots diagnostic window with estimates from both RC and SC methods

        """

        if t1 is None and t2 is None:
            t1 = self.split.meta["time"] + self.split.ts - 5.
            t2 = self.split.meta["time"] + self.split.ts + 25.

        def rot3D(inc, baz):
            """rot3D: 
            
            Defines rotation matrix from incidence angle and back-azimuth

            Returns: 
                M:  Rotation matrix

            """

            # Angles in radians
            inc = inc/180.*np.pi
            baz = baz/180.*np.pi

            # Define rotation matrix
            M = [[np.cos(inc), -np.sin(inc)*np.sin(baz), -np.sin(inc)*np.cos(baz)],
                    [np.sin(inc), np.cos(inc)*np.sin(baz), np.cos(inc)*np.cos(baz)],
                    [0., -np.cos(baz), np.sin(baz)]]

            return M

        # Copy traces to avoid overridding 
        trL_tmp = self.split.data["trL"].copy()
        trL_tmp.trim(t1, t2)
        trQ_tmp = self.split.data["trQ"].copy()
        trQ_tmp.trim(t1, t2)
        trT_tmp = self.split.data["trT"].copy()
        trT_tmp.trim(t1, t2)
        trE_tmp = self.split.data["trE"].copy()
        trE_tmp.trim(t1, t2)
        trN_tmp = self.split.data["trN"].copy()
        trN_tmp.trim(t1, t2)

        # Rotate seismograms for plots
        M = rot3D(self.split.meta["inc"], self.split.meta["baz"])
        ZEN_RC = np.dot(np.transpose(M),
                [trL_tmp.data, self.split.RC_res["trQ_c"].data, self.split.RC_res["trT_c"].data])
        E_RC = ZEN_RC[1,:]
        N_RC = ZEN_RC[2,:]
        ZEN_SC = np.dot(np.transpose(M),
                [trL_tmp.data, self.split.SC_res["trQ_c"].data, self.split.SC_res["trT_c"].data])
        E_SC = ZEN_SC[1,:]
        N_SC = ZEN_SC[2,:]

        # First plot - original time series
        taxis = np.arange(trQ_tmp.stats.npts)/trQ_tmp.stats.sampling_rate
        max1 = np.abs(trQ_tmp.data).max()
        max2 = np.abs(trT_tmp.data).max()
        max = np.amax([max1, max2])

        self.fd[1].plot(taxis, trQ_tmp.data/max, 'b--')
        self.fd[1].plot(taxis, trT_tmp.data/max, 'r')

        # Text box
        self.fd[2].text(0.5, 0.9, 'Event: '+self.split.meta["time"].ctime()+'     '+\
            str(self.split.meta["lat"])+'N  '+str(self.split.meta["lon"])+'E   '+\
            str(np.int(self.split.meta["dep"]/1000.))+'km   '+'Mw='+\
            str(self.split.meta["mag"]), horizontalalignment='center')
        self.fd[2].text(0.5, 0.7, 'Station: '+self.split.sta.station+'   Backazimuth: '+\
            str("{:.2f}").format(self.split.meta["baz"])+'   Distance: '+\
            str("{:.2f}").format(self.split.meta["gac"]), horizontalalignment='center')
        self.fd[2].text(0.5, 0.5, r'Best fit RC values: $\phi$='+\
            str(np.int(self.split.RC_res["phi"]))+r'$\pm$'+\
            str("{:.2f}").format(self.split.RC_res["ephi"])+r'   $\delta t$='+\
            str(self.split.RC_res["dtt"])+r'$\pm$'+\
            str("{:.2f}").format(self.split.RC_res["edtt"])+'s', horizontalalignment='center')
        self.fd[2].text(0.5, 0.3, r'Best fit SC values: $\phi$='+\
            str(np.int(self.split.SC_res["phi"]))+r'$\pm$'+\
            str("{:.2f}").format(self.split.SC_res["ephi"])+r'   $\delta t$='+\
            str(self.split.SC_res["dtt"])+r'$\pm$'+\
            str("{:.2f}").format(self.split.SC_res["edtt"])+'s', horizontalalignment='center')
        self.fd[2].text(0.5, 0.1, 'Is Null? '+str(self.split.null)+'    Quality? '+\
            str(self.split.quality), horizontalalignment='center')

        # Rotation-correlation
        # Second plot - corrected Fast and Slow
        sum1 = np.sum(np.abs(self.split.RC_res["trFast"].data - self.split.RC_res["trSlow"].data))
        sum2 = np.sum(np.abs(-self.split.RC_res["trFast"].data - self.split.RC_res["trSlow"].data))
        if sum1 < sum2: 
            sig = 1. 
        else: 
            sig = -1.
        taxis = np.arange(self.split.RC_res["trFast"].stats.npts)/self.split.RC_res["trFast"].stats.sampling_rate
        max1 = np.abs(self.split.RC_res["trFast"].data).max()
        max2 = np.abs(self.split.RC_res["trSlow"].data).max()
        max = np.amax([max1, max2])

        self.fd[3].plot(taxis, self.split.RC_res["trFast"].data/max, 'b--')
        self.fd[3].plot(taxis, sig*self.split.RC_res["trSlow"].data/max, 'r')

        # Third plot - Corrected Q and T
        self.fd[4].plot(taxis, self.split.RC_res["trQ_c"].data/max, 'b--')
        self.fd[4].plot(taxis, self.split.RC_res["trT_c"].data/max, 'r')

        # Fourth plot - particle motion
        self.fd[5].plot(trE_tmp.data/max, trN_tmp.data/max, 'b--')
        self.fd[5].plot(E_RC/max, N_RC/max, 'r')

        # Fifth plot - map of energy
        plt.sca(self.fd[6])
        dt = np.arange(0., cf.maxdt, cf.ddt)
        phi = np.arange(-90., 90., cf.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X,Y = np.meshgrid(dt, phi)

        E2 = np.roll(self.split.RC_res["Emat"], np.int(self.split.RC_res["phi"] - self.split.RC_res["phi_min"]), axis=0)

        Emin = self.split.RC_res["Emat"].min()
        Emax = self.split.RC_res["Emat"].max()
        dE = (Emax - Emin)/16.
        levels = np.arange(Emin, Emax, dE)
        cmap = plt.cm.RdYlBu_r
        cset1 = plt.contour(X, Y, E2, levels, \
                cmap=plt.cm.get_cmap(cmap, len(levels)))

        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        errc = self.split.RC_res["errc"]
        ecset = plt.contour(X, Y, E2, (errc, ), colors='magenta',\
                linewidths=2)

        self.fd[6].axvline(self.split.RC_res["dtt"])
        self.fd[6].axhline(self.split.RC_res["phi"])

        # Silver-Chan
        # Second plot - corrected Fast and Slow
        sum1 = np.sum(np.abs(self.split.SC_res["trFast"].data - self.split.SC_res["trSlow"].data))
        sum2 = np.sum(np.abs(-self.split.SC_res["trFast"].data - self.split.SC_res["trSlow"].data))
        if sum1 < sum2: 
            sig = 1. 
        else: 
            sig = -1.
        taxis = np.arange(self.split.SC_res["trFast"].stats.npts)/self.split.SC_res["trFast"].stats.sampling_rate
        max1 = np.abs(self.split.SC_res["trFast"].data).max()
        max2 = np.abs(self.split.SC_res["trSlow"].data).max()
        max = np.amax([max1, max2])

        self.fd[7].plot(taxis, self.split.SC_res["trFast"].data/max, 'b--')
        self.fd[7].plot(taxis, sig*self.split.SC_res["trSlow"].data/max, 'r')

        # Third plot - Corrected Q and T
        self.fd[8].plot(taxis, self.split.SC_res["trQ_c"].data/max, 'b--')
        self.fd[8].plot(taxis, self.split.SC_res["trT_c"].data/max, 'r')

        # Fourth plot - particle motion
        self.fd[9].plot(trE_tmp.data/max, trN_tmp.data/max, 'b--')
        self.fd[9].plot(E_SC/max, N_SC/max, 'r')

        # Fifth plot - map of energy
        plt.sca(self.fd[10])
        dt = np.arange(0., cf.maxdt, cf.ddt)
        phi = np.arange(-90., 90., cf.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X,Y = np.meshgrid(dt, phi)

        E2 = np.roll(self.split.SC_res["Emat"], np.int(self.split.SC_res["phi"]-self.split.SC_res["phi_min"]),axis=0)

        Emin = self.split.SC_res["Emat"].min()
        Emax = self.split.SC_res["Emat"].max()
        dE = (Emax - Emin)/16.
        levels = np.arange(Emin, Emax, dE)
        cmap = plt.cm.RdYlBu_r
        cset1 = plt.contour(X, Y, E2, levels, \
                cmap=plt.cm.get_cmap(cmap, len(levels)))

        errc = self.split.SC_res["errc"]
        ecset = plt.contour(X, Y, E2, (errc,), colors='magenta',\
                linewidths=2)
        #ecset = plt.contourf(X, Y, E2, (errc,), cmap=plt.cm.gray)

        self.fd[10].axvline(self.split.SC_res["dtt"])
        self.fd[10].axhline(self.split.SC_res["phi"])

        self.fd[0].canvas.draw()

    def save(self, file):

        self.fd[0].savefig(file)


def _calc_snr(tr, t1, dt):

    # Copy Z trace to signal and noise traces
    trSig = tr.copy()
    trNze = tr.copy()

    # Filter between 0.02 and 0.5 (dominant S wave frequencies)
    trSig.filter('bandpass',freqmin=0.02,freqmax=0.5,corners=2,zerophase=True)
    trNze.filter('bandpass',freqmin=0.02,freqmax=0.5,corners=2,zerophase=True)

    # Trim twin seconds around P-wave arrival
    trSig.trim(t1, t1 + dt)
    trNze.trim(t1 - dt, t1)

    # Calculate root mean square (RMS)
    srms = np.sqrt(np.mean(np.square(trSig.data)))
    nrms = np.sqrt(np.mean(np.square(trNze.data)))

    # Calculate signal/noise ratio in dB
    return 10*np.log10(srms*srms/nrms/nrms)

