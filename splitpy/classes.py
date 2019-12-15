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
"""

:mod:`~splitpy` defines the following base classes:

- :class:`~splitpy.classes.Split`
- :class:`~splitpy.classes.PickPlot`
- :class:`~splitpy.classes.DiagPlot`

The class :class:`~splitpy.classes.Split` contains attributes
and methods for the analysis of teleseismic shear-wave splitting 
from three-component seismograms. 

The class :class:`~splitpy.classes.PickPlot` contains figure handles 
for a picking window showing the seismograms and the predicted teleseismic
shear-wave phase arrivals. This figure is interactive and new picks can
be generated to refine the analysis.

The class :class:`~splitpy.classes.DiagPlot` contains figure handles
for a diagnostic figure showing a summary of the splitting results. It can
be called after each application of the `split.analyze` method to show 
the summary of the analysis as a figure. This figure can also be saved as
a .png file.

"""

# -*- coding: utf-8 -*-
import numpy as np
import splitpy
from obspy.core import Trace
from obspy.geodetics.base import gps2dist_azimuth as epi
from obspy.geodetics import kilometer2degrees as k2d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec


class Meta(object):
    """
    A Result object contains attributes associated with the result
    of a single splitting analysis. These are equally applicable
    to the RC or SC method - see :func:`~splitpy.classes.analyze`.

    Attributes
    ----------
    time : :class:`~obspy.core.UTCDateTime`
        Origin time of earthquake
    dep : float
        Depth of hypocenter (km)
    lon : float
        Longitude coordinate of epicenter
    lat : float
        Latitude coordinate of epicenter
    mag : float
        Magnitude of earthquake
    gac : float
        Great arc circle between station and epicenter (degrees)
    epi_dist : float
        Epicentral distance between station and epicenter (km)
    baz : float
        Back-azimuth - pointing to earthquake from station (degrees)
    az : float
        Azimuth - pointing to station from earthquake (degrees)
    slow : float
        Horizontal slowness of phase
    inc : float
        Incidence angle of phase at surface

    """

    def __init__(self, time, dep, lon, lat, mag, gac, epi_dist, baz, az):

        self.time = time
        self.dep = dep
        self.lon = lon
        self.lat = lat
        self.mag = mag
        self.gac = gac
        self.epi_dist = epi_dist
        self.baz = baz
        self.az = az
        self.slow = None
        self.inc = None


class Data(object):
    """
    A Data object contains three-component raw (NEZ) and rotated (LQT) 
    waveforms centered on the arrival time of interest.

    Attributes
    ----------

    trN : :class:`~obspy.core.Trace`
        Trace of North component of motion
    trE : :class:`~obspy.core.Trace`
        Trace of East component of motion
    trZ : :class:`~obspy.core.Trace` 
        Trace of Vertical component of motion
    trL : :class:`~obspy.core.Trace`
        Trace of longitudinal component of motion
    trQ : :class:`~obspy.core.Trace`
        Trace of radial component of motion
    trT : :class:`~obspy.core.Trace`
        Trace of tangential component of motion

    """

    def __init__(self, trE, trN, trZ):

        self.trE = trE
        self.trN = trN
        self.trZ = trZ
        self.trL = None
        self.trQ = None
        self.trT = None


class Result(object):
    """
    A Result object contains attributes associated with the result
    of a single splitting analysis. These are equally applicable
    to the RC or SC method - see :func:`~splitpy.classes.analyze`.

    Attributes
    ----------

    Emat: :class:`~numpy.ndarray`
        Error minimization matrix
    trQ_c: :class:`~obspy.core.Trace`
        Corrected radial (Q) component
    trT_c: :class:`~obspy.core.Trace`
        Corrected transverse (T) component
    trFast: :class:`~obspy.core.Trace`
        Corrected Fast component
    trSlow: :class:`~obspy.core.Trace`
        Corrected Slow component
    phi: float
        Azimuth of fast axis (deg)
    dtt: float
        Delay time between fast and slow axes (sec)
    phi_min: float  
        Azimuth used in plotting method
    ephi: float
        Error on azimuth of fast axis (deg)
    edtt: float
        Error on delay time between fast and slow axes (sec)
    errc: float  
        Error contours on `Emat`
    """

    def __init__(self, Emat, trQ_c, trT_c, trFast,
                 trSlow, phi, dtt, phi_min, edtt, ephi, errc):

        self.Emat = Emat
        self.trQ_c = trQ_c
        self.trT_c = trT_c
        self.trFast = trFast
        self.trSlow = trSlow
        self.phi = phi
        self.dtt = dtt
        self.phi_min = phi_min
        self.edtt = edtt
        self.ephi = ephi
        self.errc = errc


class Split(object):
    """
    A Split object contains dictionary attributes that associate
    station information with single event (i.e., earthquake) 
    metadata, corresponding raw and rotated seismograms and 
    splitting results.

    Note
    ----
    The object is initialized with the ``sta`` field only, and
    other attributes are added to the object as the analysis proceeds.

    Attributes
    ----------
    sta : object
        Object containing station information - from :mod:`~stdb` database.
    meta : :class:`~splitpy.classes.Meta`
        Object of metadata information for single event.
    data : :class:`~splitpy.classes.Data`
        Object containing trace data in :class:`~obspy.core.Trace` format
    RC_res : :class:`~splitpy.classes.Result`
        Object containing results of Rotation-Correlation metnod
    SC_res : :class:`~splitpy.classes.Result`
        Object containing results of Silver-Chan metnod
    err : bool
        Whether or not the `get_data_NEZ` successfully retrieved waveforms
    null : bool
        Whether or not estiamate is null result
    quality : str
        Quality of estimate ('good', 'fair', 'poor')
    ts : float
        Predicted travel time for SKS phase
    ph : str
        Name of SKS phase ('SKS')
    snrq : float
        Signal-to-noise ratio for radial (Q) component
    snrt : float
        Signal-to-noise ratio for tangential (T) component
    maxdt : float
        Max delay time between slow and fast axes in search
    ddt : float
        Sampling distance (in sec) for delay time search
    dphi : float
        Sampling distance (in degrees) for azimuth search 

    """

    def __init__(self, sta, maxdt, ddt, dphi):

        self.sta = sta
        self.maxdt = maxdt
        self.ddt = ddt
        self.dphi = dphi

    def add_event(self, event):
        """
        Adds event metadata to Split object. 

        Parameters
        ----------
        event : :class:`~obspy.core.event`
            Event metadata

        Attributes
        ----------
        meta :
            Object containing metadata information

        """

        time = event.origins[0].time
        dep = event.origins[0].depth
        lon = event.origins[0].longitude
        lat = event.origins[0].latitude

        # Problem with mag
        mag = event.magnitudes[0].mag
        if mag is None:
            mag = -9.

        # Calculate epicentral distance
        epi_dist, az, baz = epi(
            lat, lon, self.sta.latitude, self.sta.longitude)
        epi_dist /= 1000
        gac = k2d(epi_dist)

        # Store as object attributes
        self.meta = Meta(time, dep, lon, lat, mag, gac, epi_dist, baz, az)

    def add_phase(self, t, vp):
        """
        Adds phase information for SKS arrival from taup model

        Parameters
        ----------
        t : :class:`~obspy.taup.TaupPyModel`
            Travel-time table metadata

        Attributes
        ----------
        ts : float
            Travel time between earthquake and station (sec)
        ph : str
            Phase name ('SKS')
        meta.slow : float
            Horizontal slowness of phase
        meta.inc : float
            Incidence angle of phase at surface

        """

        # Store as attributes
        self.ts = t.time
        self.ph = t.name
        self.meta.slow = t.ray_param_sec_degree/111.
        self.meta.inc = np.arcsin(vp*self.meta.slow)*180./np.pi

    def add_NEZ(self, stream):
        """
        Adds seismograms from available stream

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Attributes
        ----------
        data : :class:`~splitpy.classes.Data`
            Object containing :class:`obspy.core.Trace` objects

        """

        self.data = Data(stream[0], stream[1], stream[2])

    def add_LQT(self, stream):
        """
        Adds seismograms from available stream

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Attributes
        ----------

        data : :class:`~splitpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        self.data.trL = stream[0]
        self.data.trQ = stream[1]
        self.data.trT = stream[2]

    def get_data_NEZ(self, client, dts, stdata, ndval):
        """
        Downloads seismograms based on event origin time and
        SKS phase arrival.

        Parameters
        ----------
        client : :class:`~obspy.client.fdsn.Client`
            Client object
        dts : float
            Time duration (?)
        stdata : :class:`stdb.classes.StDbElement`
            Station metadata
        ndval : float
            Fill in value for missing data

        Attributes
        ----------

        data : :class:`~splitpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        # Define start and end times for requests
        tstart = self.meta.time + self.ts - dts
        tend = self.meta.time + self.ts + dts

        # Get waveforms
        print("* Requesting Waveforms: ")
        print("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        err, trN, trE, trZ = splitpy.utils.get_data_NEZ(
            client=client,
            sta=self.sta, start=tstart,
            end=tend, stdata=stdata, ndval=ndval)

        # Store as attributes with traces in dictionay
        self.err = err
        self.data = Data(trN, trE, trZ)

    def rotate_ZEN_LQT(self):
        """
        Rotates 3-component seismograms from vertical (Z),
        east (E) and north (N) to longitudinal (L), 
        radial (Q) and tangential (T) components of motion

        Attributes
        ----------

        data : :class:`~splitpy.classes.Data`
            Object containing :class:`~obspy.core.Trace` objects

        """

        inc = self.meta.inc*np.pi/180.
        baz = self.meta.baz*np.pi/180.

        M = np.zeros((3, 3))
        M[0, 0] = np.cos(inc)
        M[0, 1] = -np.sin(inc) * np.sin(baz)
        M[0, 2] = -np.sin(inc) * np.cos(baz)
        M[1, 0] = np.sin(inc)
        M[1, 1] = np.cos(inc) * np.sin(baz)
        M[1, 2] = np.cos(inc) * np.cos(baz)
        M[2, 0] = 0.
        M[2, 1] = -np.cos(baz)
        M[2, 2] = np.sin(baz)

        # Perform 3-D rotation
        LQT = np.dot(np.array(M), np.array(
            [self.data.trZ.data, self.data.trE.data, self.data.trN.data]))

        # Store into traces and add as new items in attribute dictionary
        self.data.trL = Trace(data=LQT[0], header=self.data.trZ.stats)
        self.data.trQ = Trace(data=LQT[1], header=self.data.trN.stats)
        self.data.trT = Trace(data=LQT[2], header=self.data.trE.stats)

    def calc_snrq(self, t1=None, dt=30.):
        """
        Calculates signal-to-noise ration on radial (Q) component

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Predicted arrival time of phase
        dt : float
            Duration (sec)

        Attributes
        ----------
        snrq : float
            Signal-to-noise ratio  (dB)

        """

        if t1 is None:
            t1 = self.meta.time + self.ts - 5.
        self.snrq = _calc_snr(self.data.trQ, t1=t1, dt=dt)

    def calc_snrt(self, t1=None, dt=30.):
        """
        Calculates signal-to-noise ration on tangential (T) component

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Predicted arrival time of phase
        dt : float
            Duration (sec)

        Attributes
        ----------
        snrt : float
            Signal-to-noise ratio  (dB)

        """

        if t1 is None:
            t1 = self.meta.time + self.ts - 5.
        self.snrt = _calc_snr(self.data.trT, t1=t1, dt=dt)

    def analyze(self, t1=None, t2=None):
        """
        Calculates the shear-wave splitting parameters based 
        on two alternative method: the Rotation-Correlation (RC)
        method and the Silver-Chan (SC) method. Each set of results
        is stored in a Dictionary as attributes of the split object.

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Start time of picking window
        t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            End time of picking window

        Attributes
        ----------
        RC_res : :class:`~splitpy.classes.Result`
            Object containing results of Rotation-Correlation method
        SC_res : :class:`~splitpy.classes.Result`
            Object containing results of Silver-Chan method

        """

        if t1 is None and t2 is None:
            t1 = self.meta.time + self.ts - 5.
            t2 = self.meta.time + self.ts + 25.

        # Calculate Silver and Chan splitting estimate
        print("* --> Calculating Rotation-Correlation (RC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
            splitpy.calc.split_RotCorr(
                self.data.trQ, self.data.trT,
                self.meta.baz, t1, t2, self.maxdt, self.ddt, self.dphi)

        # Calculate error
        edtt, ephi, errc = splitpy.calc.split_errorRC(
            trT_c,
            t1, t2, 0.05, Emat, self.maxdt, self.ddt, self.dphi)

        # Store dictionary as attribute
        self.RC_res = Result(Emat, trQ_c, trT_c, trFast, trSlow,
                             phi, dtt, phi_min, edtt, ephi, errc)

        # Calculate Silver and Chan splitting estimate
        print("* --> Calculating Silver-Chan (SC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
            splitpy.calc.split_SilverChan(
                self.data.trQ, self.data.trT,
                self.meta.baz, t1, t2, self.maxdt, self.ddt, self.dphi)

        # Calculate errors
        edtt, ephi, errc = splitpy.calc.split_errorSC(
            trT_c,
            t1, t2, 0.05, Emat, self.maxdt, self.ddt, self.dphi)

        self.SC_res = Result(Emat, trQ_c, trT_c, trFast, trSlow,
                             phi, dtt, phi_min, edtt, ephi, errc)

    def is_null(self, snrTlim=3., ds=-1):
        """
        Determines if splitting result is a Null result

        Parameters
        ----------
        snrTlim : float
            Threshold for snr on T component
        ds : int
            Number of stars to print out to screen (verbiage)

        Attributes
        ----------
        null : bool 
            Boolean for Null result

        """

        self.null = False

        # Calculate Angular Difference for Null Measurement
        dphi = max(abs(self.RC_res.phi - self.SC_res.phi),
                   abs(self.SC_res.phi - self.RC_res.phi))
        if dphi > 90.:
            dphi = 180. - dphi

        # Summarize Null Measurement
        if ds >= 0:
            print("*" + " "*ds + "Null Classification: ")
            if self.snrt < snrTlim:
                print(
                    "*" + " "*ds + "  SNR T Fail: {0:.2f} < {1:.2f}".format(
                        self.snrt, snrTlim))
            else:
                print(
                    "*" + " "*ds + "  SNR T Pass: {0:.2f} > {1:.2f}".format(
                        self.snrt, snrTlim))
            if 22. < dphi < 68.:
                print("*" + " "*ds +
                      "  dPhi Fail: {0:.2f} within 22. < X < 68.".format(dphi))
            else:
                print(
                    "*" + " "*ds +
                    "  dPhi Pass:  {0:.2f} outside 22. < X < 68.".format(dphi))

        # Check snr on tangential component
        if self.snrt < snrTlim:
            self.null = True

        # Check error on `phi_min` estimate
        if 22. < dphi < 68.:
            self.null = True

    def get_quality(self, ds):
        """
        Determines the quality of the estimate (either Null or non-Null)
        based on ratio of delay times and difference in fast axis directions
        between Rotation-Correlation and Silver-Chan methods

        Parameters
        ----------
        ds : int
            Number of stars to print out to screen (verbiage)

        Attributes
        ----------
        quality : str
            String representing quality of estimate ('Good', 'Fair', 'Poor')

        """

        # Ratio of delay times
        rho = self.RC_res.dtt/self.SC_res.dtt

        # Test based on difference in fast axis directions
        dphi = max(abs(self.RC_res.phi - self.SC_res.phi),
                   abs(self.SC_res.phi - self.RC_res.phi))
        if dphi > 90.:
            dphi = 180. - dphi

        # If estimate is Null
        if self.null:
            if rho < 0.2 and (37. < dphi < 53):
                self.quality = 'Good'
            elif rho < 0.3 and (32 < dphi < 58):
                self.quality = 'Fair'
            else:
                self.quality = 'Poor'
            if ds >= 0:
                print("*" + " "*ds +
                      "Quality Estimate: Null -- {0:s}".format(self.quality))
                print("*" + " "*ds +
                      "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
                print("*" + " "*ds +
                      "      Good: rho < 0.2  &&  37 < dphi < 53")
                print("*" + " "*ds +
                      "      Fair: rho < 0.3  &&  32 < dphi < 58")
                print("*" + " "*ds +
                      "      Poor: rho > 0.3  &&  dphi < 32 | dphi > 58")

        # If estimate is non-Null
        else:
            if (0.8 < rho < 1.1) and dphi < 8.:
                self.quality = 'Good'
            elif (0.7 < rho < 1.2) and dphi < 15.:
                self.quality = 'Fair'
            else:
                self.quality = 'Poor'
            if ds >= 0:
                print(
                    "*" + " "*ds +
                    "Quality Estimate: Non-Null -- {0:s}".format(self.quality))
                print("*" + " "*ds +
                      "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
                print("*" + " "*ds +
                      "      Good: 0.8 < rho < 1.1  &&  dphi < 8")
                print("*" + " "*ds +
                      "      Fair: 0.7 < rho < 1.2  &&  dphi < 15")
                print("*" + " "*ds +
                      "      Poor: rho < 0.7 | rho > 1.3 &&  dphi > 15")

    def display_results(self, ds=0):
        """
        Prints out best fitting results to screen

        Parameters
        ----------
        ds : int
            Number of stars to print out to screen (verbiage)

        """

        print(" "*ds + ' ======= Best-fit splitting results ========')
        print()
        print(" "*ds + ' Best fit values: RC method')
        print(" "*ds + ' Phi = ' +
              str("{:3d}").format(int(self.RC_res.phi)) +
              ' degrees +/- ' + str("{:2d}").format(int(self.RC_res.ephi)))
        print(" "*ds + ' dt = ' + str("{:.1f}").format(self.RC_res.dtt) +
              ' seconds +/- ' + str("{:.1f}").format(self.RC_res.edtt))
        print()
        print(" "*ds + ' Best fit values: SC method')
        print(" "*ds + ' Phi = ' +
              str("{:3d}").format(int(self.SC_res.phi)) +
              ' degrees +/- ' + str("{:2d}").format(int(self.SC_res.ephi)))
        print(" "*ds + ' dt = ' + str("{:.1f}").format(self.SC_res.dtt) +
              ' seconds +/- ' + str("{:.1f}").format(self.SC_res.edtt))
        print()

    def display_meta(self,  ds=0):
        """
        Prints out content of metadata to screen

        Parameters
        ----------
        ds : int
            Number of stars to print out to screen (verbiage)

        """

        print(" "*ds + ' ======= Meta data ========')
        print()
        print(" "*ds + 'SNR (dB):            ' +
              str("{:.0f}").format(self.snrq))
        print(" "*ds + 'Station:             ' + self.sta.station)
        print(" "*ds + 'Time:                ' + str(self.meta.time))
        print(" "*ds + 'Event depth (km):    ' +
              str("{:.0f}").format(self.meta.dep/1000.))
        print(" "*ds + 'Magnitude (Mw):      ' +
              str("{:.1f}").format(self.meta.mag))
        print(" "*ds + 'Longitude (deg):     ' +
              str("{:.2f}").format(self.meta.lon))
        print(" "*ds + 'Latitude (deg):      ' +
              str("{:.2f}").format(self.meta.lat))
        print(" "*ds + 'GAC (deg):           ' +
              str("{:.2f}").format(self.meta.gac))
        print(" "*ds + 'Backazimuth deg):    ' +
              str("{:.2f}").format(self.meta.baz))
        print(" "*ds + 'Incidence(deg):      ' +
              str("{:.2f}").format(self.meta.inc))
        print()

    def display_null_quality(self, ds=0):
        """
        Prints out null and quality estimates to screen

        Parameters
        ----------
        ds : int
            Number of stars to print out to screen (verbiage)

        """

        print(" "*ds + ' ======= Nulls and quality ========')
        print()
        print(" "*ds + ' Is Null?     ', self.null)
        print(" "*ds + ' Quality:     ', self.quality)
        print()

    def save(self, file):
        """
        Saves Split object to file

        Parameters
        ----------
        file : str
            File name for split object

        """

        import pickle
        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()


class PickPlot(object):
    """
    A PickPlot object contains figure handles and method to plot 
    seismic data for picking/refining the SKS time window.
    The figure displays the LQT seismograms and predicted arrival 
    times for common shear-wave arrivals in the time window. 
    The figure can be picked to refine the time window for
    calculating the  splitting results.

    Note
    ----
    The object is initialized with a :class:`~splitpy.classes.Split` object,
    which is temporarily stored as an attribute. All the split attributes 
    are therefore available to the :class:`~splitpy.clases.PickPlot` object 
    for plotting.

    Attributes
    ----------
    split : :class:`~splitpy.classes.Split` object
        Split object containing attributes after analysis has been carried out.
    fp : List
        List of figure handles ([fig, ax1, ax2, ax3])

    """

    def __init__(self, split):

        def init_pickw(ax, title, ylab):
            """
            Sets default axis attributes and return axis handle

            Parameters
            ----------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle
            title : str
                Title of axis
            ylab : str
                Y-axis label

            Returns
            -------
            ax : :class:`matplotlib.pyplot.axis`
                Axis handle

            """

            ax.clear()
            ax.set_title(title)
            ax.set_ylabel(ylab)
            ax.grid(which='major', axis='x')
            ax.grid(which='both', axis='y')
            ax.set_ylim((-1.1, 1.1))
            return ax

        if not hasattr(split, 'RC_res'):
            raise(
                Exception("analysis has not yet been performed " +
                          "on split object. Aborting"))

        # Make sure to clear figure if it already exists at initialization
        if plt.fignum_exists(1):
            plt.figure(1).clf()

        # Store split as attribute
        self.split = split

        # Define plot as GridSpec object
        gs = gspec.GridSpec(24, 1)

        # Figure handle
        fig = plt.figure(num=1, facecolor='w')

        # L component seismogram
        ax1 = fig.add_subplot(gs[0:7])  # 0:3
        ax1 = init_pickw(ax=ax1, title=self.split.sta.station, ylab='L')

        # Q component seismogram
        ax2 = fig.add_subplot(gs[8:15])  # 3:7
        ax2 = init_pickw(ax=ax2, title='', ylab='Q')

        # T component seismogram
        ax3 = fig.add_subplot(gs[16:23])  # 8:11
        ax3 = init_pickw(ax=ax3, title='', ylab='T')
        ax3.set_xlabel('Time (sec)')

        fp = [fig, ax1, ax2, ax3]

        # Ensure Figure is open
        fp[0].show()

        # Store handles as attribute
        self.fp = fp

    def plot_LQT_phases(self, tt, dts, t1=None, t2=None):
        """
        Plots rotated three-components of motion for picking.

        Parameters
        ----------
        tt : :class:`~obspy.taup.TauPyModel` 
            Taup object containing travel time and phase info
        dts : float
            Time interval (?)
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Start time of picking window
        t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            End time of picking window

        Attributes
        ----------
        ll : List
            List of :class:`~matplotlib.pyplot.axvline` objects for plotting the start and
            end times of the picking window.

        """

        # Default start and end times
        if t1 is None and t2 is None:
            t1 = self.split.meta.time + self.split.ts - 5.
            t2 = self.split.meta.time + self.split.ts + 25.

        # Define time axis
        taxis = np.arange(self.split.data.trL.stats.npts) / \
            self.split.data.trL.stats.sampling_rate - dts
        tstart = self.split.data.trL.stats.starttime

        # Set uniform vertical scale
        maxL = np.abs(self.split.data.trL.data).max()
        maxQ = np.abs(self.split.data.trQ.data).max()
        maxT = np.abs(self.split.data.trT.data).max()
        max = np.amax([maxL, maxQ, maxT])

        # Plot traces
        self.fp[1].plot(taxis, self.split.data.trL.data/max)
        self.fp[2].plot(taxis, self.split.data.trQ.data/max)
        self.fp[3].plot(taxis, self.split.data.trT.data/max)

        # Set tight limits
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

        # Store List of lines as attribute
        self.ll = ll

        # Add vertical lines and text for various phases
        for t in tt:

            name = t.name
            time = t.time
            if not name[0] == 'S':
                continue
            # Direct S phase
            if name == 'S':
                self.fp[1].axvline(time - self.split.ts, color='k')
                self.fp[1].text(time - self.split.ts + 5., -1.,
                                'S', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k')
                self.fp[2].text(time - self.split.ts + 5., -1.,
                                'S', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k')
                self.fp[3].text(time - self.split.ts + 5., -1.,
                                'S', rotation=90, ha='center', va='bottom')
            # Transmitted core-refracted shear wave - Main phase of interest
            elif name == 'SKS':
                self.fp[1].axvline(time - self.split.ts, color='k')
                self.fp[1].text(time - self.split.ts + 5., -1.,
                                'SKS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k')
                self.fp[2].text(time - self.split.ts + 5., -1.,
                                'SKS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k')
                self.fp[3].text(time - self.split.ts + 5., -1.,
                                'SKS', rotation=90, ha='center', va='bottom')
            # Core-refracted shear wave with 'K' bounce at core-mantle boundary
            elif name == 'SKKS':
                self.fp[1].axvline(time - self.split.ts, color='k')
                self.fp[1].text(time - self.split.ts + 5., -1.,
                                'SKKS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k')
                self.fp[2].text(time - self.split.ts + 5., -1.,
                                'SKKS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k')
                self.fp[3].text(time - self.split.ts + 5., -1.,
                                'SKKS', rotation=90, ha='center', va='bottom')
            # Direct shear wave reflected at core-mantle boundary
            elif name == 'ScS':
                self.fp[1].axvline(time - self.split.ts, color='k')
                self.fp[1].text(time - self.split.ts + 5., -1.,
                                'ScS', rotation=90, ha='center', va='bottom')
                self.fp[2].axvline(time - self.split.ts, color='k')
                self.fp[2].text(time - self.split.ts + 5., -1.,
                                'ScS', rotation=90, ha='center', va='bottom')
                self.fp[3].axvline(time - self.split.ts, color='k')
                self.fp[3].text(time - self.split.ts + 5., -1.,
                                'ScS', rotation=90, ha='center', va='bottom')

        # Update plot
        self.fp[0].canvas.draw()

    def update_LQT(self, tp1, tp2):
        """
        Updates LQT figure with newly picked times

        Parameters
        ----------
        tp1 : float
            New start time of picking window
        tp2 : float
            New end time of picking window

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
    """
    A DiagPlot object contains figure handles and methods to plot
    the diagnostic figure, which displays the LQT seismograms, 
    the corrected/un-corrected seismograms, the particle motions,
    the minimization matrix and a text box with a summary of the 
    analysis - for each of the two analysis methods ('RC' or 'SC')

    Note
    ----
    The object is initialized with a :class:`~splitpy.classes.Split` object,
    which is temporarily stored as an attribute. All the split attributes 
    are therefore available to the :class:`~splitpy.clases.DiagPlot` object 
    for plotting.

    Attributes
    ----------
    split : :class:`~splitpy.classes.Split` object
        Split object containing attributes after analysis has been carried out.
    fd : List
        List of figure handles ([fig, ax0, axt, axRC1, axRC2, axRC3, axRC4,\
                axSC1, axSC2, axSC3, axSC4])

    """

    def __init__(self, split):

        def init_splitw(ax, title):
            """
            Initializes window for splitting results

            Parameters
            ----------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle
            title : str
                Title of axis

            Returns
            -------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle

            """

            ax.clear()
            ax.set_ylim((-1.2, 1.2))
            ax.set_title(title, fontsize=12)
            return ax

        def init_pmotion(ax):
            """
            Initializes window for particle motion

            Parameters
            ----------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle

            Returns
            -------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle

            """

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
            """
            Initializes window for matrix minimization

            Parameters
            ----------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle
            title : str
                Title of axis

            Returns
            -------
            ax : :class:`~matplotlib.pyplot.axis`
                Axis handle

            """

            ax.clear()
            ax.set_title(title, fontsize=12)
            ax.set_ylabel(r'$\phi$ (deg)', labelpad=-10)
            ax.set_xlabel(r'$\delta t$ (sec)', labelpad=-2)
            ax.set_xticks([0, 1, 2, 3, 4])

            return ax

        if not hasattr(split, 'RC_res'):
            raise(
                Exception("analysis has not yet been performed on " +
                          "split object. Aborting"))

        # Make sure to clear figure if it already exists at initialization
        if plt.fignum_exists(2):
            plt.figure(2).clf()

        self.split = split

        # Figure handle
        fig = plt.figure(num=2, figsize=(10, 7), facecolor='w')

        # Q, T component seismograms
        ax0 = fig.add_axes([0.05, 0.73, 0.2, 0.2])
        ax0 = init_splitw(ax=ax0, title='Q, T')

        # Text box
        axt = fig.add_axes([0.45, 0.73, 0.3, 0.2])
        axt.axis('off')

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

        fd = [fig, ax0, axt, axRC1, axRC2, axRC3, axRC4,
              axSC1, axSC2, axSC3, axSC4]

        # Make sure figure is open
        fd[0].show()

        # Store handes as attribute
        self.fd = fd

    def plot_diagnostic(self, t1=None, t2=None):
        """
        Plots diagnostic window with estimates from both 'RC' and 'SC' methods

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Start time of picking window
        t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            End time of picking window

        """

        if t1 is None and t2 is None:
            t1 = self.split.meta.time + self.split.ts - 5.
            t2 = self.split.meta.time + self.split.ts + 25.

        def rot3D(inc, baz):
            """
            Defines rotation matrix from incidence and back-azimuth angles

            Parameters
            ----------
            inc : float
                Incidence angle in degrees
            baz : float
                Back-azimuth angle in degrees

            Returns
            -------
            M : :class:`~numpy.ndarray` 
                Rotation matrix

            """

            # Angles in radians
            inc = inc/180.*np.pi
            baz = baz/180.*np.pi

            # Define rotation matrix
            M = [[np.cos(inc), -np.sin(inc)*np.sin(baz),
                  -np.sin(inc)*np.cos(baz)],
                 [np.sin(inc), np.cos(inc)*np.sin(baz),
                  np.cos(inc)*np.cos(baz)],
                 [0., -np.cos(baz), np.sin(baz)]]

            return M

        # Copy traces to avoid overridding
        trL_tmp = self.split.data.trL.copy()
        trL_tmp.trim(t1, t2)
        trQ_tmp = self.split.data.trQ.copy()
        trQ_tmp.trim(t1, t2)
        trT_tmp = self.split.data.trT.copy()
        trT_tmp.trim(t1, t2)
        trE_tmp = self.split.data.trE.copy()
        trE_tmp.trim(t1, t2)
        trN_tmp = self.split.data.trN.copy()
        trN_tmp.trim(t1, t2)

        # Rotate seismograms for plots
        M = rot3D(self.split.meta.inc, self.split.meta.baz)
        ZEN_RC = np.dot(
            np.transpose(M),
            [trL_tmp.data, self.split.RC_res.trQ_c.data,
             self.split.RC_res.trT_c.data])
        E_RC = ZEN_RC[1, :]
        N_RC = ZEN_RC[2, :]
        ZEN_SC = np.dot(
            np.transpose(M),
            [trL_tmp.data, self.split.SC_res.trQ_c.data,
             self.split.SC_res.trT_c.data])
        E_SC = ZEN_SC[1, :]
        N_SC = ZEN_SC[2, :]

        # Original time series
        taxis = np.arange(trQ_tmp.stats.npts)/trQ_tmp.stats.sampling_rate
        max1 = np.abs(trQ_tmp.data).max()
        max2 = np.abs(trT_tmp.data).max()
        max = np.amax([max1, max2])

        self.fd[1].plot(taxis, trQ_tmp.data/max, 'b--')
        self.fd[1].plot(taxis, trT_tmp.data/max, 'r')

        # Text box
        self.fd[2].text(
            0.5, 0.9, 'Event: ' + self.split.meta.time.ctime() + '     ' +
            str(self.split.meta.lat) + 'N  ' +
            str(self.split.meta.lon) + 'E   ' +
            str(np.int(self.split.meta.dep/1000.)) + 'km   ' + 'Mw=' +
            str(self.split.meta.mag), horizontalalignment='center')
        self.fd[2].text(
            0.5, 0.7, 'Station: ' + self.split.sta.station +
            '   Backazimuth: ' +
            str("{:.2f}").format(self.split.meta.baz) + '   Distance: ' +
            str("{:.2f}").format(self.split.meta.gac),
            horizontalalignment='center')
        self.fd[2].text(
            0.5, 0.5, r'Best fit RC values: $\phi$=' +
            str(np.int(self.split.RC_res.phi)) + r'$\pm$' +
            str("{:.2f}").format(self.split.RC_res.ephi) +
            r'   $\delta t$=' +
            str(self.split.RC_res.dtt) + r'$\pm$' +
            str("{:.2f}").format(self.split.RC_res.edtt) +
            's', horizontalalignment='center')
        self.fd[2].text(
            0.5, 0.3, r'Best fit SC values: $\phi$=' +
            str(np.int(self.split.SC_res.phi)) + r'$\pm$' +
            str("{:.2f}").format(self.split.SC_res.ephi) +
            r'   $\delta t$=' +
            str(self.split.SC_res.dtt) + r'$\pm$' +
            str("{:.2f}").format(self.split.SC_res.edtt) +
            's', horizontalalignment='center')
        self.fd[2].text(
            0.5, 0.1, 'Is Null? ' + str(self.split.null) + '    Quality? ' +
            str(self.split.quality), horizontalalignment='center')

        # Rotation-correlation
        # Corrected Fast and Slow
        sum1 = np.sum(np.abs(self.split.RC_res.trFast.data -
                             self.split.RC_res.trSlow.data))
        sum2 = np.sum(np.abs(-self.split.RC_res.trFast.data -
                             self.split.RC_res.trSlow.data))
        if sum1 < sum2:
            sig = 1.
        else:
            sig = -1.
        taxis = np.arange(self.split.RC_res.trFast.stats.npts) / \
            self.split.RC_res.trFast.stats.sampling_rate
        max1 = np.abs(self.split.RC_res.trFast.data).max()
        max2 = np.abs(self.split.RC_res.trSlow.data).max()
        max = np.amax([max1, max2])

        self.fd[3].plot(taxis, self.split.RC_res.trFast.data/max, 'b--')
        self.fd[3].plot(taxis, sig*self.split.RC_res.trSlow.data/max, 'r')

        # Corrected Q and T
        self.fd[4].plot(taxis, self.split.RC_res.trQ_c.data/max, 'b--')
        self.fd[4].plot(taxis, self.split.RC_res.trT_c.data/max, 'r')

        # Prticle motion
        self.fd[5].plot(trE_tmp.data/max, trN_tmp.data/max, 'b--')
        self.fd[5].plot(E_RC/max, N_RC/max, 'r')

        # Map of energy
        plt.sca(self.fd[6])
        dt = np.arange(0., self.maxdt, self.ddt)
        phi = np.arange(-90., 90., self.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X, Y = np.meshgrid(dt, phi)
        E2 = np.roll(self.split.RC_res.Emat, np.int(
            self.split.RC_res.phi - self.split.RC_res.phi_min), axis=0)

        Emin = self.split.RC_res.Emat.min()
        Emax = self.split.RC_res.Emat.max()
        dE = (Emax - Emin)/16.
        levels = np.arange(Emin, Emax, dE)
        cmap = plt.cm.RdYlBu_r
        cset1 = plt.contour(X, Y, E2, levels,
                            cmap=plt.cm.get_cmap(cmap, len(levels)))

        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        errc = self.split.RC_res.errc
        ecset = plt.contour(X, Y, E2, (errc, ), colors='magenta',
                            linewidths=2)

        self.fd[6].axvline(self.split.RC_res.dtt)
        self.fd[6].axhline(self.split.RC_res.phi)

        # Silver-Chan
        # Corrected Fast and Slow
        sum1 = np.sum(np.abs(self.split.SC_res.trFast.data -
                             self.split.SC_res.trSlow.data))
        sum2 = np.sum(np.abs(-self.split.SC_res.trFast.data -
                             self.split.SC_res.trSlow.data))
        if sum1 < sum2:
            sig = 1.
        else:
            sig = -1.
        taxis = np.arange(self.split.SC_res.trFast.stats.npts) / \
            self.split.SC_res.trFast.stats.sampling_rate
        max1 = np.abs(self.split.SC_res.trFast.data).max()
        max2 = np.abs(self.split.SC_res.trSlow.data).max()
        max = np.amax([max1, max2])

        self.fd[7].plot(taxis, self.split.SC_res.trFast.data/max, 'b--')
        self.fd[7].plot(taxis, sig*self.split.SC_res.trSlow.data/max, 'r')

        # Corrected Q and T
        self.fd[8].plot(taxis, self.split.SC_res.trQ_c.data/max, 'b--')
        self.fd[8].plot(taxis, self.split.SC_res.trT_c.data/max, 'r')

        # Particle motion
        self.fd[9].plot(trE_tmp.data/max, trN_tmp.data/max, 'b--')
        self.fd[9].plot(E_SC/max, N_SC/max, 'r')

        # Map of energy
        plt.sca(self.fd[10])
        dt = np.arange(0., self.maxdt, self.ddt)
        phi = np.arange(-90., 90., self.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X, Y = np.meshgrid(dt, phi)

        E2 = np.roll(self.split.SC_res.Emat, np.int(
            self.split.SC_res.phi-self.split.SC_res.phi_min), axis=0)

        Emin = self.split.SC_res.Emat.min()
        Emax = self.split.SC_res.Emat.max()
        dE = (Emax - Emin)/16.
        levels = np.arange(Emin, Emax, dE)
        cmap = plt.cm.RdYlBu_r
        cset1 = plt.contour(X, Y, E2, levels,
                            cmap=plt.cm.get_cmap(cmap, len(levels)))

        errc = self.split.SC_res.errc
        ecset = plt.contour(X, Y, E2, (errc,), colors='magenta',
                            linewidths=2)

        self.fd[10].axvline(self.split.SC_res.dtt)
        self.fd[10].axhline(self.split.SC_res.phi)

        self.fd[0].canvas.draw()

    def save(self, file):
        """
        Saves current figure into file

        Parameters
        ----------
        file : str
            File name for diagnostic figure

        """

        self.fd[0].savefig(file)


def _calc_snr(tr, t1, dt):
    """
    Calculates signal-to-noise ratio based on start time
    and a given duration in seconds

    Returns
    -------
    snr : float
        Signal-to-noise ration (dB)

    """

    # Copy Z trace to signal and noise traces
    trSig = tr.copy()
    trNze = tr.copy()

    # Filter between 0.02 and 0.5 (dominant S wave frequencies)
    trSig.filter('bandpass', freqmin=0.02, freqmax=0.5,
                 corners=2, zerophase=True)
    trNze.filter('bandpass', freqmin=0.02, freqmax=0.5,
                 corners=2, zerophase=True)

    # Trim twin seconds around P-wave arrival
    trSig.trim(t1, t1 + dt)
    trNze.trim(t1 - dt, t1)

    # Calculate root mean square (RMS)
    srms = np.sqrt(np.mean(np.square(trSig.data)))
    nrms = np.sqrt(np.mean(np.square(trNze.data)))

    # Calculate signal/noise ratio in dB
    return 10*np.log10(srms*srms/nrms/nrms)
