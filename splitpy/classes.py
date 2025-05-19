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

The class :class:`~splitpy.classes.PickPlot` defines figure handles
for a picking window showing the seismograms and the predicted teleseismic
shear-wave phase arrivals. This figure is interactive and new picks can
be generated to refine the analysis.

The class :class:`~splitpy.classes.DiagPlot` defines figure handles
for a diagnostic figure showing a summary of the splitting results. It can
be called after each application of the `split.analyze` method to show
the summary of the analysis as a figure. This figure can also be saved as
a .png file.

"""

# -*- coding: utf-8 -*-
from math import ceil
import numpy as np
from splitpy import utils, calc
from obspy import Trace, Stream
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec


class Meta(object):
    """
    A Meta object contains attributes associated with the station-event
    data for a single Teleseismic event.

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
    ttime : float
        Predicted arrival time (sec)
    ph : str
        Phase name
    slow : float
        Horizontal slowness of phase
    inc : float
        Incidence angle of phase at surface
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)
    align : str
        Alignment of coordinate system for rotation
        ('ZRT', 'LQT', or 'PVH')
    rotated : bool
        Whether or not data have been rotated to ``align``
        coordinate system
    zcomp : str
        Vertical Component Identifier. Should be a single character.
        This is different then 'Z' only for fully unknown component
        orientation (i.e., components are 1, 2, 3)

    """

    def __init__(self, sta, event, gacmin=85., gacmax=180., phase='SKS',
                 maxdt=4., ddt=0.1, dphi=1.):

        from obspy.geodetics.base import gps2dist_azimuth as epi
        from obspy.geodetics import kilometer2degrees as k2d
        from obspy.taup import TauPyModel

        # Extract event 4D parameters
        self.time = event.origins[0].time
        self.lon = event.origins[0].longitude
        self.lat = event.origins[0].latitude
        self.dep = event.origins[0].depth

        # Check if depth is valid type
        if self.dep is not None:
            if self.dep > 1000.:
                self.dep = self.dep/1000.
        else:
            self.dep = 10.

        # Magnitude
        self.mag = event.magnitudes[0].mag
        if self.mag is None:
            self.mag = -9.

        # Calculate epicentral distance
        self.epi_dist, self.az, self.baz = epi(
            self.lat, self.lon, sta.latitude, sta.longitude)
        self.epi_dist /= 1000
        self.gac = k2d(self.epi_dist)

        if self.gac > gacmin and self.gac < gacmax:

            # Get travel time info
            tpmodel = TauPyModel(model='iasp91')

            # Get Travel times (Careful: here dep is in meters)
            arrivals = tpmodel.get_travel_times(
                distance_in_degree=self.gac,
                source_depth_in_km=self.dep,
                phase_list=[phase])
            if len(arrivals) > 1:
                print("arrival has many entries: ", len(arrivals))
            elif len(arrivals) == 0:
                print("no arrival found")
                self.accept = False
                return

            arrival = arrivals[0]

            # Attributes from parameters
            self.ttime = arrival.time
            self.slow = arrival.ray_param_sec_degree/111.
            self.inc = arrival.incident_angle
            self.phase = phase
            self.accept = True
        else:
            self.ttime = np.nan
            self.slow = np.nan
            self.inc = np.nan
            self.phase = None
            self.accept = False

        # Attributes that get updated as analysis progresses
        self.snrq = None
        self.snrt = None
        self.maxdt = maxdt
        self.ddt = ddt
        self.dphi = dphi
        self.align = 'LQT'
        self.rotated = False


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
    dataZNE : :class:`~splitpy.classes.Data`
        Object containing raw trace data in :class:`~obspy.core.Trace` format
    dataLQT : :class:`~splitpy.classes.Data`
        Object containing rotated trace data in :class:`~obspy.core.Trace` format

    """

    def __init__(self, sta, zcomp='Z'):

        # Attributes from parameters
        self.sta = sta

        # Initialize meta and data objects as None
        self.meta = None
        self.dataZNE = None
        self.dataLQT = None
        self.zcomp = zcomp

    def add_event(self, event, gacmin=85., gacmax=120., phase='SKS',
                  returned=False):
        """
        Adds event metadata to Split object as Meta object.

        Parameters
        ----------
        event : :class:`~obspy.core.event`
            Event metadata

        """

        from obspy.geodetics.base import gps2dist_azimuth as epi
        from obspy.geodetics import kilometer2degrees as k2d
        from obspy.taup import TauPyModel
        from obspy.core.event.event import Event

        if not isinstance(event, Event):
            raise(Exception("Event has incorrect type"))

        # Store as object attributes
        self.meta = Meta(sta=self.sta, event=event,
                         gacmin=gacmin, gacmax=gacmax,
                         phase=phase)

        if returned:
            return self.meta.accept

    def add_data(self, stream, returned=False, new_sr=5.):
        """
        Adds stream as object attribute

        Parameters
        ----------
        stream : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms
        returned : bool
            Whether or not to return the ``accept`` attribute

        Attributes
        ----------
        dataZNE : :class:`~obspy.core.Stream`
            Stream container for NEZ seismograms

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        """

        if not self.meta:
            raise Exception("No meta data available - aborting")

        if not self.meta.accept:
            return

        if not isinstance(stream, Stream):
            raise Exception("Event has incorrect type")

        try:
            self.dataZNE = stream

            if not np.allclose(
                    [t.stats.npts for t in stream[1:]], stream[0].stats.npts):
                self.meta.accept = False

            # Filter Traces
            if not stream[0].stats.sampling_rate == new_sr:
                self.dataZNE.filter(
                    'lowpass',
                    freq=0.5*new_sr,
                    corners=2,
                    zerophase=True)
                self.dataZNE.resample(
                    new_sr,
                    no_filter=True)

        except Exception as e:
            print("Error: Not all channels are available")
            self.meta.accept = False

        if returned:
            return self.meta.accept

    def download_data(self, client, new_sr=5., dts=120.,
                      returned=False, verbose=False,
                      remove_response=False):
        """
        Downloads seismograms based on event origin time and
        P phase arrival and adds as object attribute.

        Parameters
        ----------
        client : :class:`~obspy.client.fdsn.Client`
            Client object
        new_sr : float
            New sampling rate (Hz)
        dts : float
            Time duration (sec)
        remove_response : bool
            Remove instrument response from seismogram and resitute to true ground
            velocity (m/s) using obspy.core.trace.Trace.remove_response()
        returned : bool
            Whether or not to return the ``accept`` attribute
        verbose : bool
            Output diagnostics to screen

        Returns
        -------
        accept : bool
            Whether or not the object is accepted for further analysis

        Attributes
        ----------
        dataZNE : :class:`~obspy.core.Stream`
            Stream containing ZNE :class:`~obspy.core.Trace` objects
        dataZ12 : :class:`~obspy.core.Stream`
            Stream containing Z12 :class:`~obspy.core.Trace` objects
            (for un-oriented data)

        """

        if self.meta is None:
            raise Exception("Requires event data as attribute - aborting")

        if not self.meta.accept:
            return

        # Define start time for request
        tstart = self.meta.time + self.meta.ttime - dts
        tend = self.meta.time + self.meta.ttime + dts

        # Get waveforms
        print("* Requesting Waveforms: ")
        print("*    Startime: " + tstart.strftime("%Y-%m-%d %H:%M:%S"))
        print("*    Endtime:  " + tend.strftime("%Y-%m-%d %H:%M:%S"))

        # Download data
        err, stream = utils.download_data(
            client=client, sta=self.sta, start=tstart, end=tend,
            new_sr=new_sr, verbose=verbose, remove_response=remove_response,
            zcomp=self.zcomp)

        # Store as attributes with traces in dictionary
        try:
            trE = stream.select(component='E')[0]
            trN = stream.select(component='N')[0]
            trZ = stream.select(component='Z')[0]

            self.dataZNE = Stream(traces=[trZ, trN, trE])

            # Filter Traces and resample
            self.dataZNE.filter('lowpass', freq=0.5*new_sr,
                                corners=2, zerophase=True)
            self.dataZNE.resample(new_sr, no_filter=False)

        # If there is no ZNE, perhaps there is Z12 (or zcomp12)?
        except Exception as e:

            try:
                tr1 = stream.select(component='1')[0]
                tr2 = stream.select(component='2')[0]
                trZ = stream.select(component=self.zcomp)[0]

                # Force channel name to 'Z' if zcomp is not 'Z''
                if not self.zcomp == 'Z':
                    trZ.stats.channel = trZ.stats.channel[:-1] + 'Z'

                self.dataZNE = Stream(traces=[trZ, tr1, tr2])

                # Filter Traces and resample
                self.dataZNE.filter('lowpass', freq=0.5*new_sr,
                                    corners=2, zerophase=True)
                self.dataZNE.resample(new_sr, no_filter=True)

                # Save Z12 components in case it's necessary for later
                self.dataZ12 = self.dataZNE.copy()

                # Rotate from Z12 to ZNE using StDb azcorr attribute
                self.rotate(align='ZNE')

            except Exception as e:
                self.meta.accept = False

        if returned:
            return self.meta.accept

    def rotate(self, align=None):
        """
        Rotates 3-component seismograms from vertical (Z),
        east (E) and north (N) to longitudinal (L),
        radial (Q) and tangential (T) components of motion.
        Note that the method 'rotate' from ``obspy.core.stream.Stream``
        is used for the rotation ``'ZNE->LQT'``.

        Can also rotate Z12 to ZNE.

        Parameters
        ----------
        align : str
            Alignment of coordinate system for rotation
            ('ZNE' or 'LQT')

        Returns
        -------
        rotated : bool
            Whether or not the object has been rotated

        """

        if not self.meta.accept:
            return

        if self.meta.rotated:
            print("Data have been rotated already - continuing")
            return

        # Use default values from meta data if arguments are not specified
        if not align:
            align = self.meta.align

        if align == 'ZNE':
            from obspy.signal.rotate import rotate2zne

            # Copy traces
            trZ = self.dataZNE.select(component='Z')[0].copy()
            trN = self.dataZNE.select(component='1')[0].copy()
            trE = self.dataZNE.select(component='2')[0].copy()

            azim = self.sta.azcorr

            # Try with left handed system
            Z, N, E = rotate2zne(trZ.data, 0., -90., trN.data,
                                 azim, 0., trE.data, azim+90., 0.)

            # Z, N, E = rotate2zne(trZ.data, 0., -90., trN.data,
            #                      azim, 0., trE.data, azim+90., 0.)
            trN.data = N
            trE.data = E

            # Update stats of streams
            trN.stats.channel = trN.stats.channel[:-1] + 'N'
            trE.stats.channel = trE.stats.channel[:-1] + 'E'

            self.dataZNE = Stream(traces=[trZ, trN, trE])

        elif align == 'LQT':
            data = self.dataZNE.copy()
            data.rotate('ZNE->LQT',
                        back_azimuth=self.meta.baz,
                        inclination=self.meta.inc)
            # for tr in data:
            #     if tr.stats.channel.endswith('Q'):
            #         tr.data = -tr.data
            self.meta.align = align
            self.meta.rotated = True
            self.dataLQT = data.copy()

        else:
            raise(Exception("incorrect 'align' argument"))

    def calc_snr(self, t1=None, dt=30., fmin=0.02, fmax=0.5):
        """
        Calculates signal-to-noise ratio on either Z, L or P component

        Parameters
        ----------
        t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
            Pick time of arrival (sec)
        dt : float
            Duration (sec)
        fmin : float
            Minimum frequency corner for SNR filter (Hz)
        fmax : float
            Maximum frequency corner for SNR filter (Hz)

        Attributes
        ----------
        snrq : float
            Signal-to-noise ratio on vertical component (dB)
        snrh : float
            Signal-to-noise ratio on radial component (dB)

        """

        if not self.meta.accept:
            return

        if t1 is None:
            t1 = self.meta.time + self.meta.ttime

        # Copy trace to signal and noise traces
        trSigQ = self.dataLQT.select(component='Q')[0].copy()
        trNzeQ = self.dataLQT.select(component='Q')[0].copy()
        trSigT = self.dataLQT.select(component='T')[0].copy()
        trNzeT = self.dataLQT.select(component='T')[0].copy()

        trSigQ.detrend().taper(max_percentage=0.05)
        trNzeQ.detrend().taper(max_percentage=0.05)
        trSigT.detrend().taper(max_percentage=0.05)
        trNzeT.detrend().taper(max_percentage=0.05)

        # Filter between 0.1 and 1.0 (dominant P wave frequencies)
        trSigQ.filter('bandpass', freqmin=fmin, freqmax=fmax,
                      corners=2, zerophase=True)
        trSigT.filter('bandpass', freqmin=fmin, freqmax=fmax,
                      corners=2, zerophase=True)
        trNzeQ.filter('bandpass', freqmin=fmin, freqmax=fmax,
                      corners=2, zerophase=True)
        trNzeT.filter('bandpass', freqmin=fmin, freqmax=fmax,
                      corners=2, zerophase=True)

        # Trim around S-wave arrival
        trSigQ.trim(t1, t1 + dt)
        trNzeQ.trim(t1 - dt, t1)
        trSigT.trim(t1, t1 + dt)
        trNzeT.trim(t1 - dt, t1)

        # Calculate root mean square (RMS) and SNR
        srms = np.sqrt(np.mean(np.square(trSigQ.data)))
        nrms = np.sqrt(np.mean(np.square(trNzeQ.data)))
        self.meta.snrq = 10*np.log10(srms*srms/nrms/nrms)

        srms = np.sqrt(np.mean(np.square(trSigT.data)))
        nrms = np.sqrt(np.mean(np.square(trNzeT.data)))
        self.meta.snrt = 10*np.log10(srms*srms/nrms/nrms)

    def analyze(self, t1=None, t2=None, maxdt=4., ddt=0.1, dphi=1., verbose=False):
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
        verbose : bool
            Output diagnostics to screen

        Attributes
        ----------
        RC_res : :class:`~splitpy.classes.Result`
            Object containing results of Rotation-Correlation method
        SC_res : :class:`~splitpy.classes.Result`
            Object containing results of Silver-Chan method

        """

        self.meta.maxdt = maxdt
        self.meta.ddt = ddt
        self.meta.dphi = dphi

        if t1 is None and t2 is None:
            t1 = self.meta.time + self.meta.ttime - 5.
            t2 = self.meta.time + self.meta.ttime + 25.

        # Define signal
        trQ = self.dataLQT.select(component='Q')[0].copy()
        trT = self.dataLQT.select(component='T')[0].copy()

        # Calculate Silver and Chan splitting estimate
        if verbose:
            print("* --> Calculating Rotation-Correlation (RC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
            calc.split_RotCorr(
                trQ, trT, self.meta.baz, t1, t2,
                maxdt, ddt, dphi)

        # Calculate error
        edtt, ephi, errc = calc.split_errorRC(
            trT_c, t1, t2, 0.05, Emat,
            maxdt, ddt, dphi)

        # Store dictionary as attribute
        self.RC_res = Result(Emat, trQ_c, trT_c, trFast, trSlow,
                             phi, dtt, phi_min, edtt, ephi, errc)

        # Calculate Silver and Chan splitting estimate
        if verbose:
            print("* --> Calculating Silver-Chan (SC) Splitting")
        Emat, trQ_c, trT_c, trFast, trSlow, phi, dtt, phi_min = \
            calc.split_SilverChan(
                trQ, trT, self.meta.baz, t1, t2,
                maxdt, ddt, dphi)

        # Calculate errors
        edtt, ephi, errc = calc.split_errorSC(
            trT_c, t1, t2, 0.05, Emat,
            maxdt, ddt, dphi)

        self.SC_res = Result(Emat, trQ_c, trT_c, trFast, trSlow,
                             phi, dtt, phi_min, edtt, ephi, errc)

    def is_null(self, snrTlim=3., verbose=False):
        """
        Determines if splitting result is a Null result

        Parameters
        ----------
        snrTlim : float
            Threshold for snr on T component
        verbose : bool
            Output diagnostics to screen

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
        if verbose:
            print("*" + " "*5 + "Null Classification: ")
            if self.meta.snrt < snrTlim:
                print(
                    "*" + " "*5 + "  SNR T Fail: {0:.2f} < {1:.2f}".format(
                        self.meta.snrt, snrTlim))
            else:
                print(
                    "*" + " "*5 + "  SNR T Pass: {0:.2f} > {1:.2f}".format(
                        self.meta.snrt, snrTlim))
            if 22. < dphi < 68.:
                print("*" + " "*5 +
                      "  dPhi Fail: {0:.2f} within 22. < X < 68.".format(dphi))
            else:
                print(
                    "*" + " "*5 +
                    "  dPhi Pass:  {0:.2f} outside 22. < X < 68.".format(dphi))

        # Check snr on tangential component
        if self.meta.snrt < snrTlim:
            self.null = True

        # Check error on `phi_min` estimate
        if 22. < dphi < 68.:
            self.null = True

    def get_quality(self, verbose=False):
        """
        Determines the quality of the estimate (either Null or non-Null)
        based on ratio of delay times and difference in fast axis directions
        between Rotation-Correlation and Silver-Chan methods

        Parameters
        ----------
        verbose : bool
            Output diagnostics to screen

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
            if verbose:
                print("*" + " "*5 +
                      "Quality Estimate: Null -- {0:s}".format(self.quality))
                print("*" + " "*5 +
                      "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
                print("*" + " "*5 +
                      "      Good: rho < 0.2  &&  37 < dphi < 53")
                print("*" + " "*5 +
                      "      Fair: rho < 0.3  &&  32 < dphi < 58")
                print("*" + " "*5 +
                      "      Poor: rho > 0.3  &&  dphi < 32 | dphi > 58")

        # If estimate is non-Null
        else:
            if (0.8 < rho < 1.1) and dphi < 8.:
                self.quality = 'Good'
            elif (0.7 < rho < 1.2) and dphi < 15.:
                self.quality = 'Fair'
            else:
                self.quality = 'Poor'
            if verbose:
                print(
                    "*" + " "*5 +
                    "Quality Estimate: Non-Null -- {0:s}".format(self.quality))
                print("*" + " "*5 +
                      "    rho: {0:.2f}; dphi: {1:.2f}".format(rho, dphi))
                print("*" + " "*5 +
                      "      Good: 0.8 < rho < 1.1  &&  dphi < 8")
                print("*" + " "*5 +
                      "      Fair: 0.7 < rho < 1.2  &&  dphi < 15")
                print("*" + " "*5 +
                      "      Poor: rho < 0.7 | rho > 1.3 &&  dphi > 15")

    def display_results(self, ds=0):
        """
        Prints out best fitting results to screen

        Parameters
        ----------
        ds : int
            Number of spaces to print out to screen (verbiage)

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
            Number of spaces to print out to screen (verbiage)

        """

        print(" "*ds + ' ======= Meta data ========')
        print()
        print(" "*ds + ' SNR (dB):            ' +
              str("{:.0f}").format(self.meta.snrq))
        print(" "*ds + ' Station:             ' + self.sta.station)
        print(" "*ds + ' Time:                ' + str(self.meta.time))
        print(" "*ds + ' Event depth (km):    ' +
              str("{:.0f}").format(self.meta.dep/1000.))
        print(" "*ds + ' Magnitude (Mw):      ' +
              str("{:.1f}").format(self.meta.mag))
        print(" "*ds + ' Longitude (deg):     ' +
              str("{:.2f}").format(self.meta.lon))
        print(" "*ds + ' Latitude (deg):      ' +
              str("{:.2f}").format(self.meta.lat))
        print(" "*ds + ' GAC (deg):           ' +
              str("{:.2f}").format(self.meta.gac))
        print(" "*ds + ' Backazimuth deg):    ' +
              str("{:.2f}").format(self.meta.baz))
        print(" "*ds + ' Incidence (deg):      ' +
              str("{:.2f}").format(self.meta.inc))
        print(" "*ds + ' SNR - Q:      ' +
              str("{:.2f}").format(self.meta.snrq))
        print(" "*ds + ' SNR - T:      ' +
              str("{:.2f}").format(self.meta.snrt))

        print()

    def display_null_quality(self, ds=0):
        """
        Prints out null and quality estimates to screen

        Parameters
        ----------
        ds : int
            Number of spaces to print out to screen (verbiage)

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

        from obspy.taup import TauPyModel

        # Store split as attribute
        self.split = split

        # Get travel time info
        tpmodel = TauPyModel(model='iasp91')
        self.phase_list = ['S', 'SKS', 'SKKS', 'PKS', 'ScS']

        # Get Travel times (Careful: here dep is in meters)
        self.arrivals = tpmodel.get_travel_times(
            distance_in_degree=self.split.meta.gac,
            source_depth_in_km=self.split.meta.dep,
            phase_list=self.phase_list)

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

        if not hasattr(self.split, 'RC_res'):
            raise(
                Exception("analysis has not yet been performed " +
                          "on split object. Aborting"))

        # Make sure to clear figure if it already exists at initialization
        if plt.fignum_exists(1):
            plt.figure(1).clf()

        # Define plot as GridSpec object
        gs = gspec.GridSpec(24, 1)

        # Figure handle
        fig = plt.figure(num=1, facecolor='w')

        # L component seismogram
        ax1 = fig.add_subplot(gs[0:7])  # 0:3
        title = self.split.sta.network+'.'+self.split.sta.station
        ax1 = init_pickw(ax=ax1, title=title, ylab='L')

        # Q component seismogram
        ax2 = fig.add_subplot(gs[8:15])  # 3:7
        ax2 = init_pickw(ax=ax2, title='', ylab='Q')

        # T component seismogram
        ax3 = fig.add_subplot(gs[16:23])  # 8:11
        ax3 = init_pickw(ax=ax3, title='', ylab='T')
        ax3.set_xlabel('Time (sec)')

        axes = [fig, ax1, ax2, ax3]

        # Ensure Figure is open
        # fp[0].canvas.draw()
        axes[0].show()

        # Store handles as attribute
        self.axes = axes

    def plot_LQT_phases(self, dts, t1=None, t2=None):
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

        if not self.split.meta.rotated:
            raise(Exception("Split object not rotated - aborting"))

        # Default start and end times
        if t1 is None and t2 is None:
            t1 = self.split.meta.time + self.split.meta.ttime - 5.
            t2 = self.split.meta.time + self.split.meta.ttime + 25.

        # Define signal and noise
        trL = self.split.dataLQT.select(component='L')[0].copy()
        trQ = self.split.dataLQT.select(component='Q')[0].copy()
        trT = self.split.dataLQT.select(component='T')[0].copy()

        # Define time axis
        taxis = np.arange(trL.stats.npts)/trL.stats.sampling_rate - dts
        tstart = trL.stats.starttime

        # Set uniform vertical scale
        maxL = np.abs(trL.data).max()
        maxQ = np.abs(trQ.data).max()
        maxT = np.abs(trT.data).max()
        mmax = np.amax([maxL, maxQ, maxT])

        # Plot traces
        self.axes[1].plot(taxis, trL.data/mmax)
        self.axes[2].plot(taxis, trQ.data/mmax)
        self.axes[3].plot(taxis, trT.data/mmax)

        # Set tight limits
        self.axes[1].set_xlim((taxis[0], taxis[-1]))
        self.axes[2].set_xlim((taxis[0], taxis[-1]))
        self.axes[3].set_xlim((taxis[0], taxis[-1]))

        # Add vertical lines for picked times
        ll = list(range(6))
        ll[0] = self.axes[1].axvline(t1 - tstart - dts, color='r')
        ll[1] = self.axes[1].axvline(t2 - tstart - dts, color='r')
        ll[2] = self.axes[2].axvline(t1 - tstart - dts, color='r')
        ll[3] = self.axes[2].axvline(t2 - tstart - dts, color='r')
        ll[4] = self.axes[3].axvline(t1 - tstart - dts, color='r')
        ll[5] = self.axes[3].axvline(t2 - tstart - dts, color='r')

        # Store List of lines as attribute
        self.ll = ll

        # Add vertical lines and text for various phases
        for t in self.arrivals:

            name = t.name
            time = t.time

            self.axes[1].axvline(time - self.split.meta.ttime, color='k')
            self.axes[1].text(time - self.split.meta.ttime + 5., -1.,
                              name, rotation=90, ha='center', va='bottom')
            self.axes[2].axvline(time - self.split.meta.ttime, color='k')
            self.axes[2].text(time - self.split.meta.ttime + 5., -1.,
                              name, rotation=90, ha='center', va='bottom')
            self.axes[3].axvline(time - self.split.meta.ttime, color='k')
            self.axes[3].text(time - self.split.meta.ttime + 5., -1.,
                              name, rotation=90, ha='center', va='bottom')

        # Update plot
        self.axes[0].canvas.draw()

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
        self.ll[0] = self.axes[1].axvline(tp1, color='r')
        self.ll[1] = self.axes[1].axvline(tp2, color='r')
        self.ll[2] = self.axes[2].axvline(tp1, color='r')
        self.ll[3] = self.axes[2].axvline(tp2, color='r')
        self.ll[4] = self.axes[3].axvline(tp1, color='r')
        self.ll[5] = self.axes[3].axvline(tp2, color='r')

        # Update figure
        self.axes[0].canvas.draw()

    def save(self, file):
        """
        Saves current figure into file

        Parameters
        ----------
        file : str
            File name for diagnostic figure

        """

        self.axes[0].savefig(file)


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

        self.split = split

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

        if not hasattr(self.split, 'RC_res'):
            raise(
                Exception("analysis has not yet been performed on " +
                          "split object. Aborting"))

        # Make sure to clear figure if it already exists at initialization
        if plt.fignum_exists(2):
            plt.figure(2).clf()

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

        axes = [fig, ax0, axt, axRC1, axRC2, axRC3, axRC4,
                axSC1, axSC2, axSC3, axSC4]

        # # Make sure figure is open
        axes[0].show()

        # Store handes as attribute
        self.axes = axes

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
        import matplotlib

        if t1 is None and t2 is None:
            t1 = self.split.meta.time + self.split.meta.ttime - 5.
            t2 = self.split.meta.time + self.split.meta.ttime + 25.

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
        trL_tmp = self.split.dataLQT.select(component='L')[0].copy()
        trL_tmp.trim(t1, t2)
        trQ_tmp = self.split.dataLQT.select(component='Q')[0].copy()
        trQ_tmp.trim(t1, t2)
        trT_tmp = self.split.dataLQT.select(component='T')[0].copy()
        trT_tmp.trim(t1, t2)
        trE_tmp = self.split.dataZNE.select(component='E')[0].copy()
        trE_tmp.trim(t1, t2)
        trN_tmp = self.split.dataZNE.select(component='N')[0].copy()
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
        mmax = np.amax([max1, max2])

        self.axes[1].plot(taxis, trQ_tmp.data/mmax, 'b--')
        self.axes[1].plot(taxis, trT_tmp.data/mmax, 'r')
        self.axes[1].text(taxis[0], 1, 'Q', verticalalignment='top',
                          horizontalalignment='left', color='b')
        self.axes[1].text(taxis[0], -1, 'T', verticalalignment='bottom',
                          horizontalalignment='left', color='r')

        # Text box
        self.axes[2].text(
            0.5, 0.9, 'Event: ' + self.split.meta.time.ctime() + '     ' +
            str(self.split.meta.lat) + 'N  ' +
            str(self.split.meta.lon) + 'E   ' +
            str(int(self.split.meta.dep)) + 'km   ' + 'Mw=' +
            str(self.split.meta.mag), horizontalalignment='center')
        self.axes[2].text(
            0.5, 0.7, 'Station: ' + self.split.sta.station +
            '   Backazimuth: ' +
            str("{:.2f}").format(self.split.meta.baz) + '   Distance: ' +
            str("{:.2f}").format(self.split.meta.gac),
            horizontalalignment='center')
        self.axes[2].text(
            0.5, 0.5, r'Best fit RC values: $\phi$=' +
            str(int(self.split.RC_res.phi)) + r'$\pm$' +
            str("{:.2f}").format(self.split.RC_res.ephi) +
            r'   $\delta t$=' +
            str(self.split.RC_res.dtt) + r'$\pm$' +
            str("{:.2f}").format(self.split.RC_res.edtt) +
            's', horizontalalignment='center')
        self.axes[2].text(
            0.5, 0.3, r'Best fit SC values: $\phi$=' +
            str(int(self.split.SC_res.phi)) + r'$\pm$' +
            str("{:.2f}").format(self.split.SC_res.ephi) +
            r'   $\delta t$=' +
            str(self.split.SC_res.dtt) + r'$\pm$' +
            str("{:.2f}").format(self.split.SC_res.edtt) +
            's', horizontalalignment='center')
        self.axes[2].text(
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
        mmax = np.amax([max1, max2])

        self.axes[3].plot(taxis, self.split.RC_res.trFast.data/mmax, 'b--')
        self.axes[3].plot(taxis, sig*self.split.RC_res.trSlow.data/mmax, 'r')
        self.axes[3].text(taxis[0], 1, 'Fast', verticalalignment='top',
                          horizontalalignment='left', color='b')
        self.axes[3].text(taxis[0], -1, 'Slow', verticalalignment='bottom',
                          horizontalalignment='left', color='r')

        # Corrected Q and T
        self.axes[4].plot(taxis, self.split.RC_res.trQ_c.data/mmax, 'b--')
        self.axes[4].plot(taxis, self.split.RC_res.trT_c.data/mmax, 'r')

        # Particle motion
        self.axes[5].plot(trE_tmp.data/mmax, trN_tmp.data/mmax, 'b--')
        self.axes[5].plot(E_RC/mmax, N_RC/mmax, 'r')
        self.axes[5].text(-1, 1, 'Raw', verticalalignment='top',
                          horizontalalignment='left', color='b')
        self.axes[5].text(-1, -1, 'RC', verticalalignment='bottom',
                          horizontalalignment='left', color='r')
        ang = 360. - self.split.meta.baz
        x1pos = -np.sin(ang*np.pi/180.)
        y1pos = np.cos(ang*np.pi/180.)
        x2pos = -x1pos
        y2pos = -y1pos
        self.axes[5].plot([x1pos, x2pos], [y1pos, y2pos], 'k:', lw=2)

        # Map of energy
        plt.sca(self.axes[6])
        dt = np.arange(0., self.split.meta.maxdt, self.split.meta.ddt)
        phi = np.arange(-90., 90., self.split.meta.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X, Y = np.meshgrid(dt, phi)
        E2 = np.roll(
            self.split.RC_res.Emat,
            int((self.split.RC_res.phi - self.split.RC_res.phi_min)/self.split.meta.dphi),
            axis=0)

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

        self.axes[6].axvline(self.split.RC_res.dtt)
        self.axes[6].axhline(self.split.RC_res.phi)

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
        mmax = np.amax([max1, max2])

        self.axes[7].plot(taxis, self.split.SC_res.trFast.data/mmax, 'b--')
        self.axes[7].plot(taxis, sig*self.split.SC_res.trSlow.data/mmax, 'r')

        # Corrected Q and T
        self.axes[8].plot(taxis, self.split.SC_res.trQ_c.data/mmax, 'b--')
        self.axes[8].plot(taxis, self.split.SC_res.trT_c.data/mmax, 'r')

        # Particle motion
        self.axes[9].plot(trE_tmp.data/mmax, trN_tmp.data/mmax, 'b--')
        self.axes[9].plot(E_SC/mmax, N_SC/mmax, 'r')
        self.axes[9].plot([x1pos, x2pos], [y1pos, y2pos], 'k:', lw=2)
        self.axes[9].text(-1, 1, 'Raw', verticalalignment='top',
                          horizontalalignment='left', color='b')
        self.axes[9].text(-1, -1, 'SC', verticalalignment='bottom',
                          horizontalalignment='left', color='r')

        # Map of energy
        plt.sca(self.axes[10])
        dt = np.arange(0., self.split.meta.maxdt, self.split.meta.ddt)
        phi = np.arange(-90., 90., self.split.meta.dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X, Y = np.meshgrid(dt, phi)

        E2 = np.roll(
            self.split.SC_res.Emat,
            int((self.split.SC_res.phi-self.split.SC_res.phi_min)/self.split.meta.dphi),
            axis=0)

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

        self.axes[10].axvline(self.split.SC_res.dtt)
        self.axes[10].axhline(self.split.SC_res.phi)

        self.axes[0].canvas.draw()
        # plt.show()

    def save(self, file):
        """
        Saves current figure into file

        Parameters
        ----------
        file : str
            File name for diagnostic figure

        """

        self.axes[0].savefig(file)
