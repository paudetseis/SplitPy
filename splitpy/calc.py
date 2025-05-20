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

# -*- coding: utf-8 -*-
import copy
import numpy as np
from scipy import stats
from obspy.core import Trace, Stream
from numpy.linalg import inv


def split_SilverChan(trQ, trT, baz, t1, t2, maxdt, ddt, dphi):
    """
    Calculates splitting based on the minimization of energy on 
    the corrected transverse component (Silver and Chan, 1990)

    Parameters
    ----------
    trQ : :class:`~obspy.core.Trace`
        Radial component seismogram
    trT : :class:`~obspy.core.Trace`
        Tangential component seismogram
    baz : float
        Back-azimuth - pointing to earthquake from station (degrees)
    t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time of picking window
    t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time of picking window
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)

    Returns
    -------
    Ematrix : :class:`~numpy.ndarray`
        Matrix of T component energy
    trQ_c : :class:`~obspy.core.Trace`
        Trace of corrected radial component of motion
    trT_c : :class:`~obspy.core.Trace`
        Trace of corrected tangential component of motion
    trFast : :class:`~obspy.core.Trace`
        Trace of corrected fast direction of motion
    trSlow : :class:`~obspy.core.Trace`
        Trace of corrected slow direction of motion\
    phiSC : float
        Azimuth of fast axis (deg)
    dttSC : foat
        Delay time between fast and slow axes (sec)
    phi_min : float 
        Azimuth used in plotting routine

    """

    phi = np.arange(-90.0, 90.0, dphi)*np.pi/180.
    dtt = np.arange(0., maxdt, ddt)

    M = np.zeros((2, 2, len(phi)))
    M[0, 0, :] = np.cos(phi)
    M[0, 1, :] = -np.sin(phi)
    M[1, 0, :] = np.sin(phi)
    M[1, 1, :] = np.cos(phi)

    Ematrix = np.zeros((len(phi), len(dtt)))

    trQ_tmp = trQ.copy()
    trT_tmp = trT.copy()
    trQ_tmp.trim(t1, t2)
    trT_tmp.trim(t1, t2)

    trQ_tmp.taper(max_percentage=0.1, type='hann')
    trT_tmp.taper(max_percentage=0.1, type='hann')

    # Rotation loop
    for p in range(len(phi)):

        # Test fast/slow direction
        FS_test = np.dot(np.array(M[:, :, p]), np.array(
            [trQ_tmp.data, trT_tmp.data]))

        # Compile into traces
        F0 = Trace(data=np.array(FS_test[0]), header=trQ_tmp.stats)
        F1 = Trace(data=np.array(FS_test[1]), header=trT_tmp.stats)

        # Time shift loop
        for t in range(len(dtt)):

            shift = dtt[t]
            # Shift by dtt/2. each component (+/-)
            tmpFast = tshift(F0, -shift/2.)
            tmpSlow = tshift(F1, shift/2.)

            # Rotate back to Q and T system
            corrected_QT = np.dot(
                inv(np.array(M[:, :, p])), np.array([tmpFast, tmpSlow]))
            trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
            trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)

            #plot_cmp(trQ_c, trT_c,'phi='+str(phi[p]*180./np.pi)+'dt='+str(shift))

            # Energy on transverse component (component 1)
            Ematrix[p, t] = np.sum(np.square(corrected_QT[1]))

    # Find indices of minimum value of Energy matrix
    ind = np.where(Ematrix == Ematrix.min())
    ind_phi = ind[0][0]
    ind_dtt = ind[1][0]

    # Get best-fit phi and dt
    shift = dtt[ind_dtt]
    phiSC_min = phi[ind_phi]*180./np.pi
    phiSC = np.mod((phiSC_min + baz), 180.)

    if phiSC > 90.:
        phiSC = phiSC - 180.

    FS_test = np.dot(np.array(M[:, :, ind_phi]),
                     np.array([trQ_tmp.data, trT_tmp.data]))

    F0 = Trace(data=FS_test[0], header=trQ_tmp.stats)
    F1 = Trace(data=FS_test[1], header=trT_tmp.stats)

    tmpFast = tshift(F0, -shift/2.)
    tmpSlow = tshift(F1, shift/2.)

    corrected_QT = np.dot(
        inv(np.array(M[:, :, ind_phi])), np.array([tmpFast, tmpSlow]))

    trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
    trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)

    trFast = Trace(data=tmpFast, header=trT_tmp.stats)
    trSlow = Trace(data=tmpSlow, header=trQ_tmp.stats)

    return Ematrix, trQ_c, trT_c, trFast, trSlow, \
        phiSC, shift, phiSC_min


def split_RotCorr(trQ, trT, baz, t1, t2, maxdt, ddt, dphi):
    """
    Calculates splitting based on the maximum correlation between corrected 
    radial and tangential components of motion 

    Parameters
    ----------
    trQ : :class:`~obspy.core.Trace`
        Radial component seismogram
    trT : :class:`~obspy.core.Trace`
        Tangential component seismogram
    baz : float
        Back-azimuth - pointing to earthquake from station (degrees)
    t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time of picking window
    t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time of picking window
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)

    Returns
    -------
    Ematrix : :class:`~numpy.ndarray`
        Matrix of T component energy
    trQ_c : :class:`~obspy.core.Trace`
        Trace of corrected radial component of motion
    trT_c : :class:`~obspy.core.Trace`
        Trace of corrected tangential component of motion
    trFast : :class:`~obspy.core.Trace`
        Trace of corrected fast direction of motion
    trSlow : :class:`~obspy.core.Trace`
        Trace of corrected slow direction of motion\
    phiSC : float
        Azimuth of fast axis (deg)
    dttSC : foat
        Delay time between fast and slow axes (sec)
    phi_min : float 
        Azimuth used in plotting routine

    """

    phi = np.arange(-90.0, 90.0, dphi)*np.pi/180.
    dtt = np.arange(0., maxdt, ddt)

    M = np.zeros((2, 2, len(phi)))
    M[0, 0, :] = np.cos(phi)
    M[0, 1, :] = -np.sin(phi)
    M[1, 0, :] = np.sin(phi)
    M[1, 1, :] = np.cos(phi)

    Cmatrix_pos = np.zeros((len(phi), len(dtt)))
    Cmatrix_neg = np.zeros((len(phi), len(dtt)))

    trQ_tmp = trQ.copy()
    trT_tmp = trT.copy()
    trQ_tmp.trim(t1, t2)
    trT_tmp.trim(t1, t2)
    npts = trQ_tmp.stats.npts

    trQ_tmp.taper(max_percentage=0.1, type='hann')
    trT_tmp.taper(max_percentage=0.1, type='hann')

    # Rotation loop
    for p in range(len(phi)):

        # Test fast/slow direction
        FS_test = np.dot(np.array(M[:, :, p]), np.array(
            [trQ_tmp.data, trT_tmp.data]))

        # Compile into traces
        F0 = Trace(data=np.array(FS_test[0]), header=trQ_tmp.stats)
        F1 = Trace(data=np.array(FS_test[1]), header=trT_tmp.stats)

        # Cross-correlate Fast with Slow
        ns0 = np.sum(F0.data*F0.data)
        ns1 = np.sum(F1.data*F1.data)
        norm = np.sqrt(ns0*ns1)

        cor = Trace(data=np.fft.ifftshift(np.correlate(
            F0.data, F1.data, mode='same')/norm), header=trQ_tmp.stats)

        # Time shift loop
        for t in range(len(dtt)):

            shift = dtt[t]

            # Shift by dtt each component (+/-)
            cor_pos = tshift(cor, shift)
            cor_neg = tshift(cor, -shift)
            Cmatrix_pos[p, t] = cor_pos[0]
            Cmatrix_neg[p, t] = cor_neg[0]

    # Time shift is positive: fast axis arrives after slow axis
    if abs(Cmatrix_pos).max() > abs(Cmatrix_neg).max():

        # print 'Cmatrix_pos is max'

        ind = np.where(Cmatrix_pos == max(
            Cmatrix_pos.max(), Cmatrix_pos.min(), key=abs))
        ind_phi = ind[0][0]
        ind_dtt = ind[1][0]

        # Get best-fit phi and dt
        dtRC = dtt[ind_dtt]
        phiRC_max = phi[ind_phi]*180./np.pi
        phiRC = np.mod((phiRC_max + baz - 90.), 180.)
        Cmap = Cmatrix_pos
        theta = (phiRC_max - 90.)/180.*np.pi

        S = np.sign(max(Cmatrix_pos.max(), Cmatrix_pos.min(), key=abs))

    # Time shift is negative: fast axis arrives before slow axis
    else:

        ind = np.where(Cmatrix_neg == max(
            Cmatrix_neg.max(), Cmatrix_neg.min(), key=abs))
        ind_phi = ind[0][0]
        ind_dtt = ind[1][0]

        # Get best-fit phi and dt
        dtRC = dtt[ind_dtt]
        phiRC_max = phi[ind_phi]*180./np.pi
        phiRC = np.mod((phiRC_max + baz), 180.)
        Cmap = Cmatrix_neg
        theta = (phiRC_max)/180.*np.pi

        S = np.sign(max(Cmatrix_neg.max(), Cmatrix_neg.min(), key=abs))

    Cmap = Cmap * (-S)
    shift = dtRC
    theta = theta + np.pi/2.

    if phiRC > 90.:
        phiRC = phiRC - 180.

    M2 = np.zeros((2, 2))
    M2[0, 0] = np.cos(theta)
    M2[0, 1] = -np.sin(theta)
    M2[1, 0] = np.sin(theta)
    M2[1, 1] = np.cos(theta)

    FS_test = np.dot(np.array(M2[:, :]), np.array(
        [trQ_tmp.data, trT_tmp.data]))

    F0 = Trace(data=FS_test[0], header=trQ_tmp.stats)
    F1 = Trace(data=FS_test[1], header=trT_tmp.stats)

    tmpFast = tshift(F0, shift/2.)
    tmpSlow = tshift(F1, -shift/2.)

    trFast = Trace(data=tmpFast, header=trT_tmp.stats)
    trSlow = Trace(data=tmpSlow, header=trQ_tmp.stats)

    corrected_QT = np.dot(
        inv(np.array(M2[:, :])), np.array([tmpFast, tmpSlow]))

    trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
    trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)

    return Cmap, trQ_c, trT_c, trFast, trSlow, \
        phiRC, dtRC, phiRC_max

def tshift(trace, tt):
    """
    Shifts a :class:`~obspy.core.Trace` object

    Parameters
    ----------
    trace : :class:`~obspy.core.Trace`
        Seismogram to apply shift
    tt : float
        Lag time for shifting 

    Returns
    -------
    rtrace: :class:`~obspy.core.Trace`
        Shifted version of trace

    """

    nt = trace.stats.npts
    dt = trace.stats.delta
    freq = np.fft.fftfreq(nt, d=dt)

    ftrace = np.fft.fft(trace.data)

    for i in range(len(freq)):
        ftrace[i] = ftrace[i]*np.exp(2.*np.pi*1j*freq[i]*tt)

    rtrace = np.real(np.fft.ifft(ftrace))

    return rtrace


def split_dof(tr):
    """
    Determines the degrees of freedom to calculate the
    confidence region of the misfit function. 
    From Walsh, JGR, (2013)

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Seismogram 

    Returns
    -------
    dof : float
        Degrees of freedom

    """

    F = np.abs(np.fft.fft(tr.data)[0:int(len(tr.data)/2) + 1])

    E2 = np.sum(F**2)
    E2 -= (F[0]**2 + F[-1]**2)/2.
    E4 = (1./3.)*(F[0]**4 + F[-1]**4)
    for i in range(1, len(F) - 1):
        E4 += (4./3.)*F[i]**4

    dof = int(4.*E2**2/E4 - 2.)

    return dof


def split_errorSC(tr, t1, t2, q, Emat, maxdt, ddt, dphi):
    """
    Calculate error bars based on a F-test and 
    a given confidence interval q

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Seismogram 
    t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time of picking window
    t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time of picking window
    q : float
        Confidence level
    Emat : :class:`~numpy.ndarray`
        Energy minimization matrix
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)

    Returns
    -------
    err_dtt : float
        Error in dt estimate (sec)
    err_phi : float
        Error in phi estimate (degrees)
    err_contour : :class:`~numpy.ndarray`
        Error contour for plotting

    """

    # Bounds on search
    phi = np.arange(-90.0, 90.0, dphi)*np.pi/180.
    dtt = np.arange(0., maxdt, ddt)

    # Copy trace to avoid overriding
    tr_tmp = tr.copy()
    tr_tmp.trim(t1, t2)

    # Get degrees of freedom
    dof = split_dof(tr_tmp)
    if dof < 3:
        dof = 3
        print(
            "Degrees of freedom < 3. Fixing to DOF = 3, which may " +
            "result in inaccurate errors")
    n_par = 2

    # Error contour
    vmin = Emat.min()
    err_contour = vmin*(1. + n_par/(dof - n_par) *
                        stats.f.ppf(1. - q, n_par, dof - n_par))

    # Roll copy of Emat to center
    ind = np.where(Emat == Emat.min())
    ind_phi = ind[0][0]
    E2 = np.roll(Emat, len(phi)//2-ind_phi, axis=0)

    # Estimate uncertainty (q confidence interval)
    err = np.where(E2 < err_contour)
    if len(err) == 0:
        return False, False, False
    err_phi = max(
        0.25*(phi[max(err[0])] - phi[min(err[0])])*180./np.pi, 0.25*dphi)
    err_dtt = max(0.25*(dtt[max(err[1])] - dtt[min(err[1])]), 0.25*ddt)

    return err_dtt, err_phi, err_contour


def split_errorRC(tr, t1, t2, q, Emat, maxdt, ddt, dphi):
    """
    Calculates error bars based on a F-test and 
    a given confidence level q.

    Note
    ----
    This version uses a Fisher transformation for 
    correlation-type misfit.

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Seismogram 
    t1 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time of picking window
    t2 : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time of picking window
    q : float
        Confidence level
    Emat : :class:`~numpy.ndarray`
        Energy minimization matrix
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)

    Returns
    -------
    err_dtt : float
        Error in dt estimate (sec)
    err_phi : float
        Error in phi estimate (degrees)
    err_contour : :class:`~numpy.ndarray`
        Error contour for plotting

    """

    phi = np.arange(-90.0, 90.0, dphi)*np.pi/180.
    dtt = np.arange(0., maxdt, ddt)

    # Copy trace to avoid overriding
    tr_tmp = tr.copy()
    tr_tmp.trim(t1, t2)

    # Get degrees of freedom
    dof = split_dof(tr_tmp)
    if dof <= 3:
        dof = 3.01
        print(
            "Degrees of freedom < 3. Fixing to DOF = 3, which may " +
            "result in inaccurate errors")
    n_par = 2

    # Fisher transformation
    vmin = np.arctanh(Emat.min())

    # Error contour
    zrr_contour = vmin + (vmin*np.sign(vmin)*n_par/(dof - n_par) *
                          stats.f.ppf(1. - q, n_par, dof - n_par)) *\
        np.sqrt(1./(dof-3))

    # Back transformation
    err_contour = np.tanh(zrr_contour)

    # Roll copy of Emat to center
    ind = np.where(Emat == Emat.min())
    ind_phi = ind[0][0]
    E2 = np.roll(Emat, len(phi)//2-ind_phi, axis=0)

    # Estimate uncertainty (q confidence interval)
    err = np.where(E2 < err_contour)
    err_phi = max(
        0.25*(phi[max(err[0])] - phi[min(err[0])])*180./np.pi, 0.25*dphi)
    err_dtt = max(0.25*(dtt[max(err[1])] - dtt[min(err[1])]), 0.25*ddt)

    return err_dtt, err_phi, err_contour

def split_error_average(q, Emat, maxdt, ddt, dphi, n):
    """
    Calculate error bars based on a F-test and 
    a given confidence level q

    Parameters
    ----------
    q : float
        Confidence level
    Emat : :class:`~numpy.ndarray`
        Error surface matrix
    maxdt : float
        Maximum delay time considered in grid search (sec)
    ddt : float
        Delay time interval in grid search (sec)
    dphi : float
        Angular interval in grid search (deg)
    n : int
        Number of measurements

    Returns
    -------
    err_dtt : float
        Error in dt estimate (sec)
    err_phi : float
        Error in phi estimate (degrees)
    err_contour : :class:`~numpy.ndarray`
        Error contour for plotting

    """

    # Bounds on search
    phi = np.arange(-90.0, 90.0, dphi)*np.pi/180.
    dtt = np.arange(0., maxdt, ddt)

    # Get degrees of freedom
    dof = 10*n # To be conservative - see Frederiksen et al. (2025)
    n_par = 2

    # Error contour
    Emin = Emat.min()
    if Emin < 0:
        err_contour = Emin*(1. - n_par/(dof - n_par) *
                            stats.f.ppf(1. - q, n_par, dof - n_par))
    else:
        err_contour = Emin*(1. + n_par/(dof - n_par) *
                            stats.f.ppf(1. - q, n_par, dof - n_par))

    # Roll copy of Emat to center
    ind = np.where(Emat == Emat.min())
    ind_phi = ind[0][0]
    E2 = np.roll(Emat, len(phi)//2-ind_phi, axis=0)

    # Estimate uncertainty (q confidence interval)
    err = np.where(E2 < err_contour)
    if len(err) == 0:
        return False, False, False
    err_phi = max(
        0.25*(phi[max(err[0])] - phi[min(err[0])])*180./np.pi, 0.25*dphi)
    err_dtt = max(0.25*(dtt[max(err[1])] - dtt[min(err[1])]), 0.25*ddt)

    return err_dtt, err_phi, err_contour
