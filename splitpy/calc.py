'''
SUBMODULE calc.py

Module containing calculation functions for shear-wave splitting analysis

'''
import numpy as np
from obspy.core import Trace, Stream
from splitpy import conf as cf
from numpy.linalg import inv

def split_SilverChan(trQ, trT, baz, t1, t2):
    """split_SilverChan:

    Calculates splitting based on the minimization of energy on 
    the corrected transverse component (Silver and Chan, 1990)

    Returns:
        Ematrix:    Matrix of T component energy
        trQ_c:      Trace of corrected radial component of motion
        trT_c:      Trace of corrected tangential component of motion
        trFast:     Trace of corrected fast direction of motion
        trSlow:     Trace of corrected slow direction of motion\
        phiSC:      Azimuth of fast axis (deg)
        dttSC:      Delay time between fast and slow axes (sec)
        phi_min:    Azimuth used in plotting routine
    """

    phi = np.arange(-90.0, 90.0, cf.dphi)*np.pi/180.
    dtt = np.arange(0., cf.maxdt, cf.ddt)

    M = np.zeros((2,2,len(phi)))
    M[0,0,:] = np.cos(phi)
    M[0,1,:] = -np.sin(phi)
    M[1,0,:] = np.sin(phi)
    M[1,1,:] = np.cos(phi)

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
        FS_test = np.dot(np.array(M[:,:,p]), np.array([trQ_tmp.data, trT_tmp.data]))

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
            corrected_QT = np.dot(inv(np.array(M[:,:,p])), np.array([tmpFast, tmpSlow]))
            trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
            trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)

            #plot_cmp(trQ_c, trT_c,'phi='+str(phi[p]*180./np.pi)+'dt='+str(shift))

            # Energy on transverse component (component 1)
            Ematrix[p,t] = np.sum(np.square(corrected_QT[1]))

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

    FS_test = np.dot(np.array(M[:,:,ind_phi]), np.array([trQ_tmp.data, trT_tmp.data]))

    F0 = Trace(data=FS_test[0], header=trQ_tmp.stats)
    F1 = Trace(data=FS_test[1], header=trT_tmp.stats)

    tmpFast = tshift(F0, -shift/2.)
    tmpSlow = tshift(F1, shift/2.)

    corrected_QT = np.dot(inv(np.array(M[:,:,ind_phi])), np.array([tmpFast, tmpSlow]))

    trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
    trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)
    
    trFast = Trace(data=tmpFast, header=trT_tmp.stats)
    trSlow = Trace(data=tmpSlow, header=trQ_tmp.stats)
    
    return Ematrix, trQ_c, trT_c, trFast, trSlow, \
            phiSC, shift, phiSC_min


def split_RotCorr(trQ, trT, baz, t1, t2):
    """split_RotCorr:

    Calculates splitting based on the maximum correlation between corrected 
    radial and tangential components of motion 

    Returns:
        Ematrix:    Matrix of correlation coefficients
        trQ_c:      Trace of corrected radial component of motion
        trT_c:      Trace of corrected tangential component of motion
        trFast:     Trace of corrected fast direction of motion
        trSlow:     Trace of corrected slow direction of motion\
        phiRC:      Azimuth of fast axis (deg)
        dttRC:      Delay time between fast and slow axes (sec)
        phi_min:    Azimuth used in plotting routine
    """

    phi = np.arange(-90.0, 90.0, cf.dphi)*np.pi/180.
    dtt = np.arange(0., cf.maxdt, cf.ddt)

    M = np.zeros((2,2,len(phi)))
    M[0,0,:] = np.cos(phi)
    M[0,1,:] = -np.sin(phi)
    M[1,0,:] = np.sin(phi)
    M[1,1,:] = np.cos(phi)

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
        FS_test = np.dot(np.array(M[:,:,p]), np.array([trQ_tmp.data, trT_tmp.data]))

        # Compile into traces
        F0 = Trace(data=np.array(FS_test[0]), header=trQ_tmp.stats)
        F1 = Trace(data=np.array(FS_test[1]), header=trT_tmp.stats)

        # Cross-correlate Fast with Slow
        ns0 = np.sum(F0.data*F0.data)
        ns1 = np.sum(F1.data*F1.data)
        norm = np.sqrt(ns0*ns1)

        cor = Trace(data=np.fft.ifftshift(np.correlate(F0.data, F1.data, mode='same')/norm), header=trQ_tmp.stats)

        # Time shift loop
        for t in range(len(dtt)):

            shift = dtt[t]

            # Shift by dtt each component (+/-)
            cor_pos = tshift(cor, shift)
            cor_neg = tshift(cor, -shift)
            Cmatrix_pos[p,t] = cor_pos[0]
            Cmatrix_neg[p,t] = cor_neg[0]

    # Time shift is positive: fast axis arrives after slow axis
    if abs(Cmatrix_pos).max() > abs(Cmatrix_neg).max():

        #print 'Cmatrix_pos is max'

        ind = np.where(Cmatrix_pos == max(Cmatrix_pos.max(), Cmatrix_pos.min(), key=abs))
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

        ind = np.where(Cmatrix_neg == max(Cmatrix_neg.max(), Cmatrix_neg.min(), key=abs))
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

    M2 = np.zeros((2,2))
    M2[0,0] = np.cos(theta)
    M2[0,1] = -np.sin(theta)
    M2[1,0] = np.sin(theta)
    M2[1,1] = np.cos(theta)

    FS_test = np.dot(np.array(M2[:,:]), np.array([trQ_tmp.data, trT_tmp.data]))

    F0 = Trace(data=FS_test[0], header=trQ_tmp.stats)
    F1 = Trace(data=FS_test[1], header=trT_tmp.stats)

    tmpFast = tshift(F0, shift/2.)
    tmpSlow = tshift(F1, -shift/2.)

    trFast = Trace(data=tmpFast, header=trT_tmp.stats)
    trSlow = Trace(data=tmpSlow, header=trQ_tmp.stats)

    corrected_QT = np.dot(inv(np.array(M2[:,:])), np.array([tmpFast, tmpSlow]))

    trQ_c = Trace(data=corrected_QT[0], header=trQ_tmp.stats)
    trT_c = Trace(data=corrected_QT[1], header=trT_tmp.stats)
    
    return Cmap, trQ_c, trT_c, trFast, trSlow, \
            phiRC, dtRC, phiRC_max


def tshift(trace, tt):
    """tshift:

    Shifts a Trace object with time tt

    Returns:
        rtrace:     Shifted version of trace
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
    """split_dof:

    Determines the degrees of freedom to calculate the
    confidence region of the misfit function

    Returns:
        dof:    Degrees of freedom

    From Walsh, JGR, 2013
    
    """

    ft = np.fft.fft(tr.data)[0:int(len(tr.data)/2)+1]
    ai = np.ones(len(ft))
    ai[0] = 0.5; ai[-1] = 0.5

    E2 = np.sum(np.dot(ai, np.abs(ft)**2))
    E4 = np.sum(np.dot(ai, np.abs(ft)**4))

    dof = 2.*(2.*E2**2/E4 - 1.)

    return dof


def split_errorSC(tr, t1, t2, q, Emat):
    """split_errorSC

    Calculate error bars based on a F-test and 
    a given confidence interval q

    Returns:
        err_dtt:        Error in dt estimate
        err_phi:        Error in phi estimate
        err_contour:    Error contour for plotting

    """
    from scipy import stats

    phi = np.arange(-90.0, 90.0, cf.dphi)*np.pi/180.
    dtt = np.arange(0., cf.maxdt, cf.ddt)

    tr_tmp = tr.copy()
    tr_tmp.trim(t1, t2)
    
    dof = split_dof(tr_tmp)
    n_par = 2

    # Error contour
    vmin = Emat.min()
    vmax = Emat.max()
    err_contour = vmin*(1. + n_par/(dof - n_par)*
            stats.f.ppf(1. - q, n_par, dof - n_par))

    # Estimate uncertainty (q confidence interval)
    err = np.where(Emat<err_contour)
    if len(err) == 0:
      return False, False, False
    #print err, max(err[0]), min(err[0])
    err_phi = 0.5*(phi[max(err[0])] - phi[min(err[0])])*180./np.pi
    err_dtt = 0.5*(dtt[max(err[1])] - dtt[min(err[1])])

    return err_dtt, err_phi, err_contour


def split_errorRC(tr, t1, t2, q, Emat):
    """split_errorRC

    Calculates error bars based on a F-test and 
    a given confidence interval q.

    This version uses a Fisher transformation for 
    correlation-type misfit.

    Returns:
        err_dtt:        Error in dt estimate
        err_phi:        Error in phi estimate
        err_contour:    Error contour for plotting

    """
    from scipy import stats

    phi = np.arange(-90.0, 90.0, cf.dphi)*np.pi/180.
    dtt = np.arange(0., cf.maxdt, cf.ddt)

    tr_tmp = tr.copy()
    tr_tmp.trim(t1, t2)
    
    dof = split_dof(tr_tmp)
    if dof < 3:
      return None, None, None
    n_par = 2

    # Fisher transformation
    vmin = np.arctanh(Emat.min())

    # Error contour
    zrr_contour = vmin + (vmin*np.sign(vmin)*n_par/(dof - n_par)*\
            stats.f.ppf(1. - q, n_par, dof - n_par))*\
            np.sqrt(1./(dof-3))

    # Back transformation
    err_contour = np.tanh(zrr_contour)

    # Estimate uncertainty (q confidence interval)
    err = np.where(Emat<err_contour)
    err_phi = 0.5*(phi[max(err[0])] - phi[min(err[0])])*180./np.pi
    err_dtt = 0.5*(dtt[max(err[1])] - dtt[min(err[1])])

    return err_dtt, err_phi, err_contour
