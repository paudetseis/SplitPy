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
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Program sks_plot_results.py
---------------------------

Plots the results for a given station based on ``.pkl`` files present in the
running directory.

Usage
-----

.. code-block::

   $ sks_plot_results.py -h
    Usage: sks_plot_results.py [options] <Folder> [Folder 2]

    Script to plot the average splitting results for a given station. Loads the
    available pkl files in the specified Station Directory.

    Options:
      -h, --help            show this help message and exit
      --no-figure           Specify to prevent plots from opening during
                            processing; they are still saved to disk. [Default
                            plots open and save]

      Null Selection Settings:
        Settings associated with selecting which Null or Non-Null data is
        included

        --nulls, --Nulls    Specify this flag to include Null Values in the
                            average. [Default Non-Nulls only]
        --no-nons, --No-Nons
                            Specify this flag to exclude Non-Nulls from the
                            average [Default False]

      Quality Selection Settings:
        Settings associated with selecting the qualities to include in the
        selection.

        --No-Good, --no-good
                            Specify to exclude 'Good' measurements from the
                            average. [Default Good + Fair]
        --No-Fair, --no-fair
                            Specify to exclude 'Fair' measurements from the
                            average [Default Good + Fair]
        --Poor, --poor      Specify to include 'Poor' measurements in the average
                            [Default No Poors]

      Split Type Settings:
        Settings to Select which Split types are included in the selection.

        --RC-Only, --rc-only, --RC-only
                            Specify to only include RC splits in the average.
                            [Default RC + SC]
        --SC-Only, --sc-only, --SC-only
                            Specify to only include SC splits in the average.
                            [Default RC + SC]
"""

# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gspec
from glob import glob
from os.path import exists, join
from math import ceil
import pickle


def get_options():
    from optparse import OptionParser, OptionGroup

    parser = OptionParser(usage="Usage: %prog [options] <Folder> [Folder 2]", \
        description="Script to plot the average splitting results for a given station. " \
        "Loads the available pkl files in the specified Station Directory.")

    # General Settings
    parser.add_option("--no-figure", action="store_false", dest="showfig", default=True, \
        help="Specify to prevent plots from opening during processing; they are still " \
        "saved to disk. [Default plots open and save]")

    # Null Settings
    NullGroup = OptionGroup(parser, title="Null Selection Settings", description="Settings " \
        "associated with selecting which Null or Non-Null data is included")
    NullGroup.add_option("--nulls", "--Nulls", action="store_true", dest="nulls", default=False, \
        help="Specify this flag to include Null Values in the average. [Default Non-Nulls only]")
    NullGroup.add_option("--no-nons", "--No-Nons", action="store_false", dest="nons", default=True, \
        help="Specify this flag to exclude Non-Nulls from the average [Default False]")

    # Quality Settings
    QualGroup = OptionGroup(parser, title="Quality Selection Settings", \
        description="Settings associated with selecting the qualities to include in the selection.")
    QualGroup.add_option("--No-Good", "--no-good", action="store_false", dest="goods", default=True, \
        help="Specify to exclude 'Good' measurements from the average. [Default Good + Fair]")
    QualGroup.add_option("--No-Fair", "--no-fair", action="store_false", dest="fairs", default=True, \
        help="Specify to exclude 'Fair' measurements from the average [Default Good + Fair]")
    QualGroup.add_option("--Poor", "--poor", action="store_true", dest="poors", default=False, \
        help="Specify to include 'Poor' measurements in the average [Default No Poors]")

    # Split Type Settings
    SpTypGroup = OptionGroup(parser, title="Split Type Settings", description="Settings to Select " \
        "which Split types are included in the selection.")
    SpTypGroup.add_option("--RC-Only", "--rc-only", "--RC-only", action="store_false", dest="SCinc", \
        default=True, help="Specify to only include RC splits in the average. [Default RC + SC]")
    SpTypGroup.add_option("--SC-Only", "--sc-only", "--SC-only", action="store_false", dest="RCinc", \
        default=True, help="Specify to only include SC splits in the average. [Default RC + SC]")

    parser.add_option_group(NullGroup)
    parser.add_option_group(QualGroup)
    parser.add_option_group(SpTypGroup)
    (opts, args) = parser.parse_args()

    # Check Folders
    if len(args) < 1:
        parser.error("Must include input folder to process.")

    # Check Nulls
    if not opts.nons and not opts.nulls:
        parser.error("One of Non-Nulls or Nulls must be included.")

    # Check Quality
    if not opts.goods and not opts.fairs and not opts.poors:
        parser.error("At least one Quality must be included.")

    # Check Types
    if not opts.RCinc and not opts.SCinc:
        parser.error("At leat one Splitting Tyhpe must be included.")

    # Construct Null FileName Components
    NullName = ""
    if opts.nons:
        NullName = "_Nons"
        if opts.nulls:
            NullName = NullName + "-Nulls"
    else:
        if opts.nulls: NullName = "_Nulls"
    opts.NullName = NullName

    # Construct Quality FileName Components
    QualName = ""
    if opts.goods:
        QualName = "_G"
        if opts.fairs: QualName = QualName + "-F"
        if opts.poors: QualName = QualName + "-P"
    else:
        if opts.fairs:
            QualName = "_F"
            if opts.poors: QualName = QualName + "-P"
        else:
            if opts.poors: QualName = "_P"
    opts.QualName = QualName

    # Construct Type FileName Components
    TypeName = ""
    if opts.RCinc and opts.SCinc:
        TypeName = "_RC-SC"
    elif opts.RCinc and not opts.SCinc:
        TypeName = "_RC"
    elif not opts.RCinc and opts.SCinc:
        TypeName = "_SC"
    opts.TypeName = TypeName

    return args, opts


def angle_mean(dt, phi, ddt, dphi):
    x = dt*np.cos(2*phi*np.pi/180.0)
    y = dt*np.sin(2*phi*np.pi/180.0)
    c = x + 1j*y
    m = np.mean(c)
    
    phase = np.angle(m, deg=True)/2.
    radius = np.abs(m)
    dphase = np.sqrt(np.sum(dphi**2))/len(x)
    dradius = np.sqrt(np.sum(ddt**2))/len(x)

    return phase, dphase, radius, dradius


def process(sta,opts):

    pickles = glob(join(arg, "Split.*.pkl"))
    pickles.sort()

    # Initialize Storage
    baz = []; baz_null = []; phiRC = []; DphiRC = []; dtRC = []; DdtRC = []
    phiSC = []; DphiSC = []; dtSC = []; DdtSC = []; Qual = []; Null = []

    print("  Processing {0:d} Events...".format(len(pickles)))

    # Loop over Pickle Files and read in required data
    ic = 0
    for picklef in pickles:

        # Get Pickle File Contents
        ic += 1
        print(("  {0:d}) {1:s}".format(ic, picklef)))
        try:
            fpkl = open(picklef, 'rb')
        except:
            print("     Error Opening")
            continue

        # Read in all data
        split = pickle.load(fpkl)
        fpkl.close()

        # Determine whether to accept based on Null Value
        Naccept = False
        if opts.nons and opts.nulls:
            Naccept = True
        elif opts.nons and not opts.nulls:
            if not split.null: Naccept = True
        elif not opts.nons and opts.nulls:
            if split.null: Naccept = True

        # Determine whether to accept based on Quality Selected
        Qaccept = False
        QGaccept = False
        QFaccept = False
        QPaccept = False
        if opts.goods and split.quality == "Good": QGaccept = True
        if opts.fairs and split.quality == "Fair": QFaccept = True
        if opts.poors and split.quality == "Poor": QPaccept = True
        if QGaccept or QFaccept or QPaccept: Qaccept = True

        # Accept Event?
        accept = False
        if Qaccept and Naccept:

            if split.null:
                print("      {0} Null -> Retained".format(split.quality))
            else:
                print("      {0} Non-Null -> Retained".format(split.quality))

            # BAZ, null and quality
            baz.append(split.meta.baz)
            Qual.append(split.quality)
            Null.append(split.null)

            # RC Results
            phiRC.append(split.RC_res.phi)
            DphiRC.append(split.RC_res.ephi)
            dtRC.append(split.RC_res.dtt)
            DdtRC.append(split.RC_res.edtt)

            # SC Results
            phiSC.append(split.SC_res.phi)
            DphiSC.append(split.SC_res.ephi)
            dtSC.append(split.SC_res.dtt)
            DdtSC.append(split.SC_res.edtt)

            stlat = split.sta.latitude
            stlon = split.sta.longitude
            sta = split.sta.station

        else:

            if split.null:
                print("      {0} Null -> Skipped".format(split.quality))
            else:
                print("      {0} Non-Null -> Skipped".format(split.quality))


    if not baz: return

    # Gridspec for polar plot
    gs1 = gspec.GridSpec(1, 1)
    gs1.update(left=0.57, right=1.0, bottom=0.3, top=0.7)

    baz = np.array([float(i) for i in baz])

    # Get Max DT value
    dtmax = ceil(max([max(dtRC), max(dtSC),3]))

    # Convert RC results to floats
    phi = np.array([float(i) for i in phiRC])
    Dphi = np.array([float(i) for i in DphiRC])
    dt = np.array([float(i) for i in dtRC])
    Ddt = np.array([float(i) for i in DdtRC])

    # Calculate average and STD for RC technique
    meanphiRC, stdphiRC, meandtRC, stddtRC = angle_mean(dt, phi, Ddt, Dphi)

    # Prepare Polar Plot
    ax1 = plt.subplot(gs1[0], projection='polar', aspect=1)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)

    # # Add Bazs to Plot
    # ax1.plot(baz*np.pi/180., np.ones(len(baz))*dtmax, 'co', alpha=0.5)

    # Add RC results to plot
    if opts.RCinc:
        ax1.scatter(phi*np.pi/180., dt, c='b')
        ax1.plot(np.array([meanphiRC, meanphiRC])*np.pi/180., [0, meandtRC], 'b', linewidth=2)

    # Convert SC results to floats
    phi = np.array([float(i) for i in phiSC])
    Dphi = np.array([float(i) for i in DphiSC])
    dt = np.array([float(i) for i in dtSC])
    Ddt = np.array([float(i) for i in DdtSC])

    # Calculate average and STD for SC technique
    meanphiSC, stdphiSC, meandtSC, stddtSC = angle_mean(dt, phi, Ddt, Dphi)
    
    # Add SC results to plot
    if opts.SCinc:
        ax1.scatter(phi*np.pi/180., dt, c='coral')
        ax1.plot(np.array([meanphiSC, meanphiSC])*np.pi/180., [0, meandtSC], 'coral', linewidth=2)

    ax1.set_rmax(dtmax)

    # Gridspec for panels
    gs2 = gspec.GridSpec(2, 1)
    gs2.update(left=0.12, right=0.57, bottom=0.2, top=0.8)

    ax2 = plt.subplot(gs2[0])
    ax3 = plt.subplot(gs2[1])
    
    ##### Azimuth panel
 
    phi = np.array([float(i) for i in phiRC])
    Dphi = np.array([float(i) for i in DphiRC])
    
    # Plot shaded box with RC uncertainty
    if opts.RCinc:

        ax2.axhspan(meanphiRC + stdphiRC, meanphiRC - stdphiRC, facecolor='b', alpha=0.2)
        ax2.axhline(meanphiRC, c='b')

        # Plot individual RC results
        ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='b', label='RC')

    phi = np.array([float(i) for i in phiSC])
    Dphi = np.array([float(i) for i in DphiSC])
    
    # Plot shaded box with SC uncertainty
    if opts.SCinc:

        ax2.axhspan(meanphiSC + stdphiSC, meanphiSC - stdphiSC, facecolor='orange', alpha=0.2)
        ax2.axhline(meanphiSC, c='coral')

        # Plot individual SC results
        ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='orange', label='SC')

    ax2.set_title('Station: ' + sta)
    ax2.set_ylabel(r'Fast axis, $\phi$ (degree)')
    ax2.set_ylim(-95, 95)
    ax2.legend(loc=0, numpoints=1)
    
    ##### Delay time panel
    
    dt = np.array([float(i) for i in dtRC])
    Ddt = np.array([float(i) for i in DdtRC])
    
    # Plot shaded box with RC uncertainty
    if opts.RCinc:

        ax3.axhspan(meandtRC + stddtRC, meandtRC - stddtRC, facecolor='b', alpha=0.2)
        ax3.axhline(meandtRC, c='b')

        # Plot individual RC results
        ax3.errorbar(baz, dt, yerr=Ddt, fmt='o', c='b', label='RC')

    dt = np.array([float(i) for i in dtSC])
    Ddt = np.array([float(i) for i in DdtSC])
    
    # Plot shaded box with SC uncertainty
    if opts.SCinc:

        ax3.axhspan(meandtSC + stddtSC, meandtSC - stddtSC, facecolor='orange', alpha=0.2)
        ax3.axhline(meandtSC, c='coral')

        # Plot individual SC results
        ax3.errorbar(baz, dt, yerr=Ddt, fmt='o', c='orange', label='SC')

    ax3.set_ylabel(r'Delay time, $\delta t$ (seconds)')
    ax3.set_ylim(0, dtmax)
    
    ax3.set_xlabel('Back-azimuth (degree)')
    ax3.legend(loc=0, numpoints=1)
 
    # Add labels    
    plt.figtext(0.015, 0.80, 'A', fontweight='bold')
    plt.figtext(0.015, 0.48, 'B', fontweight='bold')
    plt.figtext(0.6, 0.7, 'C', fontweight='bold')

    # Output Names
    outname = arg + opts.TypeName + opts.NullName + opts.QualName + "_results"

    # Final estimates (average of SC and RC)
    if opts.RCinc and opts.SCinc:
        PHI = (meanphiRC + meanphiSC)/2.
        DT = (meandtRC + meandtSC)/2.
        dPHI = (stdphiRC + stdphiSC)/2.
        dDT = (stddtRC + stddtSC)/2.
        ax1.plot(np.array([PHI, PHI])*np.pi/180.,[0, DT],'r',linewidth=2,alpha=0.8)
        ax2.axhline(PHI, c='coral', alpha=0.5, linewidth=2)
        ax2.axhline(PHI - dPHI, c='coral', linestyle='--', linewidth=1, alpha=0.5)
        ax2.axhline(PHI + dPHI, c='coral', linestyle='--', linewidth=1, alpha=0.5)
        ax3.axhline(DT, c='coral', alpha=0.5, linewidth=2)
        ax3.axhline(DT - dDT, c='coral', linestyle='--', linewidth=1, alpha=0.5)
        ax3.axhline(DT + dDT, c='coral', linestyle='--', linewidth=1, alpha=0.5)

    elif not opts.RCinc and opts.SCinc:
        PHI = meanphiSC
        DT = meandtSC
        dPHI = stdphiSC
        dDT = stddtSC
    else:
        PHI = meanphiRC
        DT = meandtRC
        dPHI = stdphiRC
        dDT = stddtRC

    print("")
    print("*** Station Average from {0} measurements ***".format(len(baz)))
    print("   " + arg)
    print("   Loc: {0:8.4f}, {1:7.4f}".format(stlon, stlat))
    print("   PHI: {0:7.3f} d +- {1:.3f}".format(PHI, dPHI))
    print("   DT:    {0:5.3f} s +- {1:.3f}".format(DT, dDT))
    print("   Saved to: "+join(outname))
    print("")

    # Write out Final Results
    fid = open(join(outname) + ".dat",'w')
    fid.writelines("{0}  {1:8.4f}  {2:7.4f}   {3:7.3f} {4:7.3f}   {5:5.3f} {6:5.3f}".format(\
        arg, stlon, stlat, PHI, dPHI, DT, dDT))
    fid.close()

    # Save Plot
    plt.savefig(join(outname)+'.png')

    # Display Plot
    if opts.showfig: plt.show()


if __name__ == "__main__":

    args, opts = get_options()

    print("---------------------------")
    print("Selection Criteria ")
    print(" Null Value: ")
    print("    Non Nulls:", opts.nons)
    print("    Nulls:    ", opts.nulls)
    print(" Quality Value: ")
    print("    Goods: ", opts.goods)
    print("    Fairs: ", opts.fairs)
    print("    Poors: ", opts.poors)
    print("---------------------------")


    for arg in args:

        # Check Folder Exists
        if not exists(arg):
            print(("Warning: " + arg + " does not exist"))
            continue

        print(("Working on Station: " + arg))

        # Process the Folder
        process(arg, opts)


