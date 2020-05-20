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


# -*- coding: utf-8 -*-
import pickle
import stdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
from math import ceil
from splitpy import arguments, Split
from pathlib import Path


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


def main():

    args = arguments.get_arguments_average()

    print("---------------------------")
    print("Selection Criteria ")
    print(" Null Value: ")
    print("    Non Nulls:", args.nons)
    print("    Nulls:    ", args.nulls)
    print(" Quality Value: ")
    print("    Goods: ", args.goods)
    print("    Fairs: ", args.fairs)
    print("    Poors: ", args.poors)
    print("---------------------------")

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

    # Loop over station keys
    for stkey in list(stkeys):

        # Extract station information from dictionary
        sta = db[stkey]

        split = Split(sta)

        # Results directory
        datapath = Path("DATA") / stkey
        if not datapath.is_dir():
            print('Path to ' + str(datapath) + ' doesn`t exist - continuing')
            continue

        evs = [str(x) for x in datapath.iterdir() if x.is_dir()]
        evs.sort()
        ievs = range(0, len(evs))

        # Initialize Storage
        baz = []
        baz_null = []
        phiRC = []
        DphiRC = []
        dtRC = []
        DdtRC = []
        phiSC = []
        DphiSC = []
        dtSC = []
        DdtSC = []
        Qual = []
        Null = []

        print("  Processing {0:d} Events...".format(len(evs)))

        # Loop over Pickle Files and read in required data
        ic = 0
        for iev in ievs:

            # Event Name
            evSTR = evs[iev]

            # Event data
            metafile = Path(evSTR) / "Meta_data.pkl"
            meta = pickle.load(open(metafile, "rb"))
            split.meta = meta

            # Split results
            splitfile = Path(evSTR) / "split_results_auto.pkl"
            file = open(splitfile, "rb")
            split.SC_res = pickle.load(file)
            split.RC_res = pickle.load(file)
            split.null = pickle.load(file)
            split.quality = pickle.load(file)
            file.close()

            # Determine whether to accept based on Null Value
            Naccept = False
            if args.nons and args.nulls:
                Naccept = True
            elif args.nons and not args.nulls:
                if not split.null:
                    Naccept = True
            elif not args.nons and args.nulls:
                if split.null:
                    Naccept = True

            # Determine whether to accept based on Quality Selected
            Qaccept = False
            QGaccept = False
            QFaccept = False
            QPaccept = False
            if args.goods and split.quality == "Good":
                QGaccept = True
            if args.fairs and split.quality == "Fair":
                QFaccept = True
            if args.poors and split.quality == "Poor":
                QPaccept = True
            if QGaccept or QFaccept or QPaccept:
                Qaccept = True

            # Accept Event?
            accept = False
            if Qaccept and Naccept:

                if split.null:
                    print("      {0} Null -> Retained".format(split.quality))
                else:
                    print(
                        "      {0} Non-Null -> Retained".format(split.quality))

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
                    print(
                        "      {0} Non-Null -> Skipped".format(split.quality))

        if not baz:
            return

        # Gridspec for polar plot
        gs1 = gspec.GridSpec(1, 1)
        gs1.update(left=0.57, right=1.0, bottom=0.3, top=0.7)

        baz = np.array([float(i) for i in baz])

        # Get Max DT value
        dtmax = ceil(max([max(dtRC), max(dtSC), 3]))

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
        if args.RCinc:
            ax1.scatter(phi*np.pi/180., dt, c='b')
            ax1.plot(np.array([meanphiRC, meanphiRC])*np.pi /
                     180., [0, meandtRC], 'b', linewidth=2)

        # Convert SC results to floats
        phi = np.array([float(i) for i in phiSC])
        Dphi = np.array([float(i) for i in DphiSC])
        dt = np.array([float(i) for i in dtSC])
        Ddt = np.array([float(i) for i in DdtSC])

        # Calculate average and STD for SC technique
        meanphiSC, stdphiSC, meandtSC, stddtSC = angle_mean(dt, phi, Ddt, Dphi)

        # Add SC results to plot
        if args.SCinc:
            ax1.scatter(phi*np.pi/180., dt, c='coral')
            ax1.plot(np.array([meanphiSC, meanphiSC])*np.pi /
                     180., [0, meandtSC], 'coral', linewidth=2)

        ax1.set_rmax(dtmax)

        # Gridspec for panels
        gs2 = gspec.GridSpec(2, 1)
        gs2.update(left=0.12, right=0.57, bottom=0.2, top=0.8)

        ax2 = plt.subplot(gs2[0])
        ax3 = plt.subplot(gs2[1])

        # Azimuth panel

        phi = np.array([float(i) for i in phiRC])
        Dphi = np.array([float(i) for i in DphiRC])

        # Plot shaded box with RC uncertainty
        if args.RCinc:

            ax2.axhspan(meanphiRC + stdphiRC, meanphiRC -
                        stdphiRC, facecolor='b', alpha=0.2)
            ax2.axhline(meanphiRC, c='b')

            # Plot individual RC results
            ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='b', label='RC')

        phi = np.array([float(i) for i in phiSC])
        Dphi = np.array([float(i) for i in DphiSC])

        # Plot shaded box with SC uncertainty
        if args.SCinc:

            ax2.axhspan(meanphiSC + stdphiSC, meanphiSC -
                        stdphiSC, facecolor='orange', alpha=0.2)
            ax2.axhline(meanphiSC, c='coral')

            # Plot individual SC results
            ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='orange', label='SC')

        ax2.set_title('Station: ' + sta)
        ax2.set_ylabel(r'Fast axis, $\phi$ (degree)')
        ax2.set_ylim(-95, 95)
        ax2.legend(loc=0, numpoints=1)

        # Delay time panel

        dt = np.array([float(i) for i in dtRC])
        Ddt = np.array([float(i) for i in DdtRC])

        # Plot shaded box with RC uncertainty
        if args.RCinc:

            ax3.axhspan(meandtRC + stddtRC, meandtRC -
                        stddtRC, facecolor='b', alpha=0.2)
            ax3.axhline(meandtRC, c='b')

            # Plot individual RC results
            ax3.errorbar(baz, dt, yerr=Ddt, fmt='o', c='b', label='RC')

        dt = np.array([float(i) for i in dtSC])
        Ddt = np.array([float(i) for i in DdtSC])

        # Plot shaded box with SC uncertainty
        if args.SCinc:

            ax3.axhspan(meandtSC + stddtSC, meandtSC -
                        stddtSC, facecolor='orange', alpha=0.2)
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

        # Save plot
        plotdir = Path("PLOTS")
        if not plotdir.is_dir():
            plotdir.mkdir(parents=True)

        # Output Names
        outdata = plotdir / (args.TypeName + args.NullName +
                             args.QualName + "_results.dat")
        outplot = plotdir / (args.TypeName + args.NullName +
                             args.QualName + "_results.png")

        # Final estimates (average of SC and RC)
        if args.RCinc and args.SCinc:
            PHI = (meanphiRC + meanphiSC)/2.
            DT = (meandtRC + meandtSC)/2.
            dPHI = (stdphiRC + stdphiSC)/2.
            dDT = (stddtRC + stddtSC)/2.
            ax1.plot(np.array([PHI, PHI])*np.pi/180.,
                     [0, DT], 'r', linewidth=2, alpha=0.8)
            ax2.axhline(PHI, c='coral', alpha=0.5, linewidth=2)
            ax2.axhline(PHI - dPHI, c='coral', linestyle='--',
                        linewidth=1, alpha=0.5)
            ax2.axhline(PHI + dPHI, c='coral', linestyle='--',
                        linewidth=1, alpha=0.5)
            ax3.axhline(DT, c='coral', alpha=0.5, linewidth=2)
            ax3.axhline(DT - dDT, c='coral', linestyle='--',
                        linewidth=1, alpha=0.5)
            ax3.axhline(DT + dDT, c='coral', linestyle='--',
                        linewidth=1, alpha=0.5)

        elif not args.RCinc and args.SCinc:
            PHI = meanphiSC
            DT = meandtSC
            dPHI = stdphiSC
            dDT = stddtSC
        else:
            PHI = meanphiRC
            DT = meandtRC
            dPHI = stdphiRC
            dDT = stddtRC

        if args.verb:
            print("")
            print("*** Station Average from {0} measurements ***".format(len(baz)))
            print("   Loc: {0:8.4f}, {1:7.4f}".format(stlon, stlat))
            print("   PHI: {0:7.3f} d +- {1:.3f}".format(PHI, dPHI))
            print("   DT:    {0:5.3f} s +- {1:.3f}".format(DT, dDT))
            print("   Saved to: "+str(outdata))
            print("")

        # Write out Final Results
        fid = open(outdata, 'w')
        fid.writelines(
            "{01:8.4f}  {1:7.4f}   {2:7.3f} {3:7.3f}   " +
            "{4:5.3f} {5:5.3f}".format(
                stlon, stlat, PHI, dPHI, DT, dDT))
        fid.close()

        # Save Plot
        plt.savefig(outplot)

        # Display Plot
        if args.showfig:
            plt.show()


if __name__ == "__main__":

    main()
