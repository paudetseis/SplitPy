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
import matplotlib
from math import ceil

from splitpy import Split

from argparse import ArgumentParser
from os.path import exists as exist
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

def get_arguments_average(argv=None):

    parser = ArgumentParser(
        usage="%(prog)s [arguments] <station database>",
        description="Script to plot the average splitting results for a " +
        "given station. Loads the available .pkl files in the specified " +
        "Station Directory.")

    # General Settings
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma separated list of station keys " +
        "for which to perform analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the " +
        "dictionary. For instance, providing IU will match " +
        "with all stations in the IU network [Default " +
        "processes all stations in the database]")
    parser.add_argument(
        "-V", "--verbose",
        action="store_true",
        dest="verb",
        default=False,
        help="Specify to increase verbosity.")
    parser.add_argument(
        "--show-fig",
        action="store_true",
        dest="showfig",
        default=False,
        help="Specify show plots during processing - " +
        "they are still saved to disk. [Default only saves]")
    parser.add_argument(
        "--auto",
        action="store_true",
        dest="auto",
        default=False,
        help="Specify to use automatically processed split results. "+
        "[Default uses refined ('manual') split results]")

    # Null Settings
    NullGroup = parser.add_argument_group(
        title="Null Selection Settings",
        description="Settings "
        "associated with selecting which Null or Non-Null data is included")
    NullGroup.add_argument(
        "--nulls",
        action="store_true",
        dest="nulls",
        default=False,
        help="Specify this flag to include Null Values in the average. " +
        "[Default Non-Nulls only]")
    NullGroup.add_argument(
        "--no-nons",
        action="store_false",
        dest="nons",
        default=True,
        help="Specify this flag to exclude Non-Nulls from the average " +
        "[Default False]")

    # Quality Settings
    QualGroup = parser.add_argument_group(
        title="Quality Selection Settings",
        description="Settings associated with selecting the qualities " +
        "to include in the selection.")
    QualGroup.add_argument(
        "--no-good",
        action="store_false",
        dest="goods",
        default=True,
        help="Specify to exclude 'Good' measurements from the average. " +
        "[Default Good + Fair]")
    QualGroup.add_argument(
        "--no-fair",
        action="store_false",
        dest="fairs",
        default=True,
        help="Specify to exclude 'Fair' measurements from the average " +
        "[Default Good + Fair]")
    QualGroup.add_argument(
        "--poor",
        action="store_true",
        dest="poors",
        default=False,
        help="Specify to include 'Poor' measurements in the average " +
        "[Default No Poors]")

    # Split Type Settings
    SpTypGroup = parser.add_argument_group(
        title="Split Type Settings",
        description="Settings to Select "
        "which Split types are included in the selection.")
    SpTypGroup.add_argument(
        "--RC-only",
        action="store_false",
        dest="SCinc",
        default=True,
        help="Specify to only include RC splits in the average. " +
        "[Default RC + SC]")
    SpTypGroup.add_argument(
        "--SC-only",
        action="store_false",
        dest="RCinc",
        default=True,
        help="Specify to only include SC splits in the average. " +
        "[Default RC + SC]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # Create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')

    # Check Nulls
    if not args.nons and not args.nulls:
        parser.error("One of Non-Nulls or Nulls must be included.")

    # Check Quality
    if not args.goods and not args.fairs and not args.poors:
        parser.error("At least one Quality must be included.")

    # Check Types
    if not args.RCinc and not args.SCinc:
        parser.error("At leat one Splitting Tyhpe must be included.")

    # Construct Null FileName Components
    NullName = ""
    if args.nons:
        NullName = "_Nons"
        if args.nulls:
            NullName = NullName + "-Nulls"
    else:
        if args.nulls:
            NullName = "_Nulls"
    args.NullName = NullName

    # Construct Quality FileName Components
    QualName = ""
    if args.goods:
        QualName = "_G"
        if args.fairs:
            QualName = QualName + "-F"
        if args.poors:
            QualName = QualName + "-P"
    else:
        if args.fairs:
            QualName = "_F"
            if args.poors:
                QualName = QualName + "-P"
        else:
            if args.poors:
                QualName = "_P"
    args.QualName = QualName

    # Construct Type FileName Components
    TypeName = ""
    if args.RCinc and args.SCinc:
        TypeName = "_RC-SC"
    elif args.RCinc and not args.SCinc:
        TypeName = "_RC"
    elif not args.RCinc and args.SCinc:
        TypeName = "_SC"
    args.TypeName = TypeName

    return args

def main(args=None):

    print()
    print("###############################################################")
    print("#            _ _ _                                            #")
    print("#  ___ _ __ | (_) |_     __ ___   _____ _ __ __ _  __ _  ___  #")
    print("# / __| '_ \| | | __|   / _` \ \ / / _ \ '__/ _` |/ _` |/ _ \ #")
    print("# \__ \ |_) | | | |_   | (_| |\ V /  __/ | | (_| | (_| |  __/ #")
    print("# |___/ .__/|_|_|\__|___\__,_| \_/ \___|_|  \__,_|\__, |\___| #")
    print("#     |_|          |_____|                        |___/       #")
    print("#                                                             #")
    print("###############################################################")
    print()

    if args is None:
        # Run Input Parser
        args = get_arguments_average()

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
    # stdb=0.1.4
    try:
        db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)

    # stdb=0.1.3
    except:
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
        phiRC_min = []
        DphiRC = []
        dtRC = []
        DdtRC = []
        EmatRC = []
        phiSC = []
        phiSC_min = []
        DphiSC = []
        dtSC = []
        DdtSC = []
        EmatSC = []
        evlon = []
        evlat = []
        evmag = []
        evtime = []
        Qual = []
        Null = []

        print("  Found {0:d} event folders...".format(len(evs)))
        if args.auto:
            print("  Checking 'auto' results...")
        else:
            print("  Checking 'manual' results...")

        # Loop over Pickle Files and read in required data
        ic = 0
        for iev in ievs:

            # Event Name
            evSTR = evs[iev]

            # Event data
            metafile = Path(evSTR) / "Meta_data.pkl"
            meta = pickle.load(open(metafile, "rb"))
            split.meta = meta
            maxdt = meta.maxdt
            ddt = meta.ddt
            dphi = meta.dphi

            # Split results
            if args.auto:
                splitfile = Path(evSTR) / "split_results_auto.pkl"
            else:
                splitfile = Path(evSTR) / "split_results_manual.pkl"
            if not splitfile.exists():
                continue
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

                if args.verb:
                    if split.null:
                        print("      {0} {1} Null -> Retained".format(
                            str(Path(evSTR).name), split.quality))
                    else:
                        print("      {0} {1} Non-Null -> Retained".format(
                            str(Path(evSTR).name), split.quality))

                # BAZ, null and quality
                baz.append(split.meta.baz)
                Qual.append(split.quality)
                Null.append(split.null)

                # RC Results
                phiRC.append(split.RC_res.phi)
                phiRC_min.append(split.RC_res.phi_min)
                DphiRC.append(split.RC_res.ephi)
                dtRC.append(split.RC_res.dtt)
                DdtRC.append(split.RC_res.edtt)
                EmatRC.append(split.RC_res.Emat)

                # SC Results
                phiSC.append(split.SC_res.phi)
                phiSC_min.append(split.SC_res.phi_min)
                DphiSC.append(split.SC_res.ephi)
                dtSC.append(split.SC_res.dtt)
                DdtSC.append(split.SC_res.edtt)
                EmatSC.append(split.SC_res.Emat)

                # Station info
                stlat = split.sta.latitude
                stlon = split.sta.longitude
                sta = split.sta.station

                # Event info
                evlat.append(split.meta.lat)
                evlon.append(split.meta.lon)
                evmag.append(split.meta.mag)
                evtime.append(split.meta.time)

            else:

                if args.verb:
                    if split.null:
                        print("      {0} {1} Null -> Skipped".format(
                            str(Path(evSTR).name), split.quality))
                    else:
                        print("      {0} {1} Non-Null -> Skipped".format(
                            str(Path(evSTR).name), split.quality))

        if len(baz) == 0:
            print("  No splitting results to average")
            return

        # Average error surfaces
        dt = np.arange(0., maxdt, ddt)
        phi = np.arange(-90., 90., dphi)

        extent = [phi.min(), phi.max(), dt.min(), dt.max()]
        X, Y = np.meshgrid(dt, phi)
        cmap = plt.cm.RdYlBu_r

        # RC analysis
        EmatRC_mean = 0.
        for E, pRC, pRC_min in zip(EmatRC, phiRC, phiRC_min):
            E2 = np.roll(E, int(pRC - pRC_min), axis=0)
            EmatRC_mean += E2

        EmatRC_mean = EmatRC_mean/len(EmatRC)

        # Find indices of minimum value of Energy matrix
        ind = np.where(EmatRC_mean == EmatRC_mean.min())
        ind_phi = ind[0][0]
        ind_dtt = ind[1][0]

        # Get best-fit phi and dt
        dtt_bestRC = dt[ind_dtt]
        phi_bestRC = phi[ind_phi]
        print('\n  Best dt RC: {0:.1f} sec; Best phi RC: {1:.1f} deg'.format(dtt_bestRC, phi_bestRC))

        Emin = EmatRC_mean.min()
        Emax = EmatRC_mean.max()
        dE = (Emax - Emin)/16.
        levelsRC = np.arange(Emin, Emax, dE)

        # SC analysis
        EmatSC_mean = 0.

        for E, pSC, pSC_min in zip(EmatRC, phiSC, phiSC_min):
            E2 = np.roll(E, int(pSC - pSC_min), axis=0)
            EmatSC_mean += E2

        EmatSC_mean = EmatSC_mean/len(EmatSC)

        # Find indices of minimum value of Energy matrix
        ind = np.where(EmatSC_mean == EmatSC_mean.min())
        ind_phi = ind[0][0]
        ind_dtt = ind[1][0]

        # Get best-fit phi and dt
        dtt_bestSC = dt[ind_dtt]
        phi_bestSC = phi[ind_phi]
        print('  Best dt SC: {0:.1f} sec; Best phi SC: {1:.1f} deg\n'.format(dtt_bestSC, phi_bestSC))

        Emin = EmatSC_mean.min()
        Emax = EmatSC_mean.max()
        dE = (Emax - Emin)/16.
        levelsSC = np.arange(Emin, Emax, dE)

        # Plot
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(6,3))
        ax1.contour(
            X, Y, EmatSC_mean, levelsSC,
            cmap=plt.cm.get_cmap(cmap, len(levelsSC)))
        ax1.axvline(dtt_bestSC)
        ax1.axhline(phi_bestSC)
        ax1.set_title('Silver-Chan')
        ax1.set_xlabel(r'$\delta t$ (sec)')
        ax1.set_ylabel(r'$\phi$ (deg)')
        ax2.contour(
            X, Y, EmatRC_mean, levelsRC,
            cmap=plt.cm.get_cmap(cmap, len(levelsRC)))
        ax2.axvline(dtt_bestRC)
        ax2.axhline(phi_bestRC)
        ax2.set_title('Rotation-Correlation')
        ax2.set_xlabel(r'$\delta t$ (sec)')
        plt.tight_layout()

        # Display Plot
        if args.showfig:
            plt.show()
        else:
            plt.close()

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
            ax1.scatter(phi*np.pi/180., dt, c='b', alpha=0.5)
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
            ax1.scatter(phi*np.pi/180., dt, c='orange', alpha=0.5)
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
            ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='b', label='RC', alpha=0.5)

        phi = np.array([float(i) for i in phiSC])
        Dphi = np.array([float(i) for i in DphiSC])

        # Plot shaded box with SC uncertainty
        if args.SCinc:

            ax2.axhspan(meanphiSC + stdphiSC, meanphiSC -
                        stdphiSC, facecolor='orange', alpha=0.2)
            ax2.axhline(meanphiSC, c='coral')

            # Plot individual SC results
            ax2.errorbar(baz, phi, yerr=Dphi, fmt='o', c='orange', label='SC', alpha=0.5)

        ax2.set_title('Station: ' + stkey)
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
            ax3.errorbar(baz, dt, yerr=Ddt, fmt='o', c='b', label='RC', alpha=0.5)

        dt = np.array([float(i) for i in dtSC])
        Ddt = np.array([float(i) for i in DdtSC])

        # Plot shaded box with SC uncertainty
        if args.SCinc:

            ax3.axhspan(meandtSC + stddtSC, meandtSC -
                        stddtSC, facecolor='orange', alpha=0.2)
            ax3.axhline(meandtSC, c='coral')

            # Plot individual SC results
            ax3.errorbar(baz, dt, yerr=Ddt, fmt='o', c='orange', label='SC', alpha=0.5)

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
        outav = plotdir / (stkey+args.TypeName + args.NullName +
                           args.QualName + "_results.dat")
        outev = plotdir / (stkey+args.TypeName + args.NullName +
                           args.QualName + "_events.dat")
        outplot = plotdir / (stkey+args.TypeName + args.NullName +
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
            print(
                "*** Station Average from {0} measurements ***".format(
                    len(baz)))
            print("   Loc: {0:8.4f}, {1:7.4f}".format(stlon, stlat))
            print("   PHI: {0:7.3f} d +- {1:.3f}".format(PHI, dPHI))
            print("   DT:    {0:5.3f} s +- {1:.3f}".format(DT, dDT))
            print("   Saved to: "+str(outav))
            print("")
            print(
                "*** Catalogue of events and results ***")
            print("   Saved to: "+str(outev))
            print("")

        # Write out Final Results
        fid = open(outav, 'w')
        fid.writelines(
            "sta_lon, sta_lat, PHI (deg), dPHI (deg), DT (s), dDT (s) \n")
        line = "{0:8.4f},{1:7.4f},{2:7.3f},{3:7.3f},{4:5.3f},{5:5.3f}\n".format(
                stlon, stlat, PHI, dPHI, DT, dDT)
        fid.writelines(line.replace(" ",""))
        fid.close()

        # Write out events
        fid = open(outev, 'w')
        fid.writelines(
            "eq_time (UTC), eq_lon, eq_lat, eq_mag, PHI_RC (deg), " +
            "dPHI_RC (deg), DT_RC (s), dDT_RC (s), " +
            "PHI_SC (deg), dPHI_SC (deg), DT_SC (s), dDT_SC (s)\n")
        for i in range(len(evlon)):
            line1 = "{0},{1:8.4f},{2:7.4f},{3:3.1f},".format(
                    evtime[i], evlon[i], evlat[i], evmag[i])
            line2 = "{0:7.3f},{1:7.3f},{2:5.3f},{3:5.3f},".format(
                    phiRC[i], DphiRC[i], dtRC[i], DdtRC[i])
            line3 = "{0:7.3f},{1:7.3f},{2:5.3f},{3:5.3f}\n".format(
                    phiSC[i], DphiSC[i], dtSC[i], DdtSC[i])
            fid.writelines(
                line1.replace(" ", "") + line2.replace(" ",
                                                       "") + line3.replace(" ", ""))

        fid.close()

        # Save Plot
        plt.savefig(outplot)

        # Display Plot
        if args.showfig:
            plt.show()
        else:
            plt.close()


if __name__ == "__main__":

    main()
