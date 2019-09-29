SplitPy
=======

SKS Splitting Toolbox based on the Matlab Tool SplitLab, but with modifications from
Wustefeld et al (2008). Additional error surface implementation has been added, however
these error surfaces have not been fully tested (as of 1.1.4).


Suggested Installation
----------------------

First Time Installing:

```bash
$ sudo pip install .
```

Installation after local or remote (gitlab/github) code updates

```bash
$ sudo pip install --upgrade .
```

Note that if the code will be worked on or updated semi-frequently, then I suggest 
modifying your installation to include the pip `-e` flag, which enables an editable 
installation. Essentially, installation is done using symbolic links back to the 
source folder, so any updates that are made in the local folder are updated in the 
path as well.


Package Contents
----------------

splitpy: Python Module
* data download and reading tools
* splitting calculation tools

Scripts: Python scripts designed to be run on the command line
* `sks_split.py`: all in one sks splitting calculations based on stations 
contained within a station database file
* `sks_plot_results.py`: plot the resulting split files for a given station
* `sks_prep.py`: download and prepare the data for computing splitting. Saves 
seismograms for computing splitting at a later time. Used in conjunction with 
`sks_OL_proc.py`.
* `sks_OL_proc.py`: compute the splitting for an archive of data that has been 
downloaded and prepared by sks_prep.py. This script can be run offline without 
an internet connection, but data must be in a specific format.


Dependencies and Requirments:
-----------------------------------------
* obspy (recommended >= 1.0.x) (plus all dependencies)
* stdb
* pickle
* dill
* scipy


Examples:
-----------------

Perform SKS splitting measurements on station EPYK in a database.

```bash
andrew@Iridium:~$ query_fdsn_stdb.py -N TA -C BH? -S EPYK new_list.pkl

andrew@Iridium:~$ sks_split.py --keys=TA.EPYK --local-data=/mnt/datadisk/DaySac/ new_list.pkl
```

This uses all default settings for window lengths, magnitude criteria, etc. 
In this example, data will be used from both IRIS as well as any local data 
on disk (defined with the `--local-data` flag). If no data exists on disk, then 
the program will search on the specific data sever (through `obspy` clients).

Based on the criteria specified (see `-h`), events will be processed. Processing 
will proceed for an event where the minimum SNR threshold is exceeded. Two Figure 
windows will pop up. Figure 1i is the three components for the waveforms, LQT, 
along with lins representing the SKS, SKKS, S and ScS arrivals. Red vertical 
lines denote the window.

Figure 2 summarizes the results of the splitting calculation. The top left "Q,T" 
frame shows the uncorrected Q and T components within the time window. The second 
row of panels correspond to the Rotation Correlation Results, and the third row of 
panels is for the Silver and Chan Results. In each case, the first Column shows 
the correctex Q and T fast and slow components, the second column the corrected 
Q and T components, the third column the before and after particle motion, and 
the fourth column the map of the error surfaces.

A message box will pop up asking whether to Re-pick the window. This can be done 
to refine the signal window in which the measurements are made. This can help 
improve the measurements.

The terminal will show a summary of the processing, including an examination for the 
Null/Non-Null classification as well as the quality of the estimate.

Once `'No'` is selected for the picking/re-picking of the window, a second box will 
pop up asking whether to keep the estimates. Click `'yes'` to save the results, 
or `'No'` to discard the measurement.

The results of processing are saved into a RESULTS folder in the current working 
directory, in a subfolder named after the station key. In this example, TA.EPYK.

Each measurement is stored in a separate `pickle` (.pkl) file named 
`Split.STN.YYYY.JJJ.HHMMSS.pkl`.

Plotting and subsequent processing of splitting results is carried out using 
sks_plot_results.py, where options are present to control selection of nulls 
and quality settings, as well as which methods are used. The path to a station 
`RESULT` folder is provided and all `Split.*.pkl` files are loaded. The final 
average splits are then saved in a text file for future use.


