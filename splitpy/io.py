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

Module containing a function used in one of the the `SplitPy` scripts
that accompany this package.

"""

def list_local_data_stn(lcldrs=list, sta=None, net=None, altnet=[]):
    """
    Function to take the list of local directories and recursively 
    find all data that matches the station name
    
    Parameters
    ----------
    lcldrs : List
        List of local directories
    sta : Dict
        Station metadata from :mod:`~StDb`
    net : str
        Network name
    altnet : List
        List of alternative networks

    Returns
    -------
    fpathmatch : List
        Sorted list of matched directories

    """
    from fnmatch import filter
    from os import walk
    from os.path import join
    
    if sta is None:
        return []
    else:
        if net is None:
            sstrings = ['*.{0:s}.*.SAC'.format(sta)]
        else:
            sstrings = ['*.{0:s}.{1:s}.*.SAC'.format(net, sta)]
            if len(altnet) > 0:
                for anet in altnet:
                    sstrings.append('*.{0:s}.{1:s}.*.SAC'.format(anet, sta))
    
    fpathmatch = []
    # Loop over all local data directories
    for lcldr in lcldrs:
        # Recursiely walk through directory
        for root,dirnames, filenames in walk(lcldr):
            # Keep paths only for those matching the station
            for sstring in sstrings:
                for filename in filter(filenames, sstring):
                    fpathmatch.append(join(root, filename))
            
    fpathmatch.sort()
    
    return fpathmatch