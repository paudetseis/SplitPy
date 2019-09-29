'''
SUBMODULE io.py

Module containing the input/output subroutines

'''


################################################################################
#-- List Available Local Data
################################################################################
def list_local_data_stn(lcldrs=list, sta=None, net=None, altnet=[]):
    '''
    Function to take the list of local directories and recursively find all data that matches the station name
    '''
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
    #-- Loop over all local data directories
    for lcldr in lcldrs:
        #-- Recursiely walk through directory
        for root,dirnames,filenames in walk(lcldr):
            #-- Keep paths only for those matching the station
            for sstring in sstrings:
                for filename in filter(filenames, sstring):
                    fpathmatch.append(join(root, filename))
            
    fpathmatch.sort()
    
    return fpathmatch