import stdb
import numpy as np
from pkg_resources import resource_filename
from obspy.clients.fdsn import Client


def get_stdb():
    dbfile = resource_filename('splitpy',
                               'examples/data/MMPY.pkl')
    db = stdb.io.load_db(dbfile)
    return db


