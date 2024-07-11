#%%

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
import swiftbat
import batanalysis as ba
import swifttools.swift_too as swt

# David Palmer machines only
try:
    from dmptools import proxyfix
    proxyfix()
except:
    pass

ba.datadir("~/data/mirror/swift", makepersistent=True, tdrss="tdrss")
mspattern="bat/rate/sw*brtms.lc.gz"
binedges = np.array([15,25,50,100,350])


sourcename = "Swift J1727.8-1613"
# Just the original huge peak (without its tail)
timerange = [swiftbat.string2datetime(t) for t in ("MJD60186", "MJD60200")]
minexposure = 1000     # cm^2 after cos adjust


# %%

sourceloc = swiftbat.simbadlocation(sourcename)
source_batsource = swiftbat.source(ra=sourceloc[0], dec=sourceloc[1], name=sourcename)
queryargs = dict(time=f"{timerange[0]:%Y-%m-%d} .. {timerange[1]:%Y-%m-%d}", fields='All', resultmax=0)
table_everything = ba.from_heasarc(**queryargs)

exposures = np.array([source_batsource.exposure(ra=row['RA'], dec=row['DEC'], roll=row['ROLL_ANGLE'])[0] for row in table_everything])
table_exposed = table_everything[exposures > minexposure]

# get only the 64ms data files 
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    result = ba.download_swiftdata(table_exposed, auxil=True, match=['*brtms*', '*sao.fits*'], quiet=True)
result

# %%

# Barycenter
# Side-effect of the following is it updates the "/tmp/latest_swift_tle.gz" file
orbit=swiftbat.orbit()
tlefile = Path(swiftbat.tlefile)
assert tlefile.exists()


# %%
from swifttools.swift_too import Clock, ObsQuery, VisQuery
from datetime import datetime,timedelta
cc = Clock(swifttime='2022-01-01 00:00:00')
cc

 
table_stoo = swt.ObsQuery(begin=timerange[0],end=timerange[1])
print(f"{len(table_stoo)} entries starting with")
for e in table_stoo[:5]:
    print(e)
