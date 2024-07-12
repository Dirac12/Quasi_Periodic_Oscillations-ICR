import batanalysis as ba
import swiftbat
import swifttools
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from pathlib import Path
from astropy.io import fits
import datetime
import os

sourcename = "Swift J1727.8-1613"
source = swiftbat.source(sourcename)    # Can look up positions from Simbad, and can calculate exposure for a given pointing
timerange = [swiftbat.string2datetime(t) for t in ("MJD60183", "MJD60186")]
min_exposure_area = 1000     # cm^2 after cos adjust

table_stoo = swifttools.swift_too.ObsQuery(begin=timerange[0],end=timerange[1])
print(f"{len(table_stoo)} entries starting with")
for e in table_stoo[:5]:
    print(e)

allobsids = []
obsids = []
for entry in table_stoo:
    if (entry.obsid not in allobsids):
        allobsids.append(entry.obsid)
        if source.exposure(ra=entry.ra, dec=entry.dec, roll=entry.roll)[0] > min_exposure_area:
            obsids.append(entry.obsid)
print(f"Out of {len(allobsids)} distinct OBSIDs, {len(obsids)} in FOV  to download.")
download_multi = ba.download_swiftdata(table_stoo, match=['*brtms*'], quiet=True)

lcsegments = []
# rate is rate over first 2 energy bins 15-50 keV
slice_ebins=slice(0,2)
timebin = 0.064
# Norm is mean-subtracted, stddev-scaled within a pointing
segdtype = np.dtype([('time', np.float64),('rate', np.int16),('norm', np.float32)])

for obsid, entry in download_multi.items():
    if not entry['success']:
        continue
    # if len(entry['data']) != 1:
    #     print(f"OBSID {obsid} has {len(entry['data'])} files")
    #     print(entry['data'])
    #     continue
    datafile = entry['data'][0].localpath
    obsdata = fits.getdata(datafile)
    # Split the data into arrays with no more than a second's gap
    splitlocs = np.argwhere(np.diff(obsdata['time']) > 1.1).ravel() + 1
    for segmentdata in np.split(obsdata, splitlocs):
        segment = np.empty(len(segmentdata), dtype=segdtype)
        segment['time'] = segmentdata['time']
        rate = np.sum(segmentdata['COUNTS'][:,slice_ebins], axis=1)
        segment['rate'] = rate
        segment['norm'] = (rate - np.mean(rate))/(0.001 + np.std(rate))
        lcsegments.append(segment)

# Sort by segment start time
lcsegments = sorted(lcsegments, key = lambda x:x['time'][0])

# Make sure the timebin is right
assert (0.9 * timebin) < np.median(np.diff(lcsegments[0]['time'])) < (1.1 * timebin)

# Use the segments to populate an array
t0 = lcsegments[0]['time'][0]
tmax = lcsegments[-1]['time'][-1]
ntimes = sp.fft.next_fast_len(int((tmax - t0)/timebin  + 10))    # 10 bins of slop, then round up to an FFT-friendly length
lcfull = np.zeros(ntimes)

for segment in lcsegments:
    n = len(segment)
    i0 = int((segment['time'][0] - t0)/timebin)
    lcfull[i0:i0+n] = segment['norm']

# Do the Fourier transform, and get the corresponding frequencies and powers
fnorm = sp.fft.rfft(lcfull, norm='forward')
freqs = sp.fft.rfftfreq(ntimes, timebin)
fpower = np.abs(fnorm)**2
# Ignore periods below 10 minutes when looking for peak
zerof_ignore = int(ntimes/(600/0.064))
imax = zerof_ignore + np.argmax(fpower[zerof_ignore:])
freqmax = freqs[imax]
powermax = fpower[imax]
# Don't plot all points because that takes a long time
grasslevel = np.median(fpower) * 10
wplot = np.argwhere(fpower > grasslevel).ravel()

print(f"{freqmax = } Hz, {powermax = }, {grasslevel = }")

fig,axes = plt.subplots(nrows=2, ncols=1)
axes[0].plot(freqs[wplot], fpower[wplot])
axes[0].set(yscale='log', ylim=[grasslevel, powermax*1.3], xlim=[0,1], ylabel="Power (logscale arbitrary units)")
for harmonic in range(1,5):
    irange = ((imax + np.asarray([-200,201])) * harmonic).astype(int)
    axes[1].plot(freqs[irange[0]:irange[1]]/harmonic, fpower[irange[0]:irange[1]], label=f"n={harmonic}")
axes[1].legend() #
axes[1].set(title=f"harmonic of {freqmax:.5f} Hz", xlabel="Frequency and harmonic-adjusted frequency (Hz)")
fig.tight_layout()

