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
    splitlocs = np.argwhere(np.diff(obsdata['time']) > 1.5*timebin).ravel() + 1
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
assert (0.9 * timebin) < np.median(np.diff(lcsegments[0]['time'])) < (1.5 * timebin)

# Use the segments to populate an array
t0 = lcsegments[0]['time'][0]
tmax = lcsegments[-1]['time'][-1]
ntimes = sp.fft.next_fast_len(int((tmax - t0)/timebin  + 10))    # 10 bins of slop, then round up to an FFT-friendly length
lcfull = np.zeros(ntimes)

for segment in lcsegments:
    n = len(segment)
    i0 = int((segment['time'][0] - t0)/timebin)
    lcfull[i0:i0+n] = segment['norm']
    
for datasegment in np.split(obsdata, splitlocs):
    # time of spacecraft, not always accurate because of clock error
    starttime = swiftbat.met2datetime(datasegment['TIME'][0])
    duration = datasegment['TIME'].ptp()
    print(f"{starttime:%Y-%m-%dT%H:%M:%S} + {duration:5.10f} seconds to the end of the block")
    if duration > 1300:
        longdatasegment = datasegment      


# get the highest number of datapoints for FFT
def prev_fast_FFT_len(n):
    ntry = abs(n)
    nfft = sp.fft.next_fast_len(ntry)
    while nfft > n and ntry > 1:
        ntry = int(ntry * 0.99) - 1
        nfft = sp.fft.next_fast_len(ntry)
    return nfft

# refind the amount of datapoints afte cutting off the first minute (BAT is still adjusting the first min)
n = prev_fast_FFT_len(len(datasegment[0]) - int(60 / timebin))

# Trim to the last n valuesfits.getdata(filename,  header=True)
datasegment = datasegment[-n:]
duration = datasegment['TIME'].ptp()
print(f"{duration:.3f} seconds after trimming")

# Do the Fourier transform, and get the corresponding frequencies and powers
fnorm = sp.fft.rfft(lcfull, norm='forward')
freqs = sp.fft.rfftfreq(ntimes, timebin)
fpower = np.abs(fnorm)**2
# Ignore periods below 10 minutes when looking for peak
zerof_ignore = int(ntimes/(600))
imax = zerof_ignore + np.argmax(fpower[zerof_ignore:])
freqmax = freqs[imax]
powermax = fpower[imax]
# Don't plot all points because that takes a long time
grasslevel = np.median(fpower) * 10
wplot = np.argwhere(fpower > grasslevel).ravel()

# Ignore periods below 20 minutes when looking for peak
zerof_ignore_2 = int(ntimes/(60))
imax_2 = zerof_ignore_2 + np.argmax(fpower[zerof_ignore_2:])
freqmax_2 = freqs[imax_2]
powermax_2 = fpower[imax_2]
# Don't plot all points because that takes a long time
grasslevel_2 = np.median(fpower) * 10
wplot_2 = np.argwhere(fpower > grasslevel_2).ravel()

print(f"{freqmax = } Hz, {powermax = }, {grasslevel = }")


fig, axes = plt.subplots(nrows=3, ncols=1)
axes[0].plot(freqs[wplot], fpower[wplot])
axes[0].set(yscale='log', ylim=[grasslevel, powermax*1.3], xlim=[0,1], ylabel="Power (logscale arbitrary units)")
for harmonic in range(1,5):
    irange = ((imax + np.asarray([-200,201])) * harmonic).astype(int)
    axes[1].plot(freqs[irange[0]:irange[1]]/harmonic, fpower[irange[0]:irange[1]], label=f"n={harmonic}")
axes[1].legend() #
axes[1].set(title=f"harmonic of {freqmax:f} Hz", xlabel="Frequency and harmonic-adjusted frequency (Hz)")
fig.tight_layout()

for harmonic in range(1,5):
    irange_2 = ((imax_2 + np.asarray([-200,201])) * harmonic).astype(int)
    axes[2].plot(freqs[irange_2[0]:irange_2[1]]/harmonic, fpower[irange_2[0]:irange_2[1]], label=f"n={harmonic}")
axes[2].legend() #
axes[2].set(title=f"harmonic of {freqmax_2:f} Hz", xlabel="Frequency and harmonic-adjusted frequency (Hz)")
fig.tight_layout()
plt.show()

plt.figure(1)


tzero = swiftbat.string2met('2023-08-27T00:00:00')
fig, axes = plt.subplots(nrows = 11, ncols = 1, sharex = True)

# plotting just the 1300s segment
rate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/timebin
print(len(datasegment))

segpieces = 4
# Break the segment into 4 pieces
pointsperplot = 2 * int(1/(0.064 * 0.8))
# For each segment, plot 2 cycles of data

# 2 cycles of data for all 11 data segments
for datasegment, ax in zip(np.split(obsdata, splitlocs), axes):
    n = prev_fast_FFT_len(len(datasegment) - int(60 / timebin))
    datasegment = datasegment[-n:]
    duration = datasegment['TIME'].ptp()
    segrate = np.sum(datasegment['COUNTS'][:,0:2], axis=-1)/timebin
    segcycle, segphase = np.divmod((datasegment['TIME']-tzero) * 0.8, 1)
    ax.plot(segphase[0:3675], segrate[0:3675], ".")
    
    
fig.tight_layout()
plt.show()
plt.figure(2)
