"""
Analysis of bioluminescence image by pixel. See README for details.
"""
from __future__ import division

# imports
import numpy  as np
import scipy as sp
import pandas as pd
from matplotlib import gridspec
import matplotlib.pyplot as plt

# local functions to import
from LocalImports import PlotOptions as plo
from LocalImports import Bioluminescence as blu
from LocalImports import DecayingSinusoid as dsin
from LocalImports import CellularRecording as cr

#inputs TIFF
INPUT_DIR = 'multi_frame_stacked_image.tif' # edit this
SAMPLE_INTERVAL = 1. # sample interval in h
BINNING = 5 # bin the data into larger pixes so it doesn't take a lifetime

#
# ALL PROCESSING
#

# Read TIFF
ts, locs, valid, names, traces = cr.read_tiff(INPUT_DIR, binning=BINNING, tstep=SAMPLE_INTERVAL)

# HP Detrend
detrended_times, detrended_data, trendlines = cr.hp_detrend(ts, traces)

# Eigendecomposition
denoised_times, denoised_data, eigenvalues = cr.eigensmooth(detrended_times, detrended_data, ev_threshold=0.05, dim=40)

# LS periodogram test for rhythm
lspers, pgram_data, circadian_peaks, lspeak_periods, rhythmic_or_not = cr.LS_pgram(detrended_times, detrended_data)

# Cosine fit
sinefit = cr.sinusoidal_fitting(denoised_times, denoised_data, rhythmic_or_not, forced_periods=lspeak_periods)

# Assemble Excel summary sheet
sum_df = pd.DataFrame()
sum_df['Names'] = names
sum_df['X'] = locs[:,0]
sum_df['Y'] = locs[:,1]
sum_df['Valid'] = valid
sum_df['Lomb-Scargle Peak'] = circadian_peaks
sum_df['Lomb-Scargle Rhythmic'] = rhythmic_or_not
sum_df['Cosine Phase Angle'] = sinefit['phi0all']
sum_df['Cosine Period'] = sinefit['periods']
sum_df['Cosine Amplitude'] = sinefit['amplitudes']
sum_df['Cosine Decay'] = sinefit['decays']
sum_df['Cosine R2'] = sinefit['r2s']

# Assemble sheet of raw data
raw_df = pd.DataFrame()
raw_df['Time (h)'] = ts
for idx, name in enumerate(names):
    raw_df[str(name)] = traces[:,idx]

# Assemble sheet of detrended data
detrend_df = pd.DataFrame()
detrend_df['Time (h)'] = detrended_times
for idx, name in enumerate(names):
    detrend_df[str(name)] = detrended_data[:,idx]

# Assemble sheet of denoised data
denoise_df = pd.DataFrame()
denoise_df['Time (h)'] = denoised_times
for idx, name in enumerate(names):
    denoise_df[str(name)] = denoised_data[:,idx]

# Assemble sheet of phases
phase_df = pd.DataFrame()
phase_df['Time (h)'] = sinefit['ts']
for idx, name in enumerate(names):
    phase_df[str(name)] = sinefit['phase_timeseries'][:,idx]

# Assemble sheet of cosines
cos_df = pd.DataFrame()
cos_df['Time (h)'] = sinefit['ts']
for idx, name in enumerate(names):
    cos_df[str(name)] = sinefit['cosines'][:,idx]

# Save it all
writer = pd.ExcelWriter(INPUT_DIR.split('.')[0]+'.xlsx', engine='xlsxwriter')
sum_df.to_excel(writer, sheet_name="Analysis Results", index=False)
raw_df.to_excel(writer, sheet_name="Unprocessed Data", index=False)
detrend_df.to_excel(writer, sheet_name="Detrended Data", index=False)
denoise_df.to_excel(writer, sheet_name="Denoised Data", index=False)
phase_df.to_excel(writer, sheet_name="Phase Data", index=False)
cos_df.to_excel(writer, sheet_name="Cosine Data", index=False)
writer.save()

# Plot LSPgram figure
cr.plot_lsp(sum_df, 
            title=INPUT_DIR.split('.')[0].split('/')[-1])
plt.savefig(INPUT_DIR.split('.')[0]+'_LSPgram.png')
plt.close()

# Plot phase angle figure
time_of_calc = 48
frame_of_calc = int(time_of_calc/SAMPLE_INTERVAL)
cr.plot_phaseangle(sum_df, frame_of_calc,
            title=INPUT_DIR.split('.')[0].split('/')[-1])
plt.savefig(INPUT_DIR.split('.')[0]+'_Phase_'+time_of_calc'+h.png')
plt.close()