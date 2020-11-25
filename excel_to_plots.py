"""
This file enables plotting using the excel file saved during the analysis.
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

#input excel
EXCEL_FILE = 'output_file.xlsx'
SAMPLE_INTERVAL = 1. # in hours

# load the existing file
reader = pd.ExcelFile(EXCEL_FILE)
sum_df = reader.parse(sheet_name="Analysis Results")

# Plot LSPgram figure
cr.plot_lsp(sum_df, 
            title=EXCEL_FILE.split('.')[0].split('/')[-1])
plt.savefig(EXCEL_FILE.split('.')[0]+'_LSPgram.png')
plt.close()

# Plot phase angle figure
time_of_calc = 48
frame_of_calc = int(time_of_calc/SAMPLE_INTERVAL)
cr.plot_phaseangle(sum_df, frame_of_calc,
            title=INPUT_DIR.split('.')[0].split('/')[-1])
plt.savefig(INPUT_DIR.split('.')[0]+'_Phase_'+time_of_calc'+h.png')
plt.close()