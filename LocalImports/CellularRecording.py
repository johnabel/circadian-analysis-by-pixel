'''
File containing functions useful for data analysis.
'''


from __future__ import division

# imports
import os
import numpy  as np
import scipy as sp
import pandas as pd
import spectrum
from matplotlib import gridspec
import matplotlib.pyplot as plt
from PIL import Image

from LocalImports import PlotOptions as plo
from LocalImports import Bioluminescence as blu
from LocalImports import DecayingSinusoid as dsin
from LocalImports.tol_colors import tol_cmap
from LocalImports.twilight_cmap import cmaps

# import function
def read_tiff(path, n_frames=1000, tstep=1., minthresh=10, binning=1):
    """Reads a multi-frame tiff image into a set of trajectories.
    
    Parameters
    ----------
    path : str
        Path to .tif or .tiff file.
    n_images : int, optional
        Maximum number of frames in tiff file we are reading, defaults to 1000
    tstep : float, optional
        Timestep for a frame, by default 0.25
    minthresh : int, optional
        Threshold for determining if a cell exists in a location. If all timepoints are below minthresh, valid is set to false. By default 40
    binning : int, optional
        How many adjectent pixels should be binned into one. There is in some ways good reason to not bin, notably, that we don't want to blend neuronal signals. On the other hand, a binning of 5 speeds computation 25-fold. Defaults to 1, no binning.
    
    Returns
    -------
    [type]
        [description]
    """
    img = Image.open(path)
    images = []
    for i in range(n_frames):
        try:
            img.seek(i)
            slice_ = np.zeros((img.height, img.width))
            for j in range(slice_.shape[0]):
                for k in range(slice_.shape[1]):
                    slice_[j,k] = img.getpixel((k, j))

            images.append(slice_)

        except EOFError:
                # Not enough frames in img
                break
    
    bin_images = []
    new_height = images[0].shape[0] //binning
    new_width = images[0].shape[1] //binning
    for image in images:
        # truncate so it can be reshaped
        image_trunc = image[:new_height*binning, :new_width*binning]
        new_image = rebin(image_trunc, [new_height, new_width])
        bin_images.append(new_image)
    
    images = np.array(bin_images)
    traces = []
    iloc = []
    for j in range(images[0].shape[0]):
        for k in range(images[0].shape[1]):
            traces.append(images[:,j,k])
            iloc.append([j,k])

    # get times, traces, locations
    ts = np.arange(len(traces[0]))*tstep
    traces = np.array(traces)
    locs = np.array(iloc)
    names = np.arange(len(traces))

    valid = [any(trace > minthresh) for trace in traces]

    return ts, locs, valid, names, traces.T

def alignment(original, denoised, d=40, dstart=0):
    """
    The eigensmoothing as written truncates some front-back data, as the
    input data input data is extended. This function realigns the data

    The sign is sometimes flipped on the reconstructed signal.
    This function figures out +/- alignment as well.

    dstart (default=0) tells where to start checking for alignment. Some
    of the data has an artifact, so later datapoints are more reliable
    for finding an alignment (otherwise the artifact will dominate).
    """
    original = original[dstart:]
    denoised = denoised[dstart:]
    errs = np.zeros(d-1)
    for idx in range(d-1):
        errs[idx] = np.linalg.norm(original-denoised[idx:-(d-idx-1)])
    errs_neg = np.zeros(d-1)
    for idx in range(d-1):
        errs_neg[idx] = np.linalg.norm(original+denoised[idx:-(d-idx-1)])
    pos_min = np.min(errs)
    neg_min = np.min(errs_neg)
    if neg_min < pos_min:
        return np.argmin(errs_neg), -1
    else:
        return np.argmin(errs), 1

def truncate_and_interpolate(times, data, locations, truncate_t=12,
                                outlier_std=5):
    """
    Removes the first "truncate_t" h artifact and interpolates
    any missing values and rejects outliers.
    """
    # truncate the data
    start = np.argmax(times-times[0]>=truncate_t)
    times = times[start:]
    data = data[start:]
    locations = locations[start:]

    # first, find where nans are in the data
    h1, _ = blu.nan_helper(data)

    # create output data
    outdata = np.nan*np.ones(data.shape)

    # interpolate the intermediate points
    for i in range(len(data[0,:])):
        try:
            init_goodval = np.min(np.where(~h1[:,i])[0])
            end_goodval  = np.max(np.where(~h1[:,i])[0])
            # only interpolate where we have good values
            y = np.copy(data[init_goodval:end_goodval, i])

            # find outliers - ignore where there's nans
            #  and ignore the warning here
            with np.errstate(invalid='ignore'):
                outliers = np.abs(y - np.nanmean(y)) > outlier_std * np.nanstd(y)
            outlier_idx = np.where(outliers)[0]
            # if any exist, replace them
            if len(outlier_idx)>0:
                # put in nans
                for oi in outlier_idx:
                    y[oi] = np.nan

            # find nans
            nans, x= blu.nan_helper(y)
            # interpolate
            y[nans]= np.interp(x(nans), x(~nans), y[~nans])
            # replace old values with new
            outdata[init_goodval:end_goodval, i] = y
        except ValueError:
            # if there are no nans
            pass

    # should we trim the data here as well? we can just leave it
    return times, outdata, locations

def hp_detrend(times, interpolated_data):
    """
    Does a Hodrick-Prescott detrending always using an estimated period
    of 24h.
    """
    timer = plo.laptimer()
    print "Detrending data... time: ",
    detrended_data = np.zeros(interpolated_data.shape)
    trendlines = np.zeros(interpolated_data.shape)
    for idx,idata in enumerate(interpolated_data.T):
        # copy data for editing
        cell = np.copy(idata)
        model = blu.Bioluminescence(times, cell, period_guess=24)
        # all we care about is the HP filter, no need to fit everything
        model.filter()
        model.detrend()
        baseline = model.yvals['mean']

        # find where recording starts and there are no nans
        valid = ~np.isnan(cell)
        # subtrect baseline from recording
        cell[valid] = cell[valid] - baseline
        # put the detrended cell in the detrended_data matrix,
        # and same for baseline
        detrended_data[:,idx] = cell
        trendlines[valid,idx] = baseline
    print str(np.round(timer(),1))+"s"
    return times, detrended_data, trendlines

def butterworth_lowpass(times, data, cutoff_period=4):
    """
    Does a Butterworth low pass filtering of the data.
    """
    timer = plo.laptimer()
    print "Butterworth filter... time: ",
    denoised_data = np.zeros(data.shape)
    for idx in range(len(data.T)):
        # copy data for editing
        idata = np.copy(data[:,idx])
        valid = ~np.isnan(idata)
        _, filtdata = blu.lowpass_filter(times[valid], idata[valid],
                                        cutoff_period=cutoff_period)

        # subtrect baseline from recording
        denoised_data[valid,idx] = filtdata

    print str(np.round(timer(),1))+"s"
    return times, denoised_data

def LS_pgram(times, ls_data, circ_low=18, circ_high=30, alpha=0.05):
    """Calculates a LS periodogram for each data sequence,
    and returns the p-values for each peak. If the largest significant
    peak is in the circadian range as specified by the args, it is
    rhythmic."""
    timer = plo.laptimer()
    print "Lomb-Scargle Periodogram... time: ",
    rhythmic_or_not = np.zeros(len(ls_data.T))
    pgram_data = np.zeros((300,len(ls_data.T)))
    circadian_peaks = np.zeros(len(ls_data.T))
    circadian_peak_periods = np.zeros(len(ls_data.T))
    np.seterr(divide='ignore', invalid='ignore')
    # pgram
    for data_idx, d1 in enumerate(ls_data.T):
        # remove nans
        t1 = np.copy(times[~np.isnan(d1)])
        d1 = np.copy(d1[~np.isnan(d1)])
        pers, pgram, sig = blu.periodogram(t1, d1, period_low=1,
                        period_high=60, res=300)
        peak = np.argmax(pgram)

        if (pers[peak] >= circ_low and pers[peak] <= circ_high and
                                                    sig[peak] <=alpha):
            rhythmic_or_not[data_idx] = 1
            #
            circadian_peaks[data_idx]= pgram[peak]
            circadian_peak_periods[data_idx]= pers[peak]
        else:
            minpeak = np.argmax(pers>=circ_low)
            maxpeak = np.argmin(pers<circ_high)
            circadian_peaks[data_idx] = np.max(pgram[minpeak:maxpeak])
            circadian_peak_periods[data_idx] =\
                pers[minpeak:maxpeak][np.argmax(pgram[minpeak:maxpeak])]

        # return either normed or un-normed data
        pgram_data[:,data_idx] = pgram
    print str(np.round(timer(),1))+"s"
    return pers, pgram_data, circadian_peaks, circadian_peak_periods, rhythmic_or_not

def eigensmooth(times, detrended_data, ev_threshold=0.05, dim=40, min_ev=2):
    """
    Uses an eigendecomposition to keep only elements with >threshold of the
    data. Then it returns the denoised data.

    Notes: This should take the covariance matrix of the data using fwd-backward method. Then it
    eigendecomposes it, then it finds the biggest (2+) eigenvalues and returns only the
    components associated with them.
    For an intuitive explanation of how this works:
    http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    """
    timer = plo.laptimer()
    print "Eigendecomposition... time: ",

    #set up results - by default set to NaNs
    denoised_data = np.nan*np.ones(detrended_data.shape)
    eigenvalues_list = []

    # denoise each piece of data
    for data_idx, d0 in enumerate(detrended_data.T):
        # remove nans from times, d1. keep nans in d0
        t1 = times[~np.isnan(d0)]
        d1 = np.copy(d0[~np.isnan(d0)])

        # using spectrum to get the covariance matrix
        X = spectrum.linalg.corrmtx(d1, dim-1, method='autocorrelation')
        # the embedding matrix
        X = (1/np.sqrt(len(X.T)))*np.array(X)
        XT = np.transpose(X)
        C = XT.dot(X)

        # now eigendecompose
        evals, Q = np.linalg.eig(C)

        # find evals that matter, use a minimum of 2
        eval_goodness = np.max([2,
                        np.sum(evals/(1E-12+np.sum(evals)) >= ev_threshold)])
        QT = np.transpose(Q)

        # and return the reconstruction
        P = QT.dot(XT)
        denoised = np.sum(P[:eval_goodness],0)

        # find alignment - for some reason the signal can be flipped.
        # this catches it
        # truncate the first 24h during alignment
        align, atype = alignment(d1, denoised, d=dim, dstart=96)

        # fix alignment if leading nans
        nanshift=0
        if np.isnan(d0[0]):
            nanshift=np.argmin(np.isnan(d0))

        # get the correctly-shaped denoised data
        denoised = denoised[align:align+len(d1)]*atype

        #align denoised data in matrix
        denoised_data[nanshift:len(denoised)+nanshift, data_idx] = denoised
        eigenvalues_list.append(evals/np.sum(evals+1E-12))

    print str(np.round(timer(),1))+"s"
    return times, denoised_data, eigenvalues_list

def sinusoidal_fitting(times, data, rhythmic_or_not, fit_times=None,
                       forced_periods=None):
    """
    Takes detrended and denoised data and times, and fits a damped sinusoid to
    the data. Returns: times, sine_data, phase_data, periods, amplitudes,
    decay parameters, and pseudo rsq values.

    If forced_periods are supplied, the sine fit will be within 1h of the forced period.
    """
    timer = plo.laptimer()
    tot_time = 0
    print "\nSinusoidal fitting.\n----------"
    ncells = len(data.T)

    if fit_times is None:
        fit_times = np.copy(times)

    # these will be the outputs
    sine_data = np.nan*np.ones((len(fit_times), len(data[0])))
    phase_data = np.nan*np.ones((len(fit_times), len(data[0])))
    phases = np.nan*np.ones(len(data[0]))
    meaningful_phases = np.nan*np.ones(len(data[0]))
    periods = np.nan*np.ones(len(data[0]))
    amplitudes = np.nan*np.ones(len(data[0]))
    decays = np.nan*np.ones(len(data[0]))
    pseudo_r2s = np.nan*np.ones(len(data[0]))

    for idx,idata in enumerate(data.T):
        if np.std(idata) > 1E-12:
            # try a sinusoidal fit but only if data is not perfectly flat
            if idx%100==0:
                print str(np.round(100*idx/ncells,2))+" pct complete... time: ",
                tot_time+=timer()
                print str(np.round(tot_time, 1))+"s"
            # copy data for editing
            cell = np.copy(idata)
            model = dsin.DecayingSinusoid(times, cell)
            # fit the data with the polynomial+sinusoid
            # note we are not using model averaging, we are
            # just taking the single best model by AICc

            # if no estimate is given just do the normal fitting
            # if an estimate is given, bound the period within two
            model.run()
            params = model.best_model.result.params

            if forced_periods is not None:
                # only use force if necessary
                if np.abs(params['period'].value-forced_periods[idx]) >1:
                    # force is necessary
                    model._estimate_parameters()
                    model._fit_models(period_force = forced_periods[idx])
                    model._calculate_averaged_parameters()
                    params = model.best_model.result.params
        
            phases[idx] = (fit_times[0]*2*np.pi/params['period']+params['phase'])%(2*np.pi)

            if rhythmic_or_not[idx]==1:
                phase_data[:,idx] = (fit_times*2*np.pi/params['period']+params['phase'])%(2*np.pi)
                sine_data[:,idx] = dsin.sinusoid_component(params, fit_times)
                # summary stats
                periods[idx] = model.best_model.result.params['period']
                amplitudes[idx] = model.best_model.result.params['amplitude']
                decays[idx] = model.best_model.result.params['decay']
                pseudo_r2s[idx] = model.best_model._calc_r2()
                meaningful_phases[idx] = (fit_times[0]*2*np.pi/params['period']+params['phase'])%(2*np.pi)

    sinefit = {'ts': fit_times, 'cosines':sine_data, 'phase_timeseries': phase_data, 'phi0all':phases, 'periods':periods, 'amplitudes':amplitudes, 'decays':decays, 'r2s':pseudo_r2s, 'phi0r':meaningful_phases}
    return sinefit

def rebin(arr, new_shape):
    """ Utility for rebinning a numpy arra """
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
                new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def plot_scn_heatmap(ax, X, Y, values, cmap, **kwargs):
    """Plots a heatmap of some SCN values.
    
    Parameters
    ----------
    X : [type]
        [description]
    Y : [type]
        [description]
    values : [type]
        [description]
    """
    x = np.sort(np.unique(X))
    y = np.sort(np.unique(Y))
    v = values.reshape(len(x), len(y))
    cbar = ax.pcolormesh(y, x, v, cmap=cmap, **kwargs)
    ax.invert_yaxis()
    return cbar

def plot_lsp(sum_df, title=''):
    """ Utility to plot LSPgram value:"""
    fig = plt.figure()
    ax = plt.subplot()
    X = sum_df['X'].values
    Y = sum_df['Y'].values
    values = sum_df['Lomb-Scargle Peak'].values
    valid = sum_df['Valid'].values

    # make invalid regions nans
    palette = tol_cmap(colormap='sunset')
    values[~valid] = np.nan
    cbar = plot_scn_heatmap(ax, X, Y, values, palette, vmin=0, vmax=1)
    plt.colorbar(cbar, ax=ax, label='Lomb-Scargle Value')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(title)

def plot_phaseangle(sum_df, time, title=''):
    """ Utility to plot phase angle at time (h) past start."""
    fig = plt.figure()
    ax = plt.subplot()
    X = sum_df['X'].values
    Y = sum_df['Y'].values
    values = np.copy(sum_df['Cosine Phase Angle'].values)
    periods = np.copy(sum_df['Cosine Period'].values)
    phases_at_times = (values+time*2*np.pi/periods)%(2*np.pi)
    phases_h = phases_at_times*periods/(2*np.pi)
    valid = np.copy(sum_df['Lomb-Scargle Rhythmic'].values)

    # make invalid regions nans
    values[~valid] = np.nan
    cbar = plot_scn_heatmap(ax, X, Y, phases_h, cmaps['twilight'], vmin=0, vmax=2*np.pi)
    plt.colorbar(cbar, ax=ax, label='Phase Angle, '+str(time)+'h')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(title)
