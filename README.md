# Code for per2py data analysis
This repository contains code for processing Per2Luc data on a by-pixe3l basis. It uses the scientific Python stack to identify and process circadian oscillatory data in a reproducible manner. This analysis is nearly identical to that in Shan, Abel, ... Takahashi, Neuron 2020 but applied to each pixel in a recording.

For any errors, bugs, or questions, please use the Issues board on this Github repository. Generally, this code should be adapted for the individual system employed for experimentation. If adaptation is desired, please contact John Abel at jhabel01(at)gmail(dot)com.

# Table of Contents
* [Installation](#Installation)
* [Usage](#Usage)
* [Analysis of cellular bioluminescence recordings](#Analysis-of-cellular-bioluminescence-recordings)
* [Analysis of whole-body bioluminescence recordings](#Analysis-of-whole-body-bioluminescence-recordings)
* [License](#License)
* [Authors](#Authors)
* [Funding](#Funding)

# Installation

## Using Anaconda
To set up an Anaconda environement for these scripts, run in a terminal:
```
conda env create -f ConfigFiles/per2py.yml
```
then 
```
conda activate per2py
```

## Without Anaconda
If not using Anaconda, it is sufficient to ensure that the following packages are installed in a Python 2.7 environment:

| Package | Version | Website |
|---------|---------|---------------------------|
|`jupyter`|1.0.0|https://jupyter.org/|
|`numpy`|>1.6.0|http://scipy.org |
|`scipy`|>0.13.0|http://scipy.org |
|`matplotlib`|>1.3.0|http://scipy.org |
|`future`|>0.16.0|http://python-future.org/quickstart.html |
|`spectrum`|0.7.5|https://pyspectrum.readthedocs.io/en/latest/|
|`pandas`|0.23.4|https://pandas.pydata.org/|
|`pywavelets`|1.0.3|https://pywavelets.readthedocs.io/en/latest/|
|`lmfit`|0.9.11|https://lmfit.github.io/lmfit-py/|
|`prettyplotlib`|0.1.7|https://github.com/olgabot/prettyplotlib|
|`Pillow`|6.2.1|https://pillow.readthedocs.io/|

To check for pip:
```
python -m ensurepip
```

If pip is installed, the associated packages may be installed with:
```
pip install -U [packagename]
```

# Usage
Two files are provided for analysis. The first `pixel_analysis.py` performs by-pixel analysis of a series of bioluminescence images. Results are saved in the same directory as the input image: a `.tif` image where each layer is one frame in a sequence of bioluminescence images. The second `excel_to_plots.py` takes the outputs of the first and generates two plots (although others could be made).

## Interpreting Results
The data produced during this analysis includes detrended signal, detrended and denoised signal, Lomb-Scargle periodogram results, a sinusoidal model fit to the data, the phases of the sinusoidal approximation, and parameters describing the sinusoid for each trajectory provided. These results are exported as an `.xlsx` file.

<br><br>

## Step-by-step details for analysis of single-cell data
1. **Import data.**
    Details on the data import are provided in the instructions section.
<br><br>

2. **Detrend via Hodrick-Prescott filter.**
    A Hodrick-Prescott filter is applied to detrend the data. We (JH Abel and Y Shan) found that this performs better than polynomial detrending in that it does not overfit the trend and begin to fit the oscillatory components. Parameter selection for the filter is performed as in Ravn, Uhlig 2004, and implemented as in St. John and Doyle PLOS Comp Biol 2015.<br><br>
    
3. **Eigendecompose and reconstruct the signal to denoise data.**
    An eigendecomposition of the autocorrelation matrix is performed. Any eigenvectors with eigenvalue >5% of the total sum of all eigenvalues is kept, and the signal is reconstructed from the corresponding eigenvectors. The process is explained in detail here: http://www.fil.ion.ucl.ac.uk/~wpenny/course/subspace.pdf<br><br>

4. **Apply a Lomb-Scargle periodogram to determine which timeseries.**
    The Lomb-Scargle periodogram is used to identify statistically significant circadian oscillation. Here, we have defined this as having a dominant peak with P<0.05 between periods of 18 and 30 h. The method for periodogram calculation and corresponding P-values is from WH Press, Numerical Recipes, 1994, p576. The final periodogram is normalized using the standard normalization in astropy or scipy. The detrended data are used for this analysis.<br><br>
    
5. **Fit a damped cosine to the data to yield sine fit, phase, amplitude, damping rate, goodness of fit.**
    All rhythmic cells are now fit by a damped sinusoid plus a polynomial (to fit any non-oscillatory trend) using nonlinear least squares. Only the sinusoid data is returned. Portions of this method are used from St. John, Taylor, Abel, Doyle III, Biophys J 2014.<br><br>
    
6. **Export data from the analyses.** 
    The data produced at most steps of this process is then saved in the output folder, as delineated in the [Usage](#Usage). <br><br>

7.  **Generate plots for error-checking.**
    Summaries of the Lomb-Scargle rhythmicity and the phase angles of oscillation are saved in the same folder as the TIF files.

# License

This code is licensed under the GNU General Public Licesnce version 3. Please see LICENSE for more information.

# Authors
This code was written by John Abel and Yongli Shan, with portions of code resused from other publications where noted.

# Funding
Research leading to this code was funded by NIH/NIA F32 AG064886 and NIH/NHLBI T32 HL007901 (JHA).