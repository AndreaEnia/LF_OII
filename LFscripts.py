import os, sys, subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import numpy as np
from collections import Counter
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from astropy.cosmology import z_at_value
from astropy.stats import poisson_conf_interval
from scipy import interpolate
from scipy.optimize import curve_fit
from tqdm import tqdm 
from scipy.stats import norm

def LogError(value, error):
    rel_err = (np.array(error)/np.array(value))/np.log(10)
    return rel_err, rel_err, rel_err

def partition_redshift(values, bin_extremes):
    bins = [np.array([bin_min, bin_max]) for bin_min, bin_max in zip(bin_extremes, bin_extremes[1:])]
    bins_digitized = np.digitize(values, bin_extremes, right = True)
    bin_centers = [values[bins_digitized == i].mean() for i in range(1, len(bins)+1)]
    bin_medians = np.round([np.median(values[np.where(np.logical_and(values > bin_.min(), values <= bin_.max()))[0]]) for bin_ in bins], 2)
    return [bins, bins_digitized, bin_centers, bin_medians]

def partition_luminosities(logL, number_of_bins, redshift_list, redshift_bins, bin_kind):
    bin_labels, luminosity_bins = [], []
    for idx, z_bin in enumerate(redshift_bins):
        ok = np.where(np.logical_and(redshift_list > z_bin.min(), redshift_list <= z_bin.max()))[0]
        if bin_kind == 'sspb': # Binning with same number of sources per bi
            _, bin_extremes = pd.qcut(logL[ok], number_of_bins, labels = False, precision = 5, retbins = True)
        elif bin_kind == 'eqbin': # Binning with equispatiate bins.n
            _, bin_extremes = pd.cut(logL[ok], number_of_bins, labels = False, precision = 5, retbins = True)
        bin_labels.append(_)
        luminosity_bins.append([np.array([bin_min, bin_max]) for bin_min, bin_max in zip(bin_extremes, bin_extremes[1:])])   
    return pd.concat(bin_labels, axis = 0, join = 'inner').sort_index(), luminosity_bins