# Each row in my CSV contains all observations of a single object stored as semicolon-separated strings.
# This script converts them into numeric arrays and calculates statistical and variability features for each object.

import pandas as pd
import numpy as np
from scipy.stats import kurtosis, skew, anderson
from astropy.stats import median_absolute_deviation
from scipy.optimize import curve_fit

# -----------------------
# Load CSV file containing my lightcurves
# I check the shape to ensure I loaded all objects correctly.
# -----------------------
lc_all = pd.read_csv("lightcurves_grouped_r_clean.xlsx")  # replace with your CSV path
print(f"Loaded {lc_all.shape[0]} objects")  # confirms I have all objects from the CSV

# -----------------------
# Convert semicolon-separated strings into numeric arrays
# This is critical because each cell in mag, hjd, and magerr columns contains all measurements as strings
# -----------------------
def str_to_float_array(s):
    try:
        # Convert string like "1.2;1.3;1.4" into [1.2, 1.3, 1.4]
        return np.array([float(x) for x in s.split(';') if x])
    except:
        return np.array([])  # if conversion fails, return empty array

for col in ['mag', 'hjd', 'magerr']:
    lc_all[col] = lc_all[col].apply(str_to_float_array)

# If magerr is missing or empty, I assume equal errors for all points
lc_all['magerr'] = lc_all['magerr'].apply(lambda x: x if len(x) > 0 else np.ones_like(lc_all['mag'].iloc[0]))

# -----------------------
# Define utility functions for feature extraction
# -----------------------
def weighted_mean(values, errors):
    # Standard weighted mean: gives more weight to more precise measurements
    return np.sum(values / errors**2) / np.sum(1 / errors**2)

def linear(x, a, b):
    # Simple linear function for trend fitting
    return a*x + b

def stetson_K(mags):
    # Stetson K is a robust measure of variability, useful for variable star classification
    N = len(mags)
    if N < 2:
        return np.nan
    mean_mag = np.mean(mags)
    delta = (mags - mean_mag) / np.std(mags, ddof=1)
    K = np.sum(np.abs(delta)) / (np.sqrt(np.sum(delta**2)) * np.sqrt(N))
    return K

# -----------------------
# Main loop: extract features for each object
# -----------------------
features_list = []

for idx, row in lc_all.iterrows():
    # Each row represents one object
    time = row['hjd']
    mag = row['mag']
    error = row['magerr']
    
    if len(mag) == 0 or len(time) == 0:
        continue  # skip empty objects

    # Initialize dictionary to store features for this object
    features = {}
    features['object'] = row['object']
    features['oid'] = row['oid']
    
    # Basic statistics
    features['mean'] = np.mean(mag)
    features['weighted_mean'] = weighted_mean(mag, error)
    features['standard_deviation'] = np.std(mag, ddof=1)
    features['median'] = np.median(mag)
    features['amplitude'] = np.max(mag) - np.min(mag)
    
    # Fraction of points beyond 1-sigma from mean
    features['beyond_1_std'] = np.sum(np.abs(mag - np.mean(mag)) > np.std(mag)) / len(mag)
    
    # Cumulative sum (Cusum) as a variability indicator
    features['cusum'] = np.max(np.cumsum(mag - np.mean(mag))) - np.min(np.cumsum(mag - np.mean(mag)))
    
    # Percentile-based measures
    features['inter_percentile_range_10'] = np.percentile(mag, 90) - np.percentile(mag, 10)
    
    # Shape descriptors
    features['kurtosis'] = kurtosis(mag)
    features['skew'] = skew(mag)
    
    # -----------------------
    # Linear trend fitting
    # I fit a line to the lightcurve to quantify overall trends
    # -----------------------
    try:
        popt, pcov = curve_fit(linear, time, mag)
        slope, intercept = popt
        slope_err = np.sqrt(np.diag(pcov))[0]
        residuals = mag - linear(time, *popt)
        reduced_chi2 = np.sum((residuals/error)**2) / (len(mag) - 2)

        features['linear_trend'] = slope * (time[-1] - time[0])
        features['linear_trend_sigma'] = slope_err
        features['linear_trend_noise'] = np.std(residuals)
        features['linear_fit_slope'] = slope
        features['linear_fit_slope_sigma'] = slope_err
        features['linear_fit_reduced_chi2'] = reduced_chi2
    except:
        # If fitting fails, store NaN
        features['linear_trend'] = np.nan
        features['linear_trend_sigma'] = np.nan
        features['linear_trend_noise'] = np.nan
        features['linear_fit_slope'] = np.nan
        features['linear_fit_slope_sigma'] = np.nan
        features['linear_fit_reduced_chi2'] = np.nan

    # Magnitude percentile ratios as variability measures
    sorted_mag = np.sort(mag)
    features['magnitude_percentage_ratio_40_5'] = np.percentile(sorted_mag, 40) / np.percentile(sorted_mag, 5)
    features['magnitude_percentage_ratio_20_10'] = np.percentile(sorted_mag, 20) / np.percentile(sorted_mag, 10)

    # Maximum slope between consecutive points
    if len(mag) > 1:
        slopes = np.diff(mag) / np.diff(time)
        features['maximum_slope'] = np.max(np.abs(slopes))
    else:
        features['maximum_slope'] = np.nan

    # Robust variability measures
    features['median_absolute_deviation'] = median_absolute_deviation(mag)
    med = np.median(mag)
    features['median_buffer_range_percentage_10'] = np.sum(np.abs(mag - med) < 0.1*np.std(mag)) / len(mag)
    
    # Relative variability
    features['percent_amplitude'] = (np.max(mag) - np.min(mag)) / np.mean(mag)
    features['mean_variance'] = np.var(mag) / np.mean(mag)

    # Anderson-Darling statistic to test normality of magnitudes
    try:
        ad_result = anderson(mag, dist='norm')
        features['anderson_darling_normal'] = ad_result.statistic
    except:
        features['anderson_darling_normal'] = np.nan

    # Chi-square variability measure
    features['chi2'] = np.sum((mag - np.mean(mag))**2 / error**2)

    # Stetson K index
    features['stetson_K'] = stetson_K(mag)

    # Append features of this object to list
    features_list.append(features)

# -----------------------
# Save all extracted features to CSV
# This CSV can now be used for machine learning or variability analysis
# -----------------------
features_df = pd.DataFrame(features_list)
features_df.to_csv("all_lightcurve_features.csv", index=False)
print("Feature extraction complete. Saved to all_lightcurve_features.csv")
