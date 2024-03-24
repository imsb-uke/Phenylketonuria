"""
    Gaussian Modelling and feature extraction
    Version 1.0.0 (needs Utilities.py)
    Authors: 
        - Behnam Yousefi (behnm.yousefi@zmnh.uni-hamburg.de; yousefi.bme@gmail.com)
        - Robin Khatri (robin.khatri@zmnh.uni-hamburg.de)
    Scientific Contributor: 
        - Polina Gundorova (p.gundorova@uke.de)
    UKE, Hamburg, Germany - March 2024

    To to plot the landscapes, the current version uses some codes from a Python script written by
    Italo Balestra (italo.balestra@olt-dss.com) and Andrea Turcati (andrea.turcati@olt-dss.com)
    OmegaLambdaTec GmbH - 2023
"""

# Load the packages
import os
from tqdm import tqdm
import argparse

import pandas as pd
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.optimize import curve_fit

from Utilities import *

# Main
if __name__ == "__main__":

    # Read inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True, help='Input directory; either xlsm file name or a directory')
    parser.add_argument('--save_image_dir', default="gm_output/images/", help='Directory to save images')
    parser.add_argument('--save_feature_dir', default="gm_output/features/", help='Directory to save feature csv file')
    parser.add_argument('--save_plot', action='store_true', help='If specified, the plots will be saved; this may significantly increase the running time')
    parser.add_argument('--tag', default="", help='Name tag for the feature csv file')
    parser.add_argument('--qc_thr_rmse', default=[0.2, 0.25], help='QC threshold for RMSE')
    parser.add_argument('--qc_thr_n_peaks', default=[5, 8], help='QC threshold for N Peacks')
    parser.add_argument('--qc_thr_variation', default=[0.1, 0.25], help='QC threshold for Variations')

    # Set parameters
    args = parser.parse_args()
    dir = args.input_dir
    save_image_dir = args.save_image_dir
    save_feature_dir = args.save_feature_dir
    save_plot = args.save_plot
    tag = args.tag
    qc_thr_rmse = args.qc_thr_rmse
    qc_thr_n_peaks = args.qc_thr_n_peaks
    qc_thr_variation = args.qc_thr_variation

    # mkdir the save_*_dir if not existing
    if not os.path.exists(save_image_dir) and save_plot:
        os.makedirs(save_image_dir)
    if not os.path.exists(save_feature_dir):
        os.makedirs(save_feature_dir)

    # Define the empty feature tables
    feature = pd.DataFrame(columns=['genotype', 'experiment', 'Max', 'Max_x', 'Max_y', 's_x', 's_y', 
                                    'rmse', 'n_peaks', 'variation', 'qc_result'])

    # Read the input files as a dict
    print('Read the input data ...')
    dict_data = read_experiments(dir)

    # Loop over all experiments
    print('Extract features ...')
    for exp in tqdm(dict_data.keys()):
        variants = dict_data[exp]
    
        ## Obtain the WT maximum value for normalization
        WT_av = variants['WT'][3]
        max_wt = np.around(WT_av.max().max(), decimals=2) # maximum value of WT to be used for rescaling
    
        ## Loop over variants in a single experiment
        for var in variants.keys():
            ### For each varianr, there exists 3 replicates and the last of (inx=3) contains the median one
            ### We use the median experiment for the analysis, and use the 3 replicates to obtain the QC-variations 
            e1, e2, e3, data = variants[var]
            name = exp + "_" + var
            
            #### 1. QC: replicate variations
            e1, e2, e3 = reshape(e1/max_wt), reshape(e2/max_wt), reshape(e3/max_wt)
            exp_vals = np.concatenate([e1, e2, e3], axis=1)
            median = np.median(exp_vals, axis=1)
            range = np.ptp(exp_vals, axis=1)
            variation = range / (median + 1)
            variation = variation.max()
    
            #### 2. QC: Count the number of peacks
            z_np = data.to_numpy(dtype = "float")
            filtered_z = maximum_filter(z_np, size=3)
            n_peaks = (z_np == filtered_z).sum().sum()
    
            #### 3. Gaussian modelling
            ##### Get x, y, z
            x, y, z = transform_df(data, max_wt, rescale=True)
            if save_plot:
                plot_landscape(x, y, z, name = name, 
                               show = False, save = True, save_dir = save_image_dir)
            x, y, z = np.log(x + eps), np.log(y + eps), z
    
            ##### Curve fitting
            initial_guess = (1.0, 5, 5, 2, 2)
            bounds = ([0, eps, eps, 0.5, 0.5],                       # Lower bounds
                      [120, np.log(3000), np.log(300), 1000, 1000])  # Upper bounds 1500, 150
            popt, pcov = curve_fit(gaussian_2d, (x, y), z, p0=initial_guess, bounds=bounds)
            
            a, mx, my, sx, sy = tuple(popt)
            z_hat = gaussian_2d((x, y), a, mx, my, sx, sy)
            
            #### 4. Calculate RMSE
            mse = np.mean((z/z.max() - z_hat/z_hat.max())**2)
            rmse = np.sqrt(mse)
            if save_plot:
                plot_landscape(np.exp(x), np.exp(y), z_hat, name = name + "_model",
                              show = False, save = True, save_dir = save_image_dir)
    
    
            #### 5. QC check
            if (rmse <= qc_thr_rmse[0]) and (n_peaks <= qc_thr_n_peaks[0]) and (variation <= qc_thr_variation[0]):
                qc_result = 'Pass'
            elif (rmse >= qc_thr_rmse[1]) or (n_peaks >= qc_thr_n_peaks[1]) or (variation >= qc_thr_variation[1]):
                qc_result = 'Fail'
            else:
                qc_result = 'ToCheck'
            
            #### Save features
            feature.loc[len(feature)] = [var, exp, a, mx, my, sx, sy, rmse, n_peaks, variation, qc_result]
    # End of the loop

    # Save features table
    if len(dict_data.keys()) > 1:
        n = ""
    else:
        n = "_" + name
    
    if tag == "":
        t = ""
    else:
        t = "_" + tag
    
    Path_to_save = os.path.join(save_feature_dir, "extracted_features" + n + t + ".csv")
    feature.to_csv(Path_to_save, index=False)
    print("Results are saved as:")
    print(Path_to_save)

    # End
