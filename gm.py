"""
    Gaussian Modelling and feature extraction
    Version 1.6.0
    Authors:
        - Behnam Yousefi (behnm.yousefi@zmnh.uni-hamburg.de; yousefi.bme@gmail.com)
        - Robin Khatri (robin.khatri@zmnh.uni-hamburg.de)
    Scientific Contributor:
        - Polina Gundorova (p.gundorova@uke.de)
    UKE, Hamburg, Germany - March 2024

    To plot the landscapes, the current version uses some codes from a Python script written by
    Italo Balestra (italo.balestra@olt-dss.com) and Andrea Turcati (andrea.turcati@olt-dss.com)
    OmegaLambdaTec GmbH - 2023
"""

# Load the packages
import os
from tqdm import tqdm
import argparse
import json

import pandas as pd
import numpy as np

from scipy import ndimage
from scipy.ndimage import maximum_filter
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator

import matplotlib.pyplot as plt
import matplotlib.colors as colors

eps = 0.000001

# ======================== Helper functions (utilities) ==============================

## read experiment files
def read_experiments(exp):

    cols_1 = ["BLANK", 25, 50, 75, 100, 125, 150, 200, 500, 750, 1250, 2500]
    cols_1 = [str(col) for col in cols_1]
    cols_2 = ["0.0","24.1","48.2","72.3","96.4","120.5","144.6","192.8","482.1","723.1","1205.2","2410.3"]
    cols_2 = [str(col) for col in cols_2]
    cols_final = [0, 25, 50, 75, 100, 125, 150, 200, 500, 750, 1250, 2500]

    rows = [0, 10, 25, 50, 75, 100, 150, 250]
    rows = [int(row) for row in rows]

    # the output file is a dict.
    dict_data = {}

    # if exp is a directory, list all *.xlsm files
    # if exp is a file, it consider only that one
    if os.path.isdir(exp):
        fps = [os.path.join(exp, file) for file in os.listdir(exp) if ".xlsm" in file]
    else:
        if ".xlsm" in exp:
            fps = [exp]
        else:
            raise ValueError("The file should be *.xlsm")

    for fp in tqdm(fps):
        # fp = os.path.join(exp, wb)
        xl = pd.ExcelFile(fp)
        sheets = [sheet for sheet in xl.sheet_names if "-av" in sheet]
        exp_name = fp.split("/")[-1].split("_")[1].split(".")[0]
        dict_data[exp_name] = {}
        for sheet in sheets:
            df = pd.read_excel(fp, sheet_name=sheet, skiprows=4)
            mask = np.column_stack([df[col].astype(str).str.contains("BLANK", na=False) for col in df])
            if mask.any():
                cols = cols_1
            else:
                cols = cols_2
            df.columns = [str(col) for col in df.columns]
            df = df[cols]
            df = df.dropna()
            exps = [df.iloc[0:8], df.iloc[9:17], df.iloc[18:26], df.iloc[27:35]]
            for i in range(len(exps)):
                exps[i].index = rows
                exps[i].columns = cols_final
            dict_data[exp_name][sheet.split("-av")[0]] = exps

    return dict_data


## Transform data frame to z, y, z format
def transform_df(df, max_wt, rescale=True):
    avgs = df.unstack().reset_index(level=[0, 1])
    avgs.columns = ["Phe", "BH4", "Enzyme Activity"]
    if rescale:
        avgs["Enzyme Activity"] = (avgs["Enzyme Activity"] / max_wt) * 100
    avgs[avgs < 0] = 0

    x = avgs.Phe.to_numpy()
    y = avgs.BH4.to_numpy()
    z = np.float64(avgs["Enzyme Activity"].to_numpy())
    z = np.round(z, decimals=2)
    return x, y, z


## Plot the landscape
def plot_landscape(
    x,
    y,
    z,
    name=" ",
    method="regular_grid",
    show=True,
    save=False,
    save_dir = "",
    nbins=1000,
    rescale=True,
    max_val_scale=None,
    legend=[True, True, True, True],
):
    x_dense = np.linspace(0, x.max(), nbins)
    y_dense = np.linspace(0, y.max(), nbins)
    B1, B2 = np.meshgrid(x_dense, y_dense, indexing="xy")
    dense_points = np.stack([B1.ravel(), B2.ravel()], -1)  # shape (N, 2) in 2d

    try:
        if method == "linear_ndi":
            scattered_points = np.stack(
                [x.ravel(), y.ravel()], -1
            )  # shape (N, 2) in 2d
            smooth_z = LinearNDInterpolator(
                scattered_points,
                z.ravel(),
                rescale=True,
                fill_value=0.0,
            )
            z_smoothed = smooth_z(dense_points).reshape(B1.shape)
            z_smoothed = ndimage.gaussian_filter(z_smoothed, sigma=10)
            z_smoothed[z_smoothed < 0] = 0

        if method == "regular_grid":
            Z = z.reshape(len(np.unique(x)), len(np.unique(y)))
            rgi = RegularGridInterpolator(
                (np.unique(x), np.unique(y)),
                Z,
                method="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            Z_rgi = rgi(np.array([B1.flatten(), B2.flatten()]).T).reshape(B1.shape)
            z_smoothed = ndimage.gaussian_filter(Z_rgi, sigma=7)
            z_smoothed[z_smoothed < 0] = 0
            # z_smoothed[z_smoothed > 100] = 100

        elif method not in ["linear_ndi", "regular_grid"]:
            raise NotImplementedError(
                "Smoothing method not implemented. Choose between 'regular_grid' or 'linear_ndi'."
            )

    except NotImplementedError as e:
        print("Error: ", e)
        sys.exit()

    fig, ax = plt.subplots(figsize=(6, 5))
    if rescale:
        im = ax.pcolormesh(
            B1,
            B2,
            z_smoothed,
            cmap="Spectral_r",
            norm=colors.PowerNorm(vmin=0, vmax=100, gamma=0.6),
        )
    else:
        if max_val_scale is None:
            max_val_scale = np.max(z_smoothed)

        im = ax.pcolormesh(
            B1,
            B2,
            z_smoothed,
            cmap="Spectral_r",
            norm=colors.PowerNorm(vmin=0, vmax=max_val_scale, gamma=0.6),
        )
    ax.set_xlim([np.min(x_dense), np.max(x_dense)])
    ax.set_ylim([np.min(y_dense), np.max(y_dense)])
    ax.set_xlabel("Phe [uM]")
    ax.set_ylabel("BH4 [uM]")
    # ax.set_title(f"{name}")

    CS = ax.contour(
        B1, B2, z_smoothed, 5, colors=("lightgrey"), linewidths=1, origin="lower"
    )
    ax.clabel(CS, fmt="%.0f", colors="lightgrey", fontsize=9)
    cbar = fig.colorbar(im, shrink=0.7, ax=ax)
    if rescale:
        cbar.set_label("Enzyme Activity [%]")
    else:
        cbar.set_label("Enzyme Activity")

    df = pd.DataFrame(z_smoothed, columns=x_dense, index=y_dense)

    mx = np.around(df.max().max(), decimals=2)
    bh4 = df.stack().idxmax()[0]
    phe = df.stack().idxmax()[1]

    df_tmp = df.loc[bh4, :]
    phe_min = df_tmp[df_tmp >= mx * 0.5].index.min()
    phe_max = df_tmp[df_tmp >= mx * 0.5].index.max()
    df_tmp2 = df.loc[:, phe]
    bh4_min = df_tmp2[df_tmp2 >= mx * 0.5].index.min()
    bh4_max = df_tmp2[df_tmp2 >= mx * 0.5].index.max()

    output_vals = {
        "Max": mx,
        "Phe": phe,
        "BH4": bh4,
        "50% max BH4 min": bh4_min,
        "50% max BH4 max": bh4_max,
        "50% max Phe min": phe_min,
        "50% max Phe max": phe_max,
    }

    ax.plot(
        phe,
        bh4,
        marker="x",
        c="m",
        markersize=8,
        markeredgecolor="m",
        markeredgewidth=3,
    )

    C50 = ax.contour(
        B1,
        B2,
        z_smoothed,
        levels=[mx * 0.5],
        colors=("magenta",),
        linewidths=1,
        origin="lower",
    )
    ax.clabel(C50, fmt="%.0f", colors="magenta", fontsize=9)

    info_box, max_val, peak_coords, fifty_coords = legend
    if info_box:
        pro = dict(boxstyle="round", facecolor="w", alpha=0.5)
        textstr = ""
        if max_val:
            textstr += f"Max: {mx:.0f}\n"
        if peak_coords:
            textstr += f"Phe: {phe:.0f}\nBH4: {bh4:.0f}\n"
        if fifty_coords:
            textstr += f"50% max BH4:\n{bh4_min:.0f}-{bh4_max:.0f}\n50% max Phe:\n{phe_min:.0f}-{phe_max:.0f}"
        ax.text(
            0.97,
            0.97,
            textstr,
            transform=im.axes.transAxes,
            fontsize=7,
            verticalalignment="top",
            bbox=pro,
            ha="right",
            color="k",
        )
    # im.clim(0, 100)
    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(save_dir, f"landscape_{name}.png"))
    if show:
        plt.show()
    plt.close()

    return output_vals


## Plot the 3D map
def plot_3d_map(
    x,
    y,
    z,
    name,
    method="regular_grid",
    show=False,
    save=True,
    save_dir = "",
    nbins=1000,
    elev=30,
    azim=-120,
    rescale=True,
    max_val_scale=None,
    legend=[True, True, True, True],       # legend here is not being used. Used in plot_3dsurf()
):
    x_dense = np.linspace(0, x.max(), nbins)
    y_dense = np.linspace(0, y.max(), nbins)
    B1, B2 = np.meshgrid(x_dense, y_dense, indexing="xy")
    dense_points = np.stack([B1.ravel(), B2.ravel()], -1)  # shape (N, 2) in 2d

    try:
        if method == "linear_ndi":
            scattered_points = np.stack(
                [x.ravel(), y.ravel()], -1
            )  # shape (N, 2) in 2d
            smooth_z = LinearNDInterpolator(
                scattered_points,
                z.ravel(),
                rescale=True,
                fill_value=0.0,
            )
            z_smoothed = smooth_z(dense_points).reshape(B1.shape)
            z_smoothed = ndimage.gaussian_filter(z_smoothed, sigma=10)

        if method == "regular_grid":
            Z = z.reshape(len(np.unique(x)), len(np.unique(y)))
            rgi = RegularGridInterpolator(
                (np.unique(x), np.unique(y)),
                Z,
                method="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            Z_rgi = rgi(np.array([B1.flatten(), B2.flatten()]).T).reshape(B1.shape)
            z_smoothed = ndimage.gaussian_filter(Z_rgi, sigma=7)
            z_smoothed[z_smoothed < 0] = 0
            # z_smoothed[z_smoothed > 100] = 100

        elif method not in ["linear_ndi", "regular_grid"]:
            raise NotImplementedError(
                "Smoothing method not implemented. Choose between 'regular_grid' or 'linear_ndi'."
            )

    except NotImplementedError as e:
        print("Error: ", e)
        sys.exit()

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")

    if rescale:
        im = ax.plot_surface(
            B1,
            B2,
            z_smoothed,
            cmap="Spectral_r",
            vmin=0,
            vmax=100,
            norm=colors.PowerNorm(vmin=0, vmax=100, gamma=0.6),
        )
        ax.view_init(elev=elev, azim=azim)
    else:
        if max_val_scale is None:
            max_val_scale = np.max(z)

        im = ax.plot_surface(
            B1,
            B2,
            z_smoothed,
            cmap="Spectral_r",
            vmin=0,
            vmax=np.max(z),
            norm=colors.PowerNorm(vmin=0, vmax=max_val_scale, gamma=0.6),
            rstride=10,
            cstride=10,
        )
        ax.view_init(elev=elev, azim=azim)

    cbar = plt.colorbar(im, shrink=0.4, ax=ax)
    if rescale:
        cbar.set_label("Enzyme Activity [%]")
    else:
        cbar.set_label("Enzyme Activity")

    # ax.set_title(name)
    ax.set_xlabel("Phe [uM]")
    ax.set_ylabel("BH4 [uM]")

    if z_smoothed.max() < 5:
        ax.set_zlim([0, 5])

    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(save_dir, f"3dplot_{name}.png"))
    if show:
        plt.show()
    plt.close()


## 2D Gaussian model
def gaussian_2d(xy, a, mx, my, sx, sy):
    x, y = xy
    z = a * np.exp( - 0.5 * ( ((x-mx)**2 / (sx**2)) + ((y-my)**2 / (sy**2)) ) )
    return z

## Reshape experiment table
def reshape(a):
    a = np.array(a, dtype=float)
    a[a<0]=0
    return np.log(a+1).reshape(-1,1)

# =============================== Main ================================
if __name__ == "__main__":

    # Read inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True, help='Input directory; either xlsm file name or a directory.')
    parser.add_argument('--param_dir', default="gm_parameters.json", help='Directory to parameters json file.')
    parser.add_argument('--no_plot', action='store_true', help='No plot will be created if apecified.')

    # Set inputs
    args = parser.parse_args()
    dir = args.input_dir
    param_dir = args.param_dir
    no_plot = args.no_plot

    # set parameters
    ## Read parameters json file and set parameters
    with open(param_dir, 'r') as file:
        parameters = json.load(file)
    ### Save directory
    save_feature_dir = parameters['save_feature_dir']
    save_image2d_dir = parameters['save_image2d_dir']
    save_image3d_dir = parameters['save_image3d_dir']
    tag = parameters['tag']
    ### QC thresholds
    qc_thr_rmse = parameters['qc_thr_rmse']
    qc_thr_n_peaks = parameters['qc_thr_n_peaks']
    qc_thr_variation = parameters['qc_thr_variation']
    ### 3D plot parameters
    sm_method = parameters['sm_method']
    rescale = parameters['rescale']
    max_val_scale = parameters['max_val_scale']
    elev = parameters['elev']
    azim = parameters['azim']
    nbins = parameters['nbins']
    info_box = parameters['info_box']
    max_val = parameters['max_val']
    peak_coords = parameters['peak_coords']
    fifty_coords = parameters['fifty_coords']

    if no_plot:
        save_image2d_dir = ""
        save_image3d_dir = ""

    save_plot2d = False if save_image2d_dir == "" else True
    save_plot3d = False if save_image3d_dir == "" else True

    # mkdir the save_*_dir if not existing
    if not os.path.exists(save_feature_dir):
        os.makedirs(save_feature_dir)
    if not os.path.exists(save_image2d_dir) and save_plot2d:
        os.makedirs(save_image2d_dir)
    if not os.path.exists(save_image3d_dir) and save_plot3d:
        os.makedirs(save_image3d_dir)

    # Define the empty feature tables
    feature = pd.DataFrame(columns=['genotype', 'experiment', 'Max', 'A', 'Max_x', 'Max_y', 's_x', 's_y',
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
            ##### 2D and 3D Plot of the row landscapes
            if save_plot2d:
                plot_landscape(x, y, z, name = name,
                               show = False, save = True, save_dir = save_image2d_dir,
                               method=sm_method, rescale=rescale, max_val_scale=max_val_scale,
                               legend=[info_box, max_val, peak_coords, fifty_coords])
            if save_plot3d:
                plot_3d_map(x, y, z, name = name,
                            show=False, save=True, save_dir = save_image3d_dir,
                            method=sm_method, rescale=rescale, max_val_scale=max_val_scale,
                            elev=elev, azim=azim, nbins=nbins, legend=[info_box, max_val, peak_coords, fifty_coords])

            x, y, z = np.log(x + eps), np.log(y + eps), z

            ##### Curve fitting
            initial_guess = (1.0, 5, 5, 2, 2)
            bounds = ([0, eps, eps, 0.5, 0.5],                       # Lower bounds
                      [120, np.log(1500), np.log(150), 1000, 1000])  # Upper bounds
            popt, pcov = curve_fit(gaussian_2d, (x, y), z, p0=initial_guess, bounds=bounds)

            a, mx, my, sx, sy = tuple(popt)
            z_hat = gaussian_2d((x, y), a, mx, my, sx, sy)

            #### 4. Calculate RMSE
            mse = np.mean((z/z.max() - z_hat/z_hat.max())**2)
            rmse = np.sqrt(mse)
            ##### 2D and 3D Plot of the modeled landscapes
            if save_plot2d:
                plot_landscape(np.exp(x), np.exp(y), z_hat, name = name + "_model",
                               show = False, save = True, save_dir = save_image2d_dir,
                               method=sm_method, rescale=rescale, max_val_scale=max_val_scale,
                               legend=[info_box, max_val, peak_coords, fifty_coords])
            if save_plot3d:
                plot_3d_map(np.exp(x), np.exp(y), z_hat, name = name + "_model",
                            show=False, save=True, save_dir = save_image3d_dir,
                            method=sm_method, rescale=rescale, max_val_scale=max_val_scale,
                            elev=elev, azim=azim, nbins=nbins, legend=[info_box, max_val, peak_coords, fifty_coords])

            #### 5. QC check
            if (rmse <= qc_thr_rmse[0]) and (n_peaks <= qc_thr_n_peaks[0]) and (variation <= qc_thr_variation[0]):
                qc_result = 'Pass'
            elif (rmse >= qc_thr_rmse[1]) or (n_peaks >= qc_thr_n_peaks[1]) or (variation >= qc_thr_variation[1]):
                qc_result = 'Fail'
            else:
                qc_result = 'ToCheck'

            #### Save features
            feature.loc[len(feature)] = [var, exp, np.max(z_hat), a, np.exp(mx), np.exp(my), np.exp(sx), np.exp(sy), 
                                         rmse, n_peaks, variation, qc_result]
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
