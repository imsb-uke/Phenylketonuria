import os
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from tqdm import tqdm
import re


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
        # sheets = [sheet for sheet in xl.sheet_names if "-av" in sheet]
        sheets = [sheet for sheet in xl.sheet_names if len(re.findall("-av$", sheet))]
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


## Transform df z, y, z
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
    ax.set_title(f"{name}")

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

eps = 0.000001