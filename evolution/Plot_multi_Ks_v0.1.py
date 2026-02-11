#!/usr/bin/env python3

"""
Plot_multi_Ks.py

Description:
    Plot Ks distribution for multiple species with Nature/Science style palette.
    - Reads Ks values from multiple files (one per species).
    - Automatically determines max Ks to plot based on quantiles.
    - Plots histogram + KDE curve for each species.
    - Detects and labels WGD peaks on the KDE curve.

Usage:
    python Plot_multi_Ks.py
        -i species1_Ks.csv species2_Ks.csv species3_Ks.csv
        -l "Species 1" "Species 2" "Species 3" 
        --column Ks 
        --out_prefix Ks_multi_species_nature_style 
        --save_peak_table

Note:
    1.  kde bandwidth adjustment (--bw_adjust) affects peak detection. Larger values yield smoother curves with fewer peaks.

    2.  peak detection relays on "--peak_prominence" and "--peak_distance" parameters. 
        Adjust them if needed (also need to consider ks_max and kde_points).
        
        peak_prominence: 
            How much a peak stands out due to its height. Larger value means fewer peaks detected.
            Default 0.02 is a good starting point, but you can increase it to 0.03 or higher for more stringent detection.
        
        peak_distance: 
            Minimum distance (in number of kde points) between peaks. Larger value means peaks must be more separated.
            Default 60 is a good starting point for 1200 KDE points and ks_max=10, which means peaks must be at least 0.5 Ks apart.
        
        you can calculate the point distance corresponding to a given Ks distance as:

              point_distance = (ks_distance / ks_max) * total_kde_points

        Suggestions:
          --peak_distance: <40 for closely spaced peaks, 40-80 for moderate spacing, >100 for widely spaced peaks.
          --peak_prominence: 0.01 for sensitive peak detection, 0.02 for balanced, 0.03+ for stringent detection.

    3.  kde_points affects the resolution of the KDE curve. More points yield smoother curves but increase computation time.
        Adjust according to your data size and desired detail level. Default 1200 is a good balance for most datasets. (800-2000 is typical)
    
    4.  quantile and buffer for auto ks_max:
        The auto ks_max is determined by the specified quantile of the Ks values plus a buffer. 
        Default quantile 0.98 means ks_max will be set to the value below which 98% of Ks values fall, plus the buffer (default 0.3). 
        This helps to focus on the main distribution and avoid long tails. 
        Adjust quantile (e.g., 0.95 for more inclusive, 0.99 for more exclusive) and buffer (e.g., 0.2 for tighter, 0.5 for looser) as needed based on your data distribution.
    
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import math


# Nature/Science常用配色 (Okabe-Ito palette)
OKABE_ITO = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#0072B2",  # blue
    "#000000",  # black
    "#F0E442",  # yellow
]


def get_color(idx):
    return OKABE_ITO[idx % len(OKABE_ITO)]


def read_ks_file(file, column=None):
    """
    Read Ks values from file.
    - If column is None: assume file contains only one Ks column
    - If column provided: use that column name
    """
    try:
        df = pd.read_csv(file, sep=None, engine="python")
    except Exception:
        df = pd.read_csv(file, header=None)

    if column is None:
        if df.shape[1] == 1:
            ks = df.iloc[:, 0] # first column
        else:
            raise ValueError(f"[ERROR] {file} has multiple columns, please specify --column")
    else:
        if column not in df.columns:
            raise ValueError(f"[ERROR] Column '{column}' not found in {file}. Available: {list(df.columns)}")
        ks = df[column]

    ks = pd.to_numeric(ks, errors="coerce") # convert to numeric, non-convertible become NaN
    ks = ks.replace([np.inf, -np.inf], np.nan).dropna() # remove NaN and inf
    ks = ks[ks >= 0] # remove negative Ks

    return ks.values


def auto_ks_max(ks, quantile=0.98, buffer=0.3, hard_max=10.0):
    if len(ks) == 0:
        return 5.0

    q = np.quantile(ks, quantile)
    ks_max = q + buffer
    ks_max = max(ks_max, 1.0)
    ks_max = min(ks_max, hard_max)

    return ks_max


def compute_kde(ks_values, ks_max, kde_points=1200, bw_adjust=1.2):
    """
    Compute KDE curve.
    bw_adjust: >1 means smoother, <1 means sharper.
    """
    kde = gaussian_kde(ks_values)
    kde.set_bandwidth(kde.factor * bw_adjust)

    x = np.linspace(0, ks_max, kde_points)
    y = kde(x)
    return x, y


def detect_peaks(x, y, prominence=0.02, distance=60, max_peaks=5):
    """
    Detect peaks in the KDE curve.
    """
    peaks, props = find_peaks(y, prominence=prominence, distance=distance)

    if len(peaks) == 0:
        return []

    peak_list = [(x[p], y[p], props["prominences"][i]) for i, p in enumerate(peaks)]
    peak_list = sorted(peak_list, key=lambda z: z[2], reverse=True)

    return peak_list[:max_peaks]


def setup_style():
    """
    Nature like style settings for matplotlib.
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 9
    })


def plot_multi_species(files, labels, column,
                       ks_max, auto_quantile, auto_buffer, hard_max,
                       bins, kde_points, bw_adjust,
                       peak_prominence, peak_distance_ks, max_peaks, add_line,
                       fig_width, fig_height, n_cols,
                       out_prefix, dpi, save_peak_table):

    setup_style()

    all_species_ks = []
    for file in files:
        ks = read_ks_file(file, column=column)
        all_species_ks.append(ks)

    if ks_max is None:
        merged = np.concatenate(all_species_ks)
        ks_max = auto_ks_max(
            merged,
            quantile=auto_quantile,
            buffer=auto_buffer,
            hard_max=hard_max
        )
        print(f"[INFO] Auto ks_max determined: {ks_max:.2f}")

    n = len(files)
    # for more square layout
    #ncols = math.ceil(math.sqrt(n)) 
    ncols = n_cols
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width * ncols, fig_height * nrows))
    axes = np.array(axes).reshape(-1)

    peak_records = []

    for idx, (ks, label) in enumerate(zip(all_species_ks, labels)):
        ax = axes[idx]
        color = get_color(idx)

        ks = ks[ks <= ks_max]

        if len(ks) < 20:
            ax.text(0.5, 0.5, f"{label}\nToo few Ks values",
                    ha="center", va="center", fontsize=11)
            ax.set_axis_off()
            continue

        # histogram
        ax.hist(
            ks,
            bins=bins,
            density=True,
            alpha=0.30,
            edgecolor="none",
            color=color
        )

        # KDE curve
        x, y = compute_kde(ks, ks_max=ks_max, kde_points=kde_points, bw_adjust=bw_adjust)
        ax.plot(x, y, linewidth=2.2, color=color)

        # change peak_distance from Ks units to number of points
        dx = (ks_max - 0.0) / float(kde_points - 1)

        if dx <= 0:
            peak_distance_points = 1
        else:
            peak_distance_points = int(round(peak_distance_ks / dx))
            peak_distance_points = max(1, peak_distance_points)
        
        print(f"[INFO] peak_distance_ks={peak_distance_ks} -> peak_distance_points={peak_distance_points} (dx={dx:.5f})")
        
        # detect peaks on KDE curve
        peaks = detect_peaks(
            x, y,
            prominence=peak_prominence,
            distance=peak_distance_points,
            max_peaks=max_peaks
        )

        # Sort peaks by Ks position (small to large) -> Peak1/Peak2...
        peaks_sorted = sorted(peaks, key=lambda z: z[0])

        for peak_idx, (pk_x, pk_y, pk_prom) in enumerate(peaks_sorted, start=1):
            peak_name = f"Peak{peak_idx}"

            # decide draw line or not
            if add_line:
                ax.axvline(pk_x, linestyle="--", linewidth=1.4, color=color, alpha=0.9)

                ax.text(pk_x, pk_y,
                        f"{peak_name}\nKs={pk_x:.2f}",
                        fontsize=9,
                        ha="center",
                        va="bottom",
                        color=color)

            peak_records.append({
                "Species": label,
                "Peak": peak_name,
                "Peak_Ks": pk_x,
                "Peak_density": pk_y,
                "Prominence": pk_prom,
                "N_Ks": len(ks)
            })

        #ax.set_title(f"{label} (n={len(ks)})", style='italic', fontsize=14)
        ax.set_title(f"{label}", style='italic', fontsize=14)
        ax.set_xlim(0, ks_max)
        ax.set_xlabel("Ks")
        ax.set_ylabel("Density")

    # remove unused subplots
    for j in range(n, len(axes)):
        axes[j].set_axis_off()

    #fig.suptitle("Ks Distribution (KDE + Peak Labels)", fontsize=16, y=1.02)
    plt.tight_layout()

    pdf_file = f"{out_prefix}.pdf"
    png_file = f"{out_prefix}.png"

    fig.savefig(pdf_file, dpi=dpi, bbox_inches="tight")
    fig.savefig(png_file, dpi=dpi, bbox_inches="tight")

    print(f"[Done] Saved PDF: {pdf_file}")
    print(f"[Done] Saved PNG: {png_file}")

    if save_peak_table:
        peak_df = pd.DataFrame(peak_records)
        peak_out = f"{out_prefix}.peaks.tsv"
        peak_df.to_csv(peak_out, sep="\t", index=False)
        print(f"[Done] Peak table saved: {peak_out}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Input Ks files (one per species).")
    parser.add_argument("-l", "--labels", nargs="+", required=True,
                        help="Species labels (same order as input).")
    parser.add_argument("-c", "--column", default=None,
                        help="Column name containing Ks values (optional).")

    parser.add_argument("--ks_max", type=float, default=None,
                        help="Max Ks to plot (default: auto).")

    parser.add_argument("--auto_quantile", type=float, default=0.98,
                        help="Quantile for auto ks_max (default: 0.98).")
    parser.add_argument("--auto_buffer", type=float, default=0.3,
                        help="Buffer added to quantile (default: 0.3).")
    parser.add_argument("--hard_max", type=float, default=10.0,
                        help="Hard cap for ks_max (default: 10.0).")

    parser.add_argument("--bins", type=int, default=80,
                        help="Histogram bins (default: 80).")
    parser.add_argument("--kde_points", type=int, default=1200,
                        help="Number of KDE points (default: 1200).")
    parser.add_argument("--bw_adjust", type=float, default=1.2,
                        help="KDE bandwidth adjust (default: 1.2). Larger=more smooth.")

    parser.add_argument("--peak_prominence", type=float, default=0.02,
                        help="Peak prominence threshold (default: 0.02).")
    parser.add_argument("--peak_distance_ks", type=float, default=0.25,
                        help="Min distance between peaks in Ks units (default: 0.25).")
    parser.add_argument("--max_peaks", type=int, default=5,
                        help="Max peaks per species (default: 5).")
    parser.add_argument("--add_line", action="store_true",
                        help="Add vertical line at peak positions.")
    
    parser.add_argument("--fig_width", type=float, default=4.5,
                        help="Figure width per subplot (default: 4.5).")
    parser.add_argument("--fig_height", type=float, default=3.6,
                        help="Figure height per subplot (default: 3.6).")
    parser.add_argument("--n_cols", type=int, default=1,
                        help="Number of columns in the figure (default: 1).")
    
    parser.add_argument("-o", "--out_prefix", default="Ks_distribution_NatureStyle",
                        help="Output prefix (default: Ks_distribution_NatureStyle).")
    parser.add_argument("--dpi", type=int, default=600,
                        help="Output resolution (default: 600).")

    parser.add_argument("--save_peak_table", action="store_true",
                        help="Save peak table (.tsv).")

    args = parser.parse_args()

    if len(args.input) != len(args.labels):
        raise ValueError("[ERROR] Number of input files must match number of labels!")

    plot_multi_species(
        files=args.input,
        labels=args.labels,
        column=args.column,
        ks_max=args.ks_max,
        auto_quantile=args.auto_quantile,
        auto_buffer=args.auto_buffer,
        hard_max=args.hard_max,
        bins=args.bins,
        kde_points=args.kde_points,
        bw_adjust=args.bw_adjust,
        peak_prominence=args.peak_prominence,
        peak_distance_ks=args.peak_distance_ks,
        max_peaks=args.max_peaks,
        add_line=args.add_line,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
        n_cols=args.n_cols,
        out_prefix=args.out_prefix,
        dpi=args.dpi,
        save_peak_table=args.save_peak_table
    )

if __name__ == "__main__":
    main()
