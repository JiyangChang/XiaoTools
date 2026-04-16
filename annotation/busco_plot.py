#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
A script to plot the BUSCO result

Make sure you have used the *SAME* BUSCO database for all the samples

Usage:
    python3 Busco_plot.py -i <input_file> -o <output_dir> -w <width> -p <prefix>

    input_file: A file containing the path to every BUSCO report, one per line
    output_dir: The directory to save the plot  (default: ./)
    width: A float number represent the width of the bar (default: 0.5)
    prefix: Prefix of the output files (default: "all")

Output:
    complete_busco_num.pdf: The plot of the BUSCO result in number
    complete_busco_percent.pdf: The plot of the BUSCO result in percentage

Author: Chang, Jiyang

"""

import os
import sys
import re
import argparse

import numpy as np
import matplotlib.pyplot as plt


__author__ = "Chang, Jiyang"
__mail__ = "jiyang.chang@psb.ugent.be"

_RE_SUMMARY = re.compile(r'\tC:\d+\.\d+')
_RE_FIELDS = re.compile(r'[SDFM]:(\d+\.\d+)')
_RE_NUM = re.compile(r'\d+')

_LABELS = [
    "Complete and Single BUSCO",
    "Complete and Duplicated BUSCO",
    "Fragmented BUSCO",
    "Missing BUSCO",
]
_COLORS = ["#4fcbe9", "#16519f", "#f8dd2e", "#f07e74"]


def _parse_report(report):
    basename = os.path.basename(report).split(".")[-2]
    single = duplicated = fragmented = missing = total = 0
    single_p = dup_p = frag_p = miss_p = 0.0

    with open(report, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            if "Results:" in line or "Complete BUSCOs" in line:
                continue
            if _RE_SUMMARY.search(line):
                vals = _RE_FIELDS.findall(line.split("\t")[1])
                single_p, dup_p, frag_p, miss_p = (float(v) for v in vals)
            elif "Complete and single-copy BUSCOs" in line:
                single = int(_RE_NUM.search(line).group())
            elif "Complete and duplicated BUSCOs" in line:
                duplicated = int(_RE_NUM.search(line).group())
            elif "Fragmented BUSCOs" in line:
                fragmented = int(_RE_NUM.search(line).group())
            elif "Missing BUSCOs" in line:
                missing = int(_RE_NUM.search(line).group())
            elif "Total BUSCO groups searched" in line:
                total = int(_RE_NUM.search(line).group())
            elif "Dependencies and versions:" in line:
                break
            else:
                print("Error: Unrecognized BUSCO line, please check the BUSCO report!")
                sys.exit(1)

    return basename, single, duplicated, fragmented, missing, total, single_p, dup_p, frag_p, miss_p


def _draw_bars(ax, samples, arrays, width, xlim, vline_x, annotations, annot_x):
    ax.invert_yaxis()
    left = np.zeros(len(samples))
    for values, label, color in zip(arrays, _LABELS, _COLORS):
        ax.barh(samples, values, width, label=label, left=left, color=color)
        left += values
    ax.axvline(vline_x, linestyle="--", color="grey")
    ax.set_xlim(0, xlim)
    ax.set_title("BUSCO Assessment Results", loc="center", y=1.12)
    ax.legend(ncol=2, bbox_to_anchor=(0.5, 1), loc="lower center", fontsize="small")
    for i, text in enumerate(annotations):
        ax.annotate(text, (annot_x, i), fontsize="x-small")
    plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment="center")


def plot_busco(input_file, output_dir, width, prefix):
    with open(input_file, "r", encoding="utf-8") as f:
        reports = [line.strip() for line in f if line.strip()]

    samples, singles, dups, frags, missings, totals = [], [], [], [], [], []
    singles_p, dups_p, frags_p, missings_p = [], [], [], []

    for report in reports:
        name, s, d, f, m, t, sp, dp, fp, mp = _parse_report(report)
        samples.append(name)
        singles.append(s); dups.append(d); frags.append(f); missings.append(m); totals.append(t)
        singles_p.append(sp); dups_p.append(dp); frags_p.append(fp); missings_p.append(mp)

    samples = np.array(samples)
    singles, dups, frags, missings = map(np.array, [singles, dups, frags, missings])
    singles_p, dups_p, frags_p, missings_p = map(np.array, [singles_p, dups_p, frags_p, missings_p])

    num_annotations = [
        f"C:[S:{singles[i]}; D:{dups[i]}]; F:{frags[i]}; M:{missings[i]}"
        for i in range(len(samples))
    ]
    fig1, ax1 = plt.subplots(figsize=(6, 6))
    _draw_bars(ax1, samples, [singles, dups, frags, missings], width,
               max(totals), round(0.9 * totals[0]), num_annotations, annot_x=10)
    fig1.savefig(os.path.join(output_dir, f"{prefix}_busco_summary_num.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig1)

    pct_annotations = [
        f"C:[S:{singles_p[i]}%; D:{dups_p[i]}%]; F:{frags_p[i]}%; M:{missings_p[i]}%"
        for i in range(len(samples))
    ]
    fig2, ax2 = plt.subplots(figsize=(6, 6))
    _draw_bars(ax2, samples, [singles_p, dups_p, frags_p, missings_p], width,
               100, 90, pct_annotations, annot_x=1)
    fig2.savefig(os.path.join(output_dir, f"{prefix}_busco_summary_percent.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig2)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True, help="The input file")
    parser.add_argument("-o", "--output", type=str, help="The output directory", default="./")
    parser.add_argument("-w", "--width", type=float, help="Width of the bar", default=0.5)
    parser.add_argument("-p", "--prefix", type=str, help="Prefix of the output files", default="all")
    args = parser.parse_args()

    plot_busco(args.input, args.output, args.width, args.prefix)


if __name__ == "__main__":
    main()
