#!/usr/bin/env python3
#-*- coding:utf-8 -*-

"""
A script to plot the BUSCO result
----------------------------------
Make sure you have used the *SAME* BUSCO database for all the samples

Usage:
    python3 Busco_plot.py -i <input_dir> -o <output_dir> -w <width>

    input_dir: A file containing the path to every BUSCO report, one per line
    output_dir: The directory to save the plot  (default: ./)
    width: A float number represent the width of the bar (default: 0.5)

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

def plot_busco(input_file, output_dir, width):
    """
    This function is used to plot the BUSCO result
    Inlcuding the results from TransDecoder, Braker, GeMoMa and the final EVM
    """
    # Get the BUSCO report from input file
    reports = []
    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            reports.append(line.strip())

    sample = np.array([])
    single_busco_num = np.array([])
    duplicated_busco_num = np.array([])
    fragmented_busco_num = np.array([])
    missing_busco_num = np.array([])
    total_busco_num = np.array([])

    single_busco_percent = np.array([])
    duplicated_busco_percent = np.array([])
    fragmented_busco_percent = np.array([])
    missing_busco_percent = np.array([])


    for report in reports:
        # Obtain the sample name
        basename = os.path.basename(report).split(".")[-2]
        with open(report, "r", encoding="utf-8") as f:
            # Set default values
            single_busco = 0
            duplicated_busco = 0
            fragmented_busco = 0
            missing_busco = 0
            total_busco = 0

            single_busco_per = 0
            duplicated_busco_per = 0
            fragmented_busco_per = 0
            missing_busco_per = 0

            for line in f:
                # Skip the unuseful lines
                if line.startswith("#"):
                    continue
                if re.search(r'^\s*$', line):
                    continue
                if re.search(r'Results:', line):
                    continue
                if re.search(r'Complete BUSCOs', line):
                    continue
                # Get the BUSCO result
                if re.search(r'\tC:\d+\.\d+', line):
                    percent_result = line.split("\t")[1]
                    single_busco_per = float(re.search(r"S:\d+\.\d+", percent_result).group().split(":")[1])
                    duplicated_busco_per = float(re.search(r"D:\d+\.\d+", percent_result).group().split(":")[1])
                    fragmented_busco_per = float(re.search(r"F:\d+\.\d+", percent_result).group().split(":")[1])
                    missing_busco_per = float(re.search(r"M:\d+\.\d+", percent_result).group().split(":")[1])
                elif re.search(r'Complete and single-copy BUSCOs', line):
                    single_busco = int(re.search(r'\d+', line).group())
                elif re.search(r'Complete and duplicated BUSCOs', line):
                    duplicated_busco = int(re.search(r'\d+', line).group())
                elif re.search(r'Fragmented BUSCOs', line):
                    fragmented_busco = int(re.search(r'\d+', line).group())
                elif re.search(r'Missing BUSCOs', line):
                    missing_busco = int(re.search(r'\d+', line).group())
                elif re.search(r'Total BUSCO groups searched', line):
                    total_busco = int(re.search(r'\d+', line).group())
                else:
                    print("Error: Unrecognized BUSCO line, please check the BUSCO report!")
                    sys.exit(1)

            # Append the BUSCO result to the numpy array
            sample = np.append(sample, basename)
            single_busco_num = np.append(single_busco_num, single_busco)
            duplicated_busco_num = np.append(duplicated_busco_num, duplicated_busco)
            fragmented_busco_num = np.append(fragmented_busco_num, fragmented_busco)
            missing_busco_num = np.append(missing_busco_num, missing_busco)
            total_busco_num = np.append(total_busco_num, total_busco)

            single_busco_percent = np.append(single_busco_percent, single_busco_per)
            duplicated_busco_percent = np.append(duplicated_busco_percent, duplicated_busco_per)
            fragmented_busco_percent = np.append(fragmented_busco_percent, fragmented_busco_per)
            missing_busco_percent = np.append(missing_busco_percent, missing_busco_per)

    # Plot the BUSCO result
    category_colors = ["#4fcbe9", "#16519f", "#f8dd2e", "#f07e74",]
    left = np.zeros(sample.size)

    # Plot the BUSCO number result
    fig1, ax1 = plt.subplots(figsize=(6, 6))

    ax1.invert_yaxis()
    ax1.barh(sample, single_busco_num, width, label="Complete and Single BUSCO", left=left, color=category_colors[0])
    left += single_busco_num
    ax1.barh(sample, duplicated_busco_num, width, label="Complete and Duplicated BUSCO", left=left, color=category_colors[1])
    left += duplicated_busco_num
    ax1.barh(sample, fragmented_busco_num, width, label="Fragmented BUSCO", left=left, color=category_colors[2])
    left += fragmented_busco_num
    ax1.barh(sample, missing_busco_num, width, label="Missing BUSCO", left=left, color=category_colors[3])
    ax1.axvline(1453, linestyle="--", color="grey")

    ax1.set_xlim(0, max(total_busco_num))
    ax1.set_title('BUSCO Assessment Results', loc="center", y=1.12)
    ax1.legend(ncol=2, bbox_to_anchor=(0.5, 1), loc='lower center', fontsize='small')
    for i in range(sample.size):
        result = f"C:[S:{single_busco_num[i]}; D:{duplicated_busco_num[i]}]; F:{fragmented_busco_num[i]}; M:{missing_busco_num[i]}"
        ax1.annotate(result, (10,i), fontsize='x-small')

    xlabels = ax1.get_xticklabels()
    plt.setp(xlabels, rotation=45, horizontalalignment='center')
    fig1.savefig(os.path.join(output_dir, 'all_busco_summary_num.pdf'), dpi=300, bbox_inches="tight")

    # Plot the BUSCO percentage result
    left = np.zeros(sample.size)
    fig2, ax2 = plt.subplots(figsize=(6, 6))

    ax2.invert_yaxis()
    ax2.barh(sample, single_busco_percent, width, label="Complete and Single BUSCO", left=left, color=category_colors[0])
    left += single_busco_percent
    ax2.barh(sample, duplicated_busco_percent, width, label="Complete and Duplicated BUSCO", left=left, color=category_colors[1])
    left += duplicated_busco_percent
    ax2.barh(sample, fragmented_busco_percent, width, label="Fragmented BUSCO", left=left, color=category_colors[2])
    left += fragmented_busco_percent
    ax2.barh(sample, missing_busco_percent, width, label="Missing BUSCO", left=left, color=category_colors[3])
    ax2.axvline(90, linestyle="--", color="grey")

    ax2.set_xlim(0, 100)
    ax2.set_title('BUSCO Assessment Results', loc="center", y=1.12)
    ax2.legend(ncol=2, bbox_to_anchor=(0.5, 1), loc='lower center', fontsize='small')
    for i in range(sample.size):
        result = f"C:[S:{single_busco_percent[i]}%; D:{duplicated_busco_percent[i]}%]; F:{fragmented_busco_percent[i]}%; M:{missing_busco_percent[i]}%"
        ax2.annotate(result, (1, i), fontsize='x-small')
    
    xlabels = ax2.get_xticklabels()
    plt.setp(xlabels, rotation=45, horizontalalignment='center')
    fig2.savefig(os.path.join(output_dir, 'all_busco_summary_percent.pdf'), dpi=300, bbox_inches="tight")

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True, help="The input file")
    parser.add_argument("-o", "--output", type=str, help="The output file", default="./")
    parser.add_argument("-w", "--width", type=float, help="Width of the bar", default=0.5)
    args = parser.parse_args()

    plot_busco(args.input, args.output, args.width)

if __name__ == "__main__":
    main()
