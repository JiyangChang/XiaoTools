#!/usr/bin/env python3

"""
TE Overlap Filtering for gene annotation

This script takes a GFF file and a TE annotation file, identifies overlaps between TEs and gene features (CDS, exon, intron), 
and classifies genes based on the extent of TE overlap. 

It outputs statistics for different thresholds and can optionally produce filtered GFF files.

The main steps include:
1. Parsing the input GFF to extract gene, CDS, exon, and intron information
2. Converting TE annotations to BED format
3. Using bedtools to find overlaps between TEs and gene features
4. Classifying genes based on TE overlap ratios
5. Outputting statistics and optionally filtered GFF files

Input:
- GFF file with gene annotations
- TE annotation file (can be GFF or BED format)

Usage:
    python TE_scan_filter.py --gff input.gff3 --te te_annotation.gff3 --prefix output_prefix

Parameters:
- --gff: Input GFF file with gene annotations
- --te: Input TE annotation file (GFF or BED)
- --prefix: Prefix for output files, default is 'tmp'
- --scan_output_gff: If set, outputs filtered GFF files for each threshold combination

"""

import argparse
import pandas as pd
import subprocess
from collections import defaultdict


# -----------------------------
# Parse GFF
# -----------------------------

def parse_gff(gff):
    """Parse GFF file to extract gene, CDS, exon, and intron information."""
    genes = {}
    mrna_to_gene = {}
    cds = defaultdict(list)
    exon = defaultdict(list)

    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, _, _, attr = parts
            start, end = int(start), int(end)
            attrs = dict(x.split('=') for x in attr.split(';') if '=' in x)

            if feature == 'gene':
                genes[attrs.get('ID')] = (chrom, start, end)

            elif feature in ('mRNA', 'transcript'):
                mrna_to_gene[attrs.get('ID')] = attrs.get('Parent')

            elif feature == 'CDS':
                gid = mrna_to_gene.get(attrs.get('Parent'))
                if gid:
                    cds[gid].append((start, end))

            elif feature == 'exon':
                gid = mrna_to_gene.get(attrs.get('Parent'))
                if gid:
                    exon[gid].append((start, end))

    return genes, cds, exon


# -----------------------------
# Build intron
# -----------------------------

def build_introns(exon):
    """Build intron regions from exon annotations."""
    intron = defaultdict(list)
    for gid, regions in exon.items():
        regions = sorted(regions)
        for i in range(len(regions) - 1):
            intron[gid].append((regions[i][1], regions[i+1][0]))
    return intron


# -----------------------------
# TE conversion
# -----------------------------

def convert_te_to_bed(te_file, prefix):
    """Convert TE annotation file to BED format if it's in GFF format."""
    out = f"{prefix}.te.bed"
    with open(te_file) as fin, open(out, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 9:
                chrom, _, _, start, end = parts[:5]
            else:
                chrom, start, end = parts[:3]
            fout.write(f"{chrom}\t{start}\t{end}\n")
    return out


# -----------------------------
# BED writers
# -----------------------------

def write_bed(regions_dict, genes, outfile):
    """Write BED file for given regions (CDS, exon, intron) with gene IDs."""
    with open(outfile, 'w') as f:
        for gid, regions in regions_dict.items():
            chrom = genes[gid][0]
            for s, e in regions:
                f.write(f"{chrom}\t{s}\t{e}\t{gid}\n")


def write_gene_bed(genes, outfile):
    """Write BED file for gene regions."""
    with open(outfile, 'w') as f:
        for gid, (chrom, start, end) in genes.items():
            f.write(f"{chrom}\t{start}\t{end}\t{gid}\n")


# -----------------------------
# Run command
# -----------------------------

def run(cmd):
    """Run a shell command."""
    subprocess.run(cmd, shell=True, check=True)


# -----------------------------
# Parse overlap
# -----------------------------

def parse_overlap(file):
    """Parse bedtools overlap output to calculate total TE overlap for each gene."""
    res = defaultdict(int)
    with open(file) as f:
        for line in f:
            parts = line.strip().split()
            gid = parts[3]
            overlap = int(parts[-1])
            res[gid] += overlap
    return res


# -----------------------------
# Classification
# -----------------------------

def classify_gene(c_ratio, e_ratio, i_ratio):
    """Classify gene based on TE overlap ratios."""
    if c_ratio > 0.3:
        return "TE_CDS"
    elif e_ratio > 0.3:
        return "TE_exon"
    elif i_ratio > 0.3:
        return "TE_intron"
    else:
        return "clean"


# -----------------------------
# Unified filtering
# -----------------------------

def filter_genes(genes, cds, exon, intron, overlaps, cds_high, gene_high):
    """Filter genes based on TE overlap thresholds."""
    keep = set()
    stats = []

    for gid, (chrom, start, end) in genes.items():
        gene_len = end - start
        cds_len = sum(e - s for s, e in cds.get(gid, []))
        exon_len = sum(e - s for s, e in exon.get(gid, []))
        intron_len = sum(e - s for s, e in intron.get(gid, []))

        # Calculate TE overlap ratios
        g_ratio = overlaps['gene'].get(gid, 0) / gene_len if gene_len else 0
        c_ratio = overlaps['cds'].get(gid, 0) / cds_len if cds_len else 0
        e_ratio = overlaps['exon'].get(gid, 0) / exon_len if exon_len else 0
        i_ratio = overlaps['intron'].get(gid, 0) / intron_len if intron_len else 0

        # Combined score (weighted average of CDS and exon ratios)
        score = 0.6 * c_ratio + 0.4 * e_ratio

        remove = False
        # Apply thresholds, prioritizing CDS ratio, then combined score, then gene ratio
        if c_ratio > cds_high:
            remove = True
        elif score > cds_high:
            remove = True
        elif g_ratio > gene_high:
            remove = True

        gene_type = classify_gene(c_ratio, e_ratio, i_ratio)

        stats.append([gid, g_ratio, c_ratio, e_ratio, i_ratio, score, gene_type, remove])

        if not remove:
            keep.add(gid)

    return keep, stats


# -----------------------------
# Write GFF
# -----------------------------

def write_gff(input_gff, output_gff, keep_set):
    """Write filtered GFF file containing only genes in keep_set."""
    with open(input_gff) as fin, open(output_gff, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            if '\tgene\t' in line:
                gid = line.split('ID=')[1].split(';')[0]
                if gid in keep_set:
                    fout.write(line)
            else:
                fout.write(line)


# -----------------------------
# Main
# -----------------------------

def main():
    """Main function to execute TE overlap filtering."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gff', required=True, help="Input GFF file with gene annotations")
    ap.add_argument('--te', required=True, help="Input TE annotation file")
    ap.add_argument('--prefix', default='tmp', help="Prefix for temporary files, default is 'tmp'")

    ap.add_argument('--scan_output_gff', action='store_true', help="If set, outputs filtered GFF files for each threshold combination")

    args = ap.parse_args()

    genes, cds, exon = parse_gff(args.gff)
    intron = build_introns(exon)
    te_bed = convert_te_to_bed(args.te, args.prefix)

    prefix = args.prefix

    # BED
    write_gene_bed(genes, f"{prefix}.gene.bed")
    write_bed(cds, genes, f"{prefix}.cds.bed")
    write_bed(exon, genes, f"{prefix}.exon.bed")
    write_bed(intron, genes, f"{prefix}.intron.bed")

    # overlap
    run(f"bedtools intersect -a {prefix}.gene.bed -b {te_bed} -wo > {prefix}.gene.ov")
    run(f"bedtools intersect -a {prefix}.cds.bed -b {te_bed} -wo > {prefix}.cds.ov")
    run(f"bedtools intersect -a {prefix}.exon.bed -b {te_bed} -wo > {prefix}.exon.ov")
    run(f"bedtools intersect -a {prefix}.intron.bed -b {te_bed} -wo > {prefix}.intron.ov")

    overlaps = {
        'gene': parse_overlap(f"{prefix}.gene.ov"),
        'cds': parse_overlap(f"{prefix}.cds.ov"),
        'exon': parse_overlap(f"{prefix}.exon.ov"),
        'intron': parse_overlap(f"{prefix}.intron.ov")
    }

    cds_range = [round(i * 0.1, 1) for i in range(3, 8)]
    gene_range = [round(i * 0.1, 1) for i in range(4, 9)]

    results = []

    for cds_t in cds_range:
        for gene_t in gene_range:
            keep, stats = filter_genes(genes, cds, exon, intron, overlaps, cds_t, gene_t)

            df = pd.DataFrame(stats, columns=[
                'gene', 'gene_TE_ratio', 'cds_TE_ratio', 'exon_TE_ratio', 'intron_TE_ratio', 'combined_score', 'gene_type', 'removed'
            ])

            stats_file = f"{prefix}.cds{cds_t}.gene{gene_t}.stats.tsv"
            df.to_csv(stats_file, sep='\t', index=False)

            #mean_cds = df['cds_TE_ratio'].mean()
            #mean_intron = df['intron_TE_ratio'].mean()

            results.append({
                'cds_threshold': cds_t,
                'gene_threshold': gene_t,
                'kept_genes': len(keep),
                #'mean_cds_TE': mean_cds,
                #'mean_intron_TE': mean_intron
            })

            if args.scan_output_gff:
                write_gff(args.gff, f"{prefix}.cds{cds_t}.gene{gene_t}.gff3", keep)

    pd.DataFrame(results).to_csv(f"{prefix}.scan.summary.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()
