'''
JCVI work shell generator

    -i, --input, prefix list, use space to seprate, filename should like prefix.cds.fa and prefix.gff
    -t, --type, synteny plot type, [ Macro | Micro ]
    -l, --layout, layout file name
    -m, --minspan, '--minspan' for jcvi, default is 30.
    -n, --minsize, '--min_size' for jcvi,default is 4.
    -s, --cscore,  '--cscore' for jcvi, default is .8, a C-score cutoff of .99 effectively filters the LAST hit to contain reciprocal best hit (RBH).
    -b, --iter, '--iter' for jcvi, default is 1. Set --iter to 2 will be useful to plot regions resulted from genome duplications.

'''
import os
import sys
import argparse

def prepare_micro_jcvi(name_list, layout, minspan, cscore, min_size):
    total_numbers = len(name_list)
    subject_name = name_list[0]
    query1_name = name_list[1]
    # Check file existence
    for i in name_list:
        cds_file = "{0}.cds.fa".format(i)
        gff_file = "{0}.gff".format(i)
        if not os.path.isfile(cds_file):
            print(cds_file, ": Not Exist!")
            sys.exit()
        if os.path.isfile(gff_file):
            os.system("sed -i 's/\t_\t/\t-\t/g' {0}".format(gff_file))
        else:
            print(gff_file, ": Not Exist!")
            sys.exit()

    # write command into shell
    with open("./Micro_jcvi.sh", "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("echo 'Formatting GFF files'\n")
        for j1 in name_list:
            f.write("python -m jcvi.formats.gff bed --type=mRNA --key=ID {0}.gff -o {0}.bed\n".format(j1))
        f.write("\necho 'Formatting fasta files'\n")
        for j2 in name_list:
            f.write("python -m jcvi.formats.fasta format {0}.cds.fa {0}.cds\n".format(j2))
        f.write("\necho 'Starting comparison...'\n")
        for j3 in name_list[1:]:
            f.write("python -m jcvi.compara.catalog ortholog {0} {1} --cscore={2} -n {3} --no_strip_names\n".format(subject_name, j3, cscore, min_size))
            f.write("python -m jcvi.compara.synteny screen --minspan={0} --simple {1}.{2}.anchors {1}.{2}.anchors.new\n".format(minspan, subject_name, j3))
        f.write("\necho 'Generating block files...'\n")
        for j4 in name_list[1:]:
            f.write("python -m jcvi.compara.synteny mcscan {0}.bed {0}.{1}.lifted.anchors --iter=1 -o {0}.{1}.i1.blocks\n".format(subject_name, j4))
        # Combining block files, Note if there are only two species, this step can be ignored.
        if total_numbers >= 3:
            f.write("\necho 'Combining blocks files...'\n")
            block_list = ["{0}.{1}.i1.blocks".format(subject_name, query) for query in name_list[1:]]
            f.write("python -m jcvi.formats.base join {0} --noheader | cut -f 1,{1} > {2}.blocks\n".format(" ".join(block_list), ",".join([str(col) for col in range(2, 2*total_numbers, 2)]),subject_name))
        else:
            print("There are only {0} species, no need to join blocks".format(total_numbers))
        # Combining all bed files
        f.write("\necho 'Combining bed files...'\n")
        bed_list = ["{0}.bed".format(bed) for bed in name_list]
        f.write("python -m jcvi.formats.bed merge {0} -o All.bed\n".format(" ".join(bed_list)))
        # Modify block file in case of DEBUG information
        if total_numbers >= 3:
            f.write("\necho '{0}.blocks is the original blocks file and may contain [DEBUG] information.'\n".format(subject_name))
            f.write("\ngrep -v \"^\[\" {0}.blocks | grep -v \"^\s\" > {0}.blocks.new\n".format(subject_name))
        # Plotting section
        f.write("\necho 'Plotting...'\n")
        if total_numbers >= 3:
            f.write("python -m jcvi.graphics.synteny {0}.blocks.new All.bed {1} --glyphcolor=orthogroup\n".format(subject_name, layout))
        else:
            f.write("python -m jcvi.graphics.synteny {0}.{1}.i1.blocks All.bed {2} --glyphcolor=orthogroup\n".format(subject_name, query1_name, layout))
        # Show some recommand parameters for synteny plot
        f.write("\n# Additional parameters you may intrest:\n")
        f.write("#\t--glyphstyle=arrow\tGene in arrow style.\n#\t--glyphcolor=orthogroup\tMatch the colors of the genes with respect to their homology.\n")
        f.write("#\t--genelabelsize 6\tAdd gene label.\n")
    return

def prepare_macro_jcvi(name_list, layout, minspan, cscore, min_size):
    subject_name = name_list[0]
    # Check file existence
    for i in name_list:
        cds_file = "{0}.cds.fa".format(i)
        gff_file = "{0}.gff".format(i)
        if not os.path.isfile(cds_file):
            print(cds_file, ": Not Exist!")
            sys.exit()
        if not os.path.isfile(gff_file):
            print(gff_file, ": Not Exist!")
            sys.exit()
    with open("./Macro_jcvi.sh", "w") as f:
        f.write("echo 'Formatting GFF files'\n")
        for j1 in name_list:
            f.write("python -m jcvi.formats.gff bed --type=mRNA --key=ID {0}.gff -o {0}.bed\n".format(j1))
        f.write("\necho 'Formatting fasta files'\n")
        for j2 in name_list:
            f.write("python -m jcvi.formats.fasta format {0}.cds.fa {0}.cds\n".format(j2))
        f.write("\necho 'Starting comparison and drawing dotplot...'\n")
        for j3 in name_list[1:]:
            f.write("python -m jcvi.compara.catalog ortholog {0} {1} --cscore={2} -n {3} --no_strip_names\n".format(subject_name, j3, cscore, min_size))
            f.write("python -m jcvi.compara.synteny screen --minspan={0} --simple {1}.{2}.anchors {1}.{2}.anchors.new\n".format(minspan, subject_name, j3))
            f.write("python -m jcvi.graphics.dotplot {0}.{1}.anchors --skipempty --style=white --theme=8\n".format(subject_name, j3))
            f.write("#python -m jcvi.compara.synteny depth --histogram {0}.{1}.anchors\n".format(subject_name, j3))
        f.write("\n# echo 'Plotting...'\n")
        f.write("# python -m jcvi.graphics.karyotype seqids {0}".format(layout))
    print(" Attention: The order of the genomes in [seqids] need to match the order in [layout]")


def main():
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, epilog="I think we can put our differences behind us, for science! --GLaDOS")
    parser.add_argument("-i", "--input", help="Prefix of cds and gff files, always put the subject in first place", dest="inp", required=True, nargs="*")
    parser.add_argument("-t", "--type", help="[ Macro | Micro ] Synteny plot type", dest="tp", required=True)
    parser.add_argument("-l", "--layout", help="layout file", dest="lo", required=True)
    parser.add_argument("-m", "--minspan", help="--minspan for jcvi", dest="msp", default=30)
    parser.add_argument("-s", "--cscore", help="--cscore for jcvi", dest="csc", default=.8)
    parser.add_argument("-n", "--min_size", help="--min_size for jcvi", dest="msz", default=4)
    args=parser.parse_args()

    if args.tp == "Micro":
        prepare_micro_jcvi(args.inp, args.lo, args.msp, args.csc, args.msz)
    elif args.tp == "Macro":
        prepare_macro_jcvi(args.inp, args.lo, args.msp, args.csc, args.msz)
    else:
        print("Wrong plot type, exit.")
        sys.exit()

if __name__ == "__main__":
    main()
