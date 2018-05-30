################## SUBSET GENEPOP FILE WITH LIST OF LOCI ######################
# 20180529 Natalie Lowell
# PURPOSE: to subset a Genepop file with a list of loci to keep
# INPUTS: managed by argparse below, include:
# - genepop file to filter
# - list of loci to keep
# - output filepath
# - format of header
# OUTPUT: genepop file with only desired loci
###############################################################################

# organize parameter inputs with argparse
import argparse
parser = argparse.ArgumentParser(description="Filter Genepop file to include only provided locus names")
parser.add_argument("-i", "--infile", help="input Genepop file for filtering", type=str, required=True)
parser.add_argument("-f", "--format", help="input Genepop header format, answer 1 if locus names are each on a line and 2 if locus names are in one line separated by commas", type=float, required=True)
parser.add_argument("-l", "--loci", help="text file of loci to keep, with each locus name on its own line", type=str, required=True)
parser.add_argument("-o", "--outfile", help="filtered Genepop file", type=str, required=True)
args = parser.parse_args()

# get locus names to keep from input file
locus_names_file = open(args.loci, "r")
locus_names_to_keep = []
for line in locus_names_file:
    if line.strip() != "":
        locus_names_to_keep.append(line.strip())
locus_names_file.close()

# read in genepop file that needs filtering
genepop = open(args.infile, "r")
genepop_str = genepop.read()

# get locus names out of genepop file
all_locus_names = []
if args.format == 1:
    loci = genepop_str.split("Pop")[0].split("\n")
    for locus in loci:
        if locus != "":
            all_locus_names.append(locus.strip())
elif args.format == 2:
    loci = genepop_str.split("Pop")[0].split(",")
    for locus in loci:
        if locus != "":
            all_locus_names.append(locus.strip())
else:
    print("You did not specify an appropriate Genepop header format value. Use 1 for loci names each on their own line and 2 for loci names in one line separated by commas")
genepop.close()

# get indeces of loci to keep in genepop file
indeces_loci_to_keep = []
for locus in locus_names_to_keep:
    indeces_loci_to_keep.append(all_locus_names.index(locus))

### write filtered genepop file

# start with new line
filtered_genepop = open(args.outfile, "w")

# add locus names, either one on each line or in header
if args.format == 1:
    filtered_genepop.write("\n")
    for locus in locus_names_to_keep:
        filtered_genepop.write(locus + "\n")
else:
    filtered_genepop.write("\n")
    header = ""
    for locus in locus_names_to_keep:
        header += locus + ","
    header = header [:-1]
    filtered_genepop.write(header + "\n")
    

# add pop and genotype lines
gp_chunks = genepop_str.split("Pop")[1:]
for chunk in gp_chunks:
    filtered_genepop.write("Pop\n")
    lines = chunk.split("\n")
    for line in lines:
        if line != "":
            linelist = line.split(",")
            ind_name = linelist[0]
            genotypes = linelist[1].strip().split()
            filtered_genepop.write(ind_name + ", ")
            for index in indeces_loci_to_keep:
                filtered_genepop.write(" " + genotypes[index])
            filtered_genepop.write("\n")
filtered_genepop.close()
