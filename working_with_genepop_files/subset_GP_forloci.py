################## SUBSET GENEPOP FILE WITH LIST OF LOCI ######################
# 20180711 Natalie Lowell
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
parser = argparse.ArgumentParser(description="Filter Genepop file to include only provided locus names.")
parser.add_argument("-i", "--infile", help="input Genepop file for filtering", type=str, required=True)
parser.add_argument("-f", "--format", help="input Genepop header format, answer 1 if locus names are each on a line and 2 if locus names are in one line separated by commas", type=str, required=True)
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
genepop_lines = genepop_str.split("\n")

# get locus names out of genepop file
all_locus_names = []
gp_chunks = []
gp_chunk_lines = []
for line in genepop_lines:
    if line.strip() != "pop" and line != "Pop":
        gp_chunk_lines.append(line)
    elif line.strip() == "pop" or line.strip() == "Pop":
        gp_chunks.append(gp_chunk_lines)
        gp_chunk_lines = []
gp_chunks.append(gp_chunk_lines)

if args.format == "1":
    for line in gp_chunks[0]:
        if line != "" and line[0] != "#":
            all_locus_names.append(line.strip())
if args.format == "2":
    for line in gp_chunks[0]:
        if line != "" and line[0] != "#":
            locus_names = line.split(",")
            for locus in locus_names:
                all_locus_names.append(locus.strip())        

# get indeces of loci to keep in genepop file
indeces_loci_to_keep = []
for locus in locus_names_to_keep:
    indeces_loci_to_keep.append(all_locus_names.index(locus))

### write filtered genepop file

# make new header with old header
if gp_chunks[0][0][0] == "#":
    old_header = gp_chunks[0][0]
    text_to_keep = gp_chunks[0][0][1:]
    new_header = "# Genepop made with subset_GP_forloci.py, old header: " + text_to_keep
    
# make new genepop file and add header
filtered_genepop = open(args.outfile, "w")
filtered_genepop.write(new_header+ "\n")

# add locus names, either one on each line or in header
if args.format == "1":
    for locus in locus_names_to_keep:
        filtered_genepop.write(locus + "\n")
elif args.format == "2":
    header = ""
    for locus in locus_names_to_keep:
        header += locus + ","
    header = header [:-1]
    filtered_genepop.write(header + "\n")
    
# add pop and genotype lines
for gp_chunk in gp_chunks[1:]:
    filtered_genepop.write("pop\n")
    for line in gp_chunk:
        if line != "":
            linelist = line.split(",")
            ind_name = linelist[0]
            genotypes = linelist[1].strip().split()
            filtered_genepop.write(ind_name + ", ")
            for index in indeces_loci_to_keep:
                filtered_genepop.write(" " + genotypes[index])
            filtered_genepop.write("\n")
filtered_genepop.close()

