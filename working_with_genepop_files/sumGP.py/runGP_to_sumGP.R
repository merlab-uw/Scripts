#############################################################

# 20180627 NL
# Compute basic pop gen functions in Genepop for parsing into 
# computer readible formats and a summary file

#############################################################

# get arguments from command line
args <- commandArgs(TRUE)

# make objects out of arguments
gp_file = args[1]
HWE_outfile = args[2]
Fst_outfile = args[3]
Fst_outfile_2 = args[4]
Diff_outfile = args[5]
BasicInfo_outfile = args[6]

# import library
library("genepop")

# make a file object out of genepop file
infile <- system.file(gp_file, package="genepop")
locinfile <- gp_file
check <- file.copy(infile,locinfile,overwrite = TRUE)

# run HWE, Fst, and Genetic Diff tests in genepop
test_HW(locinfile, which = 'Proba', outputFile = HWE_outfile)
Fst(locinfile, pairs = TRUE, outputFile = Fst_outfile)
Fst(locinfile, pairs = FALSE, outputFile = Fst_outfile_2)
test_diff(locinfile, genic = TRUE, pairs = TRUE, outputFile = Diff_outfile)
basic_info(locinfile, outputFile = BasicInfo_outfile)
