################################# sumGP.py ####################################
# 20180713 Natalie Lowell
# PURPOSE: to provide clear, parsed versions of Genepop output files and
# produce basic summary plots and output files of common Genepop-related 
# tasks
# INPUTS: managed by argparse below, include:
# - absolute path to genepop file with .gen extension
# OUTPUT: directory with output files, including plots and clear, parsed
# versions of Genepop output files
# DEPENDENCIES: installed R, installed 'genepop' package in R, saved the
# runGP_to_sumGP.R script in same directory as this script & sample names
# do not start with a '-'
###############################################################################

### import all necessary modules
import sys
import os
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import argparse

### organize parameter inputs with argparse
parser = argparse.ArgumentParser(description="This script generates clear, parsed versions of Genepop output files and produces basic summary plots and output files of common Genepop-related tasks. Dependences: installed R, installed 'genepop' package in R, saved runGP_to_sumGP.R script in same directory as this script, and file names do not start with '-'.")
parser.add_argument("-i", "--infile", help="absolute path of Genepop file with .gen extension", type=str, required=True)
args = parser.parse_args()

###### ---Call R script to perform functions in Genepop-- ######

# get genepop filename (gp_filename) and directory (gp_dir) from absolute path to genepop file

gp_path_list = args.infile.split("/")
gp_dir= ""
for x in gp_path_list[:-1]:
    gp_dir += x + "/"
gp_filename = gp_path_list[-1]
if gp_dir == "":
    gp_dir = "./"

# get base working directory where current script is stored
basedir = os.path.dirname(os.path.realpath(sys.argv[0]))

# generate names of genepop output files, using file name of genepop without file extension
gp_filename_list = gp_filename.split(".")
corename = ""
for x in gp_filename_list[:-1]:
    corename += x + "."
corename = corename[:-1]
HWE_out_filename = corename + "_HWE_out.txt"
Fst_out_filename = corename + "_Fst_out.txt"
Fst_out_filename2 = corename + "_Fst_global_out.txt"
Diff_out_filename = corename + "_Diff_out.txt"
BasicInfo_out_filename = corename + "_BI_out.txt"

# make directory for genepop summary files, if doesn't already exist
os.chdir(gp_dir)
if not os.path.exists(corename + "_GPsumfiles"):
    os.makedirs(corename + "_GPsumfiles")
os.chdir(corename + "_GPsumfiles")

# call gp_sum.R to run genepop functions
call_R_string = "Rscript " + basedir + "/runGP_to_sumGP.R " 
call_R_string += "../" + gp_filename
call_R_string += " " + HWE_out_filename
call_R_string += " " + Fst_out_filename
call_R_string += " " + Fst_out_filename2
call_R_string += " " + Diff_out_filename
call_R_string += " " + BasicInfo_out_filename
sp.call([call_R_string],shell=True)

### ----------------------- Parse HWE output file --------------------------###

# get relevant portion of HWE output file for parsing
HWE_file = open(HWE_out_filename, "r")
HWE_file_text = HWE_file.read()
HWE_file_text_noheader = HWE_file_text.split("Results by")[1]
HWE_file_lines = HWE_file_text_noheader.split("\n")

# iterate through lines and feed dictionary, with key locus and value a list of Pvalues for each population
pval_dict = {}
popnames = [] 
header_count = 2 
last_locus = "" 
ordered_loci = []
for line in HWE_file_lines:
    if line.startswith("Locus ") and header_count == 2:
        header_count = 0
        linelist = line.strip().split()
        locus_name = linelist[1].replace('"','')
        last_locus = locus_name
        pval_dict[locus_name] = []
        ordered_loci.append(locus_name)
    elif header_count == 2 and len(line.strip().split()) > 1 and line.startswith("All") == False and line.startswith(" ") == False:
        linelist = line.strip().split()
        popname = linelist[0]
        if popname not in popnames:
            popnames.append(popname)
        if linelist[1] == "-":
            pval_dict[locus_name].append("NA")
        else:
            pval_dict[locus_name].append(float(linelist[1]))
    elif line.startswith("-"):
        header_count += 1
HWE_file.close()

# count how often population is out of HWE for a locus & store in list
list_of_pval_lists = []
for locus in ordered_loci:
    this_locus_list = []
    for pval in pval_dict[locus]:
        this_locus_list.append(pval)
    list_of_pval_lists.append(this_locus_list)
sig_counts = []
for locus_list in list_of_pval_lists:
    count = 0
    for x in locus_list:
        if x != "NA" and x < 0.05:
            count += 1
    sig_counts.append(count) 
num_pops = len(this_locus_list)
num_out_in_x_pops_list = []
for x in range(num_pops):
    num_out_in_x_pops = sig_counts.count(x + 1)
    num_out_in_x_pops_list.append(num_out_in_x_pops)   

# write pvals array file that can be more easily used to analyze HWE data
pvals_array_file = open(corename + "_HWE_pval_array.csv", "w")
pvals_array_file_header = "locus_name" + ","
for pop in range(num_pops):
    pvals_array_file_header += "pop_" + str(pop + 1) + ","
pvals_array_file_header = pvals_array_file_header[:-1]
pvals_array_file.write(pvals_array_file_header + "\n")
for locus in ordered_loci:
    locus_vals = pval_dict[locus]
    locus_line = locus + ","
    for val in locus_vals:
        locus_line += str(val) + ","
    locus_line = locus_line[:-1]
    pvals_array_file.write(locus_line + "\n")
pvals_array_file.close()

# get xlabels into list for plot
num_pops_xlabel_list = []
for num in range(num_pops):
    num_pops_xlabel = str(num + 1)
    num_pops_xlabel_list.append(num_pops_xlabel)

# Make bar plot with counts of loci out of HWE by number of pops
fig = plt.figure()
plt.bar(num_pops_xlabel_list,num_out_in_x_pops_list, width = .25)
plt.ylabel("Number of loci")
plt.xlabel("Out of HWE in X populations")
plt.title("Number of loci out of HWE (p < 0.05) in each number of populations")
plt.savefig(corename + "_HWE_hist.png", dpi=fig.dpi)
    
### ----------------- Parse Fst and Diff output files --------------------- ###

# make a key out of population index and population names provided in genepop output file
Fst_out = open(Fst_out_filename,"r")
Fst_out_text = Fst_out.read()
ugly_key_list = Fst_out_text.split("Indices for populations:")[1].split("Estimates for each locus:")[0].split("\n")
clean_key_list = []
for line in ugly_key_list:
    if line != "" and line[0] != "-":
        clean_key_list.append(line)
pop_key_names = {}
pop_key_nums = {}
for pop_row in clean_key_list:
    pop_row_list = pop_row.split()
    pop_key_names[pop_row_list[1]] = pop_row_list[0]
    pop_key_nums[pop_row_list[0]] = pop_row_list[1]    
    
# build a dictionary that holds Fst estimates for each pair of populations
Fst_pairwise_table_lines = Fst_out_text.split("Estimates for all loci (diploid):")[1].split("\n")[3:-4]
Fst_dict = {}
for line in Fst_pairwise_table_lines:
    Fst_pairwise_linelist = line.strip().split()
    popy = Fst_pairwise_linelist[0]
    for i in range(1, len(Fst_pairwise_linelist)):
        popx = i
        pval = Fst_pairwise_linelist[popx]
        combo = "&".join([popy,str(popx)]) 
        reverse_combo = "&".join([str(popx),popy]) 
        if combo not in Fst_dict:
            Fst_dict[combo] = pval
        if reverse_combo not in Fst_dict:
            Fst_dict[reverse_combo] = pval
Fst_out.close()

# build a dictionary that holds pvalues for genic differentiation tests for each pair of populations
diff_file = open(Diff_out_filename, "r")
diff_file_text = diff_file.read()
diff_file_pairs = diff_file_text.split("P-value for each population pair across all loci\n(Fisher's method)")[1]
diff_file_pairs_lines = diff_file_pairs.split("\n")
pairwise_pvals = {}
for line in diff_file_pairs_lines[4:-2]:
    if line != "" and line[0] != "-" and line[0] != " ":
        pairs_linelist = line.strip().split()
        pop_pos1 = pop_key_names[pairs_linelist[0][0:8]]
        pop_pos2 = pop_key_names[pairs_linelist[2][0:8]]
        combo = "&".join([pop_pos1, pop_pos2])
        reverse_combo = "&".join([pop_pos2,pop_pos1])
        if combo not in pairwise_pvals:
            pairwise_pvals[combo] = pairs_linelist[-1]
        if reverse_combo not in pairwise_pvals:
            pairwise_pvals[reverse_combo] = pairs_linelist[-1]
diff_file.close()

# write a csv that is a pairwise comparison table, with Fst estimates in the lower left & pvalues from the genic differentiation tests in the upper right
Fst_diffpval_array = open(corename + "_Fst_gendiff_array.csv","w")
Fst_diffpval_array_header = ""
for popnum in range(num_pops):
    Fst_diffpval_array_header += "," + pop_key_nums[str(popnum + 1)]
Fst_diffpval_array.write(Fst_diffpval_array_header + "\n")
for popy in range(num_pops):
    this_line = pop_key_nums[str(popy + 1)]
    for popx in range(num_pops):
        if popx == popy:
            new_cell = "NA"
        elif int(popx) < int(popy):
            new_cell = Fst_dict["&".join([str(popx + 1), str(popy + 1)])]
        elif int(popx) > int(popy):
            new_cell = pairwise_pvals["&".join([str(popx + 1), str(popy + 1)])]
            if new_cell == "sign.":
                new_cell = "Highly sign."
        this_line += "," + str(new_cell)
    Fst_diffpval_array.write(this_line + "\n")
Fst_diffpval_array.close()

### ----------------- Parse Basic Info output file ------------------------- ###

# get relevant portion of basic info output file for parsing
BI_file = open(BasicInfo_out_filename, "r")
BI_file_text = BI_file.read()
BI_file_text_tables = BI_file_text.split("Tables of allelic frequencies for each locus:\n\n ")[1]
BI_file_locus_tables = BI_file_text_tables.split("Locus: ")

# get allele frequencies from genepop output file
af_dict = {}
for locus_table in BI_file_locus_tables[1:]:
    locus_table_lines = locus_table.split("\n")
    table_locus = locus_table_lines[0].strip()
    af_dict[table_locus] = []
    for line in locus_table_lines[1:]:
        if line != "" and line[0] != " " and line[0] != "-" and line[0:14] != "Normal ending.":
            locus_table_linelist = line.strip().split()
            locus_pop_afs = []
            for value in locus_table_linelist[1:-1]:
                if value != "-":
                    locus_pop_afs.append(value)
            af_dict[locus_table_lines[0].strip()].append(locus_pop_afs)

# get maximum number of alleles per locus for iterating through in next code block
num_alleles_per_locus_list = []
for locus in ordered_loci:
    num_alleles_per_locus = len(af_dict[locus][0])
    num_alleles_per_locus_list.append(num_alleles_per_locus)
max_alleles_per_locus = max(num_alleles_per_locus_list)

# write clean array file of allele frequencies per population 
afs_array_file = open(corename + "_AFs_array.csv", "w")
afs_array_file_header = "pop" + "," + "locus_name"
for x in range(max_alleles_per_locus):
    afs_array_file_header += "," + "allele_" + str(x + 1)
afs_array_file.write(afs_array_file_header + "\n")
for locus in ordered_loci:
    for pop in range(num_pops):
        this_row = "pop_" + str(pop + 1) + "," + str(locus)
        locus_afs = af_dict[locus]
        this_pop_afs = locus_afs[pop]
        afs_string = ""
        for af in this_pop_afs:
            afs_string += "," + af
        if len(this_pop_afs) < max_alleles_per_locus:
            diff = max_alleles_per_locus - len(this_pop_afs)
            for i in range(diff):
                afs_string += ","+ "NA"
        this_row += afs_string
        afs_array_file.write(this_row + "\n")
afs_array_file.close()

# make dictionary with locus as key and value is list with number of pops locus observed in,
# number of pops locus polymorphic in, number of pops locus not polymorphic, and proportion
# of pops with observed locus polymoprhic
polymorph_d = {}
num_pops_polymorph_counts = []
for locus in ordered_loci:
    pop_polymorphic_count = 0
    pop_npolymorphic_count = 0
    locus_afs = af_dict[locus]
    for pop in range(num_pops):
        this_pop_afs = locus_afs[pop]
        if len(this_pop_afs) != 0:
            first_af = this_pop_afs[0]
            if first_af != "-":
                if float(first_af) < 1 and float(first_af) > 0:
                    pop_polymorphic_count += 1
                elif float(first_af) == 0 or float(first_af) == 1:
                    pop_npolymorphic_count += 1
    num_pops_polymorph_counts.append(pop_polymorphic_count)
    pop_observed_count = pop_polymorphic_count + pop_npolymorphic_count
    polymorph_d[locus] = [pop_observed_count, pop_polymorphic_count, pop_npolymorphic_count]

# Make array with locus name and percent of populations polymorphic
polymorph_array = open(corename + "_polymorphism_array.csv", "w")
polymorph_array.write("locus_name,npops_observed,npops_polymorphic,npops_monomorphic,proportion_observed_pops_polymorphic" + "\n")
for locus in ordered_loci:
    obs = str(polymorph_d[locus][0])
    poly = str(polymorph_d[locus][1])
    mono = str(polymorph_d[locus][2])
    prop = str(float(polymorph_d[locus][1])/float(polymorph_d[locus][0]))
    polymorph_array.write(locus + "," + obs + "," + poly + "," + mono + "," + prop + "\n")
polymorph_array.close()

# make a histogram of population counts for which a locus is polymorphic
fig = plt.figure()
plt.hist(num_pops_polymorph_counts, bins = np.arange(0,num_pops+2,.5)-.25)
plt.title("Number of populations for which a locus is polymorphic")
plt.ylabel("Frequency")
plt.xlabel("Number of populations")
plt.xticks(range(0,num_pops + 1))
plt.savefig(corename + "_polymorph_counts_hist.png", dpi=fig.dpi)

### ----------------- Parse Fst global output file ------------------------- ###

# get relevant portion of Fst global output file for parsing
globalFst = open(Fst_out_filename2,"r")
globalFst_text = globalFst.read()
globalFst_table_lines = globalFst_text.split("Multilocus estimates for diploid data")[1].split("\n")

# store global Fis, Fst, and Fit data in lists
global_Fis = []
global_Fst = []
global_Fit = []
for line in globalFst_table_lines[3:-4]:
    linelist = line.strip().split()
    if linelist[1] != "-":
        global_Fis.append(float(linelist[1]))
    if linelist[2] != "-":
        global_Fst.append(float(linelist[2]))
    if linelist[3] != "-":
        global_Fit.append(float(linelist[3]))

# plot histograms of global Fis, Fst, and Fit
fig = plt.figure()
plt.hist(global_Fis, bins = np.arange(-1, 1, .05)-0.025)
plt.title("Global Fis per locus")
plt.ylabel("Frequency")
plt.xlabel("Fis")
plt.savefig(corename + "_globalFis_hist.png", dpi=fig.dpi)    
fig = plt.figure()
plt.hist(global_Fst, bins = np.arange(-1, 1, .05)-0.025)
plt.title("Global Fst per locus")
plt.ylabel("Frequency")
plt.xlabel("Fst")
plt.savefig(corename + "_globalFst_hist.png", dpi=fig.dpi)
fig = plt.figure()
plt.hist(global_Fit, bins = np.arange(-1, 1, .05)-0.025)
plt.title("Global Fit per locus")
plt.ylabel("Frequency")
plt.xlabel("Fit")
plt.savefig(corename + "_globalFit_hist.png", dpi=fig.dpi)

####################################-END-######################################