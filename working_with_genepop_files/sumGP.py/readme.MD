# sumGP.py

### Purpose

The purpose of the `sumGP.py` script is to automate the process of running Genepop, get an overview of the data, and parse output files into arrays that are computer readable. `sumGP.py` calls the Genepop program through the `genepop` package in R using the R script `runGP_to_sumGP.R`. It runs Hardy Weinberg exact tests, genic differentiation tests, estimates pairwise Fst, estimates global Fst, Fis, and Fit per locus, and generates allele frequencies through the basic info function. Then, `sumGP.py` makes a folder in the same directory as the Genepop file and stores the following:
 - Genepop output files
 - computer readable CSVs of Genepop output files for allele frequencies & HWE pvalues
 - histograms of global Fis, Fst, and Fit as well as populations where loci are out of HWE
 - summary output files, including a pairwise table comparison of Fst and pvalues for genetic differentiation tests and a file with percent of populations out of HWE per locus

### More about the output files

The script names the output file directory and output files based on the input file name. For example, if your input file is named `inputfile.gen`, then the output file directory will be named `inputfile_GPsumfiles`. And, 

- `inputfile_HWE_out.txt` will the the Genepop output file of the Hardy Weinberg exact test
- `inputfile_Fst_out.txt` and `inputfile_Fst_out.txt.MIG` wlil be the Genepop output files of the Fst estimation, between all groups
- `inputfile_Fst_global_out.txt` will be the Genepop output file of the global Fst estimation
- `inputfile_Diff_out.txt` will be the Genepop output file of the genic differentiation test 


### Dependencies

This script calls the `genepop` package in R, so you have to have installed both R and this package. Additionally, the R script `runGP_to_sumGP.R` should be saved in the same directory as the `sumGP.py` script. The parsing of Genepop output files depends on filenames that do not start with "-". The script uses the nonstandard modules `matplotlib` and `numpy`, so make sure they are installed.
### Help

For more information, call the script with the help flag, `-h`:

``$ python sumGP.py -h``

Natalie Lowell, June 28, 2018