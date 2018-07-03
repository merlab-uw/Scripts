## sumGP.py

The purpose of the `sumGP.py` script is to automate the process of running Genepop, get an overview of the data, and parse output files into arrays that are computer readable. `sumGP.py` calls the Genepop program through the `genepop` package in R using the R script `runGP_to_sumGP.R`. It runs Hardy Weinberg exact tests, genic differentiation tests, estimates pairwise Fst, estimates global Fst per locus, and generates allele frequencies through the basic info function. Then, `sumGP.py` makes a folder in the same directory as the Genepop file and stores the following:
 - Genepop output files
 - computer readable CSVs of Genepop output files for allele frequencies & HWE pvalues
 - histograms of global Fis, Fst, and Fit as well as populations where loci are out of HWE
 - summary output files, including a pairwise table comparison of Fst and pvalues for genetic differentiation tests and a file with percent of populations out of HWE per locus

 For more information, call the script with the help flag, `-h`:

``$ python sumGP.py -h``

Natalie Lowell, June 28, 2018