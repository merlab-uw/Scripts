# sumGP.py

### Purpose

The purpose of the `sumGP.py` script is to automate the process of running Genepop, get an overview of the data, and parse output files into arrays that are more user friendly. `sumGP.py` calls the Genepop program through the `genepop` package in R using the R script `runGP_to_sumGP.R`. It runs Hardy Weinberg exact tests, genic differentiation tests, estimates pairwise Fst, estimates global Fst, Fis, and Fit per locus, and generates allele frequencies through the basic info function. Then, `sumGP.py` makes a folder in the same directory as the Genepop file and stores the following:
 - Genepop output files
 - computer readable CSVs of Genepop output files for allele frequencies & HWE pvalues
 - histograms of global Fis, Fst, and Fit as well as populations where loci are out of HWE
 - summary output files, including a pairwise table comparison of Fst and pvalues for genetic differentiation tests and a file with percent of populations out of HWE per locus

### More about the output files

The script names the output file directory and output files based on the input file name. For example, if your input file is named `inputfile.gen`, then the output file directory will be named `inputfile_GPsumfiles`. To learn more about the Genepop output files, check out the [Genepop manual](http://kimura.univ-montp2.fr/~rousset/Genepop4.7.pdf). Below are the names for the Genepop output files.

- `inputfile_HWE_out.txt` will the the Genepop output file of the Hardy Weinberg exact test
- `inputfile_Fst_out.txt` and `inputfile_Fst_out.txt.MIG` wlil be the Genepop output files of the Fst estimation, between all groups
- `inputfile_Fst_global_out.txt` will be the Genepop output file of the global Fst estimation
- `inputfile_Diff_out.txt` will be the Genepop output file of the genic differentiation test
- `inputfile_BI_out.txt` will be the Genepop output file for the basic info function

In addition to the Genepop output files, `sumGP.py` produces additional plots and CSV files, listed below with descriptions.

 - `inputfile_AFs_array.csv` is an array with allele frequencies of every allele, of every locus, within every population. The script works as well for multiallelic loci. NAs exist where an allele was not observed.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/afs_for_sumgpreadme2.png?raw=true)

- `inputfile_Fst_gendiff_array.csv` is a pairwise population comparison table where values in the upper right are p- values from the genic differentiation tests, and values in the lower left are Fst estimates. NAs will the cells where populations are compared to themselves.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/diff_fst_for_readme.png?raw=true)

- `intputfile_pval_array.csv` is an array made from parsing the Hardy Weinberg exact test output into a more computer readable CSV file. It includes the p-value from the Hardy Weinberg exact test by population.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/pval_for_readme.png?raw=true)

- `inputfile_polymorphism_array.csv` is an array that contains information on the number of populations a locus was observed, the number of populations in which a locus was polymorphic, the number of populations in which a locus was monomorphic, and the proportion of the number of populations in which a locus was polymorphic compared to the total number of populations in which the locus was observed.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/polymorph_array_for_readme.png?raw=true)

- `inputfile_HWE_hist.png` is an image file of a histogram of the number of populations for which a locus was out of Hardy Weinberg.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/hwe_hist_for_readme.png?raw=true)

- `inputfile_polymorph_counts_hist.png` is an image file of a histogram of counts of populations for which a locus was polymorphic.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/polymorph_hist_for_readme.png?raw=true)

- `inputfile_globalFis_hist.png` is an image file of a histogram of global Fis per locus.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/global_fis_for_readme.png?raw=true)

- `inputfile_globalFit_hist.png` is an image file of a histogram of global Fit per locus.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/global_fit_for_readme.png?raw=true)

- `inputfile_globalFst_hist.png` is an image file of a histogram of global Fst per locus.

![img](https://github.com/nclowell/SeaCukes/blob/master/Imgs_for_Notebooks/global_fst_for_readme.png?raw=true)

### Dependencies


This script calls the `genepop` package in R, so you have to have installed both R and this package. Additionally, the R script `runGP_to_sumGP.R` should be saved in the same directory as the `sumGP.py` script. The parsing of Genepop output files depends on filenames that do not start with "-". The script uses the nonstandard modules `matplotlib` and `numpy`, so make sure they are installed.
### Help

For help running the script, call the script with the help flag, `-h`:

``$ python sumGP.py -h``

Natalie Lowell, July 13, 2018