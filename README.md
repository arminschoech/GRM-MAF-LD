# GRM-MAF-LD
Tool for fast MAF and LD dependent GRM calculation


## Introduction
This is a software tool to help estimating how the sizes of SNP effects on a given polygenic human trait or disease depend on their minor allele frequency (MAF) in the population and their level of linkage disequilibrium (LD). Details about the method can be found in the corresponding bioRxiv paper (https://www.biorxiv.org/content/early/2017/09/13/188086). In short, we use a variance component estimation method, calculating variance components with different MAF and LD weightings of SNPs and find the one with the highest profile likelihood. For currenlty available sample sizes, calculating the genetic variance component, often called genetic relatedness matrix (GRM), is computationally intensive. This software tool allows fast GRM calculation with variable frequency and LD weighting, from imputed genotype data in binary GEN file format (BGEN). The Genome-wide Complex Trait Analysis (GCTA) toolkit can then be used on the resulting output GRM files ("out.grm.bin", "out.grm.N.bin", "out.grm.id") to calculate the likelihood for each MAF and LD weighting. 


## Requirements and installation
The following instructions pertain to Linux environments. Library requirements: Gnu Science Library, Intel Math Kernel Library, zlib, C standard library, POSIX threads, all of which are freely available. Before compiling the code, the locations of these libraries have to be updated in the makefile to the appropriate location. Then use ‘make’ do compile the executable ‘grm_maf_ld’. 


## Command line arguments
Required arguments:
   
   --alpha: set the alpha parameter that determines the weighting of low vs. high frequency SNPs in the GRM (see above bioRxiv link for details)

   --bgen_file: path to compressed genotype data file in BGEN format

   --freq_file: list of SNPs to be included in plink '.frq' file format (see 'https://www.cog-genomics.org/plink2/formats#frq'), stating their frequency in the imputation panel

   --id_file: IDs of all individuals to be analyzed in a single column text file

   --bgen_pos_file: position of individuals in the bgen_file in a single column text file; start counting at 0, i.e. if the individual with ID 5th in id_file is 10th in bgen_file, the 5th row of bgen_pos_file should be 9 

Optional arguments: 
   
   --out: set name of output files; default is "out"
   
   --snp_partition: only process input data from a subset of SNPs; flag should be followed by two numbers, e.g. "2 6" means that only the second of 6 equally sized subsets of SNP data gets processed
   
   --individual_partition: same as above, but now only data from a subset of individuals gets processed
   
   --tau: set tau parameter that determines the LD based weighting of SNPs in the GRM; default is 0 (see above bioRxiv link for details)
   
   --sampled_genotypes: if this flag is used, genotype values are sampled randomly from the respective probabilities in the BGEN files
   
   --sim_pheno: if this flag is used, instead of calculating the GRM, continuous trait values using the specified frequency and LD weightings are simulated
   
   --seed: sets random seed
   
   --fraction_causal: sets random fraction of SNPs being causal in trait simulations (not recommended to use for GRM calculation)
   
   --h2: sets the heritability for trait simulations
 

## Computation time and memory needs
Time for calculating a GRM for N individuals and M SNPs scales with N^2*M. Memory requirement is about 2*N^2 bytes, i.e. the size of the output GRM. 


## Contact
Author: Armin Schoech. Please email arminschoech@g.harvard.edu for comments and questions. 

## License
GRM-MAF-LD is free software under the GNU General Public License v3.0 (GPLv3). 
