SweepR
======

Overview
--------
-

R implementation of extended haplotype heterozygosity with relative versions



## Rob's implementation of Extended Haplotype Homozygosity tools
# schae234@umn.edu 
# September 2012
#
# The MIT License (MIT)
# Copyright (c) 2012 Robert J Schaefer
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
#
### Raw Class Methods
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
Function
--------
#
#### print.sweep - prints out a sweep object in hudson MS style
## Arguments:
##  Sweep - the sweep object to print out
##  file  - an output file, NULL will print to console
##  append - flag to append to file or create new one
## Returns:
##  

Function
--------
  filter_by_index - Filters a Sweep object by index of individuals
## Arguments:
##  Sweep - sweep object to filter
##  index_list - list of indices to filter Sweep object by
## Returns:
##  A new Sweep object containing only indices specified in argument

Function
--------
 filter_by_individual - Filters a Sweep object by individuals
## Arguments:
##  Sweep - Sweep object to filter
##  individual list - list of individual ids to filter sweep object by
## Returns:
##  A new sweep object which only contains specified individuals

Function
--------
 permute_sweep - permutes so that a specified allele appears at a certain frequency.
##                    useful for when you have more 'affected' individuals than appear in the
##                    wild and you need to create a sampling which accurately represents what you
##                    would get in a random sampling of a population.
## Arguments:
##  Sweep - raw sweep object containing raw data
##  target_allele_freq : [0-.99] the frequency of the target allele, i.e. 'affected' allele
##  target_index : the SNP index of the target, 'affected' allele
##  target_allele: 1234 corresponds to ACGT. this is the target allele in which we are filtering
##                 the Sweep object by. i.e. we have an oversampling of affected=1 alleles and we
##                 are filtering so we have an expected population allele frequency.
##  N : the number of permuted individuals you would like in your new, permuted Sweep object
## Returns:
##  A new sweep object containing your targeted allele at your targeted index to have a target allele
##  frequency which is NORMALLY DISTRUBUTED AROUND YOUR TARGET ALLELE FREQUENCY. This permutation does
##  not aim to specifically hit your target frequency, but to rather create a set of permutations which 
##  are a statistical representation of what you would measure in a population.


Function
--------
  read.many - reads a .many file (as specified in the Sweep documentation) so the 
##                 individual genotype and SNP files can be read in as sweep objects
## Arguments:
##  many_file - the path of the Sweep many file                 
## Returns:
##  a data.frame containing geno and snp file paths
#### read.sweep - generates a sweep object from the raw sweep data read from genotype and SNP file
## Arguments:
##  geno_filename - path of the sweep genotype file
##  snp_filename - path of the sweep SNP filename
##  sep - maybe sometimes fields in the files are separeted by something other than tabs
## Returns:
##  a new sweep object containing the information in the input files
Function
--------
 distance_index -- returns the index of the SNP which is the target for a certain distance away
## Arguments:
##  Sweep - the sweep object you want to calculate distance index on
##  start_index - the index of the SNP which you want to start from
##  distance - the distance in base pairs from which you would like to find the index
##             of from the start_index. 
## Returns:
##  the index which is closes to the distance from the start snp specified
       
Function
--------
 closest_index - Returns the index which is closest to the specified distance from the core
##                    (which is defined as between index1 and index2).
## Arguments:
##  Sweep - the sweep object to calculate on
##  index1 - the first index of the core
##  index2 - the second index of the core
##  distance - the distance in which you would like the index from the core
## Returns:
##  the index which is closest to the distance specified 
	else{
Function
--------
 count - Return the number of times a haplotype appears in a core region
## Arguments:
##  Sweep - The sweep object to calculate from
##  haplotype - the haplotype which you would like to know how many times it appears
##              must be in list form, i.e. c(3,2,4,2)
##  core_begin - the start index of the core which you are searching
##  core_end   - the end index of the core you are searching
##  who  -  it true, return the indices of the haplotypes which match
## Returns:
##  either the number of times it found your core haplotype or the the indices
##  in which your haplotype appeared (depending on arg: who)
Function
--------
 freq - Return the frequency of a haplotype in a region (percentage)
## Arguments:
##  Sweep - as usual, the sweep object you are calculating on
##  haplotype - the haplotype which you want frequency of
##  core_being - index of core begin
##  core_end   - index of core end
## Returns:
##  frequency in percentage of core haplotype
Function
--------
 EHH - Calculates the Extended Haplotype Heterozygosity of a haplotype at a certain distance
## Arguments:
##  Sweep - sweep object which you are calculating on
##  haplotype - the core haplotype which you would like EHH information on
##  core_begin - the start index of the core 
##  core_end - the end index of the core
##  distance - the distance in base pairs which you would like EHH calculated from
## Returns:
##  returns the raw EHH of the haplotype at the specified distance


Function
--------
 REHH - Calculates the Relative Extended Haplotype Heterozygosity of a haplotype at a certain distance
## Arguments:
##  Sweep - sweep object which you are calculating on
##  haplotype - the core haplotype which you would like EHH information on
##  core_begin - the start index of the core 
##  core_end - the end index of the core
##  distance - the distance in base pairs which you would like EHH calculated from
## Returns:
##  returns the raw Relative EHH of the haplotype at the specified distance

Function
--------
 EHH_bar - Calculates the Extended Haplotype Heterozygosity BAR of a haplotype at a certain distance
## Arguments:
##  Sweep - sweep object which you are calculating on
##  haplotype - the core haplotype which you would like EHH information on
##  core_begin - the start index of the core 
##  core_end - the end index of the core
##  distance - the distance in base pairs which you would like EHH calculated from
## Returns:
##  returns the raw EHH BAR of the haplotype at the specified distance (better defined in Sabeti et al)
##  Basically, this is used in order to quantify a relative EHH and is the denominator in the term.
##  It calculates the EHH of all the other haplotypes at a core by summing the choose 2 combinations of all
##  other EXTENDED haplotypes at a core divided by the number of ways to choose 2 combinations of core haplotypes





Function
--------
 num_samples_with_haplotype -- returns the number of individual with a certain haplotype
## Argumnets:
##  genotypes - a matrix of genotypes, could be a core or simply just genotypes
##  haplotype - the haplotype you would like to find out how many individuals have
## Returns:
##  the number of individuals with the haplotype (as per the genotype matrix)
Function
--------
 homozygosity - calculates how unique a haplotype is among other genotypes in a region
## Arguments:
##  genotypes - a matrix of genotypes (could be SNPs) which contain all haplotypes in that region
##  haplotype - the haplotype of interest in which homozygosity will be calculated
## Returns:
##  the homozygosity of the haplotype within the specified genotypes
 
 
Function
--------
 haplotypes - returns the unique haplotypes from a frame of genotypes
## Arguments:
##  genotypes - a matrix of genotypes which contain haplotypes
## Returns:
##  the unique haplotypes within the genotype matrix
Function
--------
 find_core_at_core_h - find a core which contains homozygosity at a certain frequency
##                          We used this to find cores with the sme homozygosity as the empirical
##                          ones so the start of the core is always going to be index 1
## Arguments:
##  Sweep - the sweep object you want to calculat with (simulated data)
##  freq  - the target frequency to find core at
## Returns:
##   the end index with specified core homozygosity

	
Function
--------
 REHH_at_distance - Calculates REHH of all core haplotypes at a distance
## Arguments:
##  Sweep - usual sweep opject we're working on
##  core_start - index of core start
##  core_end   - index of core end
##  distance - distance from core in bp
## Returns:
##  a data frame containing all info for cores at specified distance
Function
--------
 core_homozygosity - calculates how homogenous alleles are at a certain core
## Arguments:
##  genotype - a matrix of genotypes (SNPs)
## Returns:
##  the homozygosity of the cores (as defined thoroughly in sabeti et al)

 
 
 
Function
--------
 plot_permuted_EHH_over_distance - generate a plot showing permuted EHH for all haplotypes over computed distances.
## Arguments:
##  REHH_Table - Table of values generated by function: EHH_at_all_distances_permute or another function which outputs similar data
##  title - title of the plot
##  color_order - list of haplotypes in order of appearance in the figure e.g. c(4332,2332,4132,4312)
## Returns:
##  a plot of EHH vs Distance

   
   
   
   
   
   
Function
--------
 EHH_at_all_distances_many - calculates EHH at all distances using a .many file
##                                 warning! this takes a while.. i.e. FOREVER! I dont use it
## Arguments:
##  many_file - path the sweep many file
##  core_start - start index of core of interest
##  core_end   -  end index of core of interest
## Returns:
##  a table containing all information produced by EHH_at_all_distances for each Sweep object in .many file

Function
--------
 EHH_at_all_distances_permute - given a 'raw' sweep object, first permute it according to target allele frequencies
## Arguments:
##  Raw - a 'raw' sweep object which we will permute
##  core_start - start index of core of interest
##  core_end  - end index of core of interest
##  sim - a flag indicating if this is simulated (one-sided data) and if we need to calculate the core end index (not fully implemented)
##  target_index - which index contains the allele which we need to permute frequencies on
##  target_allele - code of the allele we are permuting frequencies on
##  target_allele_freq - the final frequency we wish to see our permuted allele Sweep object to contain
##  N - the numer of individuals in the permutation
##  M - the number of permutations to perform
## Returns:
##  A huge table containing all the EHH information for all the distances possible in our data with N individuals and M permutations of the data
Function
--------
 EHH_at_all_distances - calculates EHH for each possible distance from a core
## Arguments:
##  Sweep - the object we are using to make calculations
##  core_start - the start index of the core of interest
##  core_end  -  the end index of the core of interest
##  min_core_count - filter out small core counts, which tend to be inflated for EHH
##  sim - not implemented
## Returns:
##  a huge table containing all EHH data for all cores at all distances
	
       
                        ## Takes a long time
					
                        ## Takes a long time
					
Function
--------
 REHH_Permute - Calculates a raw REHH and freq via permutation for each haplotype
## Arguments:
##  Raw - a 'raw' sweep object which we will permute
##  haplotype - depreceiated
##  core_start - start index of core of interest
##  core_end  - end index of core of interest
##  sim - a flag indicating if this is simulated (one-sided data) and if we need to calculate the core end index (not fully implemented)
##  target_index - which index contains the allele which we need to permute frequencies on
##  target_allele - code of the allele we are permuting frequencies on
##  target_allele_freq - the final frequency we wish to see our permuted allele Sweep object to contain
##  N - the numer of individuals in the permutation
##  M - the number of permutations to perform
## Returns:
##  a table containing REHH, frequency for each haplotye. Permuted over M runs.
       
       
Function
--------
 EHH_Permute - Calculates a raw EHH via permutation of the raw data for a specific haplotype
## Arguments:
##  Raw - a 'raw' sweep object which we will permute
##  haplotype - the haplotype to permute on
##  core_start - start index of core of interest
##  core_end  - end index of core of interest
##  sim - a flag indicating if this is simulated (one-sided data) and if we need to calculate the core end index (not fully implemented)
##  target_index - which index contains the allele which we need to permute frequencies on
##  target_allele - code of the allele we are permuting frequencies on
##  target_allele_freq - the final frequency we wish to see our permuted allele Sweep object to contain
##  N - the numer of individuals in the permutation
##  M - the number of permutations to perform
## Returns:
##  a list of EHHs over N inidivduals and M permutations
Function
--------
 REHH_many - calculates REHH of a certain core (target) combo at a certain distance
## Arguments:
##   many_file: this file contains the locations of the data and SNP files
##   core_start : this is the index in which the core starts, in the case of a simulation usually 1
##   target    : If this is not a simulation, the target in the index of the core_end. If it is a simulation, 
##               the target is the goal core homozygosity which will be calculated by the function
##   sim      : this flag tells the function whether or not this is a simulated file or not
## Returns:
##  the REHH_at_distance for each Sweep file in the .many file
   
       
       
       
   
Function
--------
 standard_error - calculates standard error for a vector
## Arguments: 
##  v - vector of numbers
## Returns:
##  99th percentile of standard error
Function
--------
 plot_emp_vs_sim_REHH - returns scatter plot showing simulated REHH data against permuted REHH data. Cloud plot!
## Arguments:
##  Raw - a 'raw' sweep object which we will permute
##  Sim - Simulated REHH data
##  core_start - start index of core of interest
##  core_end  - end index of core of interest
##  sim - a flag indicating if this is simulated (one-sided data) and if we need to calculate the core end index (not fully implemented)
##  target_index - which index contains the allele which we need to permute frequencies on
##  target_allele - code of the allele we are permuting frequencies on
##  target_allele_freq - the final frequency we wish to see our permuted allele Sweep object to contain
##  N - the numer of individuals in the permutation
##  M - the number of permutations to perform
##  png - flag indicating whether to write png to file or to X11
##  title - title of the saved file
## Return
##  Cloud plot

       
	
       
       
       
   
   
Function
--------
 r - reloads source file
### Arguments: none
### Returns: none


### 
# Phew! That's it!
###
