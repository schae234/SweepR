# Rob's implementation of Sweep tools COOL!
# schae234@umn.edu 
# September 2012
#

# The MIT License (MIT)
#Copyright (c) <year> <copyright holders>
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


DEBUG <- FALSE

##################################################################
### Raw Class Methods
##################################################################

# Filters a Sweep object by individuals
filter_by_individual <- function(Sweep, individual_list){
	# remember we are dealing with diploids, so each ind has 2 chromos
	individual_indices <- which(Sweep$id %in% individual_list)
	rtn <- list(
		"filename" = paste(Sweep$filename," filtered ",sep=''),
		'posit' = Sweep$posit,
		'id' = Sweep$id[individual_indices],
		'chromo' = Sweep$chromo[individual_indices],
		'geno' = Sweep$geno[individual_indices,]
	)
	class(rtn) <- 'Sweep'
	return(rtn)
}

permute_sweep <- function(Sweep, target_allele_freq, target_index=21, target_allele=1, N=50){
	# Determine who is the target and who is the non-target	
	affected_chromos<-sample(which(Sweep$geno[,target_index] == target_allele))
	unaffected_chromos<-sample(which(Sweep$geno[,target_index] != target_allele))
	# grab N accounts of random individuals at the expected frequency
    num_affected <- length(which(runif(n=N)<=target_allele_freq))
    inds <- c(affected_chromos[1:num_affected],unaffected_chromos[1:(N-num_affected)])
  rtn <- list(
    "filename" = paste(Sys.time(),"_",Sweep$filename,"_","permutation",sep=""),
    "posit" = Sweep$posit,
    "id" = Sweep$id[inds],
    "chromo" = Sweep$chromo[inds],
    "geno" = Sweep$geno[inds,]
  )  
	class(rtn) <- "Sweep"
	return(rtn)
}



###################################################################
### Sweep Class Methods
###################################################################

read.many <- function(many_file){
  files <- read.table(many_file,col.names=c("geno","snp"),colClasses=c("character","character"),)
	dir <- dirname(many_file)
	files$geno <- paste(dir,"/",files$geno,sep="")
	files$snp <- paste(dir,"/",files$snp,sep="")
	return(files)
}

read.sweep <- function(geno_filename, snp_filename, sep="\t"){
  tryCatch(
    raw_file <- read.table(geno_filename,stringsAsFactors=FALSE),
    error = function(){stop(paste("Could not read genotype file: ",geno_filename))}
  )
  tryCatch(
    snp_file <- read.table(snp_filename,header=TRUE),
    error = function(){stop(paste("Could not read snp file: ",snp_filename))}
  )
  rtn <- list(
    "filename" = geno_filename,
    "posit" = snp_file[,3],
    "id" = raw_file$V1,
    "chromo" = raw_file$V2,
    "geno" = raw_file[,seq(3,ncol(raw_file))]
  )  
  class(rtn) <- "Sweep"
  return(rtn)
}

# Distance Index -- returns the Sweep index which is the target for a certain distance away
distance_index <- function(Sweep, start_snp, distance){
		target_posit <- Sweep$posit[start_snp] + distance 
		return(which.min(abs(Sweep$posit - target_posit)))
}
# Returns the index which is closest to the specified distance from the core
closest_index <- function(Sweep, index1, index2, distance){
	pos1 <- Sweep$posit[index1]
	pos2 <- Sweep$posit[index2]
	if(distance > 0){
		if(pos1 > pos2)
			return(index1)
		else
			return(index2)
	}
	else{ # distance is negative
		if(pos1 < pos2)
			return(index1)
		else
			return(index2)
	}
}

# Return the times a haplotype appears in a core region
count <- function(Sweep, haplotype, core_begin, core_end,who=FALSE){
	indices <- which(
				apply(Sweep$geno[,seq(core_begin,core_end)],1,function(row){all(haplotype == row)})
	)
    if(who==TRUE)
        return(indices)
    else
        return(length(indices))
}
# Return the frequency of a haplotype in a region (percentage)
freq <- function(Sweep, haplotype, core_begin, core_end){
		return(count(Sweep,haplotype,core_begin,core_end)/nrow(Sweep$geno))
}

# Calculates the Extended Haplotype Heterozygosity of a haplotype at a certain distance
EHH <- function(Sweep, haplotype, core_begin, core_end, distance){
	if(DEBUG)
		cat("haplotpye ",haplotype, "\ncore begin:",core_begin,"\ncore end: ", core_end,"\ndistance ",distance,"\n" )
	if(distance == 0) 
		return(1)
	# calculate the distance from the beginning and the end
	start_index <- closest_index(Sweep,core_begin,core_end,distance)
  target_index <- distance_index(Sweep,start_index,distance)
	#calculate the core homozygosity for each core haplotype
	EHH <- core_homozygosity(Sweep$geno[
				c(which(apply(Sweep$geno[,seq(core_begin,core_end)],1,function(row){all(haplotype == row)}))),
				seq(start_index,target_index)
			])
	return(EHH)
}
REHH <- function(Sweep, haplotype, core_begin, core_end, distance){
	# Calculate the REHH 
	REHH <- (EHH(Sweep,haplotype,core_begin,core_end,distance)
											/EHH_bar(Sweep,haplotype,core_begin,core_end,distance))
	return(REHH)
}

EHH_bar <- function(Sweep,haplotype,core_begin,core_end,distance){
	# We want the start index to be the closest to sthe target distance
	start_index <- closest_index(Sweep,core_begin,core_end,distance)
	target_index  <- distance_index(Sweep,start_index,distance)
	# which cores are not our one of interest?
	other_cores <- Sweep$geno[
						c(which(apply(Sweep$geno[,seq(core_begin,core_end)],1,function(row){!all(haplotype==row)}))),
						seq(core_begin,core_end)]
	# Extract the haplotypes from the repeated cores
	other_haps <- haplotypes(other_cores)
	# remove the rare haplotypes	
	other_hap_counts <- apply(other_haps,1,function(hap){return(count(Sweep,hap,core_begin,core_end))})
	other_haps <- other_haps[other_hap_counts > 2,]
	# calculate the numerator from Sabeti et al 2002
	numer <- sum(apply(other_haps,1,function(hap){
		extended_cores <- Sweep$geno[
					c(which(apply(Sweep$geno[,seq(core_begin,core_end)],1,function(row){all(hap==row)}))),
					seq(start_index,target_index)
		]
		extended_haps <- haplotypes(extended_cores)
		sum(apply(extended_haps,1,function(hap){
			choose(num_samples_with_haplotype(extended_cores,hap),2)
		}))
	}))
	denom <- sum(apply(other_haps,1,function(hap){
		choose(num_samples_with_haplotype(other_cores,hap),2)
	}))
	return(numer/denom)
}

num_samples_with_haplotype <- function(genotypes,haplotype){
	return(length(which(apply(genotypes,1,function(row){all(haplotype == row)}))))
}

homozygosity <- function(genotypes, haplotype){
  # convert to numeric
  genotypes <- data.matrix(genotypes)
  haplotype <- data.matrix(haplotype)
  # how often do we choose the haplotype at the locus?
	choose(num_samples_with_haplotye(genotypes,haplotype),2)/choose(nrow(genotypes),2)
}

haplotypes <- function(genotypes){
  return(data.matrix(unique(genotypes)))
}

find_core_at_core_h <- function(Sweep, freq){
	# holy shit! recursion?! that CS degree wasn't worthless after all!
	index_at_freq <- function(Sweep, index, freq){
		# base case (we'd rather this not take forever)
		if(index >= ncol(Sweep$geno))
			return(index)
		else if(core_homozygosity(Sweep$geno[,seq(1,index)]) < freq){
			return(
				c((index-1),index)[
						which.min(c(abs(core_homozygosity(Sweep$geno[,seq(1,index-1)])-freq),
							abs(core_homozygosity(Sweep$geno[,seq(1,index)])-freq))
						)
					]
			)
		}
		else
			return(index_at_freq(Sweep,index+1,freq))
	}	
	return(index_at_freq(Sweep,4,freq))
}

plot_EHH_vs_Freq <- function(Sweep,core_start,core_end,distance,plot=TRUE){
	df<-NULL
	apply(haplotypes(Sweep$geno[,seq(core_start,core_end)]),1,
		function(core_hap){
				rbind(df,data.frame(
						'EHH' = EHH(Sweep,core_hap,core_start,core_end,distance),
						'marker_index' = distance_index(Sweep,core_start,distance),
		  			'freq' = freq(Sweep,core_hap,core_start,core_end),
						'hap' = paste(core_hap,collapse=''),
						'core_start' = core_start,
						'core_end' = core_end)
				) ->> df
		}
	)
	if(plot){
		# filter out EHH under 3 percent
		plot(c(df$freq,df$freq),c(df$up_EHH,df$down_EHH),xlim=c(0,1),xlab="Frequency",ylab="EHH")
	}
	return(df[order(df$freq),])
}
plot_REHH_vs_Freq <- function(Sweep,core_start,core_end,distance,plot=TRUE){
	df <- NULL
	cores<-haplotypes(Sweep$geno[,seq(core_start,core_end)])
	core_counts <- apply(cores,1,function(c){return(count(Sweep,c,core_start,core_end))})
	df<-do.call("rbind",apply(cores[core_counts > 2,],1,
			function(core_hap){
				data.frame(
						'REHH' = REHH(Sweep,core_hap,core_start,core_end,distance),
						'EHH' = EHH(Sweep,core_hap,core_start,core_end,distance),
						'not_EHH' = EHH_bar(Sweep,core_hap,core_start,core_end,distance),
						'marker_index' = distance_index(Sweep,core_start,distance),
		  			'freq' = freq(Sweep,core_hap,core_start,core_end),
		  			'count' = count(Sweep,core_hap,core_start,core_end),
						'hap' = paste(core_hap,collapse=''),
						'core_start' = core_start,
						'core_end' = core_end,
						'core_h' = core_homozygosity(Sweep$geno[,seq(core_start,core_end)]),
						'filename' = Sweep$filename
				)
			}
	))
	if(plot){
  	plot(c(df$freq,df$freq),c(df$up_REHH,df$down_REHH),xlim=c(0,1),xlab="Frequency",ylab="REHH")
	}
  return(df[order(df$freq),])
}

core_homozygosity <- function(genotypes){
	# if we only have 1 row, return zero
	if(is.null(nrow(genotypes)) | (nrow(genotypes) < 2)){return(0)}
  # coherce into numeric form
  genotypes<-data.matrix(genotypes)
  # Extract the haplotypes
  haplotypes <- unique(genotypes)
  #possible chromosomes choices
  sum(
    apply(haplotypes,1,
      function(hap){
        choose(length(which(apply(genotypes,1,function(x){all(x == hap)}))),2)
      }  
    )  
  )/choose(nrow(genotypes),2)
}

average_by_distance_over_sampling <- function(emp_all){
	require(ggplot2)
	# split into individual hap groups
	by_hap <- split(emp_all,emp_all$hap)
	
	by_hap<-lapply(by_hap,function(hap_results){
		tmp   <- split(hap_results$EHH,hap_results$distance)
		dis   <- split(hap_results$distance,hap_results$distance)
		dis   <- sapply(dis,mean)
		means <- sapply(tmp, mean)
	    stdev <- sqrt(sapply(tmp, var))
  	    n     <- sapply(tmp,length)
  	    ciw   <- qt(0.99, n) * stdev / sqrt(n)
		hap   <- unique(hap_results$hap)
		return(data.frame('distance'=dis,'mean'=means,'stdev'=stdev,'n'=n,'ciw'=ciw,'hap'=hap))
	})
	a<-do.call("rbind",by_hap)
	qplot(distance, mean, data=a, color=factor(hap), geom=c("line","point"))+
		geom_errorbar(aes(ymin=mean-ciw,ymax=mean+ciw))+
		geom_line(aes(x=as.numeric(factor(distance))))
	#sapply(by_hap,function(x){mean(x$freq)})
}

# Permute a Raw data s
permute_at_all_distance <- function(Raw,core_start,core_end,sim=FALSE){
	
}

# warning! this takes a while..
EHH_at_all_distances_many <- function(many_file,core_start,core_end,sim=FALSE){
	many <- read.many(many_file)
	# calculate all the relevant distances 
	many_EHH <- apply(many,1,function(sweep_files){
		Sweep <- read.sweep(sweep_files[1],sweep_files[2])
		EHH_at_all_distances(Sweep,core_start,core_end,sim)	
	})
	do.call("rbind",many_EHH)
}
EHH_at_all_distances_permute <- function(Raw,core_start,core_end,sim=FALSE,target_index,target_allele,target_allele_freq,N=100,M=1000){
	many_EHH<-lapply(1:M,function(m){
			cat('one permutation ',m,"\n",sep='')
			Sweep <- permute_sweep(Raw,target_index=target_index,target_allele=target_allele,target_allele_freq=target_allele_freq,N=N)
			EHH_at_all_distances(Sweep,core_start,core_end,sim=FALSE)
		}		
	)
    return(do.call("rbind",many_EHH))
}
EHH_at_all_distances <- function(Sweep,core_start,core_end,sim=FALSE){
		cat(paste("On file: ",Sweep$filename,"\n",sep=""))
		distances <- c(
			(Sweep$posit[seq(1,core_start-1)] - Sweep$posit[core_start]),
			(Sweep$posit[seq(core_end+1,length(Sweep$posit))] - Sweep$posit[core_end])
		)
		# calculate for all distances
		cores<-haplotypes(Sweep$geno[,seq(core_start,core_end)])
		core_counts <- apply(cores,1,function(c){return(count(Sweep,c,core_start,core_end))})
		many_EHH<-apply(cores[core_counts > 3,],1,
			function(core_hap){
				lapply(distances,function(distance){
					data.frame(
						#'REHH' = REHH(Sweep,core_hap,core_start,core_end,distance),
						'EHH' = EHH(Sweep,core_hap,core_start,core_end,distance),
						#'not_EHH' = EHH_bar(Sweep,core_hap,core_start,core_end,distance),
						'marker_index' = distance_index(Sweep,core_start,distance),
		 				'core_freq' = freq(Sweep,core_hap,core_start,core_end),
		 				'core_count' = count(Sweep,core_hap,core_start,core_end),
						'hap' = paste(core_hap,collapse=''),
						'core_start' = core_start,
						'core_end' = core_end,
						'core_h' = core_homozygosity(Sweep$geno[,seq(core_start,core_end)]),
						'distance' = distance,
						'filename' = Sweep$filename)
				})
		})
	do.call("rbind",do.call("rbind",many_EHH))
}

REHH_many <- function(many_file,core_start,target,distance,sim=FALSE){
	many <- read.many(many_file)
	many_REHH<-NULL
	apply(many,1,function(sweep_files){
		Sweep <- read.sweep(sweep_files[1],sweep_files[2])
		cat(paste("On file: ",Sweep$filename,"\n",sep=""))
		if(sim==TRUE){
			core_end <- find_core_at_core_h(Sweep,target)
		}
		else
			core_end <- target
		rbind(many_REHH,plot_REHH_vs_Freq(Sweep,core_start,core_end,distance,plot=FALSE))->>many_REHH
		return(1)
	})
	return(many_REHH)
}
function(){
	BEL_emp_REHH <- REHH_many("Sweep/BEL_emp_0.10/Simulations.many",19,24,500000)
	BEL_emp2_REHH <- REHH_many("Sweep/BEL_emp_0.10/Simulations.many",19,24,-650000)
	BEL_100_REHH <- REHH_many("Sweep/BEL_100_0.10/Simulations.many",1,.22,500000,sim=TRUE)
	BEL_200_REHH <- REHH_many("Sweep/BEL_200_0.10/Simulations.many",1,.22,500000,sim=TRUE)
	QH_emp_REHH <- REHH_many("Sweep/QH_emp_0.10/Simulations.many",19,24,500000)
	QH_200_REHH <- REHH_many("Sweep/QH_200_0.10/Simulations.many",1,.36,500000,sim=TRUE)

	# plot the simulated points
	sim <- BEL_200_REHH
	plot(sim$freq,sim$REHH)
	emp <- BEL_emp2_REHH
	lapply(unique(emp$hap),
		function(h){
			points(mean(emp[emp$hap==h,"freq"]),
			mean(emp[emp$hap==h,"REHH"]),col="red")
		}
	)
}
