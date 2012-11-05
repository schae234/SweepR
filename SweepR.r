# Rob's implementation of Extended Haplotype Homozygosity tools
# schae234@umn.edu 
# September 2012
#

# The MIT License (MIT)
# Copyright (c) 2012 Robert J Schaefer
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


DEBUG <- FALSE

##################################################################
### Raw Class Methods
##################################################################

print.sweep <- function(Sweep,file=NULL,append=TRUE){
	# do we have a file name? or do we just want to print to the shell?
	if(!is.null(file)){
		sink(file,append=append)
	}
    # the ms tools only take in major and minor alleles, not this 1234 garbage
    Sweep$geno <- apply(Sweep$geno,2,function(col){
        alleles <- unique(col)
        if(length(alleles) < 2)
            return(rep(0,length(col)))
        if(length(which(col==alleles[[1]])) > length(which(col==alleles[[2]]))){
            col[col==alleles[[1]]] <- 0
            col[col==alleles[[2]]] <- 1
        }
        else{
            col[col==alleles[[2]]] <- 0
            col[col==alleles[[1]]] <- 1
        }
        return(col)
    })
    # print everything out in ms style
    cat("//\n")
    cat("segsites:",ncol(Sweep$geno),"\n")
    cat("positions:",paste(Sweep$posit,collapse=" "),"\n")		
    apply(Sweep$geno,1,function(row){
        cat(paste(row,collapse=""),"\n")
    })
	# sink back to the normal console
	if(!is.null(file)){
		sink()
	}
}

# Filters a Sweep object by individuals
filter_by_index <- function(Sweep, index_list){
	# remember we are dealing with diploids, so each ind has 2 chromos
	rtn <- list(
		"filename" = paste(Sweep$filename," index filtered ",sep=''),
		'posit' = Sweep$posit,
		'id' = Sweep$id[index_list],
		'chromo' = Sweep$chromo[index_list],
		'geno' = Sweep$geno[index_list,]
	)
	class(rtn) <- 'Sweep'
	return(rtn)
}

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

permute_sweep <- function(Sweep, target_allele_freq, target_index=21, target_allele=1, N=100){
	# Determine who is the target and who is the non-target	
	affected_chromos<-sample(which(Sweep$geno[,target_index] == target_allele))
	unaffected_chromos<-sample(which(Sweep$geno[,target_index] != target_allele))
	# grab N accounts of random individuals at the expected frequency
    num_affected <- length(which(runif(n=N)<=target_allele_freq))
    if(N-num_affected > length(unaffected_chromos)){
        cat("WARNING: permute_sweep: is trying to select more chromosomes that possible!")
    }
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

#distances <- function(Sweep,core_start,core_end){
#    c(Sweep$posit[1:core_start]-Sweep$posit[core_start],seq(0,(core_end-core_start-1)),Sweep$posit[core_end]
#}

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
  	    ciw   <- qnorm(0.99) * stdev / sqrt(n)
		hap   <- unique(hap_results$hap)
		return(data.frame('distance'=dis,'mean'=means,'stdev'=stdev,'n'=n,'ciw'=ciw,'haplotype'=hap))
	})
	a<-do.call("rbind",by_hap)
    
    # For QH's
    #a$haplotype2 <- factor(a$haplotype, c(4332,2332,4132,4312)) 
    # For BEL's
    #a$haplotype2 <- factor(a$haplotype,c(4312,4334,4332,4132,2332))
    
	qplot(distance, mean, data=a, color=haplotype2, geom=c("line","point"))+
		geom_errorbar(aes(ymin=mean-ciw,ymax=mean+ciw))+
		geom_line(aes(x=as.numeric(factor(distance))))+
        scale_x_continuous("Distance (Base Pairs)")+
        scale_y_continuous("EHH")+
        ggtitle("BEL Empirical EHH")
	#sapply(by_hap,function(x){mean(x$freq)})
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

REHH_permute <- function(Raw,haplotype,target_allele_freq,distance,target_allele=1,target_index=21,core_start=20,core_end=23,sim=FALSE,N=60,M=1000){
	REHHs <- lapply(1:M,function(m){
		Sweep <- permute_sweep(Raw,target_index=target_index,target_allele=target_allele,target_allele_freq=target_allele_freq,N=N)
		r<-REHH(Sweep,haplotype,core_start,core_end,distance)
        data.frame(
            'REHH' = r,
            'freq' = freq(Sweep, haplotype, core_start, core_end) 
        )
	})
	return(do.call("rbind",REHHs))
}
EHH_permute <- function(Raw,haplotype,target_allele_freq,distance,target_allele=1,target_index=21,core_start=20,core_end=23,sim=FALSE,N=60,M=1000){
	EHHs <- sapply(1:M,function(m){
		Sweep <- permute_sweep(Raw,target_index=target_index,target_allele=target_allele,target_allele_freq=target_allele_freq,N=N)
		EHH(Sweep,haplotype,core_start,core_end,distance)
	})
	return(EHHs)
}

REHH_many <- function(many_file,core_start,target,distance,sim=FALSE){
	many <- read.many(many_file)
	REHH_many<-apply(many,1,function(sweep_files){
		Sweep <- read.sweep(sweep_files[1],sweep_files[2])
		cat(paste("On file: ",Sweep$filename,"\n",sep=""))
		if(sim==TRUE){
			core_end <- find_core_at_core_h(Sweep,target)
		}
		else
			core_end <- target
		plot_REHH_vs_Freq(Sweep,core_start,core_end,distance,plot=FALSE)
	})
	return(do.call("rbind",REHH_many))
}

standard_error <- function(v){
    qnorm(.99)*sqrt(var(v))/sqrt(length(v))
}

plot_emp_vs_sim_REHH <- function(Raw,Sim,core_start,core_end,distance,target_allele_freq,haplotype,N=60,M=5000){
    require("ggplot2")
    cat("Calculating the Permuted Raw Values....")
    a<-REHH_permute(Raw=Raw,target_allele_freq=target_allele_freq,core_start=core_start,core_end=core_end,haplotype=haplotype,distance=distance,N=N,M=M)
    x<-data.frame('freq'=mean(a$freq),'REHH'=mean(a$REHH),'freqSD'=sd(a$freq),'REHHSE'=standard_error(a$REHH))
   
    cat("Plotting the data....") 
    ggplot(data=Sim,aes(x=freq,y=REHH))+
        geom_point()+
        geom_point(data=x,aes(x=freq,y=REHH),color='red')+
        geom_errorbar(data=x,aes(ymin=REHH-REHHSE,ymax=REHH+REHHSE),width=.01, color='red')+
        geom_errorbarh(data=x,aes(xmin=freq-freqSD,xmax=freq+freqSD),color='red',height=.25)+
        ggtitle("REHH Empirical versus Simulated")+
        scale_y_continuous("REHH")+
        scale_x_continuous("Core Frequency")
}

r <- function(){
    source("SweepR.r")
}

# Now witness the firepower of this fully ARMED and OPERATIONAL battle station!
main <- function(){
# Read in the REHH for empirical and Simulated Data
    source("SweepR.r")
    Raw <- read.sweep("Empirical/Raw_GYS1_SNP_Data_Final_18Jun12.txt","Empirical/Raw_GYS1_SNP_Data_Final_18Jun12_headers.txt")
    BELa <- read.table("Empirical/BEL_affected.txt")$V1
    BELu <- read.table("Empirical/BEL_unaffected.txt")$V1
    QHa <- read.table("Empirical/QH_affected.txt")$V1
    QHu <- read.table("Empirical/QH_unaffected.txt")$V1
    QHTa <- read.table("Empirical/QHT_affected.txt")$V1
    QHTu <- read.table("Empirical/QHT_unaffected.txt")$V1

    # Filter out to just the wanted individuals
    BEL_Raw<-filter_by_individual(Raw,c(BELa,BELu)) 
    QH_Raw<-filter_by_individual(Raw,c(QHa,QHu))
    QHT_Raw<-filter_by_individual(Raw,c(QHTa,QHTu))
    BEL_Permuted <- EHH_at_all_distances_permute(BEL_Raw,core_start=20,core_end=23,target_index=21,target_allele=1,target_allele_freq=.25,N=60,M=5000)
    QH_Permuted <- EHH_at_all_distances_permute(QH_Raw,core_start=20,core_end=23,target_index=21,target_allele=1,target_allele_freq=.05,N=90,M=5000)
    QHT_Permuted <- EHH_at_all_distances_permute(QHT_Raw,core_start=20,core_end=23,target_index=21,target_allele=1,target_allele_freq=.05,N=100,M=5000)
    # read in the simulation files
	BEL_100_REHH <- REHH_many("Sweep/BEL_100_0.10/Simulations.many",1,.24,275000,sim=TRUE)
	BEL_200_REHH <- REHH_many("Sweep/BEL_200_0.10/Simulations.many",1,.24,275000,sim=TRUE)
	QH_200_REHH <- REHH_many("Sweep/QH_200_0.10/Simulations.many",1,.38,450000,sim=TRUE)
	QH_100_REHH <- REHH_many("Sweep/QH_100_0.10/Simulations.many",1,.38,450000,sim=TRUE)
    # read in the stats files

    # generate the G only permuted file
    QH_G <- filter_by_index(QH_Raw,which(QH_Raw$geno[,21]==3))
    BEL_G <- filter_by_index(BEL_Raw,which(BEL_Raw$geno[,21]==3))
    lapply(1:1000,function(x){print.sweep(filter_by_index(QH_Raw,sample(1:nrow(QH_G$geno),replace=T)),file="QH_G_Permuted.txt")})
    lapply(1:1000,function(x){print.sweep(filter_by_index(BEL_Raw,sample(1:nrow(BEL_G$geno),replace=T)),file="BEL_G_Permuted.txt")})
    # Make the stats tables 
    QH_Sim_Stats <- do.call("rbind",lapply(c(100,200,400,800,1000),function(x){
         a<-read.table(paste('Sim/QH_',x,'_0.10_stats.txt',sep=''),header=T); a$msR<-x; return(a)  
    }))
    QH_Emp_Stats <- do.call("rbind",lapply(c("QH_G_Permuted_Stats.txt","QH_Permuted_Stats.txt"),function(x){
        a<-read.table(x,header=T); if(grepl('G',x)){a$msR <- "emp_g"}else{a$msR <- "emp"}; return(a);
    }))
    QH_Stats<-rbind(QH_Sim_Stats,QH_Emp_Stats)
    # plot the recombination rates
    ggplot(QH_Stats) +
     geom_boxplot(aes(x=factor(msR,c(100,200,400,800,1000,'emp','emp_g')),y=zns))+
     ggtitle("QH Simulated and Emperical Recombination Rate")+
     scale_x_discrete(name="Recombination Rate")
}
