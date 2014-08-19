
#' check what type of OS and then figure out the path to the gsi_sim executable appropriately
#' 
#' \code{gsi_simBinaryPath} checks to see if you are running in windows or Mac and then 
#' it looks for the executable file for gsi_sim in the appropriate place in the package
#' directory.  If it doesn't find it, it notes the error and stops.  If the platform
#' is not Mac or Windows, it also throws an error. 
#' @export
gsi_simBinaryPath <- function() {
  if(.Platform$OS.type=="windows") {
     gspath <- file.path(system.file(package="gpiper"), "bin/gsi_sim.exe")
  } 
  else {
  	if(.Platform$OS.type=="unix") {
    	if(Sys.info()["sysname"]=="Darwin") {
      	gspath <- file.path(system.file(package="gpiper"), "bin/gsi_sim")
    	} 
    	else {
      	stop(paste("This appears to be a non-Mac Unix architecture.  Not supported by gpiper currently. Sys.info()[\"sysname\"] =", Sys.info()["sysname"]))
    	}
  	}
  }
  if(!file.exists(gspath)) stop(paste("gsi_sim executable should be installed at", gspath,"but does not seem to be there"))
  return(gspath)
}



#' Simple unix-like interface to run gsi_sim
#'
#' This is just a wrapper.  It runs gsi_sim with the arguments in \code{arg.string}
#' in the current working directory, and it sends the
#' standard output to \code{stdout.to} and the stderr to \code{stderr.to}
#'
#' @param arg.string a text string which includes all the options to send to gsi_sim
#' @param stdout.to path to redirect gsi_sim's standard output (stdout) to
#' @param stderr.to path to redicted gsi_sim's standard error to
#' @param ...  extra variables that get passed to system2
#' @export
gsi_Run_gsi_sim <- function(arg.string, stdout.to="GSI_SIM_Dumpola.txt", stderr.to="GSI_SIM_Stderr.txt", ...) {
	
	call.str <- paste(gsi_simBinaryPath(), paste(arg.string, collapse=" "))
	if(file.exists(stdout.to)) file.remove(stdout.to)
	if(file.exists(stderr.to)) file.remove(stderr.to)
	print(paste("Sending the following command to system2:", call.str))
	print(paste("Launching at:", date()))
	system2(command=gsi_simBinaryPath(),
			args=arg.string,
			stdout=stdout.to,
			stderr=stderr.to,
			...)
	print(paste("Done at", date()))
	if(length(readLines(stderr.to))>0) {		
		stop(paste(c("Bad News! gsi_sim produced errors as follows:\n", readLines("GSI_SIM_Stderr.txt")), collapse="\n"))
	}		
}



#' convert a \code{gPipe.data.frame} into a gsi_sim input file
#' 
#' gPdf2gsi.sim takes as input a data frame in \code{gPipe.data.frame} format
#' converts it into a gsi_sim input file. If you supply a list of pop-names 
#' (preferably as a factor with levels ordered as you want the pops to appear) to the 
#' \code{pop.ize.them} variable, then it splits them into pops, ordered as in the factor
#' 
#' @param d  the data frame.  See documentation for \code{gPipe.data.frame} for format
#' @param pop.ize.them vector that should be of length \code{nrow(d)} that gives the populations in which the individuals belong.  It this is a factor, then populations are placed in the gsi_sim output file according to the order in the \code{levels} attribute of the factor.
#' @param outfile path to the output file 
#' @export
gPdf2gsi.sim <- function(d, pop.ize.them=NULL, outfile="gsi_sim_file.txt") {
  write(dim(d)/c(1,2), outfile)
  write(cbind(names(d)[seq(1, ncol(d), 2)]), outfile, append=T)
  if(is.null(pop.ize.them)) {
    write("POP  Mixture", outfile, append=T)
    write.table(d, outfile, append=T, quote=F, row.names=T, col.names=F)
  }
  else {
    idx.list <- split(1:nrow(d), pop.ize.them)
    catchit <- lapply(names(idx.list), function(x) {
      write(paste("POP", x), outfile, append=T)
      write.table(d[idx.list[[x]],] , outfile, append=T, quote=F, row.names=T, col.names=F)
    }
    )
  }	
}


#' create a reporting units file for the gsi_sim program
#' 
#' Given a character vector \code{pops} of unique population names that are to be in a gsi_sim baseline
#' file and a factor of the same length, \code{grps} that denotes the reporting
#' group to which each population belongs. This function writes that information out to
#' the \code{repufile} that can be read by gsi_sim
#' 
#'   @param pops character vector of populations in the baseline
#'   @param grps factor parallel to pops that denotes the reporting group of each pop in pops
#'   @param repufile the name of the file to write the reporting units file to  
#'   @export
gsi_WriteReportingUnitsFile <- function(pops, grps, repufile="repunits.txt") {
  if(repufile != "") {
    unlink(repufile) # make sure that it is empty when we start appending to it 
  }
  if(anyDuplicated(pops)) stop("Duplicated population names in pops argument to gsi_WriteReportingUnitsFile")
  rl <- split(pops, grps)
  
  catchit <- lapply(names(rl), function(x) {
    write(paste("REPUNIT", x), file=repufile, append=T)
    write(cbind(paste("\t", rl[[x]], sep="")), file=repufile, append=T)
  }
  ) 
}



#' extract self-assignment results from the stdout dump file that gsi_sim produces
#' 
#' \code{gsi.simSelfAss2DF} reads the output of a gsi_sim self-assignment run and extracts the results into a list of data frames,
#' one for posteriors and the other for log-likelihood ratios.  This assumes that the individual IDs in the gsi_sim
#' output are what they would be for the pipeline, so you can get the population name by stripping
#' the numbers off them.
#' 
#' @param file the path to the gsi_sim output file
#' @export
gsi.simSelfAss2DF <- function(file="self-ass-output.txt") {
  x <- readLines(file)
  x1<-x[grep("UNSORTED_SELF_ASS_LIKE_GC_CSV:", x)]  # get just the lines we want
  x2<-read.table(textConnection(x1), sep=";", stringsAsFactors=F)
  numpops <- (ncol(x2)-4)/4  # remove one column for each of ID, NumMissingLoci, NonMissingLocusNames, and MissingLocusNames
  popnames <- as.character(x2[1,seq(from=2, by=3, length.out=numpops)])  # names of all the pops in the baseline, in that order
  IDs <- gsub("UNSORTED_SELF_ASS_LIKE_GC_CSV:/", "", x2$V1)
  FromPops <- gsub("[0-9]*","", IDs) 
  FromPops <- factor(FromPops, levels=popnames)  # make them a factor that preserves the order they went into gsi sim with
  Posteriors <- x2[,seq(from=3, by=3, length.out=numpops)]/100.0
  LogLs <- x2[, seq(from=2+3*numpops, length.out=numpops)]
  NumLoci <- x2[, ncol(x2)-2]
  
  # send result back as a list of data frames, either Posteriors or LogLs
  lapply(list(Post=Posteriors, LogLs=LogLs), function(z) {
    names(z) <- popnames
    ll <- cbind(PopulationOfOrigin=FromPops, NumberOfLoci=NumLoci, z)
    rownames(ll) <- IDs
    ll
  }
  )
}


#' sum gsi_sim posteriors over populations in reporting units
#' 
#' given a data frame x that has columns represententing scores (like posterior
#' probs) for each population in the pops.str (which should be a character
#' vector of population names) which are themselves grouped into reporting
#' units according to rep.units (preferably a factor), this returns a new data frame
#' with all the non-pop columns first and then the pop columns aggregated by summing
#' them up.
#' 
#' @param x data frame like that obtained by reading in "pop_pofz_full_em_mle.txt" or a similar output file
#' from gsi_sim (with rownames=1)
#' @param pops.str names of the populations that are the column names of x
#' @param rep.units  factor vector of same length as pops.str which given the reporting units of the populations in pops.str
#' @export
gsi_aggScoresByRepUnit <- function(x, pops.str, rep.units) {
  other.cols <- x[ !(names(x) %in% pops.str) ]
  grps <- split(pops.str, rep.units)
  sapply(grps, function(z) rowSums(x[z]))
  cbind(other.cols, sapply(grps, function(z) rowSums(x[z])))
}


#' get maximum posterior population assignment
#' 
#' pass this a data frame d that has the posterior probabilities in columns
#' with indexes (or names) in variable columns, and this will return a vector
#' of the names of the column in which the max posterior occurs.  And it also 
#' returns the maximum posterior associated with that assignment
#' @param d data frame of posterior probabilities like that obtained by reading "pop_pofz_full_em_mle.txt" with rownames=1
#' @param columns  not sure what this is for but i think it lets you select only some columns from the data frame d
#' @export
gsi_simMaxColumnAndPost <- function(d, columns) {
  tmp <- as.matrix(d[,columns])
  nam <- colnames(tmp)
  maxcolumns <- nam[max.col(tmp)]  # hold this temporarily
  maxcolumns <- factor(maxcolumns, levels=nam)  # preserve the order they were in d
  data.frame(MaxColumn=maxcolumns, MaxPosterior=apply(tmp, 1, max))
}


#' create an "assignment matrix" from gsi_sim results that have been processed
#' 
#' given a factor fp giving the "FromPopulation" and a factor mc given the "MaxColumn"
#' unit assigned to, and a numeric vector mp of the maximum posterior assignment scores,
#' and a numeric vector of cutoffs for the max posterior, this returns a list indexed by
#' cutoffs of 1) The number of individuals assigned 2)The assignment table of the assigned individuals
#' 
#' @param fp factor giving the population an individual came from
#' @param mc factor giving where the individual was assigned to
#' @param mp the value of the posterior with which the individual was assigned
#' @param cutoffs vector of cutoff values to use
#' @export
gsi_simAssTableVariousCutoffs <- function(fp, mc, mp, cutoffs=seq(0,.95, by=.05)) {
  names(cutoffs) <- paste("Cutoff",cutoffs, sep="_")
  lapply(cutoffs, function(x) {
    f <- fp[mp>=x]
    m <- mc[mp>=x]
    AssTable=table(f,m)
    names(dimnames(AssTable)) <- c("From", "To")
    list(NumAssigned=sum(mp>=x), AssTable=AssTable)
  }
  )
}





#' run  gsi_sim self-assignment exercises using different sets of loci
#'
#' given the baseline in the data frame B which is in \code{\link{gPipe.data.frame}} format
#' and a list of one or more sets of loci to use, this function computes self
#' assignments of the fish in the baseline for each set of loci.  Some day I will
#' figure out how to parallelize this with mclapply, but for now I just use lapply---there
#' seemed to be some interaction between system2 and mclapply, or something...
#'
#' @param B   a genetic baseline in \code{\link{gPipe.data.frame}} format
#' @param L   a list of vectors of names of loci that are to be used.  The name of each locus 
#'  					to be used should appear just once.  Obviously these names must correspond to locus
#'						names in the data frame B.  Note that you don't add the ".1"s to the names of the loci.
#'						If L is NULL (the default) then it is assumed that B contains only genotype data and the
#'						names of the loci are inferred from that.
#' @param pops.f 	a factor vector of populations each individual in B shall belong to.  levels of the
#'							factor should be in the order that these should appear in.  pops.f must have length
#'							equal to the number of rows of B
#' @param dir.to.use the directory to do the gsi_sim stuff in.  By default it is the current
#' 							directory, but you could use a tempfile() if you wanted to.
#' @export
gsi.a.SelfAssign <- function(B, L=NULL, pops.f, dir.to.use=NULL) {
	if(is.null(L)) L <- list(names(B)[c(T,F)])  # this picks out the odd column names
	if(length(L)==0) stop("list of locus sets has zero length") 
	if(any(duplicated(names(L)))==T) stop(paste("duplicate locus set names in L:", names(L)[duplicated(names(L))]))
	if(length(pops.f) != nrow(B)) stop("Length of pops.f must be equal to nrow(B)")
	if(!is.null(dir.to.use)) {
		td = dir.to.use
		if(!file.exists(td)) dir.create(td)
		curd <- getwd()
		setwd(td)
	} else {
		td <- getwd()
	}
	
	cols <- lapply(L, function(x) paste(rep(x, each=2), c("", ".1"), sep="")) # these include the column names with the .1's
	
	# now check to make sure the requested locus names exist:
	notThereLoci <- sort(unique(unlist(lapply(cols, function(x) setdiff(x, names(B))))))
	if(length(notThereLoci>0)) stop(paste(c("Requested locus columns not found in B: ", notThereLoci), collapse=" "))
	
	# just say where we are going to run gsi_sim
	print(paste("Running gsi_sim in", td))
	
	res <- lapply(cols, function(x) {
			
			gPdf2gsi.sim(B[,x], pops.f)  # make the gsi_sim data file with just the requested loci
			gsi_Run_gsi_sim("-b gsi_sim_file.txt --self-assign") # do the self-assignment
			SA <- gsi.simSelfAss2DF(file="GSI_SIM_Dumpola.txt")$Post # get the self-assignment posteriors
			# return a list that includes the loci used and SA
			list(LociUsed=x[c(T,F)], SelfAss=SA)
		})
		if(!is.null(dir.to.use)) setwd(curd)
#		names(res) <- names(cols)
		class(res) <- "SelfAssResult"  # give this list a class so that we can check for it when we summarize these by different reporting units
		res
}


#' Does mixture simulations from a baseline given specified mixture proportions
#' 
#' The user can specify a gpiper input file and a vector of populations
#' as well as a vector of proportions of populations and a size of sample
#' and this will do all the leg-work do simulations with gsi_sim using 
#' either the fixed-Pi or the Pi variants.  Currently it only extracts the
#' mixing proportion estimates out of it.
#' 
#' @inheritParams  gsi.a.SelfAssign
#' @param props  Vector of mixing proprtions parallel to the levels of pops.
#' These are the true mixing proportions that you want to simulate from
#' @param N the sample size of each sample you want to simulate. It will be approximate due to 
#' rounding.  I could fix that eventually, but it is a bit of a hassle at the moment.
#' @param Reps  Number of simulated samples you want to simulate
#' @param fixed Logical.  True means use the --fixed-Pi version (so that the exact
#' same number of fish are simulated each time.  False means that the origins of the fish
#' in the simulated sample are drawn from the props (by using the --Pi option of gsi_sim))
#' @param base_file_name The name for baseline file that will be created and passed to gsi_sim
#' @param repu.f A factor vector that is parallel to the levels of pops.f, each entry giving the
#' reporting unit that that particular level of pops.f belongs to.  This is used to aggregate the
#' estimated mixing proportions over reporting groups.  If NULL then it doesn't get done.
#' @param samp_unit What sampling unit to resample? Corresponds to gsi_sim's --samp_unit option.
#' Allowable values are "genes", "loci", and "multilocus".
#' @export
gsi.mixtureSims <- function(B, pops.f, props, N, Reps, fixed = TRUE, base_file_name = "mix_sim_base.txt",
                            repu.f = NULL, samp_unit = "genes") {

  if(length(props) != length(levels(pops.f))) stop("Length of props and length of levels(pops.f) don't match!")
  if(!is.null(repu.f) && length(repu.f) != length(levels(pops.f))) stop("Length of repu.f and length of levels(pops.f) don't match!")
  
  # make the gsi_sim baseline file we will need:
  gPdf2gsi.sim(B, pop.ize.them = pops.f, outfile = base_file_name)
  
  # now, figure out what numbers we will use
  if(fixed == TRUE) {
    PICOMM <- "--fixed-Pi"
    pis <- round(N * props)
  } else {
    PICOMM <- "--Pi"
    pis <- props
  }
  
  # now prepare the command lines with the --Pi's or --fixed-Pi's
  pi.commands <- paste(as.vector(rbind("-n", levels(pops.f), PICOMM, pis, "--end-pop")), collapse = " ")
  gsi.comm <- paste("-b", base_file_name, "--samp-unit", samp_unit, pi.commands, "-x", Reps, N, collapse = " ")
  
  # run gsi_sim
  gsi_Run_gsi_sim(arg.string = gsi.comm, stdout.to="MIX_SIM_Dumpola.txt", stderr.to="MIX_SIM_Stderr.txt")
  
  # parse the result
  x <- readLines("MIX_SIM_Dumpola.txt")
  mles <- x[grep("MIXFISH_PI_MLES:", x)]  # these are the MLES of Pi for each rep.
  mles <- gsub(pattern = "MIXFISH_PI_MLES: *", replacement = "", mles)  # get the Text off the front of it
  sim.mles <- simplify2array(strsplit(x = mles, "  *"))  # make them columns in an array
  rownames(sim.mles) <- levels(pops.f)
  mode(sim.mles) <- "numeric"  # make sure they are not strings any longer
  
  true_pi <- pis/sum(pis)
  names(true_pi) <- levels(pops.f)
  
  # now, if we have reporting units use those to get the estimated mixing proportions by reporting unit
  ret <- list(true_pi = true_pi, estimated_pi = sim.mles, repu_true_pi = NULL, repu_estimated_pi = NULL)
  if(!is.null(repu.f)) {
    ret$repu_true_pi <- sapply(levels(repu.f), function(y) sum(true_pi[which(repu.f==y)]))
    ret$repu_estimated_pi <- t(sapply(levels(repu.f), function(y) colSums(sim.mles[which(repu.f==y), ])))
  }
  
  ret
  
}
