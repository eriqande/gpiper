% Generated by roxygen2 (4.0.1): do not edit by hand
\name{gsi_WriteReportingUnitsFile}
\alias{gsi_WriteReportingUnitsFile}
\title{create a reporting units file for the gsi_sim program}
\usage{
gsi_WriteReportingUnitsFile(pops, grps, repufile = "repunits.txt")
}
\arguments{
\item{pops}{character vector of populations in the baseline}

\item{grps}{factor parallel to pops that denotes the reporting group of each pop in pops}

\item{repufile}{the name of the file to write the reporting units file to}
}
\description{
Given a character vector \code{pops} of unique population names that are to be in a gsi_sim baseline
file and a factor of the same length, \code{grps} that denotes the reporting
group to which each population belongs. This function writes that information out to
the \code{repufile} that can be read by gsi_sim
}

