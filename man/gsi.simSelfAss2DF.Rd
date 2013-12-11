\name{gsi.simSelfAss2DF}
\alias{gsi.simSelfAss2DF}
\title{extract self-assignment results from the stdout dump file that gsi_sim produces}
\usage{
  gsi.simSelfAss2DF(file = "self-ass-output.txt")
}
\arguments{
  \item{file}{the path to the gsi_sim output file}
}
\description{
  \code{gsi.simSelfAss2DF} reads the output of a gsi_sim
  self-assignment run and extracts the results into a list
  of data frames, one for posteriors and the other for
  log-likelihood ratios.  This assumes that the individual
  IDs in the gsi_sim output are what they would be for the
  pipeline, so you can get the population name by stripping
  the numbers off them.
}

