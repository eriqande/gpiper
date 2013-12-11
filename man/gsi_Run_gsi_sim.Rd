\name{gsi_Run_gsi_sim}
\alias{gsi_Run_gsi_sim}
\title{Simple unix-like interface to run gsi_sim}
\usage{
  gsi_Run_gsi_sim(arg.string,
    stdout.to = "GSI_SIM_Dumpola.txt",
    stderr.to = "GSI_SIM_Stderr.txt")
}
\arguments{
  \item{arg.string}{a text string which includes all the
  options to send to gsi_sim}

  \item{stdout.to}{path to redirect gsi_sim's standard
  output (stdout) to}

  \item{stderr.to}{path to redicted gsi_sim's standard
  error to}
}
\description{
  This is just a wrapper.  It runs gsi_sim with the
  arguments in \code{arg.string} in the current working
  directory, and it sends the standard output to
  \code{stdout.to} and the stderr to \code{stderr.to}
}

