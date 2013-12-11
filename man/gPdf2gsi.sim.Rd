\name{gPdf2gsi.sim}
\alias{gPdf2gsi.sim}
\title{convert a \code{gPipe.data.frame} into a gsi_sim input file}
\usage{
  gPdf2gsi.sim(d, pop.ize.them = NULL,
    outfile = "gsi_sim_file.txt")
}
\arguments{
  \item{d}{the data frame.  See documentation for
  \code{gPipe.data.frame} for format}

  \item{pop.ize.them}{vector that should be of length
  \code{nrow(d)} that gives the populations in which the
  individuals belong.  It this is a factor, then
  populations are placed in the gsi_sim output file
  according to the order in the \code{levels} attribute of
  the factor.}

  \item{outfile}{path to the output file}
}
\description{
  gPdf2gsi.sim takes as input a data frame in
  \code{gPipe.data.frame} format converts it into a gsi_sim
  input file. If you supply a list of pop-names (preferably
  as a factor with levels ordered as you want the pops to
  appear) to the \code{pop.ize.them} variable, then it
  splits them into pops, ordered as in the factor
}

