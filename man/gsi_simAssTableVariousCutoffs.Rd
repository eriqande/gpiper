\name{gsi_simAssTableVariousCutoffs}
\alias{gsi_simAssTableVariousCutoffs}
\title{create an "assignment matrix" from gsi_sim results that have been processed}
\usage{
  gsi_simAssTableVariousCutoffs(fp, mc, mp,
    cutoffs = seq(0, 0.95, by = 0.05))
}
\arguments{
  \item{fp}{factor giving the population an individual came
  from}

  \item{mc}{factor giving where the individual was assigned
  to}

  \item{mp}{the value of the posterior with which the
  individual was assigned}

  \item{cutoffs}{vector of cutoff values to use}
}
\description{
  given a factor fp giving the "FromPopulation" and a
  factor mc given the "MaxColumn" unit assigned to, and a
  numeric vector mp of the maximum posterior assignment
  scores, and a numeric vector of cutoffs for the max
  posterior, this returns a list indexed by cutoffs of 1)
  The number of individuals assigned 2)The assignment table
  of the assigned individuals
}

