% Generated by roxygen2 (4.0.1): do not edit by hand
\name{gsi_simMaxColumnAndPost}
\alias{gsi_simMaxColumnAndPost}
\title{get maximum posterior population assignment}
\usage{
gsi_simMaxColumnAndPost(d, columns)
}
\arguments{
\item{d}{data frame of posterior probabilities like that obtained by reading "pop_pofz_full_em_mle.txt" with rownames=1}

\item{columns}{not sure what this is for but i think it lets you select only some columns from the data frame d}
}
\description{
pass this a data frame d that has the posterior probabilities in columns
with indexes (or names) in variable columns, and this will return a vector
of the names of the column in which the max posterior occurs.  And it also
returns the maximum posterior associated with that assignment
}

