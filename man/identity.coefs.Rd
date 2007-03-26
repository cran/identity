\name{identity.coefs}
\alias{identity.coefs}
\title{Compute condensed identify coefficients}
\description{
  Compute condensed identify coefficients.
}
\usage{
identity.coefs (samples, pedigree)
}
\arguments{
  \item{samples}{either a vector of individual IDs or a matrix of two
    columns, each row reprepresents a pair of IDs.}
  \item{pedigree}{a matrix of three columns, id, mother id, and father id,
    respectively, specifying the pedigree structure.}
}
\value{
  a matrix of 11 columns, the first two columns give the pair of IDs,
  and the rest 9 columns give the 9 identity coefficients (delta 1 to
  delta 9).
}
\details{
  Id must be integers with 0 indicating individuals not included in the
  pedigree.

  Note that each time this function is called, it will recreate the C
  representation of the pedigree structure and the hash table.  Thus it
  is very inefficient to call this function many times with the same
  pedigree.  Normally you would only need call the function once to
  calculate the identity coefficients for all pairs of individuals.
}
\references{
  This R package is based on the C program Idcoefs 2.0 beta by
  Mark Abney (abney@uchicago.edu).

  See also Lynch and Walsh (1998) Genetics and Analysis of Quantitative
  Traits.
}
\examples{
ex <- read.table (system.file ("example", "ex.pedigree",
                  package = "identity"))
identity.coefs (c(150, 151, 152), ex)
}
\author{Na Li, \email{nali@umn.edu}}
\keyword{models}
