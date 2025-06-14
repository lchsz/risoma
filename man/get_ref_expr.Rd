% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/methods.R
\docType{methods}
\name{get_ref_expr}
\alias{get_ref_expr}
\title{Extract Expression Data for Reference miRNAs}
\usage{
get_ref_expr(x)

\S4method{get_ref_expr}{IsomirDataSet}(x)
}
\arguments{
\item{x}{An \linkS4class{IsomirDataSet} object containing isomiR data.}
}
\value{
A numeric data frame where:
  \itemize{
    \item Rownames: Reference miRNA IDs.
    \item Columns: Experimental groups.
    \item Values: Expression values (TPM).
  }
}
\description{
This method extracts the expression data for reference miRNAs from an
\linkS4class{IsomirDataSet} object. It retrieves the expression values of
reference miRNAs across experimental groups and returns a data frame with
reference miRNA IDs as row names.
}
\examples{
# Load example IsomirDataSet
data(isomirs)

# Extract expression data for reference miRNAs
get_ref_expr(isomirs)

}
