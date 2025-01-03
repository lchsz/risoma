% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mark_duplicates}
\alias{mark_duplicates}
\title{Mark duplicates in the FASTQ file}
\usage{
mark_duplicates(fq_file, min_read_num)
}
\arguments{
\item{fqFile}{FASTQ file of a sample}

\item{minReadNum}{the minimum read count for an isoform}
}
\value{
a data.frame with column 'read_seq' and 'read_num'
}
\description{
To speed up the computation, we first extract sequences of each read and
then mark and count duplicates.
}
