% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risoma.R
\name{detect_one_sample}
\alias{detect_one_sample}
\title{============================== exported functions =============================
Find isoform from a sample}
\usage{
detect_one_sample(fq_file, mirnas, max_ed_5p, max_ed_3p)
}
\arguments{
\item{fq_file}{FASTQ of a samle}

\item{mirnas}{A data.frame storing mature miRNA ID (matrue_id),
mature miRNA sequence (mature_seq), mature start (mature_start),
seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
(pre_seq)}

\item{max_ed_5p}{the maximum distance between 5’ region of reads with
reference miRNA}

\item{max_ed_3p}{the maximum distance between 3' region of reads with
reference miRNA}
}
\value{
a data.frame with columns: mature_id, read_seq, read_num and dist
}
\description{
============================== exported functions =============================
Find isoform from a sample
}
