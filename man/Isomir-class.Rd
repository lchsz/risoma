% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R
\docType{class}
\name{Isomir-class}
\alias{Isomir-class}
\title{An S4 class to represent an isomiR}
\description{
An S4 class to represent an isomiR
}
\section{Slots}{

\describe{
\item{\code{mature_id:}}{Reference miRNA ID.}

\item{\code{mature_seq:}}{Reference miRNA sequence.}

\item{\code{seed_seq:}}{The seed sequence of the mature miRNA.}

\item{\code{template_seq:}}{Template sequence including flanking regions.}

\item{\code{read_seqs:}}{Character vector of read sequences.}

\item{\code{read_num:}}{Abundance of the read.}

\item{\code{indel_5p:}}{Indels at the 5' end.}

\item{\code{indel_3p:}}{Indels at the 3' end.}

\item{\code{snp_5p:}}{Mismatches at the 5' end.}

\item{\code{snp_3p:}}{Mismatches at the 3' end.}

\item{\code{snp_seed:}}{Mismatches in the seed.}

\item{\code{dist:}}{Total edit distance.}

\item{\code{cigars:}}{CIGAR string representing alignment details.}

\item{\code{ids:}}{Character vector of unique isoform IDs.}

\item{\code{aln:}}{Multiple sequence alignment (MSA) of the isoform, reference miRNA
and template sequences}
}}

