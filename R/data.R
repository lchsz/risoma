
#' An IsomirDataSet object used for analysis of grape isomiR expression
#' profiles
#'
#' The isomirs IsomirDataSet object was constructed using publicly
#' available small RNA sequencing datasets from the Gene Expression Omnibus
#' (GEO) database (accession: GSE59802). This dataset comprises 50 samples
#' across 25 critical developmental stages, with two biological replicates per
#' stage. Processed FASTQ files were stored in the directory "fastqs", with
#' sample metadata recorded in "sample_info.csv".
#'
#'  \code{
#'    mirnas <- load_mirnas("vvi")
#'    isomirs <- detect_isomirs(sample_info_file = "sample_info.csv",
#'                                fq_dir = "fastqs",
#'                                mirnas = mirnas,
#'                                min_tpm = 5)
#'  }
#'
"isomirs"
