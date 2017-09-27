\name{CMRAnnotation}
\alias{CMRAnnotation}



\title{
  CRM Annotation
}


\description{
  This function is designed to perform annotation analysis including CMR location and enrichment profiling, CMR-related motif scanning and motif discovery, and CMR-related gene function enrichment analysis.
}


\usage{
  CMRAnnotation(cmrMat, genomic = F, UTRMat = NULL, GTF = NULL, SNR = T,
                annotation = c("location", "motifScan", "motifDetect", "GO"),
                cmrSeq = NULL, RNAseq = NULL, motifPos = NULL, plot = T, ...)
}

                  
\arguments{
  \item{cmrMat}{A BED format matrix representing the CMR position.}
  \item{genomic}{Logical, whic TRUE indicates the genomic position.}
  \item{UTRMat}{A matrix with six columns representing the relative position of 5'UTR, CDS and 3'UTR.}
  \item{GTF}{A string vector of file name, which specifies the reference genome annotation in GTF format.}
  \item{SNR}{Logical, where TRUE indicates the CMR is at single nucleotide resolution.}
  \item{annotation}{A string vector specified which annotation method would be used.}
  \item{cmrSeq}{A list containing CMR-realted sequences.}
  \item{RNAseq}{A string vector representing the RNA sequence in FASTA format.}
  \item{motifPos}{If SNR = F, the motifPos will be used, which indicates the relative position of specified consensus motif.}
  \item{plot}{Logical, where TRUE indicates perform a visualization for the results.}
  \item{\dots}{Further parameters used CMR annotation.}
}



\value{
  ##For location analysis:
  A list containing:
  \item{cmrMat}{The input CMR matrix.}
  \item{position}{The location of each CMR (e.g. 5'UTR, CDS or 3'UTR)}
  \item{distribution}{The relative position of CMR in the transcript.}
  ####For motif scan                                                                      
  
  A PWM matrix and motif logo.
  
  ####For motif detect
  
  A list containing consensus motifs in PWM (position weight matrix) format.
  
  ####For GO enrichment
  
  
  A list recording enrichment results for three subcategories (MF, BP and CC):
  \item{BP}{The GO enrichment results for BP(biological process).}
  \item{MF}{The GO enrichment results for MF(molecular function).}
  \item{CC}{The GO enrichment results for CC(cellular component).}
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  
  ######location distribution analysis######################
  GTF <- '/home/malab14/Arabidopsis_tair10.gtf'
  ##Extract the UTR position information from GTF file and perform location 
  ##distribution analysis
  UTRMat <- getUTR(GTF = GTF)
  results <- CMRAnnotation(cmrMat = cmrMat, SNR = T, UTRMat = UTRMat,
                           annotation = "location")
  
  #######################Motif analysis####################
  testSeq <- "/home/malab14/test.fa"
  results.scan <- CMRAnnotation(cmrSeq = testSeq, annotation = "motifScan")
  results.detect <- CMRAnnotation(cmrSeq = testSeq, annotation = "motifDetect")
  
  ######################GO analysis#######################
  enrichment <- CMRAnnotation(cmrMat = cmrMat, GTF = GTF, annotation = "GO",
                              topNodes = 20, dataset = "athaliana_eg_gene")
}
}




\keyword{CMR Annotation}

