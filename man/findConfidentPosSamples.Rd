\name{findConfidentPosSamples}
\alias{findConfidentPosSamples}



\title{
Find confident positive samples
}


\description{
The confident positive samples are identified if only one CMR signal existed in the corresponding peak.
}


\usage{
findConfidentPosSamples(peaks, RNAseq = NULL, motifPos = NULL, ...)
}

                  
\arguments{
   \item{peaks}{A BED format matrix containing the peak regions.}
   \item{RNAseq}{A character string specified the file name of cDNA sequences.}
   \item{motifPos}{A list containing the position of the consensus motif for each transcript.}
   \item{\dots}{Further paramters used for find confident positive samples.}
}



\value{
  A list containing two elements: the positive samples and its transcript IDs.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{

  cDNA <- '/home/malab14/tair10_cDNA.fa'
  GTF <- '/home/malab14/Arabidopsis.gtf'
  #Convert genomic position to cDNA position
  peaks <- G2T(bedPos = cmrMat, GTF = GTF)
  ##Search consensus motif in cDNA sequence
  motifPos <- searchMotifPos(RNAseq = cDNA)
  posSamples <- findConfidentPosSamples(peaks = peaks, motifPos = motifPos)
  
}
}



