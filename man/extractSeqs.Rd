\name{extractSeqs}
\alias{extractSeqs}



\title{
Extract sequences according to samples
}


\description{
Extracting the sequences of fixed length according to the samples.
}


\usage{
extractSeqs(RNAseq, samples, seqLen = 101)
}

                  
\arguments{
   \item{RNAseq}{A character string specified the file name of cDNA sequences.}
   \item{samples}{A character vetor specified the samples.}
   \item{seqLen}{The length of extracted sequence.}
}



\value{
  A list containing the extracted sequences.
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
  
  posSeq <- extractSeqs(RNAseq = cDNA, samples = posSamples$positives,
                        seqLen = 101)
  
}
}



