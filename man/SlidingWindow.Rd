\name{SlidingWindow}
\alias{SlidingWindow}



\title{
Fisher exact test based sliding window
}


\description{
  Peak calling through Fisher exact test based sliding window method (see References).
}


\usage{
  SlidingWindow(input, RIP, mappedInput = NULL, 
                mappedRIP = NULL, level = 0.05, 
                ratio = -2, ...)
}

                  
\arguments{
  \item{input}{A string vector of file name, which specifies the input control samples in BAM format.}
    \item{RIP}{A string vector of file name, which specifies the IP samples in BAM format.}
    \item{mappedInput}{Integer, which represents the number of mapped reads number in the input control samples.}
    \item{mappedRIP}{Integer, which represents the number of mapped reads number in the IP samples.}
    \item{level}{Significance level for the p-value.}
    \item{ratio}{The ratio between the normalized mapped reads number in IP and input samples, respectively.}
    \item{\dots}{Further parameters used for peak calling.}
}




\value{
  A matrix in BED format containing the enriched peaks.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  #############reads mapping######################################
  input.fq <- "/home/malab14/input.fastq"
  RIP.fq <- "/home/malab14/RIP.fastq"
  referenceGenome <- "/home/malab14/tair10.fa"
  GTF <- "/home/malab14/Arabidopsis.gtf"
  input.bam <- readsMapping(alignment = "tophat", fq = input.fq,
                            refGenome = referenceGenome, paired = F)
  RIP.bam <- readsMapping(alignment = "tophat", fq = RIP.fq,
                          refGenome = referenceGenome, paired = F)
  ########################peak calling##############################
  peaks <- SlidingWindow(RIP = RIP.bam, input = input.bam,  
                         mappedInput = 20693001, 
                         mappedRIP = 22070873, refGenome = referenceGenome)
                          
}
}




\keyword{peak calling, epitranscriptome}

\references{Meyer, K.D., ea al. Comprehensive analysis of mRNA methylation reveals enrichment in 3'UTRs and near stop codons. Cell 2012; 149(7):1635-1646}
