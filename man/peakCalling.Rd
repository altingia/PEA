\name{peakCalling}
\alias{peakCalling}



\title{
Peak Calling
}


\description{
  Identify genomic regions in a genome that have been enriched with aligned reads as a consequence of performing MeRIP-seq or ChIP-seq experiment.
}


\usage{
  peakCalling(IPBAM, inputBAM, GTF, expName = NULL, 
              method = c("SlidingWindow","exomePeak", "MetPeak", 
                          "MACS2", "BayesPeak"), paired = F, ...)
}

                  
\arguments{
  \item{IPBAM}{A string vector of file name, which specifies the IP samples in BAM format.}
  \item{inputBAM}{A string vector of file name, which specifies the input control samples in BAM format.}
  \item{GTF}{A string vector of file name, which specifies the reference genome annotation in GTF format.}
  \item{expName}{A string, which specifies the name generated in the output directory which contains all the results.}
  \item{paired}{Logical, if TRUE, the paired-end parameters will be used for peak calling, otherwise, the single-end parameters would be used.}
  \item{method}{A string, which specifies the method used for peak calling, the default is "SlidingWindow".}
  \item{\dots}{Further parameters used peak calling.}
}


\value{
  A list containing two elements:
  \item{peaks}{A enriched peak region matrix in BED format.}
  \item{method}{A string represent which method was used for peak calling.}
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
  ########peak calling using Fisher exact-based sliding window#####
  peaks <- peakCalling(IPBAM = RIP.bam, inputBAM = input.bam, 
                       refGenome = referenceGenome,
                       method = 'SlindingWindow')
                          
}
}




\keyword{peak calling, epitranscriptome}

