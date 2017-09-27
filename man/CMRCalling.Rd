\name{CMRCalling}
\alias{CMRCalling}



\title{
  Identify CMRs (chemical modifications of RNA) from epitranscriptome sequencing data
}


\description{
  CMRCalling is designed to identify different types of CMRs from epitranscriptome sequencing data.
}


\usage{
  CMRCalling(CMR = c("m6A", "m6Am", "m5C", "hm5C", "pseudoU", "m1A"),
             cpus = 1, IPBAM = NULL, inputBAM = NULL, 
             GTF = NULL,  paired = FALSE, ...)
}

                  
\arguments{
  \item{CMR}{A string, which specifies the CMR type.}
  \item{cpus}{A integer number specifying the number of cpus to be used for parallel computing.}
  \item{IPBAM}{A string vector of file name, which specifies the IP samples in BAM format.}
  \item{inputBAM}{A string vector of file name, which specifies the input control samples in BAM format.}
  \item{GTF}{A string vector of file name, which specifies the reference genome annotation in GTF format.}
  \item{paired}{Logical, if TRUE, the paired-end parameters will be used for CMR calling, otherwise, the single-end parameters would be used.}
  \item{\dots}{Further parameters used CMR calling.}
}



\value{
  A matrix in BED format containing the CMRs.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  ##########reads mapping#################################
  input.fq <- "/home/malab14/input.fastq"
  RIP.fq <- "/home/malab14/RIP.fastq"
  referenceGenome <- "/home/malab14/tair10.fa"
  GTF <- "/home/malab14/Arabidopsis.gtf"
  input.bam <- readsMapping(alignment = "tophat", fq = input.fq,
                            refGenome = referenceGenome, paired = F)
  RIP.bam <- readsMapping(alignment = "tophat", fq = RIP.fq,
                          refGenome = referenceGenome, paired = F)
  cmrMat <- CMRCalling(CMR = 'm6A', method = 'SlindingWindow', 
                       IPBAM = RIP.bam, inputBAM = input.bam)
                          
}
}




\keyword{CMR calling, epitranscriptome}

