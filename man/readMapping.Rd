\name{readMapping}
\alias{readMapping}



\title{
  Reads mapping
}


\description{
  Mapping FASTQ format reads to the reference genome with user adjustable method.
}


\usage{
  readMapping(alignment = c("bowtie", "bowtie2", "tophat", 
                              "tophat2", "hisat", "hisat2"), 
                          fq, index = NULL, refGenome = NULL,
                          paired = FALSE, softwareDir = NULL, ...)
}

                  
\arguments{
  \item{alignment}{A string, which specifies the method to be used for reads mapping.}
  \item{fq}{A character vector representing file name in FASTQ format.}
  \item{index}{A character vector representing file name of index used for alignment. If not specified, it will generate index files automaticly.}
  \item{refGenome}{A character vector representing file name of reference genome in FASTA format.}
  \item{paired}{Logical, if TRUE, the paired-end parameters will be used for reads mapping, otherwise, the single-end parameters would be used.}
  \item{softwareDir}{Optional, if the corresponding alignment tools are not in the environment variables, please specify the directory of corresponding tools.}
  \item{\dots}{Further parameters used for reads mapping.}
}



\value{
  A list recording the results including:                                            \item{command}{   The command used for reads mapping.}
  \item{resDir}{    A character vector containg the directory of alignment files.}
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  ##########Single-end alignment#####################################
  fq <- '/home/malab14/test.fastq'
  referenceGenome <- '/home/malab14/tair10.fa'
  test.bam <- readsMapping(alignment = 'tophat', fq = fq, 
                           refGenome = referenceGenome, paired = F)
  
  ########Paired-end alignment######################################
  fq1 <- '/home/malab14/test1.fastq'
  fq2 <- '/home/malab14/test2.fastq'
  referenceGenome <- '/home/malab14/tair10.fa'
  test.bam <- readMapping(alignment = 'tophat', fq = c(fq1, fq2), 
                           refGenome = referenceGenome, paired = T)
  
  ################Custome parameters#################################
  fq <- '/home/malab14/test.fastq'
  referenceGenome <- '/home/malab14/tair10.fa'
  test.bam <- readMapping(alignment = 'tophat', fq = fq, 
                           refGenome = referenceGenome, 
                           paired = F, ... = c("-p 5"))
                          
}
}




\keyword{reads mapping}

