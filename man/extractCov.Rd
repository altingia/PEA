\name{extractCov}
\alias{extractCov}



\title{
  Extract the reads coverage from BAM file
}


\description{
  Extracting reads coverage at each nucleotide using bedtools or Rsamtools.
}


\usage{
  extractCov(BAM, refGenome, method = c("bedtools", "Rsamtools"))
}

                  

\arguments{
  \item{BAM}{A character vector representing file name in bam format.}
  \item{refGenome}{A character vector representing file name of reference genome in FASTA format.}
  \item{method}{A string, which specifies the method to be used. For 'bedtools', 
  users have to install the bedtools in the terminal by typing "sudo apt-get install bedtools". if not specified, "bedtools" will be used as default.}
}



\value{
  A string containing file name of coverage.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{

   BAM <- '/home/malab14/test.bam'
   refGenome <- '/home/malab14/'
   results <- extractCov(BAM = BAM, refGenome = refGenome, method = "bedtools")
}
}



\keyword{reads coverage}

