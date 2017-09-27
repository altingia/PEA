\name{featureEncoding}
\alias{featureEncoding}



\title{
Feature encoding
}


\description{
  This function contains three feature encoding scheme, binary, k-mer and PseDNC. For binary encoding scheme, a vector of 404 (4*101) features is generated through assigning 'A', 'C', 'G', 'U' and 'N' with (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1) and (0,0,0,0), respectively. Here 'N' is a gap used to ensure the fixed features of each sample, if an m6A/non- m6A site occurs near the initiation or termination of the transcript. For K-mer encoding, the composition of short sequence with different lengths was considered to encoding samples. For PseDNC (pseudo dinucleotide composition) encoding, the local and global sequence-order information along the RNA sequence was used for scoring the each sample.
}


\usage{
  featureEncoding(RNAseq, lambda = 6, w = 0.9)
}

                  
\arguments{
  \item{RNAseq}{A list containing the FASTA format sequences.}
     \item{lambda}{The lambda parameter for the PseDNC-related features, default is 6.}
  \item{w}{The weighting parameter for PseDNC-related features, default is 0.9.}
}




\value{
  A matrix with features.
}


\author{
  Jingjing Zhai, Chuang Ma, Jie Song
}


\examples{
\dontrun{
  
                          
}
}




\keyword{feature encoding}

