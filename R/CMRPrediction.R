###########  binary coding ########

###########  binary coding ########


.Binary <- function(fasta){
  marker <- str_replace_all(string = fasta, pattern = "[^AaTUtuGgCc]", replacement = " 0 0 0 0")
  marker <- str_replace_all(string = marker, pattern = "[Aa]", replacement = " 1 0 0 0")
  marker <- str_replace_all(string = marker, pattern = "[TUtu]", replacement = " 0 1 0 0")
  marker <- str_replace_all(string = marker, pattern = "[Gg]", replacement = " 0 0 1 0")
  marker <- str_replace_all(string = marker, pattern = "[Cc]", replacement = " 0 0 0 1")
  marker <- na.omit(as.numeric(unlist(strsplit(marker, " "))))
  return(marker)
}


###########  K-mers coding ########   !!!!! just for one seq
.Kmer <- function(Seq, kmer=2){
  #### 1/0 to 1-4 sequence ###
  coding14 <- c()
  for (n in 1:(length(Seq)/4)){
    coding14[n] <- sum(which(Seq[1:4+4*(n-1)]==1))
  }
  #coding14[which(is.na(coding14))] <- 0
  #### k-mer types ###
  a=1:4
  matchtypes <- switch (kmer,
                        do.call(paste0,expand.grid(a)),
                        do.call(paste0,expand.grid(a,a)),
                        do.call(paste0,expand.grid(a,a,a)),
                        do.call(paste0,expand.grid(a,a,a,a)),
                        do.call(paste0,expand.grid(a,a,a,a,a))
  )
  #### incise sequence ####
  kmerseq <- c()
  for ( i in 1:(length(coding14)-(kmer-1))){
    kmerseq[i] <- paste(as.character( coding14[1:kmer+i-1]),collapse="")
  }
  #### matching kmers ###
  kmercoding <- table(kmerseq)[matchtypes]/length(kmerseq)
  kmercoding[which(is.na(kmercoding))] <- 0
  #print(paste0(kmer,"-mers calculate over"))
  kmercoding
}



.PseDNC <- function(Seq, lambda = 6, w = 0.9){
  
  phyrna_nucleo <- matrix(c(-12.2,-13.3,-14.2,-10.2,-7.6,-6.6,-10.2,-5.7,-8.0,-10.5,-12.2,-7.6,-7.6,-8.1,-10.2,-6.6,
                            -29.7,-35.5,-34.9,-26.2,-19.2,-18.4,-26.2,-15.5,-19.4,-27.8,-29.7,-19.2,-19.2,-22.6,-26.2,-18.4,
                            -3.26,-2.35,-3.42,-2.24,-2.08,-0.93,-2.24,-1.10,-2.36,-2.11,-3.26,-2.08,-2.11,-1.33,-2.35,-0.93),
                          16,3)
  colnames(phyrna_nucleo) <- c("Enthalpy(Ka/mol)","Entropy(eU)","FreeEnergy(Ka/mol)")
  rownames(phyrna_nucleo) <- c("33","31","34","32","13","11","14","12","43","41","44","42","23","21","24","22")
  phyrna_nucleo_scale <- apply(phyrna_nucleo, 2, scale)
  rownames(phyrna_nucleo_scale) <- c("33","31","34","32","13","11","14","12","43","41","44","42","23","21","24","22")
  
  
  
  phyrna_seq <- c()
  #### 1/0 to 1-4 sequence ###
  coding14 <- c()
  for (n in 1:(length(Seq)/4)){
    coding14[n] <- sum(which(Seq[1:4+4*(n-1)]==1))
  }
  coding14 <- coding14[which(coding14!=0)]
  #coding14[which(is.na(coding14))] <- 0
  #### k-mer types ###
  a=1:4
  matchtypes <- do.call(paste0,expand.grid(a,a))
  #### incise sequence ####
  kmerseq <- c()
  kmer <- 2
  for ( i in 1:(length(coding14)-(kmer-1))){
    kmerseq[i] <- paste(as.character( coding14[1:kmer+i-1]),collapse="")
  }
  #### tier
  k_tier <- function(kmerseq=kmerseq,lambda=lambda) {
    first_tier <- c()
    for(i in 1:(length(kmerseq)-lambda)){
      first_tier[i] <- mean((phyrna_nucleo_scale[kmerseq[i],]-phyrna_nucleo_scale[kmerseq[i+lambda],])^2)
    }
    first_tier
  }
  
  mer2 <- scale(.Kmer(Seq,2))
  p <- c()
  for (j in 1:lambda){
    p[j] <-  mean(k_tier(kmerseq,lambda=j))
  }
  for(j in 1:16){
    phyrna_seq[j] <- mer2[j]/(sum(mer2)+w*sum(p))
  }
  for(j in 17:(16+lambda)){
    phyrna_seq[j] <- (w*p[j-16])/(sum(mer2)+w*sum(p))
  }
  phyrna_seq
}

#' @export
featureEncoding <- function(RNAseq, lambda = 6, w = 0.9){
  
  ## binary feature
  cat("start converting sequences to binary features......", "\n")
  featureBinary <- t(sapply(RNAseq, .Binary))
  class(featureBinary) <- "numeric"
  
  cat("start calculating kmer-based features......", "\n")
  featureMatKmer2 <- t(apply(featureBinary, 1, .Kmer, kmer = 2))
  featureMatKmer1 <- t(apply(featureBinary, 1, .Kmer, kmer = 1))
  
  cat("start calculating PseDNC-based features......", "\n")
  featureMatPC <- t(apply(featureBinary, 1, .PseDNC, lambda = 6, w = 0.9))
  
  featureMat <- cbind(featureBinary, featureMatKmer1, featureMatKmer2, 
                      featureMatPC)
  colnames(featureMat) <- paste0("Feature", 1:ncol(featureMat))
  
  featureMat
}



.findPosSamples <- function(peaks, motifPos = NULL, RNAseq = NULL, ...){
  
  if(is.null(motifPos) & is.null(RNAseq)){
    stop("The motif positions for each RNA are not provided, please provide the
         RNA sequences")
  }else if(is.null(motifPos) & !is.null(RNAseq)){
    motifPos <- searchMotifPos(RNAseq = RNAseq, ...)
  }
  
  if( !all(peaks[,1] %in% names(motifPos) == TRUE) ){
    stop("The sequence IDs in the peaks must be included by the RNAseq!")
  }
  
  
  resPosMat <- matrix(NA, nrow = 1, ncol = 2)
  for(i in 1:nrow(peaks)){
    curID <- peaks[i,1]
    curStart <- as.numeric(peaks[i,2])
    curEnd <- as.numeric(peaks[i,3])
    curMotif <- motifPos[[curID]]
    curPos <- curMotif[which(curMotif >= curStart & curMotif <= curEnd)]
    tmpMat <- matrix(NA, nrow = length(curPos), 2)
    tmpMat[,1] <- curID
    tmpMat[,2] <- curPos
    resPosMat <- rbind(resPosMat, tmpMat)
  }
  resPosMat <- resPosMat[-1, , drop = FALSE]
  colnames(resPosMat) <- c("IDs", "Position")
  resPosMat
}

#' @export
findConfidentPosSamples <- function(peaks, RNAseq = NULL, motifPos = NULL, ...){
  
  if(is.null(motifPos) & is.null(RNAseq)){
    stop("The motif positions are not provided!")
  }
  
  if(is.null(motifPos)){
    motifPos <- searchMotifPos(sequence = RNAseq, ...)
  }
  
  posSamples <- NULL
  posMat <- matrix(NA, 1, 2)
  for(i in 1:nrow(peaks)){
    curRes <- .findPosSamples(peaks = peaks[i,,drop = FALSE], motifPos = motifPos)
    if(nrow(curRes) == 1){
      posSamples <- c(posSamples, paste0(curRes[1,1], "_", curRes[1,2]))
      posMat <- rbind(posMat, curRes)
    }
  }
  posMat <- posMat[-1,]
  posMat <- unique(posMat)
  posSamples <- unique(posSamples)
  resList <- list(positives = posSamples, cDNAID = posMat[,1])
  resList
}

#' @export
findUnlabelSamples <- function(cDNAID, motifPos, posSamples, ...){
  resSamples <- NULL
  for(i in 1:length(cDNAID)){
    curMotif <- motifPos[[cDNAID[i]]]
    curSamples <- paste0(cDNAID[i], "_", curMotif)
    curSamples <- setdiff(curSamples, posSamples)
    resSamples <- c(resSamples, curSamples)
  }
  resSamples
}



###################################

.extractSeq <- function(seqVec, seqLen, RNAseq){
  seqID <- seqVec[1]
  cenPos <- as.numeric(seqVec[2])
  cDNAseq <- as.character(RNAseq[[seqID]])
  start <- cenPos - floor(seqLen/2)
  stop <- cenPos + floor(seqLen/2)
  if(start < 1){
    curSeq <- substr(cDNAseq, 1, stop)
    kk <- seqLen - stop
    for(i in 1:kk){
      curSeq <- paste0("N", curSeq)
    }
  }else if(stop > nchar(cDNAseq)){
    curSeq <- substr(cDNAseq, start, nchar(cDNAseq))
    kk <- seqLen - (nchar(cDNAseq) - start + 1)
    for(i in 1:kk){
      curSeq <- paste0( curSeq, "N")
    }
  }else{
    curSeq <- substr(cDNAseq, start, stop)
  }

  curSeq
}
#
#' @export
extractSeqs <- function(RNAseq, samples, seqLen = 101){
   sampleMat <- do.call(rbind, strsplit(samples, "_"))
   sampleID <- sampleMat[,1]
   RNAseq <- read.fasta(RNAseq, as.string = T)
   if(!all(sampleID %in% names(RNAseq) == TRUE)){
     stop("The IDs in the sample matrix(sampleMat) are not fully included in the
          RNA sequences(RNAseq)")
   }

   resSeq <- vector("list", length = length(sampleID))
   names(resSeq) <- paste0(sampleID, "_", as.numeric(sampleMat[,2]))

   for(i in 1:length(resSeq)){
     resSeq[[i]] <- .extractSeq(seqVec = sampleMat[i,], seqLen = seqLen, RNAseq = RNAseq)
   }
 resSeq
}



#' @export
PSOL <- function(featureMatrix = NULL, positives, balanced = FALSE, ratio = 10, 
                 unlabels, cpus = 1, PSOLResDic = NULL, iterator = 10, TPR = 0.995,
                 cross = 5, method = c("randomForest", "svm"), ... ){
  
    
  if(is.null(featureMatrix)){
    stop("Parameter featureMatrix is necessary for psol.", "\n")
  }
  
  cat("Warnings: psol will cost lots of time and computational resource!", "\n")
  cat("Start using PSOL to select negative samples......", "\n")
  if(is.null(PSOLResDic)){
    PSOLResDic <- getwd()
    PSOLResDic <- paste0(PSOLResDic, "/")
  }
  intialNegatives <- .PSOL_InitialNegativeSelection(featureMatrix = featureMatrix, 
                                                   positives = positives, 
                                                   unlabels = unlabels, 
                                                   cpus = cpus, 
                                                   PSOLResDic = PSOLResDic, ...)
  negativeExpand <- .PSOL_NegativeExpansion(featureMat = featureMatrix,
                                           positives = positives, 
                                           negatives = intialNegatives$negatives, 
                                           unlabels = intialNegatives$unlabels, 
                                           cpus = cpus, iterator =  iterator, 
                                           cross = cross, PSOLResDic = PSOLResDic,
                                           method = method, ...)
  
  load(paste( PSOLResDic, "PSOL_Iteration_", iterator, ".RData", sep = ""))
  finalNegatives <- iterRes$finalNegatives
  featureMat <- featureMatrix[c(positives, finalNegatives), ]
  posLen <- length(positives)
  negLen <- length(finalNegatives)
  label <- c(rep(1, posLen), rep(0, negLen))
  fmat <- data.frame(featureMat)
  model <- randomForest(x = fmat, y = factor(label))
  iterRes[['model']] <- model
  
  if(balanced){
    finalNegatives <- sample(finalNegatives, size = length(positives))
  }
  iterRes[['finalNegatives']] <- finalNegatives
  iterRes
}
