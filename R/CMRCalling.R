###########CMR (Chemical Modification of RNA) Calling
##########author: Jingjing Zhai, Chuang Ma
##########Contact: zhaijingjing603@gmail.com


.runMACS2 <- function(IPBAM, inputBAM, broad = FALSE, 
                     expName = NULL, paired = F,MACS2Dir, ...){
  
  cat("This function is used for peak calling through MACS2, did you install this 
      tool (Y/N)", "\n")
  x <- readLines(con = stdin(), 1)
  tt <- list(...)
  if(length(tt) != 0){
    callPara <- tt[[1]]
  }else{
    callPara <- NULL
  }
  
  if(x == "y" | x == "Y"){
    
    if(is.null(expName)){
      expName <- "example"
    }
    
    if(!broad){
      if(!paired){
        macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                              " -f BAM -n ", expName, " ", callPara)
      }else{
        macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                              " -f BAMPE -n ", expName, " ", callPara)
      }
      
    }else{
      
      if(!paired){
        macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                              " --broad ", callPara)
      }else{
        macsCommand <- paste0(MACS2Dir, " callpeak -t ", IPBAM, " -c ", inputBAM, 
                              " -f BAMPE --broad ", callPara)
      }
      
    }
    
    system(command = macsCommand)
    
  }else{
    
    stop("The MACS2 is not installed, Please try:", "\n", 
         "pip install MACS2\t", "in the terminal!")
    
  }
  
  peakDir <- paste0(expName, "_peaks.narrowPeak")
  peaks <- as.matrix(read.table(file = peakDir, sep = "\t", header = F, quote = ""))
  
  resList <- list(macsCommand = macsCommand, peaks = peaks)
  resList
}



.bam2bed <- function(BAM){
  
  cat("This function is used for convert BAM file to BED format through bedtools,
      did you install this tool (Y/N)", "\n")
  x <- readLines(con = stdin(), 1)
  if(x == "y" | x == "Y"){
    BAMFile <- unlist(strsplit(x = BAM, split = "/", fixed = T))
    BAMFile <- BAMFile[length(BAMFile)]
    BAMName <- unlist(strsplit(x = BAMFile, split = ".", fixed = T))[1]
    bedCommand <- paste0("bedtools bamtobed -i ", BAM, " > ", BAMName, ".bed")
    system(bedCommand)
    resBed <- paste0(BAMName, ".bed")
  }else{
    stop("Please install bedtools before running this function!")
  }
  
  resBed
}



########################################Bisulfite-seq######################
##########################meRanTK#########################################

.m5Call <- function(meRanDir, fq, paired = FALSE, refGenome, cpus, 
                   method = c("meRanGh", "meRanGs"), GTF = NULL,  
                   resDic = NULL, maxDup = 50, conversionRate = 0.99,
                   minBaseQ = 30, ...){
  
  if(length(method) != 1){
    method <- method[1]
    cat("Warnings: methods are not defined, the first meRanGh will be used!")
  }
  
  if(method == "meRanGs" & is.null(GTF)){
    stop("The GTF/GTF file is not provided, which is required for meRanGs!")
  }
  
  if(is.null(resDic)){
    resDic <- paste0(getwd(), "/")
  }
  
  
  ############Building index
  idxDir <- paste0(resDic, method, "IDX")
  if(method == "meRanGh"){
    idxCommand <- paste0(meRanDir, method, " mkbsidx -t ", cpus, " -fa ", refGenome, 
                         " -id ", idxDir)
  }else{
    idxCommand <- paste0(meRanDir, method, " mkbsidx -t ", cpus, " -fa ", refGenome, 
                         " -GTF ", GTF, " -GTFtagEPT Parent -GTFtagEPG gene",  
                         " -id ", idxDir)
  }
  
  system(command = idxCommand)
  
  #########Alignment##########################
  
  SRAID <- unlist(strsplit(x = fq[1], split = ".", fixed = TRUE))[1]
  alignDir <- paste0(resDic, method)
  dir.create(path = alignDir, showWarnings = FALSE)
  
  if(!paired){
    if(method == "meRanGh"){
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq, 
                             " -id ", idxDir, " -GTF ", GTF, " -bg -o ", alignDir,
                             " -S ", method, "_", SRAID, ".sam -MM -un")
    }else{
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq, 
                             " -id ", idxDir, " -GTF ", GTF, " -bg -o ", alignDir,
                             " -S ", method, "_", SRAID, ".sam -MM -un --star_genomeLoad")
    }
  }else{
    if(method == "meRanGh"){
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq[1], 
                             " ", fq[2], " -id ", idxDir, " -GTF ", GTF, " -bg -o ",
                             alignDir, " -S ", method, "_", SRAID, ".sam -MM -un")
    }else{
      alignCommand <- paste0(meRanDir, method, " align -t ", cpus, " -f ", fq[1], 
                             " ", fq[2], " -id ", idxDir, " -GTF ", GTF, " -bg -o ",
                             alignDir, " -S ", method, "_", SRAID, ".sam -MM -un --star_genomeLoad")
    }
  }
  system(command = alignCommand)
  
  ########m5C Calling#############
  callCommand <- paste0(meRanDir, "meRanCall -p ", cpus, " -s ", alignDir, "/", 
                        method, "_", SRAID, "_sorted.bam -f ", refGenome, " -gref -o ",
                        resDic, method, "_meRanCall_m5C.txt", " -md ", maxDup,
                        " -cr ", conversionRate, " -mBQ ", minBaseQ, " -bed63 ")
  m5CDir <- paste0(resDic, method, "_meRanCall_m5C.txt")
  m5CresMat <- as.matrix(read.table(file = m5CDir, sep = "\t", header = T, quote = ""))
  
  resList <- list(idxCommand = idxCommand, alignCommand = alignCommand, 
                  callCommand = callCommand, m5C = m5CresMat)
  
  resList
}


#' @export
extractCov <- function(BAM, refGenome, method = c("bedtools", "Rsamtools")){
  
  if(length(method) > 1){
    method <- method[1]
    cat("Warnings: two methods are provided, only the first one (bedtools)", "\n",
        "will be used!")
  }
  
  if(method == "bedtools"){
    command <- paste0("genomeCoverageBed -split -dz -ibam ", BAM, " -g ", 
                      refGenome, " > ", paste0(BAM, ".dz"))
    system(command = command)
  }else{
    refGenome <- read.fasta(file = refGenome, as.string = T)
    Chr <- names(refGenome)
    resMat <- NULL
    for(i in 1:length(Chr)){
      seqID <- Chr[i]
      resPos <- 1:(nchar(refGenome[[i]]))
      USiteIRange <- IRanges(resPos, resPos)
      names(USiteIRange) <- seqID
      tmpList <- list()
      tmpList[[seqID]] <- USiteIRange
      which <- RangesList(unlist(tmpList))
      p1 <- ScanBamParam(which = which, what = scanBamWhat())
      res1 <- scanBam(BAM, param = p1)
      number <- lapply(res1, function(x) nrow(as.data.frame(x)))
      tmpMat <- data.frame(Chr = rep(seqID, length(number)), Position = 1:(nchar(refGenome[[i]])),
                           number = unlist(number))
      resMat <- rbind(resMat, res1)
    }
    write.table(resMat, file = paste0(BAM, ".dz"), sep = "\t",
                quote = F, row.names = F, col.names = F)
  }
  res <- paste0(BAM, ".dz")
  res
}


.findSeed <- function(index){
  res1 <- NULL
  res2 <- NULL
  index <- c(index, 0)
  for(i in 1:(length(index)-1)){
    
    if(i >= 2){
      if(((index[i+1] - index[i]) == 1) & ((index[i] - index[i-1])>=2)){
        res1 <- c(res1, i)
      }
      
      if(((index[i] - index[i-1]) == 1) & ((index[i] - index[i+1])!=-1)){
        res2 <- c(res2, i)
      }
    }else{
      if((index[i] - index[i+1]) == -1){
        res1 <- c(res1, i)
      }
    }
    
    
  }
  resMat <- cbind(res1, res2)
  resMat
}

.intervalCov <- function(dzFile){
  
  baseCov <- read.table(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
  winDow <- rep(NA, nrow(baseCov))
  baseCov <- cbind(baseCov, winDow)
  colnames(baseCov) <- c("Chr", "Position", "readsNumber", "window")
  baseCov$window <- ceiling(baseCov$Position/25)
  
  Chr <- unique(baseCov$Chr)
  resList <- vector("list", length = length(Chr))
  names(resList) <- Chr
  for(i in 1:length(Chr)){
    curChr <- Chr[i]
    curMat <- subset(x = baseCov, baseCov$Chr == curChr)
    curWave <- rep(NA, nrow(curMat))
    curMat <- cbind(curMat, curWave)
    curMat$curWave <- ave(curMat$readsNumber, curMat$window, FUN=mean)
    curMat <- subset(curMat, !duplicated(curMat$window))
    curMat <- na.omit(curMat)
    resList[[curChr]] <- curMat
  }
  resList
}

#' @export
SlidingWindow <- function(input, RIP, mappedInput = NULL, 
                          mappedRIP = NULL, level = 0.05, 
                          ratio = -2, ...){
  
  input.dz <- extractCov(BAM = input, ...)
  RIP.dz <- extractCov(BAM = RIP, ...)
  input <- .intervalCov(dzFile = input.dz)
  RIP <- .intervalCov(dzFile = RIP.dz)
  
  if(!all(names(input) == names(RIP))){
    stop("The chromosomes in the input and RIP are not consistent!")
  }
  
  resList <- vector("list", length = length(RIP))
  names(resList) <- names(RIP)
  #i=1
  for(i in 1:length(input)){
    curMat <- merge(input[[i]], RIP[[i]], by = 'window', all=TRUE)
    curMat$curWave.x[is.na(curMat$curWave.x)] <- 0
    curMat$curWave.y[is.na(curMat$curWave.y)] <- 0
    curMat[, 'windowave.input'] <- round(curMat$curWave.x)
    curMat[, 'windowave.RIP'] <- round(curMat$curWave.y)
    curMat <- curMat[curMat$windowave.RIP != 0, ]
    curPvalue <- rep(NA, nrow(curMat))
    curRatio <- rep(NA, nrow(curMat))
    curFDR <- rep(NA, nrow(curMat))
    curMat <- cbind(curMat, curPvalue, curRatio, curFDR)
    cat("Chromosome: ", names(input)[i], "\n")
    pb <- txtProgressBar(min = 0, max = nrow(curMat), style = 3,width = 75)
    for(j in 1:nrow(curMat)){
      testMat <- matrix(c(curMat$windowave.input[j],
                          mappedInput,
                          curMat$windowave.RIP[j],
                          mappedRIP), nrow = 2, ncol = 2)
      curMat$curPvalue[j] <- fisher.test(x = testMat)$p
      curMat$curRatio[j] <- log2(((curMat$windowave.input[j] + 1)*mappedRIP)/((curMat$windowave.RIP[j] + 1)*mappedInput))
      setTxtProgressBar(pb, j)
     # cat(j, "\n")
    }
    close(pb)
    curMat$curFDR <- p.adjust(curMat$curPvalue, "fdr")
    resList[[i]] <- curMat
  }
  resMat <- do.call(what = rbind, args = resList)
  sig <- rep(NA, nrow(resMat))
  
  for(i in 1:(nrow(resMat)-3)){
    if(resMat$curRatio[i]< ratio & resMat$curRatio[i+1] < ratio & 
       resMat$curRatio[i+2] < ratio & resMat$curRatio[i+3] < ratio & 
       resMat$curFDR[i] < level & resMat$curFDR[i+1] < level & 
       resMat$curFDR[i+2] < level & resMat$curFDR[i+3] < level ){
      sig[i] <- 1
    }
    else {
      sig[i] <- 0
    }
  }
  resMat <- cbind(resMat, sig)
  resMat <- subset(resMat, resMat$sig == 1)
  resPeaks <- NULL
  Chr <- names(RIP)
  for(i in 1:length(Chr)){
    curPeaks <- subset(resMat, resMat$Chr.x == Chr[i])
    curSeed <- .findSeed(index = curPeaks$window)
    if(!is.null(dim(curSeed))){
      tmpPeaks <- matrix(NA, nrow = nrow(curSeed), ncol = 3)
      tmpPeaks[,1] <- i
      for(j in 1:nrow(curSeed)){
        curSeed1 <- curSeed[j,1]
        curSeed2 <- curSeed[j,2]
        curStart <- curPeaks[which(curPeaks$window == curPeaks$window[curSeed1]), 3]
        curEnd <- curPeaks[which(curPeaks$window == curPeaks$window[curSeed2]), 3]
        tmpPeaks[j, 2:3] <- c(curStart, curEnd)
      }
    }
    
    resPeaks <- rbind(resPeaks, tmpPeaks)
  }
  resPeaks
}

#' @export
peakCalling <- function(IPBAM, inputBAM, GTF, expName = NULL,
                        method = c("SlidingWindow","exomePeak", "MetPeak", 
                                   "MACS2", "BayesPeak"), 
                        paired = F, ...){
  
  
  
  if(length(method) != 1){
    method <- method[1]
    cat("Warnings: multiple peak calling methods are provided, the first one
        will be used!")
  }
  
  if(is.null(expName)){
    expName <- "example"
  }
  
  if(method == "exomePeak"){
    cat("Start peak calling using exomePeak...", "\n")
    results <- exomepeak(GENE_ANNO_GTF = GTF, IP_BAM = IPBAM,
                         INPUT_BAM = inputBAM, ...)
    peaks <- read.table(paste0(getwd(), '/exomePeak_output/peak.bed'), 
                        sep = "\t", quote = "", stringsAsFactors = F)
    peaks <- results$all_peaks
    
  }else if(method == "BayesPeak"){
    cat("Start peak calling using BayesPeak...", "\n")
    IPBed <- bam2bed(BAM = IPBAM)
    inputBed <- bam2bed(BAM = inputBAM)
    peaks <- bayespeak(treatment = IPBed, control = inputBed, ...)
    
  }else if(method == "MACS2"){
    cat("Start peak calling using MACS2...", "\n")
    MACS2Dir <- system("which macs2", intern = TRUE)
    resList <- .runMACS2(IPBAM = IPBAM, inputBAM = inputBAM, 
                         paired = paired, MACS2Dir = MACS2Dir, ...)
    peaks <- resList$peaks
   
  }else if(method == "MetPeak"){
    cat("Start peak calling using MetPeak...", "\n")
    results <- metpeak(GENE_ANNO_GTF = GTF,IP_BAM = IPBAM, INPUT_BAM = inputBAM,
                       EXPERIMENT_NAME = expName, ...)
    resDir <- paste0(getwd(), "/", expName, "/peak.bed")
    peaks <- as.matrix(read.table(file = resDir, sep = "\t", quote = "", header = F))
  }else{
    cat("Start peak calling using Fisher exact test-based sliding window...", "\n")
    peaks <- SlidingWindow(input = inputBAM, RIP = IPBAM, ...)
  }
  resList <- list(peaks = peaks, method = method)
}


.sub <- function(number, resPos, res1){
  curNumber <- as.numeric(resPos[number])
  mapped <- as.data.frame(res1[[number]])
  if(nrow(mapped) == 0){
    ratio <- 0
  }else{
    ratio <- length(which(mapped[,5] == curNumber))/nrow(mapped)
  }
}

.pseudoURatio <- function(RNAseq, inputBAM){
  seqID <- attr(RNAseq, 'name')
  resPos <- words.pos(pattern = "[Tt]", text = RNAseq)
  names(resPos) <- seqID
  USiteIRange <- IRanges(resPos, resPos)
  names(USiteIRange) <- seqID
  tmpList <- list()
  tmpList[[seqID]] <- USiteIRange
  which <- RangesList(unlist(tmpList))
  p1 <- ScanBamParam(which = which, what = scanBamWhat())
  res1 <- scanBam(inputBAM, param = p1)
  resRatio <- apply(matrix(1:length(res1), nrow = length(res1), ncol = 1), 1, .sub, 
                    resPos = resPos, res1 = res1)
  
  resMat <- matrix(NA, nrow = length(resRatio), ncol = 3)
  resMat[,1] <- seqID
  resMat[,2] <- resPos
  resMat[,3] <- resRatio
  colnames(resMat) <- c("Transcript", "U.Position", "PseudoU.Ratio")
  #save(resMat, file = paste0(resDir, seqID, ".RData"))
  resMat
}


#' @export
pseudoURatio <- function(refGenome, inputBAM, cpus = 1){
  refGenome <- read.fasta(refGenome, as.string = T)
  resList <- vector(mode = "list", length = length(refGenome))
  names(resList) <- names(refGenome)
  if(cpus == 1){
    for(i in 1:length(resList)){
      resList[[i]] <- .pseudoURatio(RNAseq = refGenome[[i]], inputBAM = inputBAM)
    }
  }else{
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".pseudoURatio")
    sfExport(".sub")
    sfLibrary("Rsamtools", character.only = TRUE)
    sfLibrary("seqinr", character.only = TRUE)
    resList <- sfLapply(refGenome, .pseudoURatio, inputBAM = inputBAM)
    sfStop()
  }
  resList
}

#' @export
CMRCalling <- function(CMR = c("m6A", "m6Am", "m5C", "hm5C", "pseudoU", "m1A"),
                       cpus = 1, IPBAM = NULL, inputBAM = NULL, 
                       GTF = NULL,  paired = FALSE, ...){
  
  if(length(CMR) != 1){
    stop("Please specify a type of CMR!")
  }
  
  if(is.element(CMR, c("m6A", "m6Am", "m1A"))){
    resMat <- peakCalling(IPBAM = IPBAM, inputBAM = inputBAM, GTF = GTF, ...)
    resMat <- resMat$peaks
  }else if(is.element(CMR, c("m5C", "hm5C"))){
    resMat <- m5Call(paired = paired, cpus = cpus, GTF = GTF, ...)
    #resMat <- resMat$peaks
  }else{
    resMat <- pseudoURatio(refGenome = refGenome, inputBAM = inputBAM, cpus = cpus)
   # resMat <- resMat$peaks
    
  }
  
  resMat
}

