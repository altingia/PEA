# **PEA**: An integrated R toolkit for <font color=red>p</font>lant <font color=red>e</font>pitranscriptome <font color=red>a</font>nalysis </br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
<br>
PEA is an integrated R toolkit that aims to facilitate the plant epitranscriptome analysis. This toolkit generates comprehensive results for CMR (chemical modifications of RNA) calling from epitranscriptome sequencing data, CMR predictions at the transcriptome scale, and CMR annotation (location distribution analysis, motif scanning and discovery, and gene functional enrichment analysis).
<br>
## Version and download <br>
* [Version 1.0](https://github.com/cma2015/EAP/blob/master/EAP_1.0.tar.gz) -First version released on september, 27th, 2017<br>
## Depends
#### R environment <br>
* [R](https://www.r-project.org/) (>= 3.3.1) <br>
* [randomForst](https://cran.r-project.org/web/packages/randomForest/index.html) (>= 0.6) <br>
* [seqinr](https://cran.rstudio.com/web/packages/seqinr/index.html) (>= 3.4-5) <br>
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (>= 1.2.0) <br>
* [snowfall](https://cran.r-project.org/web/packages/snowfall/index.html) (>= 1.84-6.1) <br>
* [bigmemory](https://cran.r-project.org/web/packages/bigmemory/index.html) (>= 4.5.19) <br>
* [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html) (>= 1.28.0) <br>
* [ROCR](http://bioconductor.org/packages/release/bioc/html/ROCR.html) (>= 1.0-7) <br>
* [pROC](http://bioconductor.org/packages/release/bioc/html/pROC.html) (>= 1.10.0) <br>
* [devtools](https://cran.r-project.org/web/packages/devtools/index.html) (>= 1.13.3) <br>
* [rGADEM](https://bioconductor.org/packages/release/bioc/html/rGADEM.html) (>= 2.24.0) <br>
#### Global software environment <br>
* [_Tophat/Tophat2:_](http://ccb.jhu.edu/software/tophat/index.shtml) Read mapping <br>
* [_Bowtie/Bowtie2:_](bowtie-bio.sourceforge.net/) Read mapping <br>
* [_Hisat/Hisat2:_](www.ccb.jhu.edu/software/hisat/) Read mapping <br>
* [_meRanTK:_](http://www.icbi.at/software/meRanTK/merantk.shtml) CMR calling for m5c<br>
#### Python environment <br>
* [_macs2:_](https://pypi.python.org/pypi/MACS2) Peak calling <br>
#### Dependency installation <br>
```R
## Install R Dependency
dependency.packages <- c("randomForest", "seqinr", "stringr", "snowfall", "bigmemory")
install.packages(dependency.packages)
biocondctor.dependency <- c("rGADEM", "Rsamtools")
source("https://bioconductor.org/biocLite.R")
biocLite(biocondctor.dependency)
```
```bash
## Install Tophat/Tophat2
sudo apt-get update
sudo apt-get install tophat or sudo apt-get install tophat2
## Install Bowtie/Bowtie2
sudo apt-get update
sudo apt-get install bowtie or sudo apt-get install bowtie2
## Install Hisat/Hisat2
sudo apt-get update
sudo apt-get install hisat or sudo apt-get install hisat2
```
```python
pip install macs2
```
[pip](https://www.saltycrane.com/blog/2010/02/how-install-pip-ubuntu/) <br>

## Installation <br>
```R
install.packages("Download path/PEA_1.0.tar.gz",repos = NULL, type = "source")
```
## Contents <br>
#### CMR calling <br>
* Read mapping <br>
* CMR calling from read-alignment files <br>
#### Transcriptome-level CMR prediction <br>
* _Arabidopsis_ m6A benchmark dataset construction <br>
* Feature encoding <br>
* m6A predictor construction using ML-based PSOL algorithm <br>
* Ten-fold cross-validation and ROC curve analysis <br>
#### CMR annotation <br>
* CMR location distribution <br>
* Motif scanning and discovery <br>
* Functional enrichment analysis of CMR corresponded genes <br>
## Quick start <br>
Here, we showcased the utility of PEA in N6-methyladenosine(m6A) sequence datasets. <br>
More details can be seen from [user manual](https://github.com/cma2015/EAP/blob/master/PEA.pdf). <br>
#### 1.CMR calling <br>
* 1.1 Read mapping <br>
```R
input.fq <- "/home/malab14/input.fastq"  
RIP.fq <- "/home/malab14/RIP.fastq"  
referenceGenome <- "/home/malab14/tair10.fa"  
GTF <- "/home/malab14/Arabidopsis.gtf"  
input.bam <- readMapping(alignment = "tophat", fq = input.fq,   
                         refGenome = referenceGenome, paired = F,
                         bowtie1 = NULL, ... = paste0(" -p 5 –G ", GTF))
RIP.bam <- readMapping(alignment = "tophat", fq = RIP.fq,   
                       refGenome = referenceGenome, paired = F,
                       bowtie1 = NULL, ... = paste0(" -p 5 –G ", GTF)) 
```
* 1.2 CMR calling from read-alignment files <br>
```R
################m6A peak calling through "SlidingWindow"##################  
################m6A peak calling through "SlidingWindow"##################  
cmrMat <- CMRCalling(CMR = "m6A", method = "SlidingWindow",  
                     IPBAM = RIP.bam, inputBAM = input.bam,
                     refGenome = referenceGenome, 
                     mappedInput = 17472368, mappedRIP = 20072602) 

```
#### 2.Transcriptome-level CMR prediction <br>
* 2.1 _Arabidopsis_ m6A benchmark dataset construction <br>
```R
cDNA <- "/home/malab14/tair10_cDNA.fa"  
GTF <- "/home/malab14/Arabidopsis.gtf"  
###Convert genomic position to cDNA position  
peaks <- G2T(bedPos = cmrMat, GTF = GTF)  
###Search consensus motif in cDNA sequence  
motifPos <- searchMotifPos(RNAseq = cDNA)  
posSamples <- findConfidentPosSamples(peaks = peaks,  
                                      motifPos = motifPos)  
unlabelSamples <- findUnlabelSamples(cDNAID = posSamples$cDNAID,   
                                     motifPos = motifPos, posSamples = posSamples$positives)
  
```
* 2.2 Feature encoding <br>
```R
#########################Extract sequence#################################  
positives <- posSamples$positives
posSeq <- extractSeq(RNAseq = cDNA, samples = positives, seqLen = 101)  
unlabelSeq <- extractSeq(RNAseq = cDNA, samples = unlabelSamples, 
                         seqLen = 101)  
#########################Feature encoding#################################  
posFeatureMat <- featureEncoding(posSeq)  
unlabelFeatureMat <- featureEncoding(unlabelSeq) 
featureMat <- rbind(posFeatureMat, unlabelFeatureMat)
```
* 2.3 m6A predictor construction using ML-based PSOL algorithm <br>
```R
###Setting the psol directory and running the PSOL-based ML classification###  
PSOLResDic <- "/home/malab14/psol/"  
psolResults <- PSOL(featureMatrix = featureMat, positives = positives,   
                    unlabels = unlabels, PSOLResDic = PSOLResDic, cpus = 5)
```
* 2.4 Ten-fold cross-validation and ROC curve analysis <br>
```
cvRes <- cross_validation(featureMat = featureMat,   
                          positives = rownames(posFeatureMat),  
                          negatives = rownames(unlabelFeatureMat),  
                          cross = 10, cpus = 1)
```
#### 3.CMR annotation <br>
PEA also provides one-command function “CMRAnnotation” for CMR annotation, including CMR location and enrichment profiling, CMR-related motif scanning and motif discovery, and CMR-related gene function enrichment analysis. The analysis results may provide further insights into spatial and functional associations of CMRs.
* 3.1 CMR location distribution <br>
```R
GTF <- "/home/malab14/Arabidopsis_tair10.gtf"  
#####Extract the UTR position information from GTF file and perform CMR location distribution analysis.  
UTRMat <- getUTR(GTF = GTF)  
results <- CMRAnnotation(cmrMat = cmrMat, SNR = T, UTRMat = UTRMat,
                         annotation = "location")  
```
* 3.2 Motif scanning and discovery <br>
```R
RNAseq <- "/home/malab14/tair10.fa" 
testSamples <- paste0(peaks[,1], "_", peaks[,2])
cmrSeq <- extractSeqs(samples = testSamples, RNAseq = RNAseq)  
## Note: this positions in the cmrMat were transcripomic
results <- CMRAnnotation(cmrSeq = cmrSeq, annotation = "motifScan")  
results <- CMRAnnotation(cmrSeq = cmrSeq, annotation = "motifDetect")
``` 
* 3.3 Functional enrichment analysis of CMR corresponded genes <br>
```R
enrichements <- CMRAnnotation(cmrMat= cmrMat, GTF = GTF,
                              annotation = "GO", topNodes = 20,
			      dataset = "athaliana_eg_gene") 
```
## Ask questions
Please use [PEA/issues](https://github.com/cma2015/PEA/issues) for how to use PEA and reporting bugs.
