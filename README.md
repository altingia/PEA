# _PEV:_ An integrated R toolkit for _p_ lant _e_ pitranscriptome _a_ nalysis <br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
<br>
PEA is an integrated R toolkit that aims to facilitate the plant epitranscriptome analysis. This toolkit generates comprehensive results for CMR (chemical modifications of RNA) calling from epitranscriptome sequencing data, CMR predictions at the transcriptome scale, and CMR annotation (location distribution analysis, motif scanning and discovery, and gene functional enrichment analysis).
<br>
## Version and download <br>
* [Version 1.0](https://github.com/cma2015/EAP/blob/master/EAP_1.0.tar.gz) -First version released on september, 26th, 2017<br>
## Dependency <br>
#### R environment <br>
* [R](https://www.r-project.org/) (>= 3.3.1) <br>
* [randomForst](https://cran.r-project.org/web/packages/randomForest/index.html) (>= 0.6) <br>
* [seqinr](https://cran.rstudio.com/web/packages/seqinr/index.html) (>= 3.4-5) <br>
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (>= 1.2.0) <br>
* [snowfall](https://cran.r-project.org/web/packages/snowfall/index.html) (>= 1.84-6.1) <br>
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) (>= 2.2.1) <br>
* [bigmemory](https://cran.r-project.org/web/packages/bigmemory/index.html) (>= 4.5.19) <br>
* [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html) (>= 1.28.0) <br>
* [motifRG](https://www.bioconductor.org/packages/devel/bioc/html/motifRG.html) (>= 1.21.0) <br>
* [devtools](https://cran.r-project.org/web/packages/devtools/index.html) (>= 1.13.3) <br>
* [rGADEM](https://bioconductor.org/packages/release/bioc/html/rGADEM.html) (>= 2.24.0) <br>
#### Global software environment <br>
* [_Tophat/Tophat2:_](http://ccb.jhu.edu/software/tophat/index.shtml)Reads mapping <br>
* [_Bowtie/Bowtie2:_](bowtie-bio.sourceforge.net/)Reads mapping <br>
* [_Hisat/Hisat2:_](www.ccb.jhu.edu/software/hisat/)Reads mapping <br>
#### Python environment <br>
* [_macs2:_](https://pypi.python.org/pypi/MACS2)Peak calling <br>
#### Dependency installation <br>
```R
## Install R Dependency
dependency.packages <- c("randomForst", "seqinr", "stringr", "snowfall",
                          "ggplot2", "bigmemory", "Rsamtools", "motifRG",
                          "devtools")
install.pacakages(dependency.packages)
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
install.package("Download path/EAP_1.0.tar.gz",repos = NULL, type = "source")
```
## Contents <br>
#### CMR calling <br>
* Arabidopsis thaliana m6A sequencing datasets <br>
* Read mapping <br>
* CMR calling from read-alignment files <br>
#### Transcriptome-level CMR prediction <br>
* Arabidopsis m6A benchmark dataset construction <br>
* Sample vectorization with three feature encoding schemes <br>
* m6A predictor construction using ML-based PSOL algorithm <br>
* Performance evaluation using the training dataset <br>
* Comparison with other m6A predictors using the independent testing dataset <br>
#### CMR annotation <br>
* CMR location distribution <br>
* Motif scanning and discovery <br>
* Functional enrichment analysis of CMR corresponded genes <br>
* [User manual](https://github.com/cma2015/EAP/blob/master/EAP.pdf)<br>

## Quick start <br>
More details please see [user manual](https://github.com/cma2015/EAP/blob/master/EAP.pdf) <br>
#### 1.CMR calling <br>
* 1.1 Arabidopsis thaliana m6A sequencing datasets <br>
```bash
## download and convert example data
# Here, supposing that ‘sratoolkit’ has been downloaded in your device in the directory: /home/malab14/, then the following command will convert the sra format to fastq format
/home/malab14/sratoolkit/bin/fastq-dump --split-3 SRR1508369.sra 
#Otherwise, you can also download the fastqc toolkit to perform quality control for fastq formatted files.
#Downloading the fastqc
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
#decompressing
unzip fastqc_v0.11.5.zip
#Running fastq
./FastQC/fastqc SRR1508369.fastq
```
* 1.2 Read mapping <br>
```R
input.fq <- "/home/malab14/input.fastq"  
RIP.fq <- "/home/malab14/RIP.fastq"  
referenceGenome <- "/home/malab14/tair10.fa"  
GTF <- "/home/malab14/Arabidopsis.gtf"  
input.bam <- readsMapping(alignment = "tophat", fq = input.fq,   
                          refGenome = referenceGenome, paired = F,
                          bowtie1 = NULL, p = 5, G = GTF)
RIP.bam <- readsMapping(alignment = "tophat", fq = RIP.fq,   
                        refGenome = referenceGenome, paired = F,
                        bowtie1 = NULL, p = 5, G = GTF)
```
* 1.3 CMR calling from read-alignment files <br>
```R
################m6A peak calling through "SlidingWindow"##################  
cmrMat <- CMRCalling(CMR = "m6A", method = "SlidingWindow",  
                     IPBAM = RIP.bam, inputBAM = input.bam) 
```
#### 2.Transcriptome-level CMR prediction <br>
* 2.1 Arabidopsis m6A benchmark dataset construction <br>
```R
cDNA <- "/home/malab14/tair10_cDNA.fa"  
GTF <- "/home/malab14/Arabidopsis.gtf"  
###Convert genomic position to cDNA position  
peaks <- G2T(bedPos = cmrMat, GTF = GTF)  
###Search consensus motif in cDNA sequence  
motifPos <- searchMotifPos(RNAseq = cDNA)  
posSamples <- findConfidentPosSamples(peaks = peaks,  
                                      motifPos = motifPos)  
unlabelSamples <- findUnlabelSamples(cDNAID = posSamples$cDNAID,   
                                     motifPos = motifPos,   
                                     posSamples = posSamples$positives)  
```
* 2.2 Sample vectorization with three feature encoding schemes <br>
```R
#########################Extract sequence#################################  
posSeq <- extractSeq(RNAseq = cDNA, sampleMat = posSamples, seqLen = 101)  
unlabelSeq <- extractSeq(RNAseq = cDNA, sampleMat = unlabelSamples, 
                         seqLen = 101)  
#########################Feature encoding#################################  
posFeatureMat <- featureEncoding(posSeq)  
unlabelFeatureMat <- featureEncoding(unlabelSeq) 
featureMat <- rbind(posFeatureMat, unlabelFeatureMat)
```
* 3.3 m6A predictor construction using ML-based PSOL algorithm <br>
```R
###Setting the psol directory and running the PSOL-based ML classification###  
PSOLResDic <- "/home/malab14/psol/"  
psolResults <- PSOL(featureMatrix = featureMat, positives = positives,   
                    unlabels = unlabels, PSOLResDic = PSOLResDic, cpus = 5) 

####Ten-fold cross-validation and ROC curve analysis.
cvRes <- cross_validation(featureMat = featureMat,   
                          positives = rownames(posFeatureMat),  
                          negatives = rownames(unlabelFeatureMat),  
                          cross = 10, cpus = 1)
```
#### 3.CMR annotation <br>
* 3.1 CMR location distribution <br>
```R
GTF <- "/home/malab14/Arabidopsis_tair10.gtf"  
#####Extract the UTR position information from GTF file and perform CMR location distribution analysis.  
UTRMat <- getUTR(GTF = GTF)  
results <- CMRAnnotation(cmrMat = cmrMat, SNR = T, UTRMat = UTRMat)  
```
* 3.2 Motif scanning and discovery <br>
```R
RNAseq <- "/home/malab14/tair10.fa"  
testSeqs <- extractSeqs(cmrMat = cmrMat, RNAseq = RNAseq)  
results <- motifScan(sequence = testSeqs, motif = "[AG][AG]AC[ACT]")
motifs <- motifDetect (sequence = testSeqs)  
``` 
* 3.3 Functional enrichment analysis of CMR corresponded genes <br>
```R
enrichements <- runTopGO(geneID = geneID,   
                         dataset = "athaliana_eg_gene",  
                         topNodes = 20) 
```
## Ask questions
Please use [EAP/issues](https://github.com/cma2015/EAP/issues) for how to use DeepGS and reporting bugs.
