## docker run -it --mount type=bind,source=/Users/wongj4/CODEX2local,target=/TEST  migbro/codex2:3.8

library(CODEX2)
library(BSgenome.Hsapiens.UCSC.hg38)

## Initialization 

dirPath=getwd()
bamFile <- list.files(dirPath, pattern = '*.bam$')
bamdir <- file.path(dirPath, bamFile)
sampname <- substr(bamFile,1,9)
bedFile <- file.path(dirPath, "Strexome_canonical_100bp_padded_GRCh38.bed")
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile,
                       sampname = sampname, projectname = "CODEX2_demo")
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

## Getting GC content and mappability

gc <- getgc(ref, genome = BSgenome.Hsapiens.UCSC.hg38)
mapp <- getmapp(ref, genome = BSgenome.Hsapiens.UCSC.hg38) 
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))

## Getting depth coverage (returns read depth matrix & read lengths across all samples; generated for each chromosome)

coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(projectname, '_coverage.csv', sep=''), quote = FALSE)
head(Y[,1:2])

##  Quality Control (take sample-wise and exon-wise quality control procedure on depth of coverage matrix)

qcObj <- qc(Y, sampname, ref, cov_thresh = c(0, 4000),
            length_thresh = c(20, 2000), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
##################################################################################################
# filters out exons that: have extremely low coverage-median read
#      depth across all samples less than 20 or greater than 4000; are
#      extremely short-less than 20 bp; are extremely hard to map-
#      mappability less than 0.9; have extreme GC content-less than 20 or
#      greater than 80
##################################################################################################
Y_qc <- qcObj$Y_qc
sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc
qcmat <- qcObj$qcmat
gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(projectname, '_qcmat', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)


## Estimate library size factor for each sample based on genome-wide read depth after QC
## important because chrom read depth can be attenuated by large copy number abberations
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))}) ## which to use?
## pseudo.sample <- apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)
## N = library size factor, in this case ~0.72
plot(N, apply(Y,2,sum), xlab='Estimated library size factor', ylab='Total sum of reads')
dev.off()



## Running CODEX2 with negative control samples
## If there are negative control samples, use normalize_codex2_ns()
## If there are negative control regions, use normalize_codex2_nr()
## for case-control scenario, normal sample index is known (samples without spike-in signals)
chr<-c(1:22, 'X', 'Y') #This can be run for one chromosome or multiple chromosomes
## If the WES is designed
    #  under case-control setting, CODEX estimates the exon-wise Poisson
    #  latent factor using only the read depths in the control cohort,
    #  and then computes the sample-wise latent factor terms for the case
    #  samples by regression.
chr.index <- which(seqnames(ref_qc)==chr)
normObj <- normalize_codex2_ns(Y_qc = Y_qc[chr.index,],
                               gc_qc = gc_qc[chr.index], 
                               K = 1:2, norm_index = c(1,2),
                               N = N)
##################################################################################################
## Need to find optimal number of K for each run
## The larger the K is, the longer the estimation takes
## K = Number of latent Poisson factors. Can be an integer if
        #   optimal solution has been chosen or a vector of integers so
        #   that AIC, BIC, and RSS are computed for choice of optimal k.
## norm_index = indices of control samples
##################################################################################################

Yhat.ns <- normObj$Yhat; fGC.hat.ns <- normObj$fGC.hat;
beta.hat.ns <- normObj$beta.hat; g.hat.ns <- normObj$g.hat; h.hat.ns <- normObj$h.hat
AIC.ns <- normObj$AIC; BIC.ns <- normObj$BIC; RSS.ns <- normObj$RSS
# Yhat : Normalized read depth matrix
# fGC.hat : Estimated GC content bias matrix
# beta.hat : Estimated exon-specific bias as a vector
#  g.hat : Estimated Poisson latent factor
#  h.hat : Estimated Poisson latent factor

## Choose number of latent Poisson factors
## BIC is default
choiceofK(AIC.ns, BIC.ns, RSS.ns, K = 1:5 , filename = "codex2_ns_choiceofK.pdf")
par(mfrow = c(1, 3))
plot(1:4, RSS, type = "b", xlab = "Number of latent variables", pch=20)
plot(1:4, AIC, type = "b", xlab = "Number of latent variables", pch=20)
plot(1:4, BIC, type = "b", xlab = "Number of latent variables", pch=20)
par(mfrow = c(1,1))

## Running segmentation (can do Poisson or HMM but recommends Poisson)
## For negative control regions, use Yhat.nr
## For germline CNV, use 'integer' mode
## for CNV detection in heterogenous sample (somatic copy number changes in bulk cancer samples), use 'fraction' mode
finalcall.CBS <- segmentCBS(Y_qc[chr.index,],  
                            Yhat.ns, optK = which.max(BIC.ns),
                            K = 1:5,
                            sampname_qc = paste('sample',1:ncol(Y_qc),sep=''),
                            ref_qc = ranges(ref_qc)[chr.index],
                            chr = chr, lmax = 400, mode = "fraction")
write.table(finalcall.CBS, file=paste(projectname, '_fractionMode_finalcallCNV', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)                           

## Post-segmentation pruning and filtering are recommended based on CNV length (filter 1), length per exon (filter2), likelihood ratio (filter3), and number of exons (filter4).
filter1 <- finalcall.CBS$length_kb<=200
filter2 <- finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
finalcall.CBS.filter <- finalcall.CBS[filter1 & filter2, ]

filter3 <- finalcall.CBS.filter$lratio>40
filter4 <- (finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
finalcall.CBS.filter=finalcall.CBS.filter[filter3|filter4,]
write.table(finalcall.CBS.filter, file=paste(projectname, '_fractionMode_finalcallCNV.filtered', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)  


