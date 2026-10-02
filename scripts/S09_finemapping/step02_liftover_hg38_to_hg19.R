#in an R (4.4) conda environment
library(rtracklayer)
library(GenomicRanges)
library(data.table)
UC<-data.table::fread("UC526KDC.03.ff.log.Pheno.glm.logistic.hybrid")
CD<-data.table::fread("/input/data/CDKDC.250120.firthfall.04.DIS.glm.logistic.hybrid")
colnames(CD)[1]<-"CHR"
colnames(UC)[1]<-"CHR"
colnames(CD)[2]<-"BP"
colnames(UC)[2]<-"BP"
CD[,ORIGIN_CHR:=CHR]
CD[,ORIGIN_BP:=BP]
UC[,ORIGIN_CHR:=CHR]
UC[,ORIGIN_BP:=BP]

CD_gr<-GRanges(seqnames=CD$CHR,ranges=IRanges(start=CD$BP,end=CD$BP),mcols=CD[,!c("CHR","BP"),with=FALSE])
UC_gr<-GRanges(seqnames=UC$CHR,ranges=IRanges(start=UC$BP,end=UC$BP),mcols=UC[,!c("CHR","BP"),with=FALSE])
chain<-import.chain("/input/reference/hg38ToHg19.over.chain")
lifted_CD<-liftOver(CD_gr,chain)
lifted_UC<-liftOver(UC_gr,chain)
lifted_CD_gr<-unlist(lifted_CD)
lifted_UC_gr<-unlist(lifted_UC)
results_CD<-as.data.table(as.data.frame(lifted_CD_gr))
results_UC<-as.data.table(as.data.frame(lifted_UC_gr))
results_CD[,BETA:=log(mcols.OR)]
results_UC[,BETA:=log(mcols.OR)]

fwrite(results_UC,"UC526KDC.liftover.txt",quote=FALSE)
fwrite(results_UC,"UC526KDC.liftover.txt")



