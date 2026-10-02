library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(dplyr)
library(tibble)
install.packages("msigdbr")
library(msigdbr)
library(edgeR)
library(susieR)
library(LDlinkR)
library(EnsDb.Hsapiens.v86)
sessionInfo()
CD_treat<-data.table::fread("/input/data/CD152.TREAT.HWE.PC30.LR.TREATMENT.glm.linear")
UC_treat<-data.table::fread("/input/data/UC124.TREAT.HWE01.PC30.LR.TREATMENT.glm.linear")
UC_hy<-NULL
UC_treat<-data.table::fread("/input/data/UC526KDC.03.ff.log.Pheno.glm.logistic.hybrid")
CD_treat<-data.table::fread("/input/data/CDKDC.250120.firthfall.04.DIS.glm.logistic.hybrid")

chain_file<-import.chain("/input/reference/hg19ToHg38.over.chain")
colnames(CD_treat)[1:2]<-c("CHR","BP")
colnames(UC_treat)[1:2]<-c("CHR","BP")
CD_treat[,oring_CHR:=CHR]
CD_treat[,oring_BP:=BP]
UC_treat[,oring_CHR:=CHR]
UC_treat[,oring_BP:=BP]
CD_treat[,CHR:=paste0("chr",CHR)]
UC_treat[,CHR:=paste0("chr",CHR)]


CD_treat_GR<-GRanges(seqnames = CD_treat$CHR,
                     ranges=IRanges(start=CD_treat$BP,end=CD_treat$BP),
                     mcols=CD_treat[,!c("CHR","BP")])
UC_treat_GR<-GRanges(seqnames = UC_treat$CHR,
                     ranges=IRanges(start=UC_treat$BP,end=UC_treat$BP),
                     mcols=UC_treat[,!c("CHR","BP")])

lifted_CD<-liftOver(CD_treat_GR,chain_file)
lifted_UC<-liftOver(UC_treat_GR,chain_file)
lifted_CD_gr<-unlist(lifted_CD)
lifted_UC_gr<-unlist(lifted_UC)
lifted_CD_gr<-as.data.table(as.data.frame(lifted_CD_gr))
lifted_UC_gr<-as.data.table(as.data.frame(lifted_UC_gr))
lifted_CD_gr$CHR<-gsub("chr","",lifted_CD_gr$CHR)
lifted_UC_gr$CHR<-gsub("chr","",lifted_UC_gr$CHR)

setnames(lifted_CD_gr,old=c("seqnames","start"),new=c("CHR","POS"))
setnames(lifted_UC_gr,old=c("seqnames","start"),new=c("CHR","POS"))

fwrite(lifted_CD_gr,"/output/251022_logistic_CD.txt",quote=FALSE,sep="\t")
fwrite(lifted_UC_gr,"/output/251022_logistic_UC.txt",quote=FALSE,sep="\t")
fwrite(lifted_CD_gr[!(lifted_CD_gr$mcols.ID==".")&(lifted_CD_gr$mcols.P<=1e-05),
                    c("mcols.ID","CHR","POS")],"/output/251022_linear_CD.for.magma.txt",quote=FALSE,sep=" ")
fwrite(lifted_UC_gr[!(lifted_UC_gr$mcols.ID==".")&(lifted_UC_gr$mcols.P<=1e-05),
                    c("mcols.ID","CHR","POS")],"/output/251022_linear_UC.for.magma.txt",quote=FALSE,sep=" ")


fwrite(lifted_CD_gr[!(lifted_CD_gr$mcols.ID==".")&(lifted_CD_gr$mcols.P<=1e-05),
                    c("mcols.ID","mcols.P")],"/output/251022_linear_CD.for.magma.rsIDP.txt",quote=FALSE,sep=" ")

fwrite(lifted_UC_gr[!(lifted_UC_gr$mcols.ID==".")&(lifted_UC_gr$mcols.P<=1e-05),
                    c("mcols.ID","mcols.P")],"/output/251022_linear_UC.for.magma.rsIDP.txt",quote=FALSE,sep=" ")



library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
attributes <- c("entrezgene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "description", "go_id", "name_1006")
CD_magma<-read.table("/output/magma/251022_magma_CD_152.treatment.genes.out.txt",header=TRUE)
UC_magma<-read.table("/output/magma/251022_magma_UC_124.treatment.genes.out.txt",header=TRUE)

gene_info_CD <- getBM(attributes = attributes,
                   filters = "entrezgene_id",
                   values = CD_magma$GENE,
                   mart = ensembl)

gene_info_UC<-getBM(attributes = attributes,
                    filters = "entrezgene_id",
                    values =UC_magma$GENE,
                    mart = ensembl)

CD_magma_results<-merge(CD_magma,gene_info_CD,by.x="GENE",by.y="entrezgene_id",all.x=TRUE)
UC_magma_results<-merge(UC_magma,gene_info_UC,by.x="GENE",by.y="entrezgene_id",all.x=TRUE)
CD_magma_results<-CD_magma_results %>% distinct(hgnc_symbol,.keep_all=TRUE)
UC_magma_results<-UC_magma_results %>% distinct(hgnc_symbol,.keep_all=TRUE)

msigdb <- msigdbr(
  species    = "Homo sapiens",
  collection = "C5"
)
unique(msigdb$gene_symbol)
unique(msigdb$gs_description)

filtered_genes_CD <- msigdb %>%
  filter(gene_symbol %in% CD_magma_results$hgnc_symbol)
check<-filtered_genes_CD %>% distinct(gene_symbol,.keep_all=TRUE)
CD_magma_msigDB.check=merge(CD_magma_results,filtered_genes_CD[,c(1,3,8,9)],by.x="hgnc_symbol",by.y="gene_symbol",all.x=TRUE)
CD_magma_msigDB.check[is.na(CD_magma_msigDB.check)]<-"_"
CD_magma_msigDB.check$description<-gsub(" ","_",CD_magma_msigDB.check$description)
CD_magma_msigDB.check$go_id<-gsub(" ","_",CD_magma_msigDB.check$go_id)
CD_magma_msigDB.check$name_1006<-gsub(" ","_",CD_magma_msigDB.check$name_1006)
CD_magma_msigDB.check$hgnc_symbol[1]<-"_"
table(CD_magma_msigDB.check$hgnc_symbol)
CD_magma_msigDB.check[order(CD_magma_msigDB.check$P_MULTI),]


library(data.table)
fwrite(CD_magma_msigDB.check[order(CD_magma_msigDB.check$P_MULTI),],"/output/251022_treatment_magma_1e05_CD152.check.+MSigdbC5.txt",quote=FALSE)



filtered_genes_UC <- msigdb %>%
  filter(gene_symbol %in% UC_magma_results$hgnc_symbol)
check<-filtered_genes_UC %>% distinct(gene_symbol,.keep_all=TRUE)
UC_magma_msigDB.check=merge(UC_magma_results,filtered_genes_UC[,c(1,3,8,9)],by.x="hgnc_symbol",by.y="gene_symbol",all.x=TRUE)
UC_magma_msigDB.check[is.na(UC_magma_msigDB.check)]<-"_"
UC_magma_msigDB.check$description<-gsub(" ","_",UC_magma_msigDB.check$description)
UC_magma_msigDB.check$go_id<-gsub(" ","_",UC_magma_msigDB.check$go_id)
UC_magma_msigDB.check$name_1006<-gsub(" ","_",UC_magma_msigDB.check$name_1006)
UC_magma_msigDB.check$hgnc_symbol[1]<-"_"
table(UC_magma_msigDB.check$hgnc_symbol)
UC_magma_msigDB.check[order(UC_magma_msigDB.check$P_MULTI),]
fwrite(UC_magma_msigDB.check[order(UC_magma_msigDB.check$P_MULTI),],"/output/251022_treatment_magma_1e05_UC124.check.+MsigdbC5.txt",quote=FALSE)

HG19<-fread("/input/reference/snp151.txt.gz")
HG19[1:10,]
HG19<-as.data.table(as.data.frame(HG19))
lifted_CD_gr[,Variant_GR19:=paste0("chr",CHR,"_",POS)]
lifted_UC_gr[,Variant_GR19:=paste0("chr",CHR,"_",POS)]
HG19[,Variant_GR19:=paste0(V2,"_",V3)]
str(HG19)
lifted_CD_gr.1<-merge(lifted_CD_gr,HG19[,c("V2","V3","V5","Variant_GR19")],by="Variant_GR19",all.x=TRUE)
lifted_UC_gr.1<-merge(lifted_UC_gr,HG19[,c("V2","V3","V5","Variant_GR19")],by="Variant_GR19",all.x=TRUE)
lifted_CD_gr.1[1:10,]
nrow(lifted_CD_gr.1[!(is.na(lifted_CD_gr.1$Variant_GR19))])
nrow(lifted_UC_gr.1[!(is.na(lifted_UC_gr.1$Variant_GR19))])
nrow(CD_treat)
nrow(UC_treat)
fwrite(lifted_CD_gr.1[!(is.na(lifted_CD_gr.1$Variant_GR19)),],"/output/251201_CD_rsid_19_check.txt",quote=FALSE)
fwrite(lifted_UC_gr.1[!(is.na(lifted_UC_gr.1$Variant_GR19)),],"/output/251201_UC_rsid_19_check.txt",quote=FALSE)


CD_data<-fread("/output/251201_CD_rsid_19_check.txt",sep=",")
