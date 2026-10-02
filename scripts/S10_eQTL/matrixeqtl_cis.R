#08-1.library setting
library(biomaRt)
library(Biobase)
library(MatrixEQTL)
library(dplyr)
library(ggplot2)
library(data.table)

snp_file[(snp_file$rsID=="rs9743663"),4:ncol(snp_file)]
snp_file[1:10,1:10]

snp_file<-read.csv("/input/data/LTB4R_1525.csv",header=TRUE)
gene_file<-read.table("/input/data/251121_22316.CDonly.logCPM.txt",header=TRUE)
#gene_file.check<-gene_file
#pca_result <- prcomp(t(gene_file), scale = TRUE)
#write.table(pca_result$x,"/output/251121_22316.CDonly_logCPM.PC79.txt",quote=FALSE)
cvrt<-read.csv("/input/data/251121_22316.CDonly_cvrt.csv")

snp_pos<-snp_file %>% dplyr::select(c("rsID","CHR","POS"))
snp_pos$CHR<-gsub("chr","",snp_pos$CHR)
#gene_pos
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
attrs <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","transcription_start_site","strand",
           "go_id","name_1006","ucsc")
ensembl_gene_id <- sub("\\.\\d+$", "", rownames(gene_file))
#gene_pos<-gene_file %>% dplyr::select(c("ensembl_gene_id", "chromosome_name", "transcription_start_site"))
#260317_추가
snp_file<-read.csv("/input/data/HLADRB5_9533.v1.csv")
gene_file<-read.table("/input/data/251121_21198.UConly.logCPM.txt",header=TRUE)
snp_pos<-snp_file %>% dplyr::select(c("rsID","CHR","POS"))
snp_pos$CHR<-gsub("chr","",snp_pos$CHR)
#gene_pos
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
attrs <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","transcription_start_site","strand",
           "go_id","name_1006","ucsc")
ensembl_gene_id <- sub("\\.\\d+$", "", rownames(gene_file))


# 타임아웃 300초로 증가
results <- getBM(
  attributes = attrs,
  filters    = "ensembl_gene_id",
  values     = ensembl_gene_id,
  mart       = ensembl
)
results <- getBM(attributes = attrs, filters = "ensembl_gene_id", values = ensembl_gene_id, mart = ensembl)
results.1 <- results %>% dplyr::select(c("ensembl_gene_id","chromosome_name","start_position","end_position"))
results.1 <- results.1 %>% distinct(ensembl_gene_id,.keep_all = TRUE)
results.2 <- results %>% dplyr::select(c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","transcription_start_site","end_position","go_id","name_1006"))

snps.3 = SlicedData$new()
snps.3$fileDelimiter = ","     # the "," in snp_file and " " in gene file
snps.3$fileOmitCharacters = "NA" # denote missing values
snps.3$fileSkipRows = 1          # one row of column labels
snps.3$fileSkipColumns = 3       # one column of row labels 
snps.3$fileSliceSize = 2000     # read file in pieces of 2,000 rows
#snps.3$LoadFile("/input/data/LTB4R_1525.v1.csv",rowNamesColumn = 3)
snps.3$LoadFile("/input/data/HLADRB5_9533.v1.csv",rowNamesColumn = 3)
ncol(snps.3)


gene.3 = SlicedData$new()
gene.3$fileDelimiter = " "      # the TAB character
gene.3$fileOmitCharacters = "NA" # denote missing values;
gene.3$fileSkipRows = 1          # one row of column labels
gene.3$fileSkipColumns = 1      # one column of row labels
gene.3$fileSliceSize = 2000      # read file in pieces of 2,000 rows
#gene.3$LoadFile("/input/data/251121_22316.CDonly.logCPM.txt")

gene.3$LoadFile("/input/data/251121_21198.UConly.logCPM.txt")
rownames(gene.3) <- sub("\\.\\d+$", "", rownames(gene.3))
ncol(gene.3)

cvrt.3 = SlicedData$new()
cvrt.3$fileDelimiter = ","      # the TAB character
cvrt.3$fileOmitCharacters = "NA" # denote missing values;
cvrt.3$fileSkipRows = 1          # one row of column labels
cvrt.3$fileSkipColumns = 1      # one column of row labels
cvrt.3$fileSliceSize = 2000      # read file in pieces of 2,000 rows
#cvrt.3$LoadFile("/input/data/251121_22316.CDonly_cvrt_v1.csv",rowNamesColumn = 1)
cvrt.3$LoadFile("/input/data/251121_21198.UConly_cvrt_v1.csv",rowNamesColumn = 1)
colnames(cvrt.3)<-colnames(gene.3)
ncol(cvrt.3)

me.5 = Matrix_eQTL_main(
  snps = snps.3,
  gene = gene.3,
  cvrt = cvrt.3,
  pvOutputThreshold = 0.05,
  pvOutputThreshold.cis = 0.05,
  snpspos = snp_pos,
  genepos = results.1,
  cisDist = 1e6,
  #output_file_name.cis = "250925_CD_316.cis.txt",
  #output_file_name = "250925_CD_316.trans_eqtl.txt",
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

CD_cis<-me.5$cis$eqtls
CD_trans<-me.5$trans$eqtls
CD_cis %>% filter(snps == "rs9743663")
CD_cis %>% filter(snps == "rs67987819")

results.2$name_1006<-gsub(" ","_",results.2$name_1006)
results.2$go_id<-gsub(" ","_",results.2$go_id)
results.3<-results.2 %>% distinct(ensembl_gene_id,.keep_all=TRUE)
results.3 <- results.3 %>%
  mutate(
    go_id = ifelse(is.na(go_id) | go_id == "", "-", go_id),
    name_1006 = ifelse(is.na(name_1006) | name_1006 == "", "-", name_1006)
  )

CD_cis.real<-merge(CD_cis,results.3,by.x="gene",by.y="ensembl_gene_id",all.x=TRUE)
CD_cis.real<-CD_cis.real[order(CD_cis.real$FDR),]
CD_trans.real<-merge(CD_trans,results.3,by.x="gene",by.y="ensembl_gene_id",all.x=TRUE)
CD_trans.real<-CD_trans.real[order(CD_trans.real$FDR),]
CD_trans.real.1 <- CD_trans.real %>% distinct(snps,.keep_all=TRUE)
write.table(CD_cis.real,"/output/260313_CD_only_ciseQTL_check.txt",quote=FALSE)
write.table(CD_trans.real,"/output/260313_CD_only_transeQTL_check.txt",quote=FALSE)
write.table(CD_cis,"/output/260313_UC_only_ciseQTL_check.txt",quote=FALSE)
write.table(CD_trans,"/output/260313_UC_only_transeQTL_check.txt",quote=FALSE)

#boxplot
Frame_making<-function(genotype,exp,snpid,geneid){
  genotype<-genotype
  genotype.1<-genotype[(genotype$rsID==snpid),]
  genotype.1<-as.numeric(genotype.1[,c(4:ncol(genotype.1))])
  rownames(exp)<-sub("\\.\\d+$", "", rownames(exp))
  exp<-exp
  exp.1<-exp[(rownames(exp)==geneid),]
  exp.1<-as.numeric(exp.1)
  
  print("importing process was done")
  
  results<-data.frame(GENO=genotype.1,
                      TMM=exp.1,
                      Sample=colnames(exp))
  
  results<-results %>% mutate(Group=case_when(startsWith(Sample,"GTEX")~"Ctrl",
                                              TRUE ~ "Case"),
                              GENO_fac=case_when(GENO=="0" ~ "z",
                                                 GENO=="1" ~ "o",
                                                 GENO=="2" ~ "t"))
  print("generating results was done")
  
  results$Group=factor(results$Group,levels=c("Ctrl","Case"))
  results$GENO_fac=factor(results$GENO_fac,levels=c("z","o","t"))
  return(results)
  print("results omitting was done")
  
}

#CD_cis.real<-read.table("/output/251106_vcor+_CDCIS.txt",sep=" ",header=TRUE)
#CD_trans.real<-read.table("/output/251106_vcor+_CDTRANS.txt",sep=" ",header=TRUE)
snp_file<-snp_file[,1:82]
ncol(snp_file)

CD_cis.real[1:10,]
CD_cis.real %>% distinct(gene,.keep_all=TRUE) %>% filter(FDR<=0.05)
CD_cis.real %>% filter(snps=="rs9743663")

CD_ggplot<-Frame_making(snp_file,gene_file,"rs9743663","ENSG00000213903")#ENSG00000259431  rs35578109
UC_ggplot<-Frame_making(snp_file,gene_file,"rs67987819","ENSG00000196126") #HLA-DRB1
UC_ggplot<-Frame_making(snp_file,gene_file,"rs67987819","ENSG00000179344") #HLA-DQB1
UC_ggplot<-Frame_making(snp_file,gene_file,"rs67987819","ENSG00000198502") #HLA-DRB5 ~ two,exised.
#CD 확인 시 ENSG00000157379- DHRS1, 지질 대사 관련, IBD와 연결점 확인 필요는 ggplot상 확인이 X
#
CD_cis
CD_cis.real
p1<-ggplot(UC_ggplot,aes(GENO_fac,TMM,fill=GENO_fac))+
  geom_boxplot()+
  #facet_wrap(~ Group)+
  theme_classic()
print(p1)
t.test(UC_ggplot[(UC_ggplot$GENO_fac=="z"),"TMM"],
       UC_ggplot[(UC_ggplot$GENO_fac=="o"),"TMM"])
t.test(UC_ggplot[(UC_ggplot$GENO_fac=="z"),"TMM"],
       UC_ggplot[(UC_ggplot$GENO_fac=="t"),"TMM"])
t.test(UC_ggplot[(UC_ggplot$GENO_fac=="o"),"TMM"],
       UC_ggplot[(UC_ggplot$GENO_fac=="t"),"TMM"])


p<-ggplot(CD_ggplot,aes(GENO_fac,TMM,fill=GENO_fac))+
  geom_boxplot()+
  #facet_wrap(~ Group)+
  theme_classic()
print(p)
t.test(CD_ggplot[(CD_ggplot$GENO_fac=="z"),"TMM"],
       CD_ggplot[(CD_ggplot$GENO_fac=="o"),"TMM"])

summary(CD_ggplot)

CD_trans.real[1:2,1:5]

read_vcor_files <- function(list, path, pattern) {
  R2_list.CD<-list
  # pattern을 변수로 사용하여 파일 목록 추출
  files_to_read <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  R2_list.CD <- lapply(files_to_read, function(file) {
    if (file.info(file)$size == 0) {
      warning(paste("빈 파일 무시:", file))
      return(NULL)
    }
    tryCatch(read.table(file, header = FALSE),
             error = function(e) {
               warning(paste("읽기 실패 파일:", file, " - ", e$message))
               return(NULL)
             })
  })
  
  R2_list.CD <- Filter(Negate(is.null), R2_list.CD)
  
  library(data.table)
  library(dplyr)
  
  R2_list.CD.1 <- rbindlist(R2_list.CD, use.names = TRUE, fill = TRUE)
  
  R2_list.CD.1.1 <- R2_list.CD.1 %>%
    group_by(V3) %>%
    summarise(vcor = mean(V7, na.rm = TRUE))
  
  R2_list.CD.1.2 <- R2_list.CD.1 %>%
    group_by(V6) %>%
    summarise(vcor = mean(V7, na.rm = TRUE))
  
  colnames(R2_list.CD.1.1) <- c("ID","vcor")
  colnames(R2_list.CD.1.2) <- c("ID","vcor")
  
  Results_R2.CD <- bind_rows(R2_list.CD.1.1, R2_list.CD.1.2) %>%
    distinct(ID, .keep_all = TRUE)
  
  return(Results_R2.CD)
}
R2_list.CD<-list()

# 함수 사용 예시
R2_list.CD <- list()
path <- "/output/finemapping/"
pattern <- "1kb.*CDKDC.*\\.vcor$"
result.1kb <- read_vcor_files(R2_list.CD,path, pattern)

R2_list.CD <- list()
path <- "/output/finemapping/"
pattern <- "50kb.*CDKDC.*\\.vcor$"
result.50kb <- read_vcor_files(R2_list.CD,path, pattern)

R2_list.CD <- list()
path <- "/output/finemapping/"
pattern <- "100kb.*CDKDC.*\\.vcor$"
result.100kb <- read_vcor_files(R2_list.CD,path, pattern)

path <- "/output/finemapping/"
pattern <- "250kb.*CDKDC.*\\.vcor$"
result.250kb <- read_vcor_files(R2_list.CD,path, pattern)

path <- "/output/finemapping/"
pattern <- "500kb.*CDKDC.*\\.vcor$"
result.500kb <- read_vcor_files(R2_list.CD,path, pattern)

colnames(result.1kb)[2]<-"1kb_vcor"
colnames(result.50kb)[2]<-"50kb_vcor"
colnames(result.100kb)[2]<-"100kb_vcor"
colnames(result.250kb)[2]<-"250kb_vcor"
colnames(result.500kb)[2]<-"500kb_vcor"

merged_dt <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), list(result.1kb,result.50kb,result.100kb,result.250kb,result.500kb))

CDCIS_merge<-merge(CD_cis.real, merged_dt, by.x="snps",by.y="ID",all.x=TRUE)
CDCIS_merge[is.na(CDCIS_merge)]<-0
CDCIS_merge<-CDCIS_merge[order(CDCIS_merge$pvalue),]
write.table(CDCIS_merge,"/output/251106_vcor+_CDCIS.txt",quote=FALSE,row.names=FALSE)

CDTRANS_merge<-merge(CD_trans.real, merged_dt, by.x="snps",by.y="ID",all.x=TRUE)
CDTRANS_merge[is.na(CDTRANS_merge)]<-0
CDTRANS_merge<-CDTRANS_merge[order(CDTRANS_merge$pvalue),]
write.table(CDTRANS_merge,"/output/251106_vcor+_CDTRANS.txt",quote=FALSE,row.names=FALSE)

CD_CIS<-read.table("/output/251106_vcor+_CDCIS.txt",sep=" ",header=TRUE)
