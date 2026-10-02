#00.import library and dataset

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(dplyr)

CD<-data.table::fread("/input/data/CDKDC.250120.firthfall.04.DIS.glm.logistic.hybrid")
CD<-as.data.frame(CD)

UC<-data.table::fread("UC526KDC.03.ff.log.Pheno.glm.logistic.hybrid")
UC<-as.data.frame(UC)

colnames(CD)[1:2]<-c("CHROM","POS")
colnames(UC)[1:2]<-c("CHROM","POS")

#01.Data Manipulation

#01-1.chromosome-position give
CD <- CD %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(CD, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot)

UC <- UC %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(UC, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot)



#01-2. axisdf definition

axisCD <- CD %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisUC <- UC %>% group_by(CHROM) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

#01-3.known, novel, finemap annotate

anno<-data.table::fread("/input/data/chr7_annovar_check.hg38_multianno.1-12.txt")
anno<- anno %>% mutate(Variant_GRCh38=paste0(anno$Chr,":",anno$Start,":",anno$Ref,":",anno$Alt))

anno.1<-data.table::fread("/input/data/IBDKDC_total.hg38_multianno.Gene.ref.check.txt")
anno.1<- anno.1 %>% mutate(Variant_GRCh38=paste0(anno.1$Chr,":",anno.1$Start,":",anno.1$Ref,":",anno.1$Alt))

anno.2<-merge(anno.1[,c("Variant_GRCh38","Func.refGene","Gene.refGene","GeneDetail.refGene")],anno[,c("Variant_GRCh38","Func.refGene","Gene.refGene","GeneDetail.refGene")],all=TRUE)

print(anno.2[1:10,])

CD <- CD %>% mutate(suggest=case_when(P<=1e-05 ~ "yes",
                                      TRUE ~ "no"))
UC <- UC %>% mutate(suggest=case_when(P<=1e-05 ~ "yes",
                                      TRUE ~ "no"))

CD<- CD %>% mutate(Variant_GRCh38=paste0("chr",CD$CHROM,":",CD$POS,":",CD$REF,":",CD$ALT))
UC<- UC %>% mutate(Variant_GRCh38=paste0("chr",UC$CHROM,":",UC$POS,":",UC$REF,":",UC$ALT))

CD.1<-merge(CD, anno.2,all=TRUE)
UC.1<-merge(UC, anno.2,all=TRUE)

top10_CD <- CD.1 %>% filter(suggest == "yes") %>% arrange(P) %>% slice_head(n = 10)
top10_UC <- UC.1 %>% filter(suggest == "yes") %>% arrange(P) %>% slice_head(n = 10)

#01-4.manhattan plot drawing
#library(ggplot2)
CD.2<-ggplot(data=CD.1, aes(x=BPcum,y=-log10(P)))+
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values=rep(c("orange","blue"),22)) + 
 # geom_label(data = top10_CD, color = "black", aes(label = Gene.refGene), size = 6)+
 # geom_point(data=subset(CD.1,suggest=="yes"), color="red",size=4)+
 # geom_label(data=subset(CD.1,suggest=="yes"), color="black",aes(label=Gene.refGene),size=4)+
 # geom_label(data=subset(CD.1,suggest="yes"),color="black",aes(label=ID),size=4)+
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

name="CD"
#tiff(paste0("250812_",name,"_DIS.manhattan.tiff"), width=12, height=8, units="in", res=300)
tiff(paste0("250910_",name,"_DIS.no_annotated.CD237KDC.HLA+.manhattan.tiff"), width=12, height=8, units="in", res=300)


print(CD.2)

dev.off()

write.table(CD.1,"250912_CD_HLA+_labelled.hybrid.",quote=FALSE, row.names=FALSE)
print("CD manhattan plot complete")

#library(ggplot2)
UC.2<-ggplot(data=UC.1, aes(x=BPcum,y=-log10(P)))+
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values=rep(c("orange","blue"),22)) + 
 # geom_label(data = top10_UC, color = "black", aes(label = Gene.refGene), size = 6)+
 # geom_point(data=subset(UC.1,suggest=="yes"), color="red",size=4)+
 # geom_label(data=subset(UC.1,suggest=="yes"), color="black",aes(label=Gene.refGene),size=4)+
 # geom_label(data=subset(UC.1,suggest="yes"),color="black",aes(label=ID),size=4)+
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

name="UC"
#tiff(paste0("250812_",name,"_DIS.manhattan.tiff"), width=12, height=8, units="in", res=300)
tiff(paste0("250910_",name,"_DIS.no_annotated.UC576.HLA+.manhattan.tiff"), width=12, height=8, units="in", res=300)

print(UC.2)

dev.off()

write.table(UC.1,"250912_UC_HLA+_labelled.hybrid.",quote=FALSE,row.names=FALSE)

print("UC manhattan plot complete")
