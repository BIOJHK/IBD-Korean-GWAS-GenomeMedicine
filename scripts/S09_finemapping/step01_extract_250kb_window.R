library(dplyr)
library(data.table)

CD_hy<-data.table::fread("250912_CD_HLA+_labelled.hybrid.")
UC_hy<-data.table::fread("250912_UC_HLA+_labelled.hybrid.")

twomilon_extract <- function(rsid, hy){
  window_size_bp <- 5e5 / 2  # ±250kb = 총 500kb 윈도우
  # 또는 명확히: window_size_bp <- 250000  # 250kb
  half_window_bp <- window_size_bp / 2  # 125kb (focal SNP 기준 앞뒤)
  
  results <- lapply(rsid, function(snp) {
    idx <- which(hy[["ID"]] == snp)
    if (length(idx) == 0) return(NULL)
    focal_bpcum <- hy[["BPcum"]][idx]
    start_bpcum <- focal_bpcum - half_window_bp
    end_bpcum <- focal_bpcum + half_window_bp
    subset(hy, BPcum >= start_bpcum & BPcum <= end_bpcum)
  })

CD_c <- CD_hy %>% filter(P<=1e-05)
UC_c <- UC_hy %>% filter(P<=1e-05)

CD_list<-list()
UC_list<-list()

for (i in 1:nrow(CD_c)){
	CD_list[[i]]<-twomilon_extract(CD_c$ID[i],CD_hy)
	write.table(CD_list,paste0(CD_c$ID,"_250kb.txt"),quote=FALSE)
}

for (i in 1:nrow(UC_c)){
	UC_list[[i]]<-twomilon_extract(UC_c$ID[i],UC_hy)
	write.table(UC_list,paste0(UC_c$ID,"_250kb.txt"),quote=FALSE)
}


#print(twomilon_extract(CD_hy.causal$ID[1]))
