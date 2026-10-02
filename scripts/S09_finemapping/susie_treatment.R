#00.importing library=======================
library(data.table)
library(susieR)
library(dplyr)


#01. importing dataset======================
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
id_cs <- sub(".*(rs[0-9]+).*", "\\1", basename(filename))
geno<-fread(filename)
pheno<-fread(args[2])
output_file<-args[3]
ld<-data.table::fread(args[4])


#geno<-fread("rs55951892.out.raw" or "rs56167332.UC.out.raw")
#pheno<-fread("CDKDC.250120.03.pca.pheno.txt" or "UC526KDC.pheno.txt")
#ld_check<-fread("rs55951892.ld.vcor" or "rs56167332.UC.out.check.vcor")

#02.QC with ld_check

snps_with_ld <- unique(c(
  ld[PHASED_R2 >= 0.1, ID_A],
  ld[PHASED_R2 >= 0.1, ID_B]
))

length(snps_with_ld)

snps <- unique(c(ld$ID_A, ld$ID_B))
R <- matrix(0, length(snps), length(snps))
rownames(R) <- colnames(R) <- snps

idxA <- match(ld$ID_A, snps)
idxB <- match(ld$ID_B, snps)

R[cbind(idxA, idxB)] <- ld$PHASED_R2
R[cbind(idxB, idxA)] <- ld$PHASED_R2
diag(R) <- 1


R_filtered <- R
R_filtered[abs(R_filtered) < 0.1] <- 0 #R<=0.1 ignored

ld_partner_count <- apply(abs(R) >= 0.1, 1, sum) - 1
summary(ld_partner_count) # if snp< nn, LD region normal // snp > nnn - long or complex LD region

cut(
  abs(R[upper.tri(R)]),
  breaks = c(0, 0.05, 0.1, 0.3, 0.6, 0.9, 1),
  include.lowest = TRUE
) |> table() # distributed 0.1-0.6 -- normal // mass would be in 0.6> ~ long or complex LD region


#03. matrix generating and phenotype making ==

X<-as.matrix(geno[,7:ncol(geno)])
colnames(X)<-colnames(geno)[7:ncol(geno)]
y<-pheno

# NA → 평균 대치
X <- apply(X, 2, function(v) {
  v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})

print(all(geno$IID == pheno$IID))
maf <- colMeans(X) / 2
X <- X[, maf > 0.01 & maf < 0.99]

#y.1<-as.numeric(ifelse(y$Pheno==2,1,0))

y <- y %>%
  mutate(TREATMENT = as.numeric(TREATMENT) - 1)

y_vec<-y$TREATMENT

#y$DIS in CD and y$Pheno in UC

# 04. susie_making with credibleset ==========

fit<-susie(X,y_vec,L = 3,                 # causal variant 수 (보통 5–15) -> sample size가 작으면 1-3
  standardize = TRUE,     # ⭐ 매우 중요
  intercept = TRUE,
  estimate_residual_variance = TRUE,
  scaled_prior_variance = 0.02,  # case-control에서는 보수적으로, 0.1-0.3사이에서 만족되면 됌 / linear regresion은 sample이 120-150, 0.01-0.05
  tol = 1e-4,
  max_iter = 300
)
print(fit)


if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {

  # ✅ CS 존재: credible set 기반 결과 생성
  cs_df <- bind_rows(lapply(seq_along(fit$sets$cs), function(l) {
    cs <- fit$sets$cs[[l]]
    tibble(
      cs_index = l,
      rsid = colnames(X)[cs],
      cs_pip = fit$pip[cs]
    )
  }))

  message("Credible sets detected")

} else {

  # ⚠️ CS 없음: PIP-only 결과
  cs_df <- tibble(
    rsid = colnames(X),
    pip = fit$pip
  )

  message("No credible sets detected; using PIP only")

}

#05.writing results=====================
write.table(cs_df, paste0(output_file, id_cs , ".genopheno.susie.tsv" ),quote=FALSE, row.names=FALSE)
if (length(fit$sets$cs) == 0) {
  write.table(
    data.frame(
      lead_id = id_cs,
      status = "no_credible_set"
    ),
    paste0(output_file, id_cs, ".genopheno.susie.tsv"),
    quote = FALSE,
    row.names = FALSE
  )
  quit(save = "no")
}

#output = CD_causal/ or UC_causal/
#id_cs = extracted in filename, rsXXXXX
