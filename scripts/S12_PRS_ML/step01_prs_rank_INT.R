library(dplyr)
library(ggplot2)
library(data.table)
library(MASS)
library(forecast)
library(bestNormalize)
library(mclust)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1] #directory : "03.train_test_split/UC526KDC/UC526KDC.grp1.5e-08.only.best"
output <- args[2] #UC526KDC.grp1.5e-08.only.best
# 데이터 불러오기
Total <- read.table(input, header=TRUE)

# Group 변수 생성
Total["Group"] <- ifelse(startsWith(Total$IID, "KDC"), "0", "1")

## 1. Box-Cox 변환 (양수만 적용 가능하므로 음수는 shift 필요할 수 있음)
#shift_value_UC <- ifelse(min(Total$PRS) <= 0, abs(min(Total$PRS)) + 1e-5, 0)
#Total$PRS_shifted <- Total$PRS + shift_value_UC
#fit_UC <- lm(PRS_shifted ~ 1, data = Total)
#boxcox_result_UC <- boxcox(fit_UC, lambda = seq(-2, 2, 0.1), plotit = FALSE)
#lambda_opt_UC <- boxcox_result_UC$x[which.max(boxcox_result_UC$y)]
#Total$PRS_boxcox <- BoxCox(Total$PRS_shifted, lambda_opt_UC)

# CDTotal Box-Cox 적용 (마찬가지 shift)
#shift_value_CD <- ifelse(min(CDTotal$PRS) <= 0, abs(min(CDTotal$PRS)) + 1e-5, 0)
#CDTotal$PRS_shifted <- CDTotal$PRS + shift_value_CD
#fit_CD <- lm(PRS_shifted ~ 1, data = CDTotal)
#boxcox_result_CD <- boxcox(fit_CD, lambda = seq(-2, 2, 0.1), plotit = FALSE)
#lambda_opt_CD <- boxcox_result_CD$x[which.max(boxcox_result_CD$y)]
#CDTotal$PRS_boxcox <- BoxCox(CDTotal$PRS_shifted, lambda_opt_CD)

# 2. Yeo-Johnson 변환 (음수도 처리 가능)
yj_obj_UC <- yeojohnson(Total$PRS)
Total$PRS_yeojohnson <- yj_obj_UC$x.t

# 3. GMM 클러스터링 (원본 PRS 기준)
gmm_UC <- Mclust(Total$PRS)
Total$cluster_gmm <- gmm_UC$classification

# 4. INT 변환 함수 선언 및 적용
prs_int_transform <- function(df, prs_col = "PRS") {
  df %>%
    mutate(!!paste0(prs_col, "_INT") := qnorm((rank(.data[[prs_col]], na.last = "keep") - 0.5) / sum(!is.na(.data[[prs_col]]))))
}
Total <- prs_int_transform(Total, "PRS")
write.table(Total,paste0(output,".best.tsv"),quote=FALSE)

# 5. t-test 함수 (PRS 컬럼명 인자 지원)
t_return <- function(genotype, prs_col="PRS") {
  con <- genotype[genotype$Group == "0", prs_col]
  cas <- genotype[genotype$Group == "1", prs_col]
  tt <- t.test(con, cas)
  print(tt)
  return(tt)
}

# 6. 그래프 함수 (PRS 컬럼명, x축 라벨 인자 지원)
GRS_ggplot <- function(data, tt, prs_col="PRS", x_label="Genetic Risk Score (GRS)") {
  con <- data[data$Group=="0", prs_col]
  cas <- data[data$Group=="1", prs_col]
  tt <- t.test(con, cas)  # 다시 수행하여 정확한 결과 삽입
  data <- as.data.frame(data)
  ggplot(data, aes_string(x=prs_col, fill="Group")) +
    geom_density(alpha=0.5) +
    scale_fill_manual(values=c("1"="red", "0"="blue")) +
    scale_color_manual(values=c("1"="red", "0"="blue")) +
    theme_classic() +
    labs(
      title="GRS Density Plot by Case/Control",
      subtitle=sprintf("t-test results: t = %.2f, df = %.2f, p-value = %.3e",
                       tt$statistic, tt$parameter, tt$p.value),
      x = x_label,
      y = "Density",
      fill = "Group", color = "Group"
    )
}

# 7. 그래프 저장 함수 (tiff 형식, 해상도 300)
GRS_ggplot_figure <- function(g, output_prefix) {
  tiff(paste0(output_prefix, "Total_PRS.tiff"), width=12, height=8, units="in", res=300)
  print(g)
  dev.off()
}

# 8. 각 변환 컬럼별로 t-test, 그래프 생성, 저장 실행 예시

prs_cols <- c("PRS", "PRS_yeojohnson", "PRS_INT")

for (prs_col in prs_cols) {
  # Total 처리
  tt_uc <- t_return(Total, prs_col)
  g_uc <- GRS_ggplot(Total, tt_uc, prs_col, paste0("Genetic Risk Score (", prs_col, ") - UC"))
  GRS_ggplot_figure(g_uc, paste0(output, prs_col, "_"))
  
}

# 9. GMM 클러스터링 결과 시각화 예 (원본 PRS 기준)
g_gmm<-ggplot(Total, aes(x=PRS, fill=factor(cluster_gmm))) +
  geom_density(alpha=0.5) +
  labs(title="GMM Clustering of PRS - UC", fill="Cluster") +
  theme_minimal()

#GRS_ggplot_figure(g_gmm, "UC526KDC.grp1_5e08.250925.best")
#GRS_ggplot_figure(g_gmm_cd, "CD238KDC.grp5_5e08.250925.best")
