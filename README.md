# Genome-wide Association and Genomic Prediction Study of Inflammatory Bowel Disease in a Large-Scale Korean Cohort

Analysis code accompanying the manuscript

> **Genome-wide Association and Genomic Prediction Study of Inflammatory Bowel Disease in a Large-Scale Korean Cohort**
> *Manuscript submitted to Genome Medicine.*

Author list, article DOI, and full citation will be added upon publication. Until then, please cite this repository using the metadata in [`CITATION.cff`](CITATION.cff) (Zenodo DOI: _to be added_).

---

## Repository structure

Directory names follow the section numbers of **Supplementary Methods** (Additional file 1), so each folder can be read alongside the corresponding parameter table.

| Directory | Supplementary Methods | Contents |
|---|---|---|
| [`scripts/S01_WGS_QC`](scripts/S01_WGS_QC) | S1. Sequence alignment, variant calling, and QC (WGS) | GATK BQSR → HaplotypeCaller (GVCF) → GenomicsDBImport → GenotypeGVCFs → hard filtering → SNP / biallelic SNP extraction → MAF / missingness filtering |
| [`scripts/S03_S04_CNV`](scripts/S03_S04_CNV) | S3. CNV calling and merging<br>S4. Genome-wide window-level CNV association | Delly / CNVkit call extraction, reciprocal-overlap merging (20–50%), depth and overdispersion diagnostics (Poisson GLM, permutation), 1-Mb window-level Firth logistic association at 30% / 40% thresholds, candidate tagging, rarity and database checks |
| [`scripts/S06_GWAS`](scripts/S06_GWAS) | S6. GWAS (susceptibility and treatment persistence) | PLINK2 genotype/sample missingness, HWE and MAF filters, PCA, Firth-fallback logistic regression, suggestive-SNP extraction and LD calculation |
| [`scripts/S08_MAGMA`](scripts/S08_MAGMA) | S8. Gene- and pathway-level analysis (MAGMA) | hg38 → hg19 liftover, SNP-to-gene annotation (`window=1,0.5`), gene-based tests (MAGMA v1.10, EAS reference), MSigDB C5 gene-set annotation |
| [`scripts/S09_finemapping`](scripts/S09_finemapping) | S9. Statistical fine-mapping | ±250 kb candidate windows, genotype and LD-matrix generation (PLINK2), SuSiE fine-mapping with individual-level data (susceptibility and treatment persistence), Manhattan plot |
| [`scripts/S10_eQTL`](scripts/S10_eQTL) | S10. Cis-eQTL analysis | MatrixEQTL cis-eQTL mapping (±1 Mb, FDR < 0.05) for CD and UC colon tissue |
| [`scripts/S12_PRS_ML`](scripts/S12_PRS_ML) | S12. Polygenic risk score and machine learning | PRSice-2 scoring with rank-based inverse normal transformation; logistic regression / random forest / SVM / gradient boosting models for the three locus sets |

Within each directory, scripts are numbered `stepNN_*` in execution order. Files suffixed `_CD`, `_UC`, `_CD_treat`, `_UC_treat` run the same step for Crohn's disease, ulcerative colitis, and the corresponding treatment-persistence phenotypes.

### Locus sets used in `S12_PRS_ML`

| Script | Locus set (Supplementary Methods S12) |
|---|---|
| `ml_{CD,UC}_IIBDGC227.py` | 227 of 320 IIBDGC loci present in the dataset (Supplementary Fig. S13A, C) |
| `ml_{CD,UC}_EAS_default.py`, `ml_{CD,UC}_EAS_tuned.py` | 22 CD / 27 UC East Asian–specific loci (Supplementary Fig. S13B, D) |
| `ml_CD_KoreanGWAS44.py`, `ml_UC_KoreanGWAS45.py` | 44 CD / 45 UC Korean GWAS-derived loci (Figure 5A, B) |

### Analyses not included in this release

Code for the following analyses described in the manuscript is not part of this release and may be added in a later version:
RNA-seq preprocessing (S2), population structure with 1000 Genomes (S5), HLA imputation with HIBAG (S7), metagenomic analysis (S11), and the deep neural network model (S12).

---

## Running the scripts

### Paths

All site-specific paths have been replaced with generic placeholders. Edit them to match your environment before running:

| Placeholder | Meaning |
|---|---|
| `/input/tools/` | Executables (e.g. `plink2`, `PRSice`, `PRSice.R`, `magma`) |
| `/input/data/` | Input data (genotypes, phenotypes, summary statistics, intermediate results) |
| `/input/reference/` | Reference resources (liftover chain files, dbSNP, 1000 Genomes EAS LD reference) |
| `/output/` | Output location |

### Cluster submission

Most shell scripts in `S01_WGS_QC`, `S06_GWAS`, and `S09_finemapping` were run on an SGE cluster via `qsub`. Queue name, notification e-mail, and slot count are read from environment variables:

```bash
QUEUE=your_cluster_queue NOTIFY_EMAIL=you@example.com SLOTS=12 \
  bash scripts/S01_WGS_QC/step00_bqsr.sh
```

If a variable is unset, a generic placeholder (`your_queue.q`, `your_email@example.com`) is used and job submission will fail until it is set.

---

## Software environment

### Languages

- **R** v4.6.0
- **Python** v3.10 — package list in [`requirements.txt`](requirements.txt)

### Command-line tools

Versions as reported in the manuscript; see Supplementary Methods (Additional file 1) for the complete list.

| Tool | Version |
|---|---|
| BWA-mem | 0.7.17 |
| samtools | 1.9 |
| GATK | 4.3.0.0 (VariantFiltration step: GATK 3.8) |
| vcftools | 0.1.17 |
| bcftools | 1.14 |
| SHAPEIT5 | 5.1.1 |
| ANNOVAR | 2020Jun08 |
| snpEff | 4.3t |
| CNVkit | 0.9.10 |
| Delly | 0.9.1 |
| SURVIVOR | 1.0.7 |
| bedtools | 2.31.1 |
| PLINK | 2.0 |
| MAGMA | 1.10 |
| PRSice-2 | — |

### R packages

| Package | Version | Used in |
|---|---|---|
| susieR | — | S09 |
| MatrixEQTL | — | S10 |
| edgeR | 4.2.2 | S10 |
| LDlinkR | 1.4.0 | S09 |
| msigdbr | — | S08 |
| rtracklayer, GenomicRanges | — | S08, S09 |
| biomaRt, Biobase, EnsDb.Hsapiens.v86 | — | S10 |
| bestNormalize, mclust, MASS, forecast | — | S12 |
| data.table, dplyr, tibble, ggplot2, ggrepel, ggpubr | — | all |

### Python packages

scikit-learn 1.6.1, pandas, numpy, scipy, statsmodels, firthlogist, intervaltree, matplotlib, seaborn, shap (see [`requirements.txt`](requirements.txt)).

---

## Data availability

Individual-level genomic and clinical data are **not** included in this repository. They are available under **controlled access** from the Korea BioData Station (K-BDS, <https://kbds.re.kr>) upon approval of a data access request:

| Data | Samples | Accession |
|---|---|---|
| Whole-genome sequencing | 3,818 | K-BDS [KAP230556](https://kbds.re.kr/service/utilization/view/KAP230556) |
| RNA-seq | 277 | K-BDS [KAP242592](https://kbds.re.kr/service/utilization/view/KAP242592) |
| Gut metagenome (reused from Kang et al., 2023) | — | NCBI SRA [PRJNA945504](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA945504); K-BDS [KAP220304](https://kbds.re.kr/service/utilization/view/KAP220304) (formerly PRJKA220304) |

Requests are subject to approval by the relevant data access committee.

---

## License

This code is released under the [MIT License](LICENSE).

## Contact

For questions about the code, please open an issue in this repository.
