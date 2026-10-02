#!/usr/bin/env bash
# MAGMA gene-based analysis commands (Windows build of MAGMA v1.10).
# See Methods: "Gene-based analysis using multiple regression adjusted for
# LD ... MAGMA (Multi-marker Analysis of Genomic Annotation v1.10) ...
# preprocessing parameters --annotate window=1, 0.5."
#
# Reconstructed from the analyst's terminal history: this previously
# contained the full captured stdout/log of each MAGMA run (gene counts,
# warnings, etc., specific to that machine) along with a leaked local
# username/hostname/path in each command's shell prompt. Only the actual
# commands are kept below, in the order they were run; rerun MAGMA to
# regenerate its own logs.

# --- CD, disease susceptibility ---
/input/tools/magma --annotate window=1,0.5 --snp-loc 251223_CD_lift_gwas.txt --gene-loc NCBI37.3/NCBI37.3.gene.loc --out 251223_CD_lift_gwas
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_CD_lift_gwas.genes.annot.txt --pval 251223_CD_lifP_gwas_P.txt N=78 --gene-model multi --out 251223_CD_lifPanno

# --- UC, disease susceptibility ---
/input/tools/magma --annotate window=1,0.5 --snp-loc 251223_UC_lift_gwas.txt --gene-loc NCBI37.3/NCBI37.3.gene.loc --out 251223_UC_liftgwas
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_UC_liftgwas.genes.annot.txt --pval 251223_UC_lifP_gwas_P.txt N=342 --gene-model multi --out 251223_UC_lifPanno

# --- CD, treatment persistence ---
/input/tools/magma --annotate window=1,0.5 --snp-loc 251223_CD_treat_lift_gwas.txt --gene-loc NCBI37.3/NCBI37.3.gene.loc --out 251223_CD_treatgwas
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_CD_treatgwas.genes.annot.txt --pval 251223_CD_treat_lifP_gwas_P.txt N=61 --gene-model multi --out 251223_CD_treat_lifPanno

# --- UC, treatment persistence (re-run twice in the original session) ---
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_UC_treatgwas.genes.annot.txt --pval 251223_UC_treat_lifP_gwas_P.txt N=124 --gene-model multi --out 251223_UC_treat_lifPanno
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_UC_treatgwas.genes.annot.txt --pval 251223_UC_treat_lifP_gwas_P.txt N=124 --gene-model multi --out 251223_UC_treat_lifPanno

# --- CD, treatment persistence (different N) ---
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_CD_treatgwas.genes.annot.txt --pval 251223_CD_treat_lifP_gwas_P.txt N=152 --gene-model multi --out 251223_CD_treat_lifPanno

# --- UC, disease susceptibility (full cohort N) ---
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_UC_liftgwas.genes.annot.txt --pval 251223_UC_lifP_gwas_P.txt N=3510 --gene-model multi --out 251223_UC_lifPanno

# --- CD, disease susceptibility (full cohort N) ---
/input/tools/magma --bfile /input/reference/g1000_eas --gene-annot 251223_CD_lift_gwas.genes.annot.txt --pval 251223_CD_lifP_gwas_P.txt N=3238 --gene-model multi --out 251223_CD_lifPanno
