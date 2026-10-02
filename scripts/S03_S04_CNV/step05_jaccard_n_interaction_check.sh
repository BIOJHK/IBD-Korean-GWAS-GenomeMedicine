#!/usr/bin/env bash
# Jaccard/interaction check between the IBD and control (KDC) CNVkit target
# BED files, to confirm both call sets cover the same genomic regions before
# CNV burden comparison.
#
# Reconstructed from the analyst's terminal history: this previously
# contained a leaked internal server hostname/IP/SSH port/account used to
# fetch the BED files, plus captured rsync/conda-install/sort output that
# isn't part of the actual analysis. Only the real commands are kept below;
# fetch step02.IBD_targets.bed and step02.KDC_targets.bed from wherever you
# store your own CNVkit target BED outputs, then run the rest as-is.
# Requires bedtools (v2.31.1, matching Methods).

# rsync -avzr --rsh="ssh -p <PORT>" <user>@<host>:<path-to-cnvkit-output>/step02*bed ./

sed -i 's/^chr//' step02.IBD_targets.bed
# step02.KDC_targets.bed had no "chr" prefix in the original session (not shown being stripped here)

bedtools jaccard \
  -a <(sort -k1,1 -k2,2n step02.IBD_targets.bed) \
  -b <(sort -k1,1 -k2,2n step02.KDC_targets.bed)
