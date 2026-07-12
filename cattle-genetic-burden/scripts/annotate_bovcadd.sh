#!/bin/bash
set -euo pipefail

module load bcftools htslib

VCF_IN=$1        # e.g. results/merged/chr1.merged.filtered.vcf.gz
BOVCADD=$2       # e.g. resources/bovcadd/bovcadd_scores.tsv.gz
VCF_OUT=$3       # e.g. results/annotated/chr1.bovcadd.vcf.gz

# Header line to declare the new INFO field
cat > /tmp/bovcadd_header.hdr << 'EOF'
##INFO=<ID=BovCADD,Number=1,Type=Float,Description="BovCADD PHRED-scaled deleteriousness score">
EOF

bcftools annotate \
    -a "$BOVCADD" \
    -h /tmp/bovcadd_header.hdr \
    -c CHROM,POS,REF,ALT,-,BovCADD \
    -Oz -o "$VCF_OUT" \
    "$VCF_IN"

bcftools index -f "$VCF_OUT"
