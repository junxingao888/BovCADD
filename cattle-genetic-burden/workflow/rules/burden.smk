rule annotate_bovcadd:
    input:
        vcf="results/merged/{chrom}.merged.filtered.vcf.gz",
        bovcadd="resources/bovcadd/bovcadd_scores.tsv.gz"
    output:
        "results/annotated/{chrom}.bovcadd.vcf.gz"
    shell:
        "bash scripts/annotate_bovcadd.sh {input.vcf} {input.bovcadd} {output}"

rule hom_burden:
    input:
        "results/annotated/{chrom}.bovcadd.vcf.gz"
    output:
        "results/burden/{chrom}.hom_burden.csv"
    shell:
        """
        python scripts/compute_hom_burden.py \
            --vcf {input} --out {output} \
            --del-threshold 20 --benign-threshold 5
        """

rule merge_burden:
    input:
        expand("results/burden/{chrom}.hom_burden.csv", chrom=CHROMS)
    output:
        "results/burden/genomewide_hom_burden.csv"
    script:
        "../../scripts/merge_chrom_burden.py"

rule population_stats:
    input:
        "results/burden/genomewide_hom_burden.csv"
    output:
        "results/burden_summary.csv",
        "results/plots/burden_by_breed.png"
    script:
        "../../scripts/population_stats.py"
