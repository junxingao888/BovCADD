"""
 Module to generate whole genome PHRED-like CADD scores.
 First all possible variants are generated, these are then annotated and scored
 using the annotate_variants and train_test_model modules of the workflow.
 Afterwards the raw model scores are sorted in descending order and are
 assigned the final PHRED-like CADD score. These are then again sorted by
 chromosome and position. Finally, the scores for each position
 are summarised by a mean, min and max value per position.

 :Author: Job van Schipstal
 :Date: 19-06-2024

 The scripts and workflow have been adopted from the work of Christian Gross.
"""

CADD_P = "../scripts/cadd_score/"

"""
Generates all possible variants for a chromosome in blocks.
The files are also directly bgzipped.
"""
rule generate_all_variants:
    input:
         reference=config["generate_variants"]["reference_genome_wildcard"],
         script=workflow.source_path(CADD_P + "create_variants.py")
    params:
         blocksize=config["parallelization"]["whole_genome_positions_per_file"]
    conda:
         "../envs/mainpython.yml"
    priority: -10
    log:
        "results/whole_genome_variants/chr{chr}/stats.txt"
    output:
         out_dir=directory("results/whole_genome_variants/chr{chr}"),
    shell:
         """
         python3 {input.script} -o {output.out_dir} \
         -s {params.blocksize} -r {input.reference} \
         -c {wildcards.chr} > {log} && \
         for file in {output.out_dir}/*.vcf
         do
           bgzip "$file"
         done
         """

"""
Combines all log files for variant generation for all chromosomes.
This rule also serves as a checkpoint, after which it is checked how 
many blocks were generated for each chromosome and thus how many runs 
of the annotation and scoring modules will be needed.
"""
checkpoint summarize_generation:
    input:
         expand("results/whole_genome_variants/chr{chr}/stats.txt",
                chr=CHROMS)
    output:
         report("results/logs/whole_genome_variants.txt", category="Logs")
    priority: -10
    shell:
         "tail -n +2 {input} > {output}"


def get_all_score_files(wildcards):
    """
    Fetches all input files and requests the annotated versions as input.
    The exact number of input files is not known before execution of the pipeline,
    hence a checkpoint is called. The output of this function, and thus the
    input of the sort_raw_scores is only determined after summarize_generation.
    :param wildcards: Any, Wildcards of rule, unused
    :return: list of str, input files
    """
    checkpoints.summarize_generation.get()
    globed = glob_wildcards(f"results/whole_genome_variants/chr{wildcards.chr}/chr{wildcards.chr}_{{part}}.vcf.gz")
    return expand(f"results/whole_genome_scores/raw_parts/All/chr{wildcards.chr}_{{part}}.csv",
                  part=globed.part)

"""
Sort raw scores based on predicted probability to be deleterious (descending).
in this step we merge sort the part files into one file for the whole chromosome.
"""
rule sort_raw_scores:
    input:
         get_all_score_files
    threads: 8
    resources:
        mem_mb=200000
    output:
         "results/whole_genome_scores/RAW_scores_chr{chr}.csv"
    shell:
        """
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
         {input} > {output}
        """

"""
Count number of lines in file. Since there is no header and it is one line per variant that is the number of variants.
Needed for PHRED-like scaling of scores.
"""
rule count_chrom_variants:
    input:
        "results/whole_genome_scores/RAW_scores_chr{chr}.csv"
    output:
        "results/whole_genome_scores/counts/chr{chr}.txt"
    shell:
        """wc -l {input} > {output}"""

"""
Sort raw scores based on predicted probability to be deleterious (descending)
in this step we merge sort the pre sorted chromosome files together
"""
rule merge_raw_scores:
    input:
         expand("results/whole_genome_scores/RAW_scores_chr{chr}.csv",
         chr=CHROMS)
    threads: 8
    resources:
        mem_mb=200000,
        tmpdir="results/whole_genome_tmp"
    output:
         "results/whole_genome_scores/full_RAW_scores.csv"
    shell:
        """
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
         {input} > {output}
        """

# Old shell
# head -1 {input[0]} >> {output} && \
# tail -n +2 --quiet {input} | \
# LC_ALL=C sort \
# -t "," \
# -k5gr \
# -S {resources.mem_mb}M \
# --parallel={threads} >> {output}

"""
Take file with all variants, sorted by raw score and phred scale them.
outputs scaled variants in a separate file for each chromosome (still sorted by score).
"""
rule assign_phred_scores:
    input:
        data="results/whole_genome_scores/full_RAW_scores.csv",
        counts=expand("results/whole_genome_scores/counts/chr{chr}.txt",
                      chr=CHROMS),
        script=workflow.source_path(CADD_P + "assign_phred_score.py")
    params:
        outmask="results/whole_genome_scores/phred/chrCHROM.tsv",
        chromosomes=CHROMS,

    output:
        expand("results/whole_genome_scores/phred/chr{chr}.tsv",
               chr=CHROMS)
    shell:
        """python3 {input.script} \
        -i {input.data} \
        -o {params.outmask} \
        --chroms {params.chromosomes} \
        --count-file {input.counts}"""

"""
Sort PHRED score variants of one chromosome by position
"""
rule sort_phred:
    input:
        "results/whole_genome_scores/phred/chr{chr}.tsv"
    threads: 4
    resources:
        mem_mb=100000
    params:
        tmpdir="results/tmp/chr{chr}"
    conda:
        "../envs/mainpython.yml"
    output:
        "results/cadd_scores/chr{chr}.tsv.gz"
    shell:
        """
        mkdir -p {params.tmpdir}
        tail -n +2 {input} | \
        LC_ALL=C sort \
            -k2n \
            -S {resources.mem_mb}M \
            --parallel={threads} | \
        bgzip -c > {output}
        rm -rf {params.tmpdir}
        """

"""
Generate table with the count of each VEP annotation class by CADD score bins (e.g. 0 to 1, 1 to ...)
"""
rule cadd_consequence_bins:
    input:
        data="results/whole_genome_variants/chr{chr}/chr{chr}_{part}_full_annotation.tsv",
        annotaton="results/whole_genome_scores/phred/chr{chr}.tsv",
        script=workflow.source_path(CADD_P + "consequence_bin.py")
    conda:
         "../envs/mainpython.yml"
    output:
        "results/consequence_bins/chr{chr}_{part}.csv"
    shell:
        """python3 {input.script} \
        -i {input.data} \
        -a {input.annotaton} \
        -o {output}"""
