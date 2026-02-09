"""
 Module that simulates variants based on the mutation rates found
 between the reference genome and the ancestral genome.
 An equal number of variants will be generated.

 :Author: Job van Schipstal
 :Date: 4-10-2023

 Based upon the work of Seyan Hu.
"""

ANNOTATE_P = "../scripts/annotate_variants/"
CONVERSION_P = "../scripts/conversion_tools/"

rule basic_annotation:
    input:
         vcf="{folder}/{file}.vcf.gz",
         ref_genomes=expand(config["generate_variants"]["reference_genome_wildcard"],
                            chr=CHROMS),
         idx=expand(config["generate_variants"]["reference_genome_wildcard"] + ".fai",
                            chr=CHROMS),
         shape_file=config["basic_annotation"]["shape_file"] \
             if config["basic_annotation"]["shape_file"] != "None" else [],
         script=workflow.source_path(ANNOTATE_P + "basic_annotation.py")
    params:
         include_masked="--include-masked" if config["basic_annotation"]["include_masked"] else " ",
         shape=config["basic_annotation"]["shape_file"],
         ref_genome=config["generate_variants"]["reference_genome_wildcard"].replace("{chr}","[CHROM]"),
    conda:
         "../envs/mainpython.yml"
    output:
          "{folder}/{file}.basic.tsv"
    shell:
         """python3 {input.script} \
         -i {input.vcf} \
         -s {params.shape} \
         {params.include_masked} \
         -r {params.ref_genome} \
         -o {output}"""


"""
Optional rule that installs the needed VEP cache using the vep_install tool
included with VEP. It is used with the -n no update flag and a set version for reproducibility.
The VEP cache and program should be from the same release, hence care should be taken to update them together.
"""
rule vep_cache:
    params:
          version_species=config["vep"]["cache"]["install_params"]
    conda:
         "vep"
    output:
          directory(config["vep"]["cache"]["directory"])
    shell:
         "vep_install -a cf -n {params.version_species} -c {output} --CONVERT"

"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""  # TODO make sift a config option
rule run_vep:
    input:
         vcf="{folder}/{file}.vcf.gz",
         script=workflow.source_path(ANNOTATE_P + "vep.sh"),
         cache=rules.vep_cache.output if
            config["vep"]["cache"]["should_install"] == "True" else []
    params:
          cache_dir=config["vep"]["cache"]["directory"],
          species_name=config["species_name"]
    conda:
         "../envs/ensembl.yml"
    # Parts are at most a few million variants, 2 threads is already fast.
    threads: 2
    output:
          "{folder}/{file}_vep_output.tsv"
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.vcf} {output} "
         "{params.cache_dir} {params.species_name} {threads} && "
         "[[ -s {output} ]]"

"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""
rule process_vep:
    input:
         vcf="{folder}/chr{chr}.vcf.gz",
         index="{folder}/chr{chr}.vcf.gz.tbi",
         vep="{folder}/chr{chr}_vep_output.tsv",
         exons=config["vep"]["exons"],
         grantham=workflow.source_path(ANNOTATE_P + "grantham_matrix.tsv"),
         script=workflow.source_path(ANNOTATE_P + "VEP_processing.py")
    params:
        multiple = " "
    conda:
         "../envs/mainpython.yml"
    output:
          "{folder}/chr{chr}.vep.tsv"
    shell:
         """python3 {input.script} \
         -v {input.vep} \
         -s {input.vcf} \
         -e {input.exons} \
         -g {input.grantham} \
         {params.multiple} \
         -o {output}"""


"""
Process whole genome variant VEP annotation with a lower priority.
If the model is trained before the full genome variants have been fully processed,
it can continuously schedule tasks without having to wait on the model training.
"""
use rule process_vep as process_genome_vep with:
    input:
         vcf="results/whole_genome_variants/{batch}/{file}.vcf.gz",
         index="results/whole_genome_variants/{batch}/{file}.vcf.gz.tbi",
         vep="results/whole_genome_variants/{batch}/{file}_vep_output.tsv",
         exons=config["vep"]["exons"],
         grantham=workflow.source_path(ANNOTATE_P + "grantham_matrix.tsv"),
         script=workflow.source_path(ANNOTATE_P + "VEP_processing.py"),
    params:
        multiple = "--multiple"
    priority: -8
    output:
          "results/whole_genome_variants/{batch}/{file}.vep.tsv"

"""
Processes VEP output for the validation sets.
Does the process for all chromosomes, sequentially since few variants are expected here
"""
use rule process_vep as validation_vep with:
    input:
         vcf="results/validation_variants/{file}.vcf.gz",
         index="results/validation_variants/{file}.vcf.gz.tbi",
         vep="results/validation_variants/{file}_vep_output.tsv",
         exons=config["vep"]["exons"],
         grantham=workflow.source_path(ANNOTATE_P + "grantham_matrix.tsv"),
         script=workflow.source_path(ANNOTATE_P + "VEP_processing.py"),
    output:
          "results/validation_variants/{file}.vep.tsv"


rule convert2bed:
    input:
         infile=lambda wildcards: config["bed_annotation"][wildcards.label][
             "file"],
         wig_script=workflow.source_path(
             CONVERSION_P + "fixStepToBedGraph.pl"),
         bigwig_script=workflow.source_path(CONVERSION_P + "bigWigToBedGraph")
    params:
          file_type=lambda wildcards: config["bed_annotation"][
              wildcards.label]["file"].replace(".gz", "").split(".")[-1]
    threads: 8
    output:
          temp("results/annotations/{label}.converted.bed")
    shell:
         """
         infile={input.infile}
         outfile={output}
         if [ {params.file_type} == "wig" ]; then
             if gzip -t $infile; then
                 gzip -dc $infile | perl {input.wig_script} > $outfile
             else
                 cat $infile | perl {input.wig_script} > $outfile 
             fi
         elif [ {params.file_type} == "bw" ]; then
             chmod +x {input.bigwig_script}
             {input.bigwig_script} $infile $outfile
         fi
         """


"""
Prepares additional annotation files for use with BCFtools annotate.
We convert to a bedgraph, which is #Chrom to from score1 (score2).
Annotations are bgzipped and indexed to enable accessing only the specific blocks needed.
"""
rule prepare_bed_annotation:
    input:
        infile=lambda wildcards: config["bed_annotation"][wildcards.label]["file"]
            if config["bed_annotation"][wildcards.label]["file"].replace(".gz", "").split(".")[-1] == "bed"
            else f"results/annotations/{wildcards.label}.converted.bed"
    params:
        loading_arg=lambda wildcards: "gzip -dc " if config["bed_annotation"][wildcards.label]["file"].endswith("bed.gz") else "cat ",
        sed_optional=lambda wildcards: "sed 's/^chr//'" if config["bed_annotation"][wildcards.label].get("strip_chr", False) else "cat",
        sort_optional=lambda wildcards: "LC_ALL=C sort -k1,1V -k2,2n -k3,3n |" if config["bed_annotation"][wildcards.label]["should_sort"] else ""
    conda:
        "../envs/mainpython.yml" # Contains bgzip
    threads: 16
    output:
        "results/annotations/{label}.bed.gz"
    shell:
        """
        {params.loading_arg} {input} | {params.sed_optional} | {params.sort_optional} bgzip -c > {output}
        """

rule index_annotation_file:
    input:
         "results/annotations/{label}.bed.gz"
#    params:
#          is_0_based=lambda wildcards: "-0 " if
#            config["bed_annotation"][wildcards.label]["coord_base"] == 0 else " "
    conda:
         "../envs/mainpython.yml" # Contains tabix
    output:
          "results/annotations/{label}.bed.gz.tbi"
    shell:
        "tabix -s1 -b2 -e3 -0 -f {input}"



"""
Use BCFtools to annotate vcf file.
The python wrapper directly extracts the found score and outputs this as a tsv file.
It also saves coverage statistics to a log file.
"""
def format_annotations(wildcards):
    formatted_annotations = ""
    annotations = config["bed_annotation"][wildcards.label]["annotations"]
    for name, row in annotations.items():
        formatted_annotations += f"{row}={name} "
    return formatted_annotations

"""
rule bed_annotation:
    input:
         vcf="{folder}/{file}.vcf.gz",
         bed="results/annotations/{label}.bed.gz",
         index="results/annotations/{label}.bed.gz.tbi",
         script=workflow.source_path(ANNOTATE_P + "bed_annotate_vcf.py"),
    params:
         annotations=format_annotations,
         is_0_based= lambda wildcards: "--zero-based " if
            config["bed_annotation"][wildcards.label]["coord_base"] == 0 else " "
    conda:
         "../envs/mainpython.yml"
    output:
        "{folder}/{file}.{label}.bed.tsv"  # TODO make temp
    shell:
         python3 {input.script} \
         -i {input.vcf} \
         -b {input.bed} \
         -a {params.annotations} \
         {params.is_0_based} \
         -o {output}
"""

"""
Use BCFtools to annotate vcf file.
The python wrapper directly extracts the found score and outputs this as a tsv file.
It also saves coverage statistics to a log file.
"""
rule vcf_annotation:
    input:
         vcf="{folder}/{file}.vcf.gz",
         vcf_index="{folder}/{file}.vcf.gz.tbi",
         bed=config["vcf_annotation"]["file"],
         script=workflow.source_path(ANNOTATE_P + "vcf_annotation_wrapper.py"),
    params:
         annotations=" ".join(config["vcf_annotation"]["labels"])
    conda:
         "../envs/mainpython.yml"
    log:
       "results/logs/{folder}/{file}.vcf_annotation.log"
    output:
       "{folder}/{file}.vcf.tsv"
    shell:
         "python3 {input.script} -i {input.vcf} -a {input.bed} "
         "-l {params.annotations} --logfile {log} -o {output}"

"""
Use pysam to determine minimum distances to transcription start and end sites.
Requires pre-treated gtf to lookup transcription sites.
Assumes chromosomes do not have a prefix, e.g. chromosome 19 is just '19'
set --prefix to what the prefix is or just remove from the file.
"""
rule transcript_annotation:
    input:
         vcf="{folder}/{file}.vcf.gz",
         vcf_index="{folder}/{file}.vcf.gz.tbi",
         gtf=config["transcript_annotation"]["transcript_file"],
         script=workflow.source_path(ANNOTATE_P + "annotate_transcription_sites.py")
    conda:
         "../envs/mainpython.yml"
    output:
        "{folder}/{file}.transcript.tsv"
    shell:
         """python3 {input.script} \
         -i {input.vcf} \
         -a {input.gtf} \
         -o {output}"""


"""
Use BCFtools to annotate vcf file.
The python wrapper directly extracts the found score and outputs this as a tsv file.
It also saves coverage statistics to a log file.
"""
rule bed_annotation:
    input:
         vcf="{folder}/{file}.vcf.gz",
         vcf_index="{folder}/{file}.vcf.gz.tbi",
         bed="results/annotations/{label}.bed.gz",
         bed_index="results/annotations/{label}.bed.gz.tbi",
         script=workflow.source_path(ANNOTATE_P + "bed_annotation_wrapper.py"),
    params:
         annotations=format_annotations
    conda:
         "../envs/mainpython.yml"
    log:
       "results/logs/{folder}/{file}.{label}.annotation.log"
    output:
       "{folder}/{file}.{label}.bed.tsv"
    shell:
         "python3 {input.script} -i {input.vcf} -a {input.bed} "
         "-l {params.annotations} --logfile {log} -o {output}"


"""
Merges all annotations for a chromosome(part) into a single file
paste adds files line by line, separated with a tab by default.
Since we make use of tab separated format this yields a final tsv.
It is assumed that all tsv files are of same length, since they originate from the same vcf.
"""
rule merge_annotations:
    input:
         basic="{folder}/{file}.basic.tsv",
         vep="{folder}/{file}.vep.tsv",
         vcf="{folder}/{file}.vcf.tsv"
            if config["vcf_annotation"]["enabled"] else [],
         trans="{folder}/{file}.transcript.tsv"
            if config["transcript_annotation"]["enabled"] else [],
         bed=expand("{{folder}}/{{file}}.{label}.bed.tsv",
            label=list(config["bed_annotation"].keys()))
    output:
          "{folder}/{file}_full_annotation.tsv"
    shell:
         """paste \
         {input.basic} \
         {input.vep} \
         {input.vcf} \
         {input.trans} \
         {input.bed} > {output}"""

"""
No longer done in a separate step but during generation of dataset
rule merge_chrom_annotations:
    input:
        lambda wildcards: expand(
        "results/simulated_variants/trimmed_snps/chr{chr}_full_annotation.tsv",
        chr=CHROMS)
        if wildcards.variant == "simulated" else expand(
        "results/derived_variants/singletons/chr{chr}_full_annotation.tsv",
        chr=CHROMS)
    output:
          "results/dataset/unprocessed/{variant}_snps.tsv"
    wildcard_constraints:
         variant="derived|simulated"
    shell:
         "head -1 {input[0]} >> {output} &&"
         "tail --quiet -n +2 {input} >> {output}"
"""

"""
Missing values are imputed based on the mean of all simulated variants 
for some annotations, e.g. GC or GERP score. This script gathers all relevant 
columns for the different chromosomes and outputs their means to a text file 
so the next steps can be performed in parallel, for derived and simulated on 
each chromosome, using the already determined means. Scaling also has to be 
done centrally, so this is delayed till after the parallel processing.
"""
rule derive_impute_means:
    input:
        tsv=lambda wildcards: expand(
        "results/simulated_variants/trimmed_snps/chr{chr}_full_annotation.tsv",
        chr=CHROMS),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(ANNOTATE_P + "derive_means.py"),
    conda:
         "../envs/mainpython.yml"
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs")
    shell:
        "python3 {input.script} -i {input.tsv} "
        "--processing-config {input.processing} -o {output}"


rule column_analysis:
    input:
         derived=expand("results/derived_variants/singletons/chr{chr}_full_annotation.tsv",
                        chr=CHROMS),
         simulated=expand("results/simulated_variants/trimmed_snps/chr{chr}_full_annotation.tsv",
                          chr=CHROMS),
         script=workflow.source_path(ANNOTATE_P + "column_analysis.py")
    conda:
         "../envs/mainpython.yml"
    params:
         out_folder="results/figures/column_analysis/"
    output:
         relevance=report("results/figures/column_analysis/relevance.tsv",
                          category="Column Analysis"),
         derived_cor=report("results/figures/column_analysis/derived_variants_corr.tsv",
                            category="Column Analysis"),
         simulated_cor=report("results/figures/column_analysis/simulated_variants_corr.tsv",
                              category="Column Analysis"),
         combined_cor=report("results/figures/column_analysis/combined_variants_corr.tsv",
                             category="Column Analysis")
    shell:
        """python3 {input.script} \
        -s {input.simulated} \
        -d {input.derived} \
        -o {params.out_folder}"""

# Helper functions for the prepare_data rule:
def get_input_variants(wildcards):
    """
    Return the right file format based on the type of variant being processed
    :param wildcards: namespace(like) at least containing wildcard file and variant
    :return: str, required input file
    """
    if wildcards.variant == "derived":
        return f"results/derived_variants/singletons/{wildcards.file}_full_annotation.tsv"
    elif wildcards.variant == "simulated":
        return f"results/simulated_variants/trimmed_snps/{wildcards.file}_full_annotation.tsv"
    elif wildcards.variant == "validation":
        return f"results/validation_variants/{wildcards.file}_full_annotation.tsv"

def get_y(wildcards):
    """
    Determine y value based on variant type.
    :param wildcards: namespace(like) at least containing wildcard variant and file
    :return: str, argument -y <value> or " " if no y is needed
    """
    if wildcards.variant == "derived":
        return "-y 0.0"
    elif wildcards.variant == "simulated":
        return "-y 1.0"
    elif wildcards.variant == "validation":
        if wildcards.file.endswith("y0"):
            return "-y 0.0"
        return "-y 1.0"
    # Not needed for whole_genome
    return " "

"""
Prepare data takes the fully annotated variants and processes 
it as defined in the processing config file.
Means for imputation are already calculated and taken as an input.
It is saved as a sparse matrix in npz format, since npz does not support
column names they are in a separate file, metadata is also stored separately.
"""
rule prepare_data:
    input:
         data=get_input_variants,
         imputaton="results/dataset/imputation_dict.txt",
         processing=config["annotation_config"]["processing"],
         interactions=config["annotation_config"]["interactions"],
         script=workflow.source_path(ANNOTATE_P + "data_preparation.py"),
    params:
         derived_variants=lambda wildcards: "-d" if wildcards.variant == "derived" else " ",
         y=lambda wildcards: get_y(wildcards)
    output:
         npz="results/dataset/{variant}_snps/{file}.npz",
         meta="results/dataset/{variant}_snps/{file}.npz.meta.csv.gz",
         cols="results/dataset/{variant}_snps/{file}.npz.columns.csv"
    wildcard_constraints:
        variant="(derived|simulated|validation)"
    conda:
         "../envs/mainpython.yml"
    priority: 10
    log:
        "results/logs/data_preparation/{variant}_{file}.log"
    shell:
         "python3 {input.script} -i {input.data} --npz {output.npz} "
         "--processing-config {input.processing} "
         "--interaction-config {input.interactions} "
         "--imputation-dict {input.imputaton} "
         "{params.derived_variants} {params.y} > {log}"

"""
Run whole_genome with lower priority.
Building the model has single-threaded elements so it is better to save these
Final computations for last, when the main tasks are bottleneck or finished.
"""
use rule prepare_data as prepare_whole_genome with:
    input:
         data="results/whole_genome_variants/chr{chr}/chr{chr}_{part}_full_annotation.tsv",
         imputaton="results/dataset/imputation_dict.txt",
         processing=config["annotation_config"]["processing"],
         interactions=config["annotation_config"]["interactions"],
         script=workflow.source_path(ANNOTATE_P + "data_preparation.py"),
    params:
         derived_variants=" ",
         y=" "
    priority: -7
    output:
         npz="results/dataset/whole_genome_snps/chr{chr}_{part}.npz",
         meta="results/dataset/whole_genome_snps/chr{chr}_{part}.npz.meta.csv.gz",
         cols="results/dataset/whole_genome_snps/chr{chr}_{part}.npz.columns.csv"
    log:
        "results/logs/data_preparation/whole_genome/{chr}_{part}.log"