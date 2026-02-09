"""
 Module that generates annotations.
 PhastCons and PhyloP are ran locally,
 gerp scores are downloaded from Ensembl.

 Other annotations are taken as an external input in bedgraph(-like) format.
 These are not handled by this module.

 :Author: Job van Schipstal
 :Date: 23-9-2023

 Based upon the work of Christian Groß.
"""
import sys

PHAST_BIN_F = "../scripts/phast_bin/"
GEN_VAR_F = "../scripts/generate_annotations/"

rule split_maf:
    input:
         maf="results/alignment/sorted/{name}/chr{chr}.maf",
         script=workflow.source_path(GEN_VAR_F + 'split_maf.sh')
    params:
          block_size=config["parallelization"]["phast_blocks_per_split"]
    output:
          folder=directory("results/phast/alignment/{name}/chr{chr}"),  #TODO mark temp
          mock=temp("results/phast/alignment/{name}/chr{chr}/finished_splitting.txt")
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.maf} {params.block_size} {output.folder} && "
         "touch {output.mock}"

checkpoint split_stats:
    input:
        expand("results/phast/alignment/{{name}}/chr{chr}/finished_splitting.txt",
               chr=CHROMS)
    output:
        "results/logs/{name}/split_log.txt"
    shell:
        "ls -alh {input} > {output}"

rule phylo_fit:
    input:
        maf=lambda wildcards: f"results/phast/alignment/{config['phast'][wildcards.name]['alignment']}/chr{wildcards.chr}/{wildcards.part}.maf",
        bin=workflow.source_path(PHAST_BIN_F + 'phyloFit')
    params:
        tree=lambda wildcards: config["phast"][wildcards.name]["tree"],
        precision=lambda wildcards: config["phast"][wildcards.name]["train_precision"],
        out="results/phast/phylo_model/{name}/chr{chr}/{part}"
    log:
        "results/logs/{name}/chr{chr}_{part}_phylo_fit_log.txt"
    output:
        "results/phast/phylo_model/{name}/chr{chr}/{part}.mod"
    shell:
        """
        echo 'Using tree: {params.tree}' >> {log}
        chmod +x {input.bin}
        {input.bin} \
            --tree "{params.tree}" \
            -p {params.precision} \
            --subst-mod REV \
            --out-root {params.out} \
            {input.maf} 2>> {log}
        """

def aggregate_model(name):
    aggregate_input = []
    alignment_name = config["phast"][name]["alignment"]
    chroms = CHROMS
    checkpoints.split_stats.get(name=alignment_name)
    for chrom in chroms:
        parts = glob_wildcards(
            f"results/phast/alignment/{alignment_name}/chr{chrom}/{{part}}.maf").part
        aggregate_input.extend(expand(
            f"results/phast/alignment/{alignment_name}/chr{chrom}/{{part}}.maf",
            part=parts))
        aggregate_input.extend(
            expand(f"results/phast/phylo_model/{name}/chr{chrom}/{{part}}.mod",
                   part=parts))
    return aggregate_input


rule merge_models:
    input:
         infiles=lambda wildcards: aggregate_model(wildcards.name),
         script=workflow.source_path(
             GEN_VAR_F + 'phyloFit_model_weighted_mean.py')
    params:
          model_p="results/phast/phylo_model/{name}/",
          chunk_p=lambda
              wildcards: f"results/phast/alignment/{config['phast'][wildcards.name]['alignment']}/"
    conda:
         "../envs/mainpython.yml"
    output:
          "results/phast/phylo_model_merged/{name}.mod"
    shell:
         "python3 {input.script} -m {params.model_p} -c {params.chunk_p} -o {output}"

rule run_phastCons:
    input:
         maf=lambda wildcards: f"results/phast/alignment/{config['phast'][wildcards.name]['alignment']}/chr{{chr}}/{{part}}.maf",
         mod=lambda wildcards:
         config["phast"][wildcards.name].get("phast_cons_model",
                                             f"results/phast/phylo_model_merged/{wildcards.name}.mod"),
         bin=workflow.source_path(PHAST_BIN_F + 'phastCons')
    params:
          species_interest=lambda wildcards: config["alignments"][
              config["phast"][wildcards.name]["alignment"]
                ]["name_species_interest"],
          phast_params=lambda wildcards: config["phast"][wildcards.name]["phast_cons_params"]
    output:
          "results/phast/phastCons/{name}/chr{chr}/{part}.wig"
    shell:
         "chmod +x {input.bin} "
         "&& {input.bin} --msa-format MAF "
         "--not-informative={params.species_interest} "
         "{params.phast_params} {input.maf} {input.mod} > {output}"
 
rule run_phyloP:
    input:
         maf=lambda wildcards: f"results/phast/alignment/{config['phast'][wildcards.name]['alignment']}/chr{{chr}}/{{part}}.maf",
         mod=lambda wildcards:
         config["phast"][wildcards.name].get("phylo_p_model",
                                             f"results/phast/phylo_model_merged/{wildcards.name}.mod"),
         bin=workflow.source_path(PHAST_BIN_F + 'phyloP')
    params:
          species_interest=lambda wildcards: config["alignments"][
              config["phast"][wildcards.name]["alignment"]
                ]["name_species_interest"],
          phylo_params=lambda wildcards: config["phast"][wildcards.name]["phylo_p_params"]
    benchmark:
        "logs/phast/phyloP/{name}/chr{chr}/{part}.tsv"
    output:
          "results/phast/phyloP/{name}/chr{chr}/{part}.wig"
    shell:
         "cat {input.maf} | sed 's/^a .*/a/' > "
         "tmp.{wildcards.chr}_{wildcards.part}.maf "
         "&& chmod +x {input.bin} "
         "&& {input.bin} --msa-format MAF "
         "--chrom {wildcards.chr} --wig-scores "
         "--not-informative={params.species_interest} "
         "{params.phylo_params} {input.mod}"
         " tmp.{wildcards.chr}_{wildcards.part}.maf > {output} "
         "&& rm tmp.{wildcards.chr}_{wildcards.part}.maf"

def get_parts(wildcards):
    alignment_name = config["phast"][wildcards.name]["alignment"]
    checkpoints.split_stats.get(name=alignment_name)
    parts = glob_wildcards(f"results/phast/alignment/{alignment_name}/chr{wildcards.chr}/{{part}}.maf").part
    # ensure ascending order of parts
    parts = sorted(parts, key=int)
    return expand("results/phast/{{type}}/{{name}}/chr{{chr}}/{part}.wig", part=parts)

"""
Since the phyloP and phastCons were generated for blocks their scores need to be merged into a single file.
The chromosome was not extracted from the data so we add it back in here with sed find/replace.
"""
rule merge_phast_chr:
    input:
        get_parts
    output:
        "results/phast/{type}/{name}/chr{chr}.wig"
    wildcard_constraints:
        type="[^/]+",
        name="[^/]+"
    shell:
        "cat {input} | sed 's/chrom=(null)/chrom={wildcards.chr}/g' > {output}"
"""
Concatenates phyloP/phastCons conservation score files for all chromosomes
"""
rule merge_phast_full:
    input:
         expand("results/phast/{{type}}/{{name}}/chr{chr}.wig", chr=CHROMS)
    output:
          "results/phast/{type}/{name}/all_chr.wig"
    wildcard_constraints:
        type="[^/]+",
        name="[^/]+"
    shell:
         "cat {input} > {output}"


rule install_gerp_api:
    input:
        workflow.source_path(GEN_VAR_F + "install_ensembl_apiclient.sh")
    output:
        directory(config["ensembl-api"]["directory"])
    shell:
         """mkdir {output} && cd {output} && \
         source {input} && chmod -R +x ."""


def get_chrom_length(wildcards):
    fasta_t = config["generate_variants"]["reference_genome_wildcard"]
    index = fasta_t.replace("{chr}",wildcards.chr) + ".fai"
    try:
        with open(index, 'r') as infile:
            first_line = infile.readline().strip().split("\t")
            if not len(first_line) == 5:
                sys.exit(f"Malformed fasta index: {first_line} in {index}")
            return first_line[1]
    except FileNotFoundError:
        sys.stderr.write(f"Missing fasta index file: {index}, "
                         f"it should be generated automatically later.\n")
        return "None" # Should only be returned on dry-runs, unneeded there.

rule get_gerp:
    input:
        install=rules.install_gerp_api.output if config["ensembl-api"]["directory"] else [],
        setup=workflow.source_path(GEN_VAR_F + "setup_api_env.sh"),
        script=workflow.source_path(GEN_VAR_F + "get_gerp_ensembl.pl"),
    params:
         dir=config["ensembl-api"]["directory"],
         is_elem=lambda wildcards: "--elem" if wildcards.type == "elem" else "",
         chrom_length=get_chrom_length,
         tmp_file="{name}.{set}.{chr}.{type}.bed"
    conda:
        "../envs/ensembl.yml"
    wildcard_constraints:
        type="elem|score"
    output:
         "results/gerp/{name}/{set}/chr{chr}_{type}.bed"
    shell:
        """pushd {params.dir} && \
        source {input.setup} . && \
        perl {input.script} \
        --chr {wildcards.chr} \
        --start 0 \
        --end {params.chrom_length} \
        --sp {wildcards.name} \
        --set {wildcards.set} \
        {params.is_elem} > {params.tmp_file} && \
        popd && \
        mv {params.dir}/{params.tmp_file} {output}"""

"""
Concatenates gerp conservation score files for all chromosomes
"""
rule merge_gerp_full:
    input:
         expand("results/gerp/{{name}}/{{set}}/chr{chr}_{{type}}.bed",
             chr=CHROMS)
    conda:
         "../envs/mainpython.yml" # Contains bgzip
    output:
          "results/gerp/{name}/{set}/{type}_all.bed.gz"
    wildcard_constraints:
        type="[^/]+",
        set="[^/]+",
        name="[^/]+"
    shell:
         "cat {input} | bgzip -c > {output}"