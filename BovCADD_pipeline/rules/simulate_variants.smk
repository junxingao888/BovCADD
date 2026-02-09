SIMULATED_P = "../scripts/simulate_variants/"

rule create_parameters:
    input:
        ancestral=f"results/ancestral_seq/{config['derive_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        reference=config["generate_variants"]["reference_genome_wildcard"],
        script=workflow.source_path(SIMULATED_P + "create_parameters.py")
    conda: "../envs/mainpython.yml"
    output: "results/simulated_variants/parameters/chr{chr}.txt"
    shell:
        "python3 {input.script} -a {input.ancestral} -r {input.reference} -c {wildcards.chr} -o {output}"

rule required_simulations:
    input:
        derived_count="results/derived_variants/singletons/total.count",
        ancestral_coverage=f"results/ancestral_seq/{config['derive_ancestor']['name_ancestor']}/coverage.txt"
    params:
        factor=config["generate_variants"]["simulate"]["overestimation_factor"]
    output: "results/simulated_variants/simulation.count"
    shell:
        "echo '{params.factor}' | paste {input.derived_count} {input.ancestral_coverage} - | awk '{{printf \"%d\", $1/$2*$3}}' > {output}"

rule process_parameters:
    input:
        parameters=expand("results/simulated_variants/parameters/chr{chr}.txt", chr=CHROMS),
        simulations="results/simulated_variants/simulation.count",
        script=workflow.source_path(SIMULATED_P + "process_parameters.py")
    conda: "../envs/mainpython.yml"
    output:
        parameters="results/simulated_variants/params.pckl",
        log=report("results/logs/process_parameters.log", category="Logs")
    shell:
        "python3 {input.script} -n $(cat {input.simulations}) -p {input.parameters} -l {output.log} -o {output.parameters}"

rule gen_simulated_snps:
    input:
        reference=config["generate_variants"]["reference_genome_wildcard"],
        params="results/simulated_variants/params.pckl",
        script=workflow.source_path(SIMULATED_P + "simulate_variants.py")
    conda: "../envs/mainpython.yml"
    output: "results/simulated_variants/raw_snps/chr{chr}.vcf"
    shell:
        "python3 {input.script} -i {input.reference} -c {wildcards.chr} -p {input.params} --snps {output}"

rule filter_for_anc_site:
    input:
        variants="results/simulated_variants/raw_{type}/chr{chr}.vcf",
        ancestral=f"results/ancestral_seq/{config['derive_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        script=workflow.source_path(SIMULATED_P + "filter_ancestor_site.py")
    conda: "../envs/mainpython.yml"
    output: "results/simulated_variants/filtered_{type}/chr{chr}.vcf"
    shell:
        "python3 {input.script} -i {input.variants} -a {input.ancestral} -o {output}"

rule merge_sim_by_chr:
    input: expand("results/simulated_variants/filtered_{{type}}/chr{chr}.vcf", chr=CHROMS)
    output: "results/simulated_variants/filtered_{type}/all_chr.vcf"
    shell:
        """
        echo '##fileformat=VCFv4.1' > {output}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="CpG context">' >> {output}
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> {output}
        grep -vh '^#' {input} >> {output}
        """
#####in middle not last
rule remove_overlaps_simulated_derived:
    input:
        sim="results/simulated_variants/filtered_snps/all_chr.vcf.gz",
        derived="results/derived_variants/singletons/merged/all_chr.vcf.gz"
    output: "results/simulated_variants/filtered_snps/all_chr.no_overlap.vcf"
    shell:
        "bcftools isec -C -w1 {input.sim} {input.derived} -o {output}"

rule trim_vcf:
    input:
        vcf="results/simulated_variants/filtered_snps/all_chr.no_overlap.vcf",
        simulated_count="results/simulated_variants/filtered_snps/all_chr.no_overlap.vcf.count",
        derived_count="results/derived_variants/singletons/total.count",
        script=workflow.source_path(SIMULATED_P + "trim_vcf.py")
    conda: "../envs/mainpython.yml"
    output: "results/simulated_variants/trimmed_snps/all_chr.clean.vcf"
    shell:
        "python3 {input.script} -i {input.vcf} -o {output} -c $(cat {input.simulated_count}) -d $(cat {input.derived_count})"

rule bgzip_index_trimmed:
    input: "results/simulated_variants/trimmed_snps/all_chr.clean.vcf"
    output:
        vcf="results/simulated_variants/trimmed_snps/all_chr.clean.vcf.gz",
        index="results/simulated_variants/trimmed_snps/all_chr.clean.vcf.gz.tbi"
    shell:
        "bgzip -c {input} > {output.vcf} && tabix -p vcf {output.vcf}"

rule split_by_chrom:
    input:
        vcf="results/simulated_variants/trimmed_snps/all_chr.clean.vcf.gz",
        index="results/simulated_variants/trimmed_snps/all_chr.clean.vcf.gz.tbi"
    output: "results/simulated_variants/trimmed_snps/chr{chr}.vcf"
    conda: "../envs/mainpython.yml"
    shell:
        "bcftools view {input.vcf} --regions {wildcards.chr} -o {output} -O v"
