"""
 Module that applies a logistic regression model to the dataset generated.
 Trains a model which is then validated using a hold-out test set
 and databases of known variants.
 Additionally, the model is used to score all variants,
 which are then used to generate the whole-genome CADD scores.

 :Author: Job van Schipstal
 :Date: 23-10-2023

 Scripts are based upon the work of Christian Gross,
 but written for scikit-learn instead of turi create.
"""

MODEL_P = "../scripts/train_test_model/"
wildcard_constraints:
    cols="[^/]+"


def get_folds(excluding = None) -> list:
    """
    Get list of numbers, one for each fold that is to be taken as input.
    :param excluding: optional int(-like), fold to exclude (def None)
    :return: List of numbers, usefull for snakemake.expand
    """
    folds = list(range(config["model"]["n_folds"]))
    if excluding:
        excluding = int(excluding)
        if excluding in folds:
            folds.remove(excluding)
    return folds

"""
Loads the different dataset chunks and merges them.
The dataset is then split into n_folds which are each written to disk
"""
rule fold_data:
    input:
         derived=expand("results/dataset/derived_snps/chr{chr}.npz",
                        chr=CHROMS),
         derived_m=expand("results/dataset/derived_snps/chr{chr}.npz.meta.csv.gz",
                          chr=CHROMS),
         derived_c=expand("results/dataset/derived_snps/chr{chr}.npz.columns.csv",
                          chr=CHROMS),
         simulated=expand("results/dataset/simulated_snps/chr{chr}.npz",
                          chr=CHROMS),
         simulated_m=expand("results/dataset/simulated_snps/chr{chr}.npz.meta.csv.gz",
                            chr=CHROMS),
         simulated_c=expand("results/dataset/simulated_snps/chr{chr}.npz.columns.csv",
                            chr=CHROMS),
         script=workflow.source_path(MODEL_P + "fold_data.py"),
         lib=workflow.source_path(MODEL_P + "data_helper.py")
    conda:
         "../envs/mainpython.yml"
    priority: 20
    threads: 4
    resources:
        mem_mb=int(config["dataset_memory_mb"] * 2)
    output:
          test=expand("results/dataset/fold_{fold}.npz",
                      fold=get_folds()),
          test_m=expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                        fold=get_folds()),
          test_c=expand("results/dataset/fold_{fold}.npz.columns.csv",
                        fold=get_folds())
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         -n {threads} \
         -i {input.derived} {input.simulated} \
         -o {output.test}"""

rule train_model:
    input:
          test="results/dataset/fold_{fold}.npz",
          test_m="results/dataset/fold_{fold}.npz.meta.csv.gz",
          test_c="results/dataset/fold_{fold}.npz.columns.csv",
          train=lambda wildcards: expand("results/dataset/fold_{fold}.npz",
                                         fold=get_folds(wildcards.fold)),
          train_m=lambda wildcards: expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                                           fold=get_folds(wildcards.fold)),
          train_c=lambda wildcards: expand("results/dataset/fold_{fold}.npz.columns.csv",
                                           fold=get_folds(wildcards.fold)),
          sel_cols= lambda wildcards: [] if wildcards.cols == "All" else \
              config["model"]["column_subsets"][wildcards.cols],
          script=workflow.source_path(MODEL_P + "train_model.py"),
          lib=workflow.source_path(MODEL_P + "data_helper.py")
    params:
        c=config["model"]["test_params"]["c"],
        max_iter=config["model"]["test_params"]["max_iter"],
        file_pattern="results/model/{cols}/fold_{fold}_[C]C_[ITER]iter.mod",
        sel_cols= lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
         "../envs/mainpython.yml"
    priority: 20
    resources:
        mem_mb=config["dataset_memory_mb"]
    threads: len(config["model"]["test_params"]["c"]) * \
             len(config["model"]["test_params"]["max_iter"])
    output:
          model=expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pickle",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
          stats=expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.stats.txt",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
          weights=expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.weights.csv",
                      c=config["model"]["test_params"]["c"],
                      iter=config["model"]["test_params"]["max_iter"]),
          probs=expand("results/model/{{cols}}/fold_{{fold}}_{c}C_{iter}iter.mod.pred.csv.gz",
                       c=config["model"]["test_params"]["c"],
                       iter=config["model"]["test_params"]["max_iter"]),
          scaler="results/model/{cols}/fold_{fold}.scaler.pickle"
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         --train {input.train} \
         --test {input.test} \
         --columns {params.sel_cols} \
         -c {params.c} \
         -i {params.max_iter} \
         --file-pattern {params.file_pattern} \
         -n {threads} \
         --save-weights \
         --save-scaler {output.scaler}"""

rule final_model:
    input:
          train=expand("results/dataset/fold_{fold}.npz",
                       fold=get_folds()),
          train_m=expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                         fold=get_folds()),
          train_c=expand("results/dataset/fold_{fold}.npz.columns.csv",
                         fold=get_folds()),
          sel_cols = lambda wildcards: [] if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols],
          script=workflow.source_path(MODEL_P + "train_model.py"),
          lib=workflow.source_path(MODEL_P + "data_helper.py")
    params:
        c=config["model"]["test_params"]["c"],
        max_iter=config["model"]["test_params"]["max_iter"],
        file_pattern="results/model/{cols}/full.mod",
        sel_cols= lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
         "../envs/mainpython.yml"
    priority: 20
    resources:
        mem_mb=config["dataset_memory_mb"]
    output:
          model="results/model/{cols}/full.mod.pickle",
          scaler="results/model/{cols}/full.scaler.pickle",
          weights="results/model/{cols}/full.mod.weights.csv"
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         --train {input.train} \
         --columns {params.sel_cols} \
         -c {params.c} \
         -i {params.max_iter} \
         --file-pattern {params.file_pattern} \
         --save-weights \
         --save-scaler {output.scaler}"""

"""
Scores the predicted probability for all possible variants to be of class 1,
(proxy) deleterious. Saved as an csv with chr, pos, ref, alt and raw score.
"""
rule score_variants:
    input:
        data="results/dataset/whole_genome_snps/{file}.npz",
        data_m="results/dataset/whole_genome_snps/{file}.npz.meta.csv.gz",
        data_c="results/dataset/whole_genome_snps/{file}.npz.columns.csv",
        scaler="results/model/{cols}/full.scaler.pickle",
        model="results/model/{cols}/full.mod.pickle",
        script=workflow.source_path(MODEL_P + "model_predict.py"),
    conda:
         "../envs/mainpython.yml"
    priority: -5
    output:
        temp("results/whole_genome_scores/raw_parts/{cols}/{file}.csv")
    shell:
        """python3 {input.script} \
        -i {input.data} \
        --model {input.model} \
        --scaler {input.scaler} \
        -o {output} \
        --sort \
        --no-header"""

"""
Compress VCF file using bgzip if needed and index using tabix, 
also moves the validation variants into the results folder
"""
rule get_validation_variants:
    input:
         vcf=lambda wildcards: config["validation"][wildcards.name][wildcards.type],
         script=workflow.source_path(MODEL_P + "sort_validation_variants.sh"),
    params:
          is_zipped=lambda wildcards: "True" if
            config["validation"][wildcards.name][wildcards.type].endswith(
              ".gz") else "False"
    conda:
         "../envs/mainpython.yml"
    priority: -1
    output:
          vcf=temp("results/validation_variants/{name}_{type}.vcf.gz"),
          index=temp("results/validation_variants/{name}_{type}.vcf.gz.tbi"),
    wildcard_constraints:
        type="(y0|y1)"  # the two variant types
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.vcf} {output.vcf} {params.is_zipped}"

"""
Validates the pre-trained model against datasets of known variants.
Calculates the ROC-AUC score and draws the ROC curve graph.
"""
rule validate_model:
    input:
         npz=expand("results/dataset/validation_snps/{{set}}_{type}.npz",
                        type=["y0", "y1"]),
         meta=expand("results/dataset/validation_snps/{{set}}_{type}.npz.meta.csv.gz",
                        type=["y0", "y1"]),
         cols=expand("results/dataset/validation_snps/{{set}}_{type}.npz.columns.csv",
                        type=["y0", "y1"]),
         sel_cols=lambda wildcards: [] if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols],
         model="results/model/{cols}/full.mod.pickle",
         scaler="results/model/{cols}/full.scaler.pickle",
         script=workflow.source_path(MODEL_P + "test_model.py"),
         lib=workflow.source_path(MODEL_P + "data_helper.py")
    params:
         sel_cols=lambda wildcards: "All" if wildcards.cols == "All" else \
            config["model"]["column_subsets"][wildcards.cols]
    conda:
         "../envs/mainpython.yml"
    output:
          log=report("results/model/{cols}/{set}.stats.txt",
                     category="Validation sets"),
          roc=report("results/figures/validation/{cols}_{set}_roc.png",
                     category="Validation sets"),
    shell:
         """python3 {input.script} \
         -i {input.npz} \
         --model {input.model} \
         --scaler {input.scaler} \
         --columns {params.sel_cols} \
         --module {input.lib} \
         --roc-plot {output.roc} > {output.log}"""

"""
Reproduces the prediction performance plot of (only) CADD by genomic region
as in the work from Gross et al.
"""
rule average_auc_by_region:
    input:
         npz=expand("results/dataset/fold_{fold}.npz",
                      fold=get_folds()),
         meta=expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                        fold=get_folds()),
         cols=expand("results/dataset/fold_{fold}.npz.columns.csv",
                        fold=get_folds()),
         preds=expand("results/model/{{cols}}/fold_{fold}_{{c}}C_{{iter}}iter.mod.pred.csv.gz",
                      fold=get_folds()),
         script=workflow.source_path(MODEL_P + "average_auc_by_region.py"),
         lib=workflow.source_path(MODEL_P + "data_helper.py")
    conda:
         "../envs/mainpython.yml"
    resources:
        mem_mb=int(config["dataset_memory_mb"] / config["model"]["n_folds"])
    output:
          log="results/model/{cols}/auc_by_region_{c}C_{iter}iter.txt",
          fig=report("results/figures/auc_by_region/{cols}_{c}C_{iter}iter.png",
                     category="ROC-AUC by region", caption="../report/auc_by_region.rst",
                     subcategory=lambda wildcards: wildcards.cols)
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         -i {input.npz} \
         -p {input.preds} \
         -o {output.fig} > {output.log}"""

"""
Reproduces part of the prediction performance plot by genomic region
as in the work from Gross et al.
"""
rule column_auc_by_region:
    input:
         npz=expand("results/dataset/fold_{fold}.npz",
                      fold=get_folds()),
         meta=expand("results/dataset/fold_{fold}.npz.meta.csv.gz",
                        fold=get_folds()),
         cols=expand("results/dataset/fold_{fold}.npz.columns.csv",
                        fold=get_folds()),
         script=workflow.source_path(MODEL_P + "plot_column_auc_by_region.py"),
         lib=workflow.source_path(MODEL_P + "data_helper.py")
    params:
         col = config["generate_roc_auc_for_cols"],
         figures = expand("results/figures/column_auc_by_region/{col}.png",
                          col=config["generate_roc_auc_for_cols"])
    conda:
         "../envs/mainpython.yml"
    resources:
        mem_mb=int(config["dataset_memory_mb"])
    output:
          log=report("results/figures/column_auc_by_region/scores.txt",
                     category="ROC-AUC by region",
                     caption="../report/auc_by_region.rst",
                     subcategory="specific columns"),
          fig=directory(report("results/figures/column_auc_by_region",
                     category="ROC-AUC by region",
                     caption="../report/auc_by_region.rst",
                     subcategory="specific columns",
                     patterns=expand("{col}.png",
                                     col=config["generate_roc_auc_for_cols"])))
    shell:
         """python3 {input.script} \
         -m {input.lib} \
         -i {input.npz} \
         -c {params.col} \
         -o {params.figures} > {output.log}"""
