
import os
import pandas as pd
from os import path

if not workflow.overwrite_configfile:
    configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)

include: "snakelib/utils.snake"

rule dump_versions:
    output:
        ver = "versions.txt"
    conda: "env.yml"
    shell:"""
    command -v conda > /dev/null && conda list > {output.ver}
    """

rule merge_tsvs:
    params:
        control = config["control_tails"],
        treated = config["treated_tails"]
    output:
        mt = "merged/all_tails.tsv"
    run:
        dfs = []
        for sample, tsv in params.control.items():
            df = pd.read_csv(tsv, sep="\t")
            df['sample'] = sample
            df['group'] = "control"
            dfs.append(df)
        for sample, tsv in params.treated.items():
            df = pd.read_csv(tsv, sep="\t")
            df['sample'] = sample
            df['group'] = "treatment"
            dfs.append(df)
        adf = pd.concat(dfs)
        adf.to_csv(output.mt, sep="\t", index=False)

rule polya_diff:
    input:
        tails = rules.merge_tsvs.output.mt
    params:
        x = config["per_transcript_plots"],
        min_cov = config["min_reads_per_transcript"],
        pdf = "report/polya_diff_report.pdf",
    output:
        pdf = "report/polya_diff_report.pdf",
        global_tsv = "report/polya_diff_global.tsv",
        tr_tsv = "report/polya_diff_per_transcript.tsv",
    conda: "env.yml"
    shell:"""
    X=""
    if [ {params.x} == "true" ];
    then
        X="-x"
    fi
    {SNAKEDIR}/scripts/polya_diff.py -i {input.tails} -g {output.global_tsv} -t {output.tr_tsv} -r {params.pdf} $X -c {params.min_cov}
    """


rule all:
    input:
        ver = rules.dump_versions.output.ver,
        report = rules.polya_diff.output.pdf,
        global_tsv = rules.polya_diff.output.global_tsv,
