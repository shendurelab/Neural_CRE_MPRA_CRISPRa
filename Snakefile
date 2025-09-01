# Parallelizing SCEPTRE CRISPRaQTL tests across guides
# @author Troy, Diego

import os, sys, pickle, re
os.system('mkdir -p slurm_files ; ')

N_partitions = 100 # set this for the number of partitions

BASE_PATH = '/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/SCEPTRE_PT/'

# this asks snakemake to generate the output file, then snakemake
# looks for the rule that does this task.
rule all:
    input:
        expand(BASE_PATH+"outputs/sceptre_parallel_part{partition}/results_run_discovery_analysis.rds",
            partition=range(1, N_partitions+1))

# Snakemake realizes this is the rule to generate the output so it 
# runs a job with this rule.
rule run_gRNA_tests:
    output: BASE_PATH+"outputs/sceptre_parallel_part{partition}/results_run_discovery_analysis.rds"
    params:
        error_out_file=BASE_PATH + "/slurm_files/tests_{partition}",
        # can set run parameters here
        # lower mem footprint means more jobs to run
        run_time="8:00:00", cores="1", memory="40", job_name="tests{partition}"
    shell:
        "module unload R ; "
        "module load R/4.3.2 ; "
        "Rscript {BASE_PATH}SCEPTRE_snake_V2.r {wildcards.partition} ; "




