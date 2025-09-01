# # snakemake job submission script that submits 100 jobs at a time
# if this is going really slow, can use '-P sage' to send jobs to sage instead
# snakemake --keep-going -j 350 --cluster "qsub -N {params.job_name} -o {params.error_out_file}.out -e {params.error_out_file}.error -l h_rt={params.run_time} -pe serial {params.cores} -l h_vmem={params.memory}G"

# new version with snakemake 8
snakemake --keep-going -j 350 --executor cluster-generic --cluster-generic-submit-cmd "qsub -N {params.job_name} -o {params.error_out_file}.out -e {params.error_out_file}.error -l h_rt={params.run_time} -pe serial {params.cores} -l h_vmem={params.memory}G"
