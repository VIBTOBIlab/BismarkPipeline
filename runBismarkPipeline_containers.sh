PipelineDir=/scratch/gent/vo/002/gvo00224/vsc44629/pipelines/BismarkPipeline
snakefile="${PipelineDir}/Bismark_pipeline_PE_containers.snakefile"
profile="${PipelineDir}/slurm_profile"
clusterTime="${PipelineDir}/clusterLong.json"

snakemake --use-singularity -s ${snakefile} --cluster-config ${clusterTime} --profile ${profile} --jobs 200 --rerun-incomplete
