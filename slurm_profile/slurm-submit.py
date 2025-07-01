#!/usr/bin/env python3
"""
Snakemake SLURM submit script.
"""
from snakemake.utils import read_job_properties
import slurm_utils
from CookieCutter import CookieCutter

##### Get SnakeMake basedir #####
## cluster_config is specific to each workflow. Hence, get config path
## relative to project root or snakemake workflow directory.
## SMKDIR variable is required BEFORE executing snakemake bash command,
## and should point to base path to workflow dir.
import os.path

## Workdir - could be different than workflow dir
workingdir = os.getcwd()

## Do not emit ant stdout in this script except the very last print statement.

# cookiecutter arguments
SBATCH_DEFAULTS = CookieCutter.SBATCH_DEFAULTS
CLUSTER = CookieCutter.get_cluster_option()
CLUSTER_CONFIG = CookieCutter.CLUSTER_CONFIG
ADVANCED_ARGUMENT_CONVERSION = CookieCutter.get_advanced_argument_conversion()

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes"),
    "ntasks-per-node": ("ntasks-per-node", "ntasks_per_node", "ntasks")
}

# parse job
jobscript = slurm_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

sbatch_options = {}
cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)

# 1) sbatch default arguments and cluster
sbatch_options.update(slurm_utils.parse_sbatch_defaults(SBATCH_DEFAULTS))
sbatch_options.update(slurm_utils.parse_sbatch_defaults(CLUSTER))

# 2) cluster_config defaults
sbatch_options.update(cluster_config["__default__"])

# 3) Convert resources (no unit conversion!) and threads
sbatch_options.update(
    slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING)
)

# 4) cluster_config for particular rule
sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 5) cluster_config options
sbatch_options.update(job_properties.get("cluster", {}))

# 6) Advanced conversion of parameters
if ADVANCED_ARGUMENT_CONVERSION:
    sbatch_options = slurm_utils.advanced_argument_conversion(sbatch_options)

# 7) Format pattern in snakemake style
#sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)


# create rule-wildcards based stdout and stderr files
# From Ben Parks @bnprks, https://github.com/bnprks/snakemake-slurm-profile/blob/c967347bbebe123af1533272ae06fa88ba8ec02e/slurm-submit.py#L41-L51
if job_properties["type"] == "single":
    sbatch_options["job-name"] = "snake_" + job_properties["rule"]
    if len(job_properties["wildcards"]) > 0:
        sbatch_options["job-name"] += "_" + "_".join([key + "_" + slurm_utils.file_escape(value) for key,value in job_properties["wildcards"].items()])
    sbatch_options["output"] = os.path.join(workingdir, "logs", "slurm", job_properties["rule"], "") + sbatch_options["job-name"] + "_%j.out"
    sbatch_options["error"] = os.path.join(workingdir, "logs", "slurm", job_properties["rule"], "") + sbatch_options["job-name"] + "_%j.err"
elif job_properties["type"] == "group":
    sbatch_options["job-name"] = "snake_" + job_properties["groupid"]
    sbatch_options["output"] = os.path.join(workingdir, "logs", "slurm", "") + job_properties["groupid"] + "_%j.out"
    sbatch_options["error"] = os.path.join(workingdir, "logs", "slurm", "") + job_properties["groupid"] + "_%j.err"
else:
    print("Error: slurm-submit.py doesn't support job type {} yet!".format(job_properties["type"]))
    sys.exit(1)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(slurm_utils.submit_job(jobscript, **sbatch_options))
