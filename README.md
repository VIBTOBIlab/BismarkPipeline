# BismarkPipeline v1.0.0
> [!WARNING]
**Deprecated:** This pipeline was used in the [Benchmarking] paper but is no longer maintained. We strongly recommend using the actively supported [MethylSeq pipeline](https://github.com/VIBTOBIlab/MethylSeq) instead of this BismarkPipeline.

## Set-up
Install the pipeline with:
```bash
git clone git@github.com:VIBTOBIlab/BismarkPipeline.git
cd BismarkPipeline/
(echo g/^PipelineDir=/s@.*@PipelineDir="$PWD"@; echo wq) | ed runBismarkPipeline_containers.sh 
```

Then, create a conda environment to run snakemake using the following [yaml file](resources/environment.yml):

```bash
# If you are in the pipeline directory
conda env create -f ./environment.yml
conda activate snakemake-6.1.0
```

## Required files
### Reduced representation (RR) reference genome
The reference genome can be converted in silico to an RR genome with [mkrrgenome](https://github.com/peterstockwell/DMAP). It needs to be run per chromosome, so first you have to convert your reference genome to individual chromosomes.

```bash
mkdir -p chromosomes

awk '/^>/{if (out){close(out)}; out="chromosomes/"substr($1,2)".fa"} {print >> out}' hg19.fa
```

Once this is done, you can run the following command, which scans the genome for MspI recognition sites (C′CGG) saving only those fragments that fall in the specified size range of 20-200 bp.

```bash
# If you are in the pipeline directory
./resources/mkrrgenome/meth_progs_dist/src/mkrrgenome -g ./chromosomes/chr -M 20,200 -m 20,200 > RRBS_regions20-200_hg19_chr.txt
```
Note that the chromosomes/ directory was previously generated to store the individual chromosome fasta files, while "chr" is the prefix of the individual chromosome fasta file (e.g. if chr1.fa, chr2.fa etc, adding "chr" will allow the script to scan files with this prefix, otherwise by default it will search for 1.fa,2.fa, etc).

Finally, you can run the following script to convert the region file to a bed format:

```python
python3 ../resources/mkrrgenome2bed.py -i RRBS_regions20-200_hg19_chr.txt -o RRBS_regions20-200_hg19_chr.bed
```

### Target file
The bed file generated above can be converted to an interval list, target file. The target file can be generated using the Picard [BedToIntervalList](https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard) function, given a dictionary file in input for the corresponding genome. The dictionary can be generated using the Picard [CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard) using the fasta genome file (e.g. hg19.fa).

To run Picard HSmetrics, this target file is necessary. In the [resources](resources/) directory you can find 2 interval list files that contain all of the MspI regions between 20 and 200 bp in hg19 and hg38 genomes.


## Initiating the pipeline
To initiate the pipeline, navigate to a folder where you want to perform the analysis. This folder should only contain the `config_bismark.yaml` with the correct paths. Finally you can run [runBismarkPipeline_containers.sh](runBismarkPipeline_containers.sh):

```bash
screen 
# In folder containing config_bismark.yaml
conda activate snakemake-6.1.0
bash <path-to-BismarkPipeline>/runBismarkPipeline_containers.sh
ctr a d
```

Optionally, you can specify a profile that you can use to optimally submit jobs in an HPC environment.

## Output
The output of the pipeline is a html report with the workflow and MultiQC.

## Error
`subprocess.CalledProcessError: Command '<path>/BismarkPipeline/slurm_profile/slurm-status.py 56792793' returned non-zero exit status 126` can be solved by command:
```bash
chmod +x <path-to-BismarkPipeline>/slurm_profile/slurm-submit.py
```

## Contact
This pipeline was built by [Ruben Van Paemel](https://github.com/rmvpaeme) and [Andries De Koker](https://github.com/AndriesDeKoker). Maintainers of the pipeline are [Edoardo Giuili](https://github.com/edogiuili) and [Sofie Van de Velde](https://github.com/sofvdvel).

## Overview
<img src="./report/workflow.png" width="100%" height="100%">
