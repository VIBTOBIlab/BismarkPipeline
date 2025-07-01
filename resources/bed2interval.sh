#!/bin/bash
#SBATCH -J BED2INTERVAL  
#SBATCH -o BED2INTERVAL.out 
#SBATCH -e BED2INTERVAL.err
#SBATCH --no-requeue
#SBATCH -c 1 #cores or cpu's
#SBATCH -t 00:10:00 ## time in hr:min:sec
#SBATCH --mem 10GB
#SBATCH --mem 27200MB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sofvdvel.vandevelde@ugent.be

## GENERAL
date  ## print start time

## MODULES
ml purge && ml picard/2.25.1-Java-11

## PARAMETERS
workDir='/scratch/gent/vo/002/gvo00224/vsc44629/pipelines/BismarkPipeline/resources'
RRBS_regions='RRBS_regions20-200_hg19_chr'

## SCRIPT 
cd $workDir

# Convert .bed file into .interval_list file 
java -jar $EBROOTPICARD/picard.jar BedToIntervalList I=${RRBS_regions}.bed  SD=hg19.dict O=${RRBS_regions}.interval_list  

echo "Script completed"
date  ## print end time

