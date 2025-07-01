#!/bin/bash
#SBATCH -J DICT  
#SBATCH -o DICT.out 
#SBATCH -e DICT.err
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

## SCRIPT 
cd $workDir

# Creates a sequence dictionary for the reference sequence. 
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary -R /scratch/gent/vo/002/gvo00224/TOBI/Resources/Genomes/hg19_uxm/hg19.fa  -O hg19.dict

echo "Script completed"
date  ## print end time
