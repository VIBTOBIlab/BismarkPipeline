#!/bin/bash
#SBATCH -J MKRRGENOME  
#SBATCH -o MKRRGENOME.out 
#SBATCH -e MKRRGENOME.err
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

## PARAMETERS
workDir='.'
RRBS_regions='RRBS_regions20-200_hg19_chr'

## SCRIPT 
cd $workDir

# Scans genome for MspI recognition sites (C′CGG) saving only those fragments that fall in the specified size range of 20-200 bp.
ml purge && ml Seaborn 
./mkrrgenome/meth_progs_dist/src/mkrrgenome -g ./ -M 20,200 -m 20,200 > ${RRBS_regions}.txt

# Convert .txt file into .bed file 
python ./Scripts/mkrrgenome2bed.py -i ${RRBS_regions}.txt -o ${RRBS_regions}.bed

# Convert .bed file into .interval_list file 
ml purge && ml picard/2.25.1-Java-11
java -jar $EBROOTPICARD/picard.jar BedToIntervalList I=${RRBS_regions}.bed  SD=GRCm38.dict O=${RRBS_regions}.interval_list  


echo "Script completed"
date  ## print end time

