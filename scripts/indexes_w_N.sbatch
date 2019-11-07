#!/usr/bin/env bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=IwN	        ### Job Name
#SBATCH --output=index_w_N.out         ### File in which to store job output
#SBATCH --error=index_w_N.err          ### File in which to store job error messages
#SBATCH --time=0-6:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task
#SBATCH --mail-user=rlancio5@uoregon.edu
#SBATCH --mail-type=ALL

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR == 2 || NR % 2 ==0' | awk ' NR % 2 ==1' | grep -c 'N'

zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR == 2 || NR % 2 ==0' | awk ' NR % 2 ==1' | grep -c 'N'

