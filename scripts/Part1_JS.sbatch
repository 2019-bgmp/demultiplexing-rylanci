#!/usr/bin/env bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=DemultiP1	        ### Job Name
#SBATCH --output=meanscores.out         ### File in which to store job output
#SBATCH --error=meanscores.err          ### File in which to store job error messages
#SBATCH --time=0-12:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task
#SBATCH --mail-user=rlancio5@uoregon.edu
#SBATCH --mail-type=ALL

conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3

File1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
File2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
File3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
File4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz


# Calling my python script for each file to produce plots

/usr/bin/time -v ./q_score_dist.py -l 101 -f $File1 -n meanscores_f1


/usr/bin/time -v ./q_score_dist.py -l 8 -f $File2 -n meanscores_f2


/usr/bin/time -v ./q_score_dist.py -l 8 -f $File3 -n meanscores_f3


/usr/bin/time -v ./q_score_dist.py -l 101 -f $File4 -n meanscores_f4
