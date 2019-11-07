#!/usr/bin/env bash

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=dem          ### Job Name
#SBATCH --output=dem.out         ### File in which to store job output
#SBATCH --error=dem.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task




# Full path to input
# /projects/bgmp/rlancio5/demultiplexing-rylanci/input_files

# load env
conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3

# To run the demultiplexing python script
/usr/bin/time -v ./demultiplexing.py -fq input_files/*fastq* -i input_files/indexes.txt
