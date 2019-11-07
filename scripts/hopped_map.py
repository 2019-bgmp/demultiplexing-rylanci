#!/usr/bin/env python3

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=map          ### Job Name
#SBATCH --output=map.out         ### File in which to store job output
#SBATCH --error=map.err          ### File in which to store job error messages
#SBATCH --time=0-02:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task

### Packages
import itertools
import statistics
from collections import defaultdict
import matplotlib.pyplot as plt


# Input File
File = ("/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/R1_hopped.fastq")

R1_indexs =['TACGCTAC', 'ATCGATCG', 'CCTTGATC', 'TCGCTGTT',
 'CATGGCTA', 'GATTACCG', 'ATCCAGAG', 'ATCCGGTA', 'TGAGCTAG',
  'GTGAAGTG', 'AGAGTAGC', 'CTGATCGT', 'GTGCCATA', 'ACGGAACA',
   'CTTAGGAC', 'CTTGTCGA', 'GTCGAAGA', 'CGCATGAT', 'ACCACGAT',
    'ACTCTCGA', 'GAATCCGA', 'GCAAGATC', 'TGGACTCT', 'GCTATCCT',]
R1_i_dict = { i : 0 for i in R1_indexs }


R2_indexs = ['GTAGCGTA','CGATCGAT','GATCAAGG', 'AACAGCGA', 'TAGCCATG',
    'CGGTAATC','CTCTGGAT','TACCGGAT', 'CTAGCTCA', 'CACTTCAC',
    'GCTACTCT','ACGATCAG', 'TATGGCAC','TGTTCCGT', 'GTCCTAAG',
    'TCGACAAG', 'TCTTCGAC','ATCATGCG','ATCGTGGT', 'TCGAGAGT',
    'TCGGATTC', 'GATCTTGC', 'AGAGTCCA', 'AGGATAGC']
R2_i_dict = { i : 0 for i in R2_indexs }

#print("hello")
    # Function from helpers.py via Jared
def fastq_records_iterator(file_pointer):
    '''
    This function will take in a file pointer
    and returns an iterator which yeilds four lines of a file
    (stripped of \n) in a list of strings.
    '''
    while True:
        record = []
        for i in range(4):
            record.append(file_pointer.readline().strip())
        if record[0] == '':
            break
        else:
            yield record
    return None

# Entering my lazy file iterators in a list with function from helpers.py
frg = fastq_records_iterator(open(File,"rt"))


for r in frg:
    line = (r[0])
    line = (line.split(":"))
    bcd1 = line[9][2:]
    bcd2 = line[10]

    if bcd1 in R1_i_dict:
#        print("True")
        R1_i_dict[bcd1] +=1
    if bcd1 in R2_i_dict:
        R2_i_dict[bcd1] +=1
    if bcd2 in R1_i_dict:
        R1_i_dict[bcd2] +=1
    if bcd2 in R2_i_dict:
        R2_i_dict[bcd2] +=1





#print(R1_i_dict.values(), R2_i_dict.values())

#R1_out_file = open("R1_IH_output.txt", "w")

for i, j, k, l in zip(R1_i_dict.keys(), R2_i_dict.values(), R2_i_dict.keys(), R2_i_dict.values()):
    print(i," ",j, " ",k, " ",l)
