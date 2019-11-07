#!/usr/bin/env python3

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=dem          ### Job Name
#SBATCH --output=dem.out         ### File in which to store job output
#SBATCH --error=dem.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task

### Packages
import itertools
import statistics
import argparse
import gzip



def get_args():
	parser = argparse.ArgumentParser(description="A program to demultiplex Illumina sequencing data")
	parser.add_argument("-fq", help="Program requires 4 gzipeed fastq files",
     required=True, type=str, nargs='+')
	parser.add_argument("-i", help="A file with list of barcodes", required=True)
	return parser.parse_args()

args = get_args()


# Barcodes
indi_file = open(args.i,"r")
indi = []
for i, line in enumerate(indi_file):
    if i > 0:
        indi.append(line.strip().split()[4])




### Functions
def is_rev_comp(Rec1,Rec2):
    '''
    Takes two recrds as input. Selects for the sequence line in the records.
    Determines seq1's reverse complement. Then asks if the seq in rec2 is
    the reverse complement of the calulated sequence. Returns boolian.
    '''
    nuc_dict={"A":"T", "C":"G","T":"A","G":"C", "N":"N"}
    LC=0
    def listToString(s):
        str1 = ""
        for ele in s:
            str1 += ele
        return str1

    for l, k  in zip(Rec1, Rec2):
        LC+=1
        if LC%4 == 2:
            temp_list=list()
            for char in l:
                temp_list.append(nuc_dict[char])
            rev_comp=listToString(temp_list)
            if rev_comp[::-1] == k:
                return True
            else:
                return False

def rev_comp_list_maker(list):
    '''
    Takes a list of strings as input. Converts the reverse complement for
    each string and returns them as a new list.
    '''
    nuc_dict={"A":"T", "C":"G","T":"A","G":"C", "N":"N"}
    LC=0
    rev_indi = []

    def listToString(s):
        str1 = ""
        for ele in s:
            str1 += ele
        return str1

    for l in list:
        LC+=1
        temp_list=['']
        for char in l:
            temp_list.append(nuc_dict[char])
            rev_comp=listToString(temp_list)
        rev_indi.append(rev_comp[::-1])
    return rev_indi




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


def convert_phred(letter):
    """Character as input. Converts a single character into a phred score"""
    return ord(letter) - 33


def mean_score(record):
    '''
    File or record as input. Calulates the mean quality score of of a
    sequence using the qs line of the fastq file. Returns mean.
    Dependent on funtion convert_phred().
    '''
    scores = []
    NLN2 = 0
    counter = 0
    for line in record:
        NLN2 +=1
        if NLN2%4 == 0:
            line = line.strip()
            scores= [0] * len(line)
            for i,char in enumerate(line):
                scores[i] += convert_phred(char)

    return statistics.mean(scores)

'''
def if_N_in_line(file):
    Interates through seq and checks for presence of char "N". Returns boolian.
    NLN=0
    for line in file:
        NLN +=1
        if NLN%4 == 2:
            if "N" in line:
                return True
            else:
                return False
'''

def barc_score(record):
    '''
    Record or file as input. Iterates through qs line of fastq
    and if a characters qs is lower than a set parameter it returns true
    Dependent on funtion convert_phred().
    '''
    scores = []
    LN = 0
    counter = 0
    for line in record:
        LN +=1
        if LN%4 == 0:
            line = line.strip()
            for char in line:
                if convert_phred(char) < 15:
                    return True
                else:
                    return False


def known_barcode(record):
    '''
    Record or file input. Checks the seq line of fastq and asks if seq in
    the list of known barcodes "full_bars" which is created from the list indicis
    and their reverse complements E.G. output of rev_comp_list_maker(). Returns boolian.
    '''
    LN=0
    for line in record:
        LN +=1
        if LN%4 == 2:
            line = line.strip()
            if line in full_bars:
                return True
            else:
                return False



# Entering my lazy file iterators in a list with function from helpers.py
frgs = [fastq_records_iterator(gzip.open(args.fq[n],"rt")) for n in range(4)]

# Create a list of our barcodes reverse complements
rev_indi= rev_comp_list_maker(indi)
full_bars = indi + rev_indi


# Output file creation. 48 unqique index files, 2 index hopped files
# and 2 undetermined files for lowquality reads/indicis and unrecognized indicis
file_dict = {}
for line in full_bars:
    index = line.strip()
    file_pointer = open(f"/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/{index}.fastq","w")
    file_dict[f"{index}"] = file_pointer

file_dict["R1_index_hopped"] = open(f"/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/R1_hopped.fastq","w")
file_dict["R2_index_hopped"] = open(f"/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/R2_hopped.fastq","w")
file_dict["R1_undetermined"] = open(f"/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/R1_undefined.fastq","w")
file_dict["R2_undetermined"] = open(f"/projects/bgmp/rlancio5/demultiplexing-rylanci/demulti_output/R2_undefined.fastq","w")

# Some bins for a summary
Index_Hopped = 0
Undetermined = 0
Gd_seqs = 0

# Main Loop
# Lets iterate through fragments and run conditions for binning
for r1_read, r2_id, r3_id, r4_read in zip(frgs[0],frgs[1],frgs[2],frgs[3]):

    #Modify read file headers to include Barcodes for output
    r1_read[0] = f"{r1_read[0]}_{r2_id[1]}:{r3_id[1]}"
    r4_read[0] = f"{r4_read[0]}_{r2_id[1]}:{r3_id[1]}"

    # To bin records that do not pass the barcode or sequence quality filters
    # Bin "undetermined" will also contain records with unidentified barcodes
    if barc_score(r2_id) or barc_score(r3_id):# or mean_score(r1_read) <= 28 or mean_score(r4_read) <= 28:
        Undetermined +=1
        file_dict["R1_undetermined"].write("\n".join(r1_read)+"\n")
        file_dict["R2_undetermined"].write("\n".join(r1_read)+"\n")

    # Records with known barcodes and reverse complements get stored in their unique index file
    elif known_barcode(r2_id) and known_barcode(r3_id) and is_rev_comp(r2_id,r3_id):
        Gd_seqs += 1
        file_dict[f"{r2_id[1]}"].write("\n".join(r1_read)+"\n")
        file_dict[f"{r3_id[1]}"].write("\n".join(r1_read)+"\n")

    # Records with known barcodes that are NOT reverse complements are considered index hopped
    elif known_barcode(r2_id) and known_barcode(r3_id) and is_rev_comp(r2_id,r3_id) == False:
        Index_Hopped +=1
        file_dict["R1_index_hopped"].write("\n".join(r1_read)+"\n")
        file_dict["R2_index_hopped"].write("\n".join(r1_read)+"\n")

    # Bins the unrecognized barcode records in "undetermined"
    else:
        Undetermined +=1
        file_dict["R1_undetermined"].write("\n".join(r1_read)+"\n")
        file_dict["R2_undetermined"].write("\n".join(r1_read)+"\n")


print("Correct Pairs")
print(Gd_seqs)

print("Undetermined/LowQuality")
print(Undetermined)

print("Index Hopped")
print(Index_Hopped)
