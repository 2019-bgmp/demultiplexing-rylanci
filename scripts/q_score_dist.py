#!/usr/bin/env python3
import numpy as np
import argparse
import matplotlib.pyplot as plt
import gzip


def get_args():
    parser = argparse.ArgumentParser(description="A program to display mean phred score per bp plots")
    parser.add_argument("-l", "--readlength", help="read length", required=True, type=int)
    parser.add_argument("-f", "--file", help="filename", required=True)
    parser.add_argument("-n", "--namefig", help="figurename", required=True)
    return parser.parse_args()

args = get_args()
FILE=args.file
rlen=args.readlength
FigName=str(args.namefig)

#TFile1= "R1_test.fq"
#TFile2= "R2_test.fq"
#TFile3= "R3_test.fq"
#TFile4= "R4_test.fq"

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33

mean_scores = []
for i in range(0,rlen):
    mean_scores.append(0.0)
#with open(FILE) as nf:
with gzip.open(FILE, "rt") as nf:
    NLN2 = 0
    counter = 0
    for line in nf:
        NLN2 +=1
        if NLN2%4 == 0:
            line = line.strip()
#            print(len(line))
            pos_counter=0
            for i,char in enumerate(line):
#                np.append(all_qscores[new_counter2], (convert_phred(char)))
                mean_scores[pos_counter] += convert_phred(char)
                pos_counter +=1

#print(mean_scores)
for i,v in enumerate(mean_scores):
    mean_scores[i] = mean_scores[i]/(NLN2/4)
#    print(i,v)




plt.bar(range(rlen), mean_scores)
plt.xlabel("Base Pos")
plt.ylabel("Mean Phed Score")
plt.title("Mean Phred Score per Base")#k-{} cc-{}.format(args.ksize,args.cov_cutoff))
plt.savefig(FigName)


#print(all_qscores)
#all_qscores = np.zeros((101,4000000), dtype=int)
#indexes=[GTAGCGTA, CGATCGAT, GATCAAGG, AACAGCGA, TAGCCATG, CGGTAATC, CTCTGGAT, TACCGGAT, CTAGCTCA, CACTTCAC, GCTACTCT, ACGATCAG, TATGGCAC, TGTTCCGT, GTCCTAAG, TCGACAAG, TCTTCGAC, ATCATGCG, ATCGTGGT, TCGAGAGT, TCGGATTC, GATCTTGC, AGAGTCCA, AGGATAGC]
#1452986940 num lines in R1
# Write func to check for flip indexes
