#################################
# author: Caitlin Ross
# Due date: 11/6/15
# CSCI 6971
# Assignment 5
#################################

import sys
import random
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO

def consensusSeq(alignment):
    con_seq = ""
    for j in range(len(alignment[0])):
        counts={}
        counts["A"] = 0
        counts["G"] = 0
        counts["C"] = 0
        counts["U"] = 0
        counts["-"] = 0
        for i in range(len(alignment)):
            counts[alignment[i][j]] += 1
        con_seq += max(counts, key=lambda k: counts[k])
    return con_seq

def Main():
    # read from clustal file
    alignment = []
    f = open("muscle-align.clwstrict")
    for record in AlignIO.read(f, "clustal"):
        alignment.append(str(record.seq))
    f.close()

    con_seq = consensusSeq(alignment)
    print("consensus sequence: " + con_seq)

if __name__ == '__main__':
    Main()
