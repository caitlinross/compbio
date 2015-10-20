#################################
# author: Caitlin Ross
# Due date: 10/23/15
# CSCI 6971
# Assignment 4
#################################


import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


def Main():
    # get seqs from fasta file
    all_seqs = []
    f = open("test.fa", "rU")
    for record in SeqIO.parse(f, "fasta"):
        all_seqs.append(str(record.seq))
    f.close()

    # initialization
    a = []
    for i in range(len(all_seqs)):
        a.append([])
        for j in range(len(all_seqs[i])):
            a[i].append(random.randint(0, len(all_seqs[i])-w))


if __name__ == '__main__':
    Main()
