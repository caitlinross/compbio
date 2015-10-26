#################################
# author: Caitlin Ross
# Due date: 11/6/15
# CSCI 6971
# Assignment 5
#################################


import random
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



if __name__ == '__main__':
    Main()
