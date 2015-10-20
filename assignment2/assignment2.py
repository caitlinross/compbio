#################################
# author: Caitlin Ross
# Due date: 9/24/15
# CSCI 6971
# Assignment 2
#################################

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

# this function returns two lists from the CpG sequence file
# one with 200 randomly chosen sequences for the + training set
# the other is the rest of the sequences from the CpG file that makes up the + test set
def pos_samp():
    all_seqs = []
    f = open("CpG", "rU")
    # get seqs from fasta file
    for record in SeqIO.parse(f, "fasta"):
        tmp = str(record.seq)
        all_seqs.append(tmp)
    f.close()

    # get 200 sequence for training set
    train_set = []
    for i in range(200):
        pos = random.randint(0, len(all_seqs) - 1)
        train_set.append(all_seqs[pos])
        del all_seqs[pos]

    # train_set is the set of 200 sequences, all_seqs contains all of the other sequences (but no longer the train_set sequences)
    return train_set, all_seqs

# this function returns a randomly chosen sequence from the chromosome sequence file
# input samp_len is the length the chosen sequence should be
# input seq is the full chromosome sequence
def random_seq(samp_len, seq):
    # randomly choose sequence, ignoring ends of chromosomes
    pos = random.randint(10100, len(seq)-10100)
    # make sure the sequence doesn't contain any 'N's
    while not seq[pos:pos+samp_len].find("N") == -1:
        pos = random.randint(10100, len(seq)-10100)
    return seq[pos:pos+samp_len].upper()

# this function calculates the transition matrix for a set of sequences
# input seq_set is a set of sequences
def calc_trans_matrix(seq_set):
    # trans_matrix holds the counts for each pair of nucleotides
    trans_matrix = {}
    trans_matrix["A"] = {}
    trans_matrix["C"] = {}
    trans_matrix["T"] = {}
    trans_matrix["G"] = {}
    for k, v in list(trans_matrix.items()):
        v["A"] = 0
        v["C"] = 0
        v["T"] = 0
        v["G"] = 0

    # get the counts for each sequence in the set
    for seq in seq_set:
        for i in range(1,len(seq)):
            trans_matrix[seq[i-1]][seq[i]] += 1

    # divide each component by the total for that row
    for k,v in list(trans_matrix.items()):
        total = v["A"] + v["G"] + v["T"] + v["C"]
        v["A"] /= float(total)
        v["C"] /= float(total)
        v["T"] /= float(total)
        v["G"] /= float(total)

    return trans_matrix
        
# this function calculates the log odds ratio for a sequence
# input seq is some sequence of which to calc log odds ratio
# pos_matrix is the positive training set's transition matrix
# neg_matrix is the neg training set's transition matrix
def calc_log_odds_ratio(seq, pos_matrix, neg_matrix):
    ratio = 0.0
    for i in range(1,len(seq)):
        ratio += math.log(pos_matrix[seq[i-1]][seq[i]]/neg_matrix[seq[i-1]][seq[i]], 2)
    # returns the log odds ratio normalized by seq length
    return ratio/len(seq)

def Main():
    # get positive training and test sets from CpG file
    pos_train, pos_test = pos_samp()

    # get chr 12 seq from fasta file
    seq = ""
    f = open("chr12.fa", "rU")
    for record in SeqIO.parse(f, "fasta"):
        seq = record.seq
    f.close()

    # create neg training and test sets
    neg_train = []
    neg_test = []
    for i in range(len(pos_train)):
        length = len(pos_train[i])
        neg_train.append(str(random_seq(length, seq)))

    for i in range(len(pos_test)):
        length = len(pos_test[i])
        neg_test.append(str(random_seq(length, seq)))

    # get pos and neg transition matrices
    pos_matrix = calc_trans_matrix(pos_train)
    neg_matrix = calc_trans_matrix(neg_train)

    # print the transition matrices
    print("transition matrix for CpG regions")
    print("\tA\tC\tG\tT")
    print("A  %.5f  %.5f  %.5f  %.5f" % (pos_matrix["A"]["A"], pos_matrix["A"]["C"], pos_matrix["A"]["G"], pos_matrix["A"]["T"]))
    print("C  %.5f  %.5f  %.5f  %.5f" % (pos_matrix["C"]["A"], pos_matrix["C"]["C"], pos_matrix["C"]["G"], pos_matrix["C"]["T"]))
    print("G  %.5f  %.5f  %.5f  %.5f" % (pos_matrix["G"]["A"], pos_matrix["G"]["C"], pos_matrix["G"]["G"], pos_matrix["G"]["T"]))
    print("T  %.5f  %.5f  %.5f  %.5f" % (pos_matrix["T"]["A"], pos_matrix["T"]["C"], pos_matrix["T"]["G"], pos_matrix["T"]["T"]))
    print("\ntransition matrix for non-CpG regions")
    print("\tA\tC\tG\tT")
    print("A  %.5f  %.5f  %.5f  %.5f" % (neg_matrix["A"]["A"], neg_matrix["A"]["C"], neg_matrix["A"]["G"], neg_matrix["A"]["T"]))
    print("C  %.5f  %.5f  %.5f  %.5f" % (neg_matrix["C"]["A"], neg_matrix["C"]["C"], neg_matrix["C"]["G"], neg_matrix["C"]["T"]))
    print("G  %.5f  %.5f  %.5f  %.5f" % (neg_matrix["G"]["A"], neg_matrix["G"]["C"], neg_matrix["G"]["G"], neg_matrix["G"]["T"]))
    print("T  %.5f  %.5f  %.5f  %.5f" % (neg_matrix["T"]["A"], neg_matrix["T"]["C"], neg_matrix["T"]["G"], neg_matrix["T"]["T"]))
    print("\n")

    # calculate the log odds ratio for each sequence
    pos_log_odds = []
    neg_log_odds = []

    for s in pos_test:
        pos_log_odds.append(calc_log_odds_ratio(s, pos_matrix, neg_matrix))
    for s in neg_test:
        neg_log_odds.append(calc_log_odds_ratio(s, pos_matrix, neg_matrix))

    # plot histogram of the log odds ratios
    plt.hist(pos_log_odds, bins=15,  color="blue", alpha=0.5, label="CpG sequences")
    plt.hist(neg_log_odds, bins=15, color="red", alpha=0.5, label="Random Chr 12 sequences")
    plt.legend(loc="upper right")
    plt.xlabel("Log Odds Ratio")
    plt.ylabel("Count")
    plt.title("Log Odds Ratio normalized by length")
    #plt.show()
    plt.savefig("histogram.pdf", format="pdf")
    
    # determine if each sequence in the test_sequences.fasta file is likely to be CpG or not
    test_seqs = {}
    f = open("test_sequences.fasta", "rU")
    for record in SeqIO.parse(f, "fasta"):
        test_seqs[record.id] = calc_log_odds_ratio(str(record.seq), pos_matrix, neg_matrix)
    
    for k,v in list(test_seqs.items()):
        prstr = k + " is probably "
        if v < 0:
            prstr += "not "
        prstr += "a CpG sequence, with a log-odds ratio of " + str(v)
        print(prstr)

if __name__ == '__main__':
    Main()
