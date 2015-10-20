###########################
# author: Caitlin Ross
# due date: 9/17/2015
# CSCI 6971
# Assignment 1
###########################

import sys
import numpy as np
import matplotlib.pyplot as plt

# reads the fasta file to get the sequence
# seq contains the DNA sequence
def readfile(filename):
    seq = ""
    f = open(filename, "r")
    for line in f:
        if not line.startswith(">"):
            seq += line.strip("\n")
    return seq

# calculates the frequencies of the window in the sequence
# count is a dictonary that stores the counts of nucleotides
# window is the window size
# start is the starting position of the window in the genome
def calcfreq(seq, window, start):
    count={}
    count["A"] = 0
    count["G"] = 0
    count["T"] = 0
    count["C"] = 0
    for i in range(start, start + window):
        count[seq[i]] += 1

    # combine counts for AT and GC content
    count["AT"] = count["A"] + count["T"]
    count["GC"] = count["G"] + count["C"]
    return count["A"]/float(window), count["G"]/float(window), count["T"]/float(window), count["C"]/float(window), count["AT"]/float(window), count["GC"]/float(window)

filename = sys.argv[1]
seq = readfile(filename)

window = [500, 1000, len(seq)//20] # window sizes
x = {} # x axis for plotting 
freqs ={} # dictionary to store frequencies for each nucleotide for each window
# for each window, calculate frequences
for i in window:
    x[i] = []
    x[i].extend(np.linspace(0, len(seq)-i, len(seq)-i))
    freqs[i] ={}
    freqs[i]["A"] = []
    freqs[i]["G"] = []
    freqs[i]["T"] = []
    freqs[i]["C"] = []
    freqs[i]["AT"] = []
    freqs[i]["GC"] = []

    for j in range(len(seq)-i):
        a,g,t,c,at,gc = calcfreq(seq, i, j)
        freqs[i]["A"].append(a)
        freqs[i]["G"].append(g)
        freqs[i]["T"].append(t)
        freqs[i]["C"].append(c)
        freqs[i]["AT"].append(at)
        freqs[i]["GC"].append(gc)

    # plot nucleotide frequencies
    plt.plot(x[i], freqs[i]["A"], label="A")
    plt.plot(x[i], freqs[i]["G"], label="G")
    plt.plot(x[i], freqs[i]["T"], label="T")
    plt.plot(x[i], freqs[i]["C"], label="C")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.title("Frequencies of nucleotides when window size is " + str(i))
    plt.legend(loc="upper right")
    #plt.show()
    plt.savefig("nuc-window-"+str(i)+".pdf",format="pdf")
    plt.close()

    # plot AT and GC content frequencies
    plt.plot(x[i], freqs[i]["AT"], label="AT")
    plt.plot(x[i], freqs[i]["GC"], label="GC")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.title("Frequencies of AT content and GC content when window size is " + str(i))
    plt.legend(loc="upper right")
    #plt.show()
    plt.savefig("dinuc-window-"+str(i)+".pdf",format="pdf")
    plt.close()
