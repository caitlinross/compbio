#################################
# author: Caitlin Ross
# Due date: 10/25/15
# CSCI 6971
# Assignment 4
#################################


import random
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

def sampler(numIterations, isBurnin):
    theta_ij = {}
    theta_ij["A"] = []
    theta_ij["C"] = []
    theta_ij["G"] = []
    theta_ij["T"] = []

    theta_Bj = {}
    theta_Bj["A"] = 0
    theta_Bj["C"] = 0
    theta_Bj["G"] = 0
    theta_Bj["T"] = 0

    counts = {}
    counts["A"] = 0
    counts["G"] = 0
    counts["C"] = 0
    counts["T"] = 0
    for itr in range(numIterations):
        for i in range(len(all_seqs)):
            a[i] = 0
            

            # calculate theta of the motif by counting nucleotides in each position
            for pos in range(w):
                for j in range(len(all_seqs)):
                    if not i == j:
                        counts[all_seqs[j][a[j]+pos]] += 1
                theta_ij["A"].append((float(counts["A"])+ 1) / (sum(counts.values()) + 4))
                theta_ij["G"].append((float(counts["G"])+ 1) / (sum(counts.values()) + 4))
                theta_ij["C"].append((float(counts["C"])+ 1) / (sum(counts.values()) + 4))
                theta_ij["T"].append((float(counts["T"])+ 1) / (sum(counts.values()) + 4))
                counts["A"] = 0
                counts["G"] = 0
                counts["C"] = 0
                counts["T"] = 0
                
            # calculate theta of background by counting nucleotides
            for pos in range(len(all_seqs[0])):
                for j in range(len(all_seqs)):
                    if not i == j:
                        counts[all_seqs[j][pos]] += 1
            theta_Bj["A"] = (float(counts["A"])+ 1) / (sum(counts.values()) + 4)
            theta_Bj["G"] = (float(counts["G"])+ 1) / (sum(counts.values()) + 4)
            theta_Bj["C"] = (float(counts["C"])+ 1) / (sum(counts.values()) + 4)
            theta_Bj["T"] = (float(counts["T"])+ 1) / (sum(counts.values()) + 4)
            counts["A"] = 0
            counts["G"] = 0
            counts["C"] = 0
            counts["T"] = 0

            r = []
            # calculate P(x) for each theta
            for pos in range(len(all_seqs[0])-w):
                pMotif = 1
                pBack = 1
                for m_pos in range(w):
                    pMotif *= theta_ij[all_seqs[i][pos+m_pos]][m_pos]
                    pBack *= theta_Bj[all_seqs[i][pos+m_pos]]
                r.append(pMotif/pBack)
            
            # normalize
            sumR = sum(r)
            for j in range(len(r)):
                r[j] = r[j]/sumR

            tmpA = np.random.choice(len(all_seqs[0])-w, 1, p=r)
            a[i] = tmpA[0]

            if not isBurnin:
                c[i][a[i]] += 1

                if i == numIterations - 1:
                    for j in range(w):
                        allTMotif["A"][j] += theta_ij["A"][j]
                        allTMotif["C"][j] += theta_ij["C"][j]
                        allTMotif["G"][j] += theta_ij["G"][j]
                        allTMotif["T"][j] += theta_ij["T"][j]
                    allTBack["A"] += theta_Bj["A"]
                    allTBack["C"] += theta_Bj["C"]
                    allTBack["G"] += theta_Bj["G"]
                    allTBack["T"] += theta_Bj["T"]


def Main():
    global a
    global w
    global c
    global allTMotif
    global allTBack
    global all_seqs
    
    # get seqs from fasta file
    all_seqs = []
    f = open("test.fa", "rU")
    for record in SeqIO.parse(f, "fasta"):
        all_seqs.append(str(record.seq))
    f.close()

    # initialization
    a = []
    w = 8
    c = []
    for i in range(10):
        c.append([])
        for j in range(len(all_seqs[0])-w):
            c[i].append(0)
    allTMotif = {}
    allTMotif["A"] = []
    allTMotif["G"] = []
    allTMotif["C"] = []
    allTMotif["T"] = []
    for i in range(w):
        allTMotif["A"].append(0)
        allTMotif["G"].append(0)
        allTMotif["C"].append(0)
        allTMotif["T"].append(0)
    allTBack = {}
    allTBack["A"] = 0
    allTBack["G"] = 0
    allTBack["C"] = 0
    allTBack["T"] = 0
    burnin = 10
    sampling = 20

    # run 5 chains
    for i in range(5):
        for i in range(len(all_seqs)):
            a.append(random.randint(0, len(all_seqs[i])-w))
        sampler(burnin, True)
        sampler(sampling, False)

    for i in range(len(c)):
        for j in range(len(c)):
            c[i][j] /= sampling*(i+1)

    print(allTMotif)
    for i in range(w):
        total = allTMotif["A"][i] + allTMotif["G"][i] + allTMotif["C"][i] + allTMotif["T"][i]
        allTMotif["A"][i] /= total
        allTMotif["C"][i] /= total
        allTMotif["G"][i] /= total
        allTMotif["T"][i] /= total

    total = allTBack["A"] + allTBack["G"] + allTBack["C"] + allTBack["T"]
    allTBack["A"][i] /= total
    allTBack["C"][i] /= total
    allTBack["G"][i] /= total
    allTBack["T"][i] /= total

    # print the results
    print("Motif Model\n")
    print(allTMotif)

    print("Background Model:\n")
    print(allTBack)

    print("\nlocations\n")
    for i in range(len(all_seqs)):
        for j in range(len(all_seqs[i] - w)):
            if (c[i][j] > .5):
                print("Sequence " + str(i+1) + "\n")
                print("pos: " + str(j+1) + "\n")
                print("seq: " + all_seqs[i][j:j+w] + "\n")
                print("samp freq: " + str(c[i][j]) + "\n")
    
    # create graphs
    x = np.arange(1, 93)
    for i in range(len(all_seqs)):
        plt.bar(x, c[i])
        plt.xlabel("Position")
        plt.ylable("Probability")
        plt.title("Probability of each position being sampled for Seq " + str(i+1))
        plt.show()


if __name__ == '__main__':
    Main()
