#################################
# author: Caitlin Ross
# Due date: 11/6/15
# CSCI 6971
# Assignment 5
#################################

import sys
import random
import math
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
    
def mutualInformation(alignment):
    matrix = [[0 for j in range(len(alignment[0]))] for i in range(len(alignment[0]))]
    counts = [[0 for j in range(len(alignment[0]))] for i in range(5)]
    
    # get the counts of each base in each column
    for i in range(len(alignment[0])):
        for j in range(len(alignment[0])):
            counts[alphabet[alignment[i][j]]][j] += 1
    
    # calc frequencies of the bases in each column        
    for i in range(len(counts)):
        for j in range(len(counts[i])):
            counts[i][j] /= len(alignment)
        
    # this is probably wrong
    # calc mutual information matrix from frequencies
    for row in range(len(alignment[0])):
        for col in range(len(alignment[0])):
            total = 0
            for i in range(len(counts)):
                for j in range(len(counts)):
                    dinuc = rev_alph[i]+rev_alph[j]
                    if dinuc in ('AU', 'UA', 'CG', 'GC', 'GU', 'UG'):
                        fib = counts[alphabet[dinuc[0]]][row]
                        fjb = counts[alphabet[dinuc[1]]][col]
                        if fib != 0 and fjb != 0:
                            total += (fib+fjb) * math.log((fib+fjb)/(fib*fjb))
            
            if total < 0:
                matrix[row][col] = 0
            else:                
                matrix[row][col] = total
                             
#    # check that it's symmetric
#    test1 = np.array(matrix)
#    test2 = test1.T
#    check = test1 == test2
#    for i in range(len(check)):
#        for j in range(len(check[i])):
#            if check[i][j] == False:
#                print(str(test1[i][j]) + " " + str(test2[i][j]))
              
    return matrix
    
def secondaryStructure(mI):
    D = [[0 for j in range(len(mI))] for i in range(len(mI))]
    for i in range(len(mI)):
        for j in range(i+2, len(mI)-1):
            tmp = [0, 0, 0, 0]
            tmp[0] = D[i+1][j]
            tmp[1] = D[i][j-1]
            tmp[2] = D[i+1][j-1] + mI[i][j]
            maxtmp = 0
            for k in range(i+1,j):
                t = D[i][k]+D[k+1][j]
                if t > maxtmp:
                    maxtmp = t
            tmp[3] = maxtmp
            D[i][j] = max(tmp)
            
    
def BackTrace(g, seq):
    """
    BackTrace - backtrace throug Nussinov matrix
    g - matrix returned by Nussinov
    seq - sequence string
    return
    pairs - pairs array; i, i+1 base paired
    M - numpy pairs array M[i,j] == 1 if i and j base paired
    """
    
    def traceback(i, j):
        """
        traceback - recursive back trace through g matrix
                    fills pairs matrix defined in outer scope
        """
        nonlocal pairs
        
        if i < j:
            if g[i,j] == g[i+1, j]:
                traceback(i+1, j)
            elif g[i,j] == g[i, j-1]:
                traceback(i, j-1)
            elif g[i,j] == g[i+1, j-1] + delta(seq[i], seq[j]):
                if abs(i-j) > 2:
                    pairs += [i,j]
                traceback(i+1, j-1)
            else:
                for k in range(i+1, j):
                    if g[i,j] == g[i, k] + g[k+1, j]:
                        traceback(i, k)
                        traceback(k+1, j)
                        break
    
    L = g.shape[0]
    pairs = []
    traceback(0, L-1)

    M = np.zeros([L, L], dtype =int)
    for i in range(0, len(pairs), 2):
        M[pairs[i], pairs[i+1]] = 1

    return pairs, M
    
def delta(si, sj):
    """
    delta - score base pairs for Nussinov algorithm
    si, si - nucleotide letters
    return
    d - delta score 1 for complement or G-U wobble pair, otherwise 0
        1 indicates possible base pairng
    """
    
    nt2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}
    nt2comp = {'A': 3, 'C': 2, 'G': 1, 'T': 0, 'U': 0}

    d = 0
    if (si != '-' and sj != '-'):
        if nt2int[si] == nt2comp[sj]:
            d = 1
        elif (si == 'G' and sj == 'U') or (si == 'U' and sj == 'G'):
            d = 1

    return d

def Main():
    global alphabet, rev_alph
    alphabet = {'A': 0, 'C':1, 'G': 2, 'U':3, '-':4}
    rev_alph = {0:'A', 1:'C', 2:'G', 3:'U', 4:'-'}
    # read from clustal file
    alignment = []
    f = open("muscle-align.clwstrict")
    for record in AlignIO.read(f, "clustal"):
        alignment.append(str(record.seq))
    f.close()

    con_seq = consensusSeq(alignment)
    print("consensus sequence: " + con_seq)
    
    miMatrix = mutualInformation(alignment)
    print("Mutual Information matrix")
    print(miMatrix)
    
    secondaryStructure(miMatrix)
    
    pairs, M = BackTrace(np.array(miMatrix), con_seq)
    print(pairs)
    print(M)

if __name__ == '__main__':
    Main()