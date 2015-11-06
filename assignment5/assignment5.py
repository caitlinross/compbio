#################################
# author: Caitlin Ross
# Due date: 11/6/15
# CSCI 6971
# Assignment 5
#################################

import math
import numpy as np
from Bio import AlignIO

# returns the consensus sequence from all of the alignment sequences
# alignment is the list of the alignments
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
    
# calculates the mutual information matrix based on the sums of frequencies of base pairs
# takes in the list of alignments
def mutualInformation(alignment):
    matrix = [[0 for j in range(len(alignment[0]))] for i in range(len(alignment[0]))]
    counts = [[0 for j in range(len(alignment[0]))] for i in range(5)]
    
    # get the counts of each base in each column
    for i in range(len(alignment)):
        for j in range(len(alignment[i])):
            counts[alphabet[alignment[i][j]]][j] += 1
    
    # calc frequencies of the bases in each column 
    freq = [[0 for j in range(len(alignment[0]))] for i in range(4)]  
    for j in range(len(counts[0])):
        total = 0
        for i in range(len(counts)):
            if i != 4:
                total += counts[i][j]
        for i in range(len(freq)):
            freq[i][j] = counts[i][j] / total
        
    # calc mutual information matrix from frequencies
    for i in range(len(alignment[0])):
        for j in range(len(alignment[0])):
            joint = {}
            joint["A"] = {}
            joint["G"] = {}
            joint["C"] = {}
            joint["U"] = {}
            for k, v in list(joint.items()):
                joint[k]["A"] = 0 
                joint[k]["G"] = 0
                joint[k]["C"] = 0
                joint[k]["U"] = 0
            total = 0
            # get counts for this i, j pair            
            if i != j:            
                for k in range(len(alignment)):
                    dinuc = alignment[k][i] + alignment[k][j]
                    if alignment[k][i] != '-' and alignment[k][j] != '-':
                        joint[alignment[k][i]][alignment[k][j]] += 1
                        total += 1
            # calculate frequences for this i, j pair
            if i != j:
                for k1, v1 in list(joint.items()):
                    for k2, v2 in list(v1.items()):
                        joint[k1][k2] /= len(alignment[0])
                        
            # calculate i,jth value for mutual information 
            tot = 0.0
            #if j >= i:
            for b in range(len(freq)):
                for bp in range(len(freq)):   
                    dinuc = rev_alph[b]+rev_alph[bp]
                    if dinuc in ('AU', 'UA', 'CG', 'GC', 'GU', 'UG'):
                        fib = freq[alphabet[dinuc[0]]][i]
                        fjb = freq[alphabet[dinuc[1]]][j]
                        if fib != 0 and fjb != 0 and joint[dinuc[0]][dinuc[1]] != 0:
                            tot += joint[dinuc[0]][dinuc[1]] * math.log(joint[dinuc[0]][dinuc[1]]/(fib*fjb), 2)
            
            # put the summed value in the appropriate matrix position
            if tot < 0:
                matrix[i][j] = 0
            else:                
                matrix[i][j] = tot
              
    return matrix
    
# takes in the mutual information matrix to calculate the maximum mutual 
# information secondary structure
def maxMI(mI):
    D = [[0 for j in range(len(mI))] for i in range(len(mI))]
    for i in range(len(mI)-1,-1,-1):
        for j in range(i+3, len(mI)):
            #if abs(i-j) > 2:
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
    return D
            
# backtrace taken from nussinov.py from course website 
def BackTrace(g, seq, mI):
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
            elif g[i,j] == g[i+1, j-1] + mI[i][j]:
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

# taken from nussinov.py from course website to calculate
# the structure based on the pairs from the backtrace 
def StructureFromPairs(pairs, L):
    struct = list('.' * L)
    for i in range(0, len(pairs), 2):
        struct[pairs[i]] = '('
        struct[pairs[i+1]] = ')'
        
    return ''.join(struct)

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

    #get consensus sequence
    con_seq = consensusSeq(alignment)
    print("consensus sequence: " + con_seq)
    
    # get mutual information matrix
    miMatrix = mutualInformation(alignment)
  
    # get maximum mutual information secondary structure
    D = maxMI(miMatrix)
    
    # get pairs using backtrace
    pairs, M = BackTrace(np.array(D), con_seq, miMatrix)
    print("pairs: ")
    print(pairs)

    # determine the structure from the pairs
    struct = StructureFromPairs(pairs, len(alignment[0]))
    print("structure: " + struct)

if __name__ == '__main__':
    Main()
