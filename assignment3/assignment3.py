#################################
# author: Caitlin Ross
# Due date: 10/9/15
# CSCI 6971
# Assignment 3
#################################

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

# use forward algorithm described in textbook
# e is emission matrix
# a is transition matrix
# f is forward probabilities being calculated
# seq is the sequence
# i is the position in the sequence
# l is the state it's being calculated for
def forward_alg(e, a, f, seq, i, l):
    sumK = 0.0
    if i > 0:
        for k in range(len(f)):
            sumK += f[k][i-1]*a[k][l]
        f[l].append(e[l][seq[i-1]]*sumK)
    return f[l][i]

# use backward algorithm described in textbook
# e is emission matrix
# a is transition matrix
# b is backward probabilities being calculated
# seq is the sequence
# i is the position in the sequence
# k is the state it's being calculated for
def backward_alg(e, a, b, seq, i, k):
    sumL = 0.0
    if i < len(seq):
        for l in range(len(b)):
            if not i==0:
                sumL +=a[k][l]*e[l][seq[i]]*b[l][i+1]
        b[k][i] = sumL
    return b[k][i]  

def Main():
    # get seqs from fasta file
    all_seqs = []
    f = open("assign3.fasta", "rU")
    for record in SeqIO.parse(f, "fasta"):
        all_seqs.append(str(record.seq))
    f.close()

    #initialize emission probabilities
    e = []
    E = []
    tmpE1 = []
    tmpE2 = []
    for i in range(2):
        e.append({})
        e[i]["A"] = random.random()
        e[i]["G"] = random.random()
        e[i]["T"] = random.random()
        e[i]["C"] = random.random()
        total = sum(e[i].values())
        e[i]["A"] /= total
        e[i]["G"] /= total
        e[i]["T"] /= total
        e[i]["C"] /= total
        E.append({})
        E[i]["A"] = 0
        E[i]["G"] = 0
        E[i]["T"] = 0
        E[i]["C"] = 0
        tmpE1.append({})
        tmpE1[i]["A"] = 0
        tmpE1[i]["G"] = 0
        tmpE1[i]["T"] = 0
        tmpE1[i]["C"] = 0      
        tmpE2.append({})
        tmpE2[i]["A"] = 0
        tmpE2[i]["G"] = 0
        tmpE2[i]["T"] = 0
        tmpE2[i]["C"] = 0
        

    # initialize transition probabilities
    a = []
    A = []
    tmpA1 = []
    tmpA2 = []

    for i in range(2):
        a.append([])
        a[i].append(.9)
        a[i].append(.1)
        A.append([])
        A[i].append(0)
        A[i].append(0)
        tmpA1.append([])
        tmpA1[i].append(0)
        tmpA1[i].append(0)
        tmpA2.append([])
        tmpA2[i].append(0)
        tmpA2[i].append(0)
        
    a[1][0] = .1
    a[1][1] = .9
    
    #define start and end transitions
    a0k = .5
    ak0 = 1
    
    # used to determine when to end BW
    old_log_sum = -9999999
    log_sum = 0
    index = 0

    #Baum-Welch
    while abs(old_log_sum - log_sum) > .000001 and index < 1000:
        old_log_sum = log_sum
        log_sum = 0
        index +=1
    
        for seq in all_seqs:
            # initialize f
            f = []
            for i in range(2):
                f.append([])
                f[i].append(0)
            f[0][0] = 1
            
            # initialize b
            b=[]
            for i in range(2):
                b.append([])
                b[i].append(0)
                
            for i in range(1, len(seq)+1):
                f[0][i] = forward_alg(e, a, f, seq, i, 0)
                f[1][i] = forward_alg(e, a, f, seq, i, 1)
                b[0].append(0)
                b[1].append(0)
            
            p_x_for = f[0][-1]*ak0 + f[1][-1]*ak0
            
            for i in range(2):
                b[i][-1] = ak0
                
            for i in range(len(seq)-1, 0, -1):
                b[0][i] = backward_alg(e, a, b, seq, i, 0)
                b[1][i] = backward_alg(e, a, b, seq, i, 1)
                
            p_x_back = a0k*e[0][seq[0]]*b[0][1] + a0k*e[1][seq[0]]*b[1][1]
#        
            log_sum += math.log(p_x_for)
            
            # add contribution of seq j to A using eqn (3.20)
            for k in range(2):
                for l in range(2):
                    tmpA1[k][l] += 1/p_x_for
                    for i in range(len(seq)):
                        tmpA2[k][l] += f[k][i]*a[k][l]*e[l][seq[i]]*b[l][i+1]
                
            # add contribution of seq j to E using 
            for k in range(2):
                for key,val in list(E[k].items()):
                    tmpE1[k][key] += 1/p_x_for
                    for i in range(len(seq)):
                        if key == seq[i]:
                            tmpE2[k][key] += f[k][i+1]*b[k][i+1]
        print(log_sum)                      
        for k in range(2):
            for l in range(2):
                A[k][l] = tmpA1[k][l]*tmpA2[k][l]
            for key, val in list(E[k].items()):
                E[k][key] = tmpE1[k][key]*tmpE2[k][key]
    
        for k in range(2):
            total = sum(A[k])
            for l in range(2):
                a[k][l] = A[k][l]/total
            total = sum(E[k].values())
            for key,val in list(E[k].items()):
                e[k][key] = E[k][key]/total

    # print final information
    print("num iterations = " + str(index))
    print("transition matrix")
    print(a)
    print("emission matrix")
    print(e)
    
if __name__ == '__main__':
    Main()
