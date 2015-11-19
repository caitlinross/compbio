#################################
# author: Caitlin Ross
# Due date: 11/20/15
# CSCI 6971
# Assignment 6
#################################

import math
from Bio import Phylo

# takes in a tree node and recurses through to find the probability
# for each nucleotide at each node
def felsenstein(root):
    P_Lk = {}
    P_Lk["A"] = 0
    P_Lk["C"] = 0
    P_Lk["G"] = 0
    P_Lk["T"] = 0 
      
    if root.is_terminal():
        P_Lk[str(root)[-1]] = 1
        print("probability at leaf " + str(root)[-1] + " " + str(P_Lk))
    else:
        [leftChild, rightChild] = root.clades
        newPlk1 = felsenstein(leftChild)
        St1 = jukesCantor(leftChild.branch_length)
        newPlk2 = felsenstein(rightChild)
        St2 = jukesCantor(rightChild.branch_length)        

        
        for i in range(4):
            tmp1 = 0
            tmp2 = 0
            for b in range(4):
                tmp1 += St1[b][i] * newPlk1[rev_alphabet[b]] 
            for c in range(4):
                tmp2 += St2[c][i] * newPlk2[rev_alphabet[c]] 
            P_Lk[rev_alphabet[i]] = tmp1 * tmp2
        print("probability at node " + str(P_Lk))
    return P_Lk
    
# finds Jukes-Cantor substitution matrix given the branch length
def jukesCantor(branchLen):
    alpha = math.pow(10,-9)
    t = branchLen * 1000000
    St = [[0.0 for j in range(4)] for i in range(4)]
    for i in range(len(St)):
        for j in range(len(St[i])):
            if i == j:
                St[i][j] = 1/4*(1+3*math.exp(-4*alpha*t))
            else:
                St[i][j] = 1/4*(1-math.exp(-4*alpha*t))
    return St

# after getting probablities from felsenstein() for each tree, compute
# the termination step using a uniform background composition 
# (i.e. q_a = .25 for all a in A,C,G,T)    
def finalProb(plk):
    total = 0
    for k,v in list(plk.items()):
        total += (v * .25)
    print("P(data|tree) = " + str(total))

def Main():
    global alphabet, rev_alphabet
    alphabet = {"A": 0, "C": 1, "G": 2, "T": 3}
    rev_alphabet = {0: "A", 1: "C", 2: "G", 3: "T"}
    
    # get trees from each file
    tree1 = Phylo.read("tree1.txt", "newick")
    tree1.rooted = True

    tree2 = Phylo.read("tree2.txt", "newick")
    tree2.rooted = True
    
    tree3 = Phylo.read("tree3.txt", "newick")
    tree3.rooted = True
    
    root1 = tree1.clade
    root2 = tree2.clade
    root3 = tree3.clade
    
    
    print("------------------ Tree 1 ------------------")
    plk = felsenstein(root1)
    finalProb(plk)
    
    print("\n\n------------------ Tree 2 ------------------")
    plk = felsenstein(root2)
    finalProb(plk)
    
    print("\n\n------------------ Tree 3 ------------------")
    plk = felsenstein(root3)
    finalProb(plk)

if __name__ == '__main__':
    Main()
