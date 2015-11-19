#################################
# author: Caitlin Ross
# Due date: 11/20/15
# CSCI 6971
# Assignment 6
#################################

import math
from Bio import Phylo

def felsenstein(n, root):
    k = 2*n-1
    P_Lk = {}
    P_Lk["A"] = 0
    P_Lk["C"] = 0
    P_Lk["G"] = 0
    P_Lk["T"] = 0 
      
    if root.is_terminal():
        P_Lk[str(root)[-1]] = 1
        print(root)
        
    else:
        Phylo.draw_ascii(root) 
        [leftChild, rightChild] = root.clades
        newPlk1 = felsenstein(n-1, leftChild)
        St1 = jukesCantor(root.distance(leftChild))
        newPlk2 = felsenstein(n-1, rightChild)
        St2 = jukesCantor(root.distance(rightChild))
        
        for i in range(4):
            tmp1 = 0
            tmp2 = 0
            for b in range(4):
                tmp1 += St1[b][i] * newPlk1[rev_alphabet[b]] 
            for c in range(4):
                tmp2 += St2[c][i] * newPlk2[rev_alphabet[c]] 
            P_Lk[rev_alphabet[i]] = tmp1 * tmp2
            print(rev_alphabet[i] + ": " + str(P_Lk[rev_alphabet[i]]))
     
    return P_Lk
    
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

def Main():
    global alphabet, rev_alphabet
    alphabet = {"A": 0, "C": 1, "G": 2, "T": 3}
    rev_alphabet = {0: "A", 1: "C", 2: "G", 3: "T"}
    
    # get trees from each file
    print("------------------ Test 1 ------------------")
    test = Phylo.read("test.txt", "newick")
    test.rooted = True
    plk = felsenstein(6, test.clade)
    total = 0
    composition = {'A': 0.2, 'C':0.3, 'G': 0.3, 'T': 0.2}
    for k,v in list(plk.items()):
        total += (v * composition[k])
    print("P(data|tree) = " + str(total))
    

    tree1 = Phylo.read("tree1.txt", "newick")
    tree1.rooted = True
    #Phylo.draw_ascii(tree1)

    tree2 = Phylo.read("tree2.txt", "newick")
    tree2.rooted = True
    #Phylo.draw_ascii(tree2)

    tree3 = Phylo.read("tree3.txt", "newick")
    tree3.rooted = True
    #Phylo.draw_ascii(tree3)
    
    root1 = tree1.clade
    root2 = tree2.clade
    root3 = tree3.clade
    
    print("------------------ Tree 1 ------------------")
    plk = felsenstein(6, root1)
    total = 0
    for k,v in list(plk.items()):
        total += (v * .25)
    print("P(data|tree) = " + str(total))
    
    print("\n\n------------------ Tree 2 ------------------")
    plk = felsenstein(6, root2)
    total = 0
    for k,v in list(plk.items()):
        total += (v * .25)
    print("P(data|tree) = " + str(total))
    
    print("\n\n------------------ Tree 3 ------------------")
    plk = felsenstein(6, root3)
    total = 0
    for k,v in list(plk.items()):
        total += (v * .25)
    print("P(data|tree) = " + str(total))

if __name__ == '__main__':
    Main()
