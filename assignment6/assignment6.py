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
        return P_Lk, str(root)[-1]
    else:
        Phylo.draw_ascii(root) 
        [leftChild, rightChild] = root.clades
        newPlk1 = felsenstein(n-1, leftChild)
        St1, b1 = jukesCantor(root.distance(leftChild))
        newPlk2 = felsenstein(n-1, rightChild)
        St2, b2 = jukesCantor(root.distance(rightChild))
        
        P_Lk["A"] = St1[alphabet[b1]][alphabet["A"]] * newPlk1[b1]
        
        return P_Lk, ""
     
    
    
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
    
    felsenstein(6, root1)
    

if __name__ == '__main__':
    Main()
