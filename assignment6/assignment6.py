#################################
# author: Caitlin Ross
# Due date: 11/20/15
# CSCI 6971
# Assignment 6
#################################


import random
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo

def Main():
    # get trees from each file
    tree1 = Phylo.read("tree1.txt", "newick")
    tree1.rooted = True

    tree2 = Phylo.read("tree2.txt", "newick")
    tree2.rooted = True

    tree3 = Phylo.read("tree3.txt", "newick")
    tree3.rooted = True

if __name__ == '__main__':
    Main()
