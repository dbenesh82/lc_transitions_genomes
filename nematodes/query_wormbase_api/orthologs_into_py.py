# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 15:37:53 2017

@author: phosp
"""
import csv

orthologs = []

with open("get_orthologs\platy_orthologs.csv", newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        orthologs.append(row[0])

orthologs = orthologs[1:] # remove first element (header)


ortho_test = orthologs[0:10]


# need to import these functions

for o in range(0, len(ortho_test)):
    gene = ortho_test[o]
    out = call_wormbase(gene)
    fasta = extract_returned_aa_seq(out)
    
    file_name = "my_fasta" + str(o) + ".txt"
    ofile = open(file_name, "w")
    ofile.write(fasta)
    ofile.close()

