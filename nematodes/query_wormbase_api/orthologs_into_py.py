# -*- coding: utf-8 -*-
"""
Script uses list of orthologues to query wormbase 
and saves fasta file for each orthologue queried
"""

import csv
# import functions to call api, extract sequences, return fasta
from call_wormbase_api import call_wormbase, extract_returned_aa_seq


# import list of orthologs acquired from biomart
orthologs = []
with open("..\\get_orthologs\\nem_orthologs.csv", newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        orthologs.append(row[0])
orthologs = orthologs[1:] # remove first element (header)



 # just take a random number to test
ortho_test = orthologs[0:10]



# query orthologs against wormbase API to aquire amino acid sequences
for o in range(0, len(ortho_test)):
    gene = ortho_test[o]
    out = call_wormbase(gene)
    fasta = extract_returned_aa_seq(out)
    
    file_name = "..\\fasta_files\\" + gene + ".txt"
    ofile = open(file_name, "w")
    ofile.write(fasta)
    ofile.close()

