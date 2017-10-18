# -*- coding: utf-8 -*-
"""
This script is playing with the wormbase api to get ortholog data
"""
import requests, sys

# create function to call API
def call_wormbase(gene_id):
    server = "https://parasite.wormbase.org"
    ext = "/rest-9/homology/id/" + gene_id + "?sequence=protein;type=orthologues"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json() # json returned as nested dictionary
    return(decoded)
    

# create another function to make fasta from returned data
def extract_returned_aa_seq(data):
    # extract desired data from json returned by wormbase
    
    d2 = data['data'][0]['homologies']
    
    i = 0
    for x in range(0, len(d2)):
        
        species, aa_seq, ortho = (d2[i]['target']['species'], 
                                  d2[i]['target']['align_seq'],
                                  d2[i]['type'])
        
        if i == 0:
            # first time through loop get source taxon used in api query
            src_sp, src_aa_seq = (d2[0]['source']['species'], 
                                  d2[0]['source']['align_seq'])
            # save results to fasta
            fasta = ">" + src_sp + "\n" + src_aa_seq + "\n"
        
        
        if ortho != "ortholog_one2one" or not is_platyhelminth(species):
            next
        else:
            fasta = fasta + ">" + species + "\n" + aa_seq + "\n"
        
        i = i + 1    
    
    return(fasta)
