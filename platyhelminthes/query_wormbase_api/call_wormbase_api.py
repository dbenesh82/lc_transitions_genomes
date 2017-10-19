# -*- coding: utf-8 -*-
"""
This script calls the wormbase api to get ortholog sequence data
"""

import requests, sys
from is_platyhelminth_func import is_platyhelminth

# function to call API; borrowed from wormbase API examples
def call_wormbase(gene_id):
    server = "https://parasite.wormbase.org"
    ext = "/rest-9/homology/id/" + gene_id + "?sequence=protein;type=orthologues"
    # request parameters are for amino acid sequence and just orthologues
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json() # json by api coerced to nested dictionary
    return(decoded)
    

# function to make fasta file from dictionary returned by wormbase
def extract_returned_aa_seq(data):
        
    d2 = data['data'][0]['homologies'] # this goes to genome level
    
    i = 0
    for x in range(0, len(d2)):
        
        species, aa_seq, ortho = (d2[i]['target']['species'], 
                                  d2[i]['target']['align_seq'],
                                  d2[i]['type'])
        
        if i == 0:
            # first time through loop get source taxon used in api query
            src_sp, src_aa_seq = (d2[0]['source']['species'], 
                                  d2[0]['source']['align_seq'])
            # save results to fasta string
            fasta = ">" + src_sp + "\n" + src_aa_seq + "\n"
        
        
        # only save entries if they are flatworms and one to one orthologues
        if ortho != "ortholog_one2one" or not is_platyhelminth(species):
            next
        else:
            fasta = fasta + ">" + species + "\n" + aa_seq + "\n"
        
        i = i + 1    
    
    return(fasta)
