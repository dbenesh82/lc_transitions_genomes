library(dplyr)
library(ggplot2)

dat <- read.csv(file = "data/combined_data_files_filledin.csv", header = TRUE, sep = ',')

# THIS LIBRARY DID NOT WORK WELL FOR ME
#
# library(seqinr)
# chooseb
# choosebank("genbank")
# query("my1stquery", "SP=Schistocephalus solidus AND M=dna")
# closebank()

# This uses the ncbi API and seems to work better
library(rentrez)
entrez_dbs() # databases to query

# test queries
search_result <- entrez_search(db = "nuccore", term = "Schistocephalus pungitii")
search_result <- entrez_search(db = "nucleotide",
                              term = "KY552799.1[ACCN]")

search_result <- entrez_search(db = "nucleotide",
                               term = "Schistocephalus[ORGN] AND 18s[WORD]")
search_result$ids

# return dna sequence as character
dna <- entrez_fetch(db="nuccore", id=search_result$ids, rettype="fasta")

# write to fasta file for aligning
write(dna, file="quick_dna.fasta")




# these are queries for each species with a genome
queries <- paste(dat$Parasite.species, "[ORGN] AND 18s[WORD]", sep = "")


for(i in seq_along(queries[1:10]) ){
  print(queries[i])
  search_result <- entrez_search(db = "nuccore", term = queries[i])
  
  if( !exists("ids_dat") ) { # if no ids exist, make them
    ids_dat <- data.frame(species = dat$Parasite.species[i], 
                          ids = search_result$ids)
  } else { # otherwise add to them
    ids_dat <- bind_rows(ids_dat, 
                         data.frame(species = dat$Parasite.species[i], 
                                    ids = search_result$ids)
                         )
  }
  
  
}



seq_length_rough <- numeric()
for(i in seq_along(ids_dat$ids[1:20]) ){
  dna <- entrez_fetch(db = "nuccore", id=ids_dat$ids[i], rettype = "fasta") # too many
  seq_length_rough <- c(seq_length_rough, nchar(dna))
  print( paste(ids_dat$species[i], " ", nchar(dna) ) )
}


qplot(x = seq_length_rough) + scale_x_log10()











library(seqinr)
library(refGenome)
library(ape)

cs <- read.fasta(file = "data/clonorchis_sinensis.PRJDA72781.WBPS9.protein.fa.gz", seqtype = "AA")

cs_att <- attr(cs[[1]], which = "Annot")
cs_att

# read GTF file into ensemblGenome object
ens <- ensemblGenome()
read.gtf(ens, "data/clonorchis_sinensis.PRJDA72781.WBPS9.canonical_geneset.gtf")

ens
my_genes <- getGenePositions(ens)

head(my_genes)
sapply(my_genes, unique)

filter(my_genes, !is.na(gene_name))

