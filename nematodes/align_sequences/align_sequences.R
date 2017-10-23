library(seqinr)
library(dplyr)
library(msa)
source("platyhelminthes/align_sequences/clean_alignment_function.R") # for cleaning alignment
options(stringsAsFactors = FALSE)

# read in ortholog names, make into character vector
orthologs <- read.csv(file = "nematodes/get_orthologs/nem_orthologs.csv")
orthologs <- orthologs$x


# open every alignment fasta and check the number of species
# skip the ones with fewer than 70 species
ortho_to_test <- character()
for(o in orthologs){
  
  # make file name string
  file_name <- paste("nematodes/fasta_files/", o, ".txt", sep = "")
  
  if(!file.exists(file_name)) { # if no file, skip it
    next
  } else {
    
    al <- read.alignment(file_name, format = 'fasta')
    if(al$nb >= 70) { # if sequences for at least 70 worms, skip
      ortho_to_test <- c(ortho_to_test, o)
    }
  }
}





# loop through files, align sequences
for(o in ortho_to_test){
  
  # check if alignment file already exists
  file_name <- paste("nematodes/align_sequences/alignments/", o, ".txt", sep = "")
  if(file.exists(file_name)) {
    next # if alignment exists, skip to next fasta file
    
  } else {
    
    # make fasta file name string to load fasta file
    file_name <- paste("nematodes/fasta_files/", o, ".txt", sep = "")
    # read in fasta in msa alignment format
    fasta <- readAAStringSet(file = file_name)
    # align using clustalW with default parameters
    alignment <- msa(fasta)
    # convert to seqinr alignment format, which is a list
    alignment <- msaConvert(alignment, type = "seqinr::alignment")
    # clean up alignment; remove gaps, non-conserved regions
    alignment <- cleanAlignment(alignment, 30, 30)
    # file name for writing alignment file
    file_name <-
      paste("nematodes/align_sequences/alignments/", o, ".txt", sep = "")
    # write aligned fasta
    write.fasta(
      sequences = as.list(alignment$seq),
      names = alignment$nam,
      file.out = file_name,
      as.string = TRUE
    )
  }
}
