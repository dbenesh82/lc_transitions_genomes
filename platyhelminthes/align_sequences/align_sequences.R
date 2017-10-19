library(seqinr)
library(ape)
library(ggplot2)
library(dplyr)
library(msa)
source("platyhelminthes/align_sequences/clean_alignment_function.R") # for cleaning alignment
options(stringsAsFactors = FALSE)

# read in ortholog names, make into character vector
orthologs <- read.csv(file = "platyhelminthes/get_orthologs/platy_orthologs.csv")
orthologs <- orthologs$x


# loop through files, align sequences
for(o in orthologs){
  
  # make file name string
  file_name <- paste("platyhelminthes/fasta_files/", o, ".txt", sep = "")
  
  # check if file exists
  if(!file.exists(file_name)) {
    next
  } else {
    # read in the fasta file
    fasta <- read.alignment(file = file_name, format = 'fasta')
    # if sequences for < 26 flatworms, skip
    if(fasta$nb < 26) {
      next
      } else {
      # read in fasta in msa alignment format
      fasta <- readAAStringSet(file = file_name)
      # align using clustalW with default parameters
      alignment <- msa(fasta)
      # convert to seqinr alignment format, which is a list
      alignment <- msaConvert(alignment, type="seqinr::alignment")
      # clean up alignment; remove gaps, non-conserved regions
      alignment <- cleanAlignment(alignment, 30, 30)
      
      file_name <- paste("platyhelminthes/align_sequences/alignments/", o, ".txt", sep = "")
      
      write.fasta(sequences = as.list(alignment$seq),
                  names = alignment$nam,
                  file.out = file_name,
                  as.string=TRUE)
      
      }
    }
}









