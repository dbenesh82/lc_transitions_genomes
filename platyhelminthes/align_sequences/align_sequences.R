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
      # make alignment into data frame
      align_df <- data.frame(name = alignment$nam, seq = alignment$seq)
      names(align_df) <- c("name", o)
      
      # output alignments into a df
      if(!exists('out_df')) { # if a combined out df does not exist create it
        out_df <- align_df
      } else {
        out_df <- full_join(out_df, align_df)  # join existing df
      }
      
      # after join, some species lack orthologues for given gene
      # replace those missing cells with "--"
      align_length <- nchar(alignment$seq[1])
      missing_str <- paste(rep("-", times = align_length), collapse = "")
      out_df[,o] <- if_else( is.na(out_df[,o]), missing_str, out_df[,o])
      
      }
    }
}



sapply(out_df, function(x)sum(is.na(x)))


cols <- names(out_df)[-1]
out_df2 <- data.frame( seq = do.call(paste, c(out_df[cols], sep="")))
for (co in cols) out_df[co] <- NULL

out_df2 <- bind_cols(out_df, out_df2)


alignx <- list(nb = length(out_df2$name),
               nam = out_df2$name,
               seq = out_df2$seq,
               com = NA)
attr(alignx, "class") <- "alignment"



d <- dist.alignment(alignx, matrix = "identity")
nj_tree <- nj(d)
plot(nj_tree)









