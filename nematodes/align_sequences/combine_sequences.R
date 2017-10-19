library(seqinr)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)



# the alignment for each gene will be combined with this df; it's column is genome ids
out_df <- data.frame(name = c('clonorchis_sinensis_prjda72781','diphyllobothrium_latum_prjeb1206',
                              'echinococcus_canadensis_prjeb8992','echinococcus_granulosus_prjeb121',
                              'echinococcus_granulosus_prjna182977','echinococcus_multilocularis_prjeb122',
                              'echinostoma_caproni_prjeb1207','fasciola_hepatica_prjeb6687',
                              'fasciola_hepatica_prjna179522','gyrodactylus_salaris_prjna244375',
                              'hydatigera_taeniaeformis_prjeb534','hymenolepis_diminuta_prjeb507',
                              'hymenolepis_microstoma_prjeb124','hymenolepis_nana_prjeb508',
                              'macrostomum_lignano_prjna284736','mesocestoides_corti_prjeb510',
                              'opisthorchis_viverrini_prjna222628','protopolystoma_xenopodis_prjeb1201',
                              'schistocephalus_solidus_prjeb527','schistosoma_curassoni_prjeb519',
                              'schistosoma_haematobium_prjna78265','schistosoma_japonicum_prjea34885',
                              'schistosoma_mansoni_prjea36577','schistosoma_margrebowiei_prjeb522',
                              'schistosoma_mattheei_prjeb523','schistosoma_rodhaini_prjeb526',
                              'schmidtea_mediterranea_prjna12585','spirometra_erinaceieuropaei_prjeb1202',
                              'taenia_asiatica_prjeb532','taenia_asiatica_prjna299871',
                              'taenia_saginata_prjna71493','taenia_solium_prjna170813',
                              'trichobilharzia_regenti_prjeb4662')
)





align_files <- list.files(path = "platyhelminthes/align_sequences/alignments/.")
ortho_to_test <- sub(align_files, pattern = ".txt", replacement = "")


# combine alignments into single data frame
for(i in seq_along(ortho_to_test)){
  
  # open alignment
  file_name <- paste("platyhelminthes/align_sequences/alignments/", align_files[i], sep="")
  alignment <- read.alignment(file_name, format = 'fasta')
  
  # make alignment into data frame
  align_df <- data.frame(name = alignment$nam, seq = unlist(alignment$seq))
  gene <- ortho_to_test[i]
  names(align_df) <- c("name", gene)
  
  # join alignment with out df
  out_df <- full_join(out_df, align_df)
  
  # after join, some species lack orthologues for given gene
  # replace those missing cells with "--"
  align_length <- nchar(alignment$seq[[1]])
  missing_str <- paste(rep("-", times = align_length), collapse = "")
  out_df[,gene] <- if_else( is.na(out_df[,gene]), missing_str, out_df[,gene])
}



# sapply(out_df, function(x)sum(is.na(x))) # check if any missing values in combined alignment


# combine alignments into single character vector
cols <- names(out_df)[-1]
out_df2 <- data.frame( seq = do.call(paste, c(out_df[cols], sep="")))
for (co in cols) out_df[co] <- NULL
out_df <- bind_cols(out_df, out_df2)
rm(out_df2)


# simpler species names for tree
sp_names <- substr(out_df$name, start = 1, stop = regexpr(pattern = "_prj", text = out_df$name) - 1)
sp_names <- sub(sp_names, pattern = "_", replacement = " ")


# return alignment to 'alignment' class 
alignx <- list(nb = length(out_df$name),
               nam = sp_names,
               seq = out_df$seq,
               com = NA)
attr(alignx, "class") <- "alignment"

# write alignment fasta
write.fasta(sequences = as.list(alignx$seq),
            names = alignx$nam,
            file.out = "platyhelminthes/align_sequences/platy_aa_alignment.fasta",
            as.string=TRUE)

