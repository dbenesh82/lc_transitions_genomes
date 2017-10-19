library(seqinr)
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)



# the alignment for each gene will be combined with this df; it's column is genome ids
out_df <- data.frame(name = c('acanthocheilonema_viteae_prjeb4306',
                              'ancylostoma_caninum_prjna72585',
                              'ancylostoma_ceylanicum_prjna231479',
                              'ancylostoma_ceylanicum_prjna72583',
                              'ancylostoma_duodenale_prjna72581',
                              'angiostrongylus_cantonensis_prjeb493',
                              'angiostrongylus_costaricensis_prjeb494',
                              'anisakis_simplex_prjeb496',
                              'ascaris_lumbricoides_prjeb4950',
                              'ascaris_suum_prjna62057',
                              'ascaris_suum_prjna80881',
                              'brugia_malayi_prjna10729',
                              'brugia_pahangi_prjeb497',
                              'brugia_timori_prjeb4663',
                              'bursaphelenchus_xylophilus_prjea64437',
                              'caenorhabditis_angaria_prjna51225',
                              'caenorhabditis_brenneri_prjna20035',
                              'caenorhabditis_briggsae_prjna10731',
                              'caenorhabditis_elegans_prjna13758',
                              'caenorhabditis_japonica_prjna12591',
                              'caenorhabditis_remanei_prjna53967',
                              'caenorhabditis_sinica_prjna194557',
                              'caenorhabditis_tropicalis_prjna53597',
                              'cylicostephanus_goldi_prjeb498',
                              'dictyocaulus_viviparus_prjeb5116',
                              'dictyocaulus_viviparus_prjna72587',
                              'dirofilaria_immitis_prjeb1797',
                              'ditylenchus_destructor_prjna312427',
                              'dracunculus_medinensis_prjeb500',
                              'elaeophora_elaphi_prjeb502',
                              'enterobius_vermicularis_prjeb503',
                              'globodera_pallida_prjeb123',
                              'globodera_rostochiensis_prjeb13504',
                              'gongylonema_pulchrum_prjeb505',
                              'haemonchus_contortus_prjeb506',
                              'haemonchus_contortus_prjna205202',
                              'haemonchus_placei_prjeb509',
                              'heligmosomoides_polygyrus_prjeb1203',
                              'heligmosomoides_polygyrus_prjeb15396',
                              'heterorhabditis_bacteriophora_prjna13977',
                              'litomosoides_sigmodontis_prjeb3075',
                              'loa_loa_prjna246086',
                              'loa_loa_prjna60051',
                              'meloidogyne_floridensis_prjeb6016',
                              'meloidogyne_hapla_prjna29083',
                              'meloidogyne_incognita_prjea28837',
                              'necator_americanus_prjna72135',
                              'nippostrongylus_brasiliensis_prjeb511',
                              'oesophagostomum_dentatum_prjna72579',
                              'onchocerca_flexuosa_prjeb512',
                              'onchocerca_ochengi_prjeb1204',
                              'onchocerca_ochengi_prjeb1809',
                              'onchocerca_volvulus_prjeb513',
                              'panagrellus_redivivus_prjna186477',
                              'parascaris_equorum_prjeb514',
                              'parastrongyloides_trichosuri_prjeb515',
                              'pristionchus_exspectatus_prjeb6009',
                              'pristionchus_pacificus_prjna12644',
                              'rhabditophanes_sp._kr3021_prjeb1297',
                              'romanomermis_culicivorax_prjeb1358',
                              'soboliphyme_baturini_prjeb516',
                              'steinernema_carpocapsae_prjna202318',
                              'steinernema_feltiae_prjna204661',
                              'steinernema_glaseri_prjna204943',
                              'steinernema_monticolum_prjna205067',
                              'steinernema_scapterisci_prjna204942',
                              'strongyloides_papillosus_prjeb525',
                              'strongyloides_ratti_prjeb125',
                              'strongyloides_stercoralis_prjeb528',
                              'strongyloides_venezuelensis_prjeb530',
                              'strongylus_vulgaris_prjeb531',
                              'syphacia_muris_prjeb524',
                              'teladorsagia_circumcincta_prjna72569',
                              'thelazia_callipaeda_prjeb1205',
                              'toxocara_canis_prjeb533',
                              'toxocara_canis_prjna248777',
                              'trichinella_britovi_prjna257433',
                              'trichinella_murrelli_prjna257433',
                              'trichinella_nativa_prjna179527',
                              'trichinella_nativa_prjna257433',
                              'trichinella_nelsoni_prjna257433',
                              'trichinella_papuae_prjna257433',
                              'trichinella_patagoniensis_prjna257433',
                              'trichinella_pseudospiralis_prjna257433',
                              'trichinella_pseudospiralis_prjna257433',
                              'trichinella_pseudospiralis_prjna257433',
                              'trichinella_pseudospiralis_prjna257433',
                              'trichinella_pseudospiralis_prjna257433',
                              'trichinella_spiralis_prjna12603',
                              'trichinella_spiralis_prjna257433',
                              'trichinella_sp._t6_prjna257433',
                              'trichinella_sp._t8_prjna257433',
                              'trichinella_sp._t9_prjna257433',
                              'trichinella_zimbabwensis_prjna257433',
                              'trichuris_muris_prjeb126',
                              'trichuris_suis_prjna179528',
                              'trichuris_suis_prjna208415',
                              'trichuris_suis_prjna208416',
                              'trichuris_trichiura_prjeb535',
                              'wuchereria_bancrofti_prjeb536',
                              'wuchereria_bancrofti_prjna275548')
)






align_files <- list.files(path = "nematodes/align_sequences/alignments/.")
ortho_to_test <- sub(align_files, pattern = ".txt", replacement = "")


# combine alignments into single data frame
for(i in seq_along(ortho_to_test)){
  
  # open alignment
  file_name <- paste("nematodes/align_sequences/alignments/", align_files[i], sep="")
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



sapply(out_df, function(x)sum(is.na(x))) # check if any missing values in combined alignment


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
            file.out = "nematodes/align_sequences/nematode_aa_alignment.fasta",
            as.string=TRUE)


