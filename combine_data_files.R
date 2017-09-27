# This script combines the list of existing helminth genomes with life cycle length data

library(dplyr)
library(ggplot2)
library(tidyr)
options(stringsAsFactors = FALSE)

# import lc database data
dataH <- read.csv(file="data/CLC_database_hosts.csv", header = TRUE, sep=",")

# reduce to just list of parasite species and their life cycle length.
lcl <- group_by(dataH, Parasite.species)%>%
  summarize(lcl = max(Host.no, na.rm =T))

# import list of worms with sequenced genomes
nem <- read.csv(file="data/species_Nematoda--___.csv", header = TRUE, sep=",")
platy <- read.csv(file="data/species_Platyhelminthes--___.csv", header = TRUE, sep=",")

# combine nematode and platyhelminth lists
genome_worms <- bind_rows( select(nem, Parasite.species = Species.Name, Clade), 
                           select(platy, Parasite.species = Species.Name, Clade) )%>%distinct()

# add life cycle length info
genome_worms <- left_join(genome_worms, lcl)

# most worms lack life cycle length info
table(!is.na(genome_worms$lcl))

# export table and fill it in
write.table(genome_worms, file = "data/combined_data_files.csv", sep = ',', row.names = F)
