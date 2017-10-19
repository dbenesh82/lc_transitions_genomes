library(biomaRt)
library(dplyr)
library(ggplot2)

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "parasite.wormbase.org")

bm_filters <- listFilters(mart)
bm_attr <- listAttributes(mart)


# get non-exhaustive list of genes for nematodes by querying biomart
# we'll choose 3 distantly related roundworms: T.spiralis, C.elegans, S.ratti

# T.spiralis, C.elegans orthologs
genes1 <- getBM(mart = mart, 
               filters = c("species_id_1010", "only_elegaprjna13758_homologue", "with_paralog"),
               value = list("spiraprjna257433", TRUE, FALSE),
               attributes = c("wbps_gene_id", "elegaprjna13758_gene"))

# T.spiralis, S.ratti orthologs
genes2 <- getBM(mart = mart, 
                filters = c("species_id_1010", "only_rattiprjeb125_homologue", "with_paralog"),
                value = list("spiraprjna257433", TRUE, FALSE),
                attributes = c("wbps_gene_id", "rattiprjeb125_gene"))

# C.elegans, S.ratti orthologs
genes3 <- getBM(mart = mart, 
                filters = c("species_id_1010", "only_rattiprjeb125_homologue", "with_paralog"),
                value = list("elegaprjna13758", TRUE, FALSE),
                attributes = c("wbps_gene_id", "rattiprjeb125_gene"))


# make column names more user-friendly
names(genes1) <- c("Tspiralis", "Celegans")
names(genes2) <- c("Tspiralis", "Sratti")
names(genes3) <- c("Celegans", "Sratti")


# join tables; creates NAs
nem_genes <- full_join(genes1, genes2)
nem_genes <- full_join(nem_genes, genes3)

nem_genes <- na.omit(nem_genes) # remove cases where orthologs in one species pair but not another

# still returns duplicates
Tspiralis_remove <- nem_genes%>%
  group_by(Tspiralis)%>%
  summarize(n = n())%>%
  filter(n != 1)

Celegans_remove <- nem_genes%>%
  group_by(Celegans)%>%
  summarize(n = n())%>%
  filter(n != 1)

Sratti_remove <- nem_genes%>%
  group_by(Sratti)%>%
  summarize(n = n())%>%
  filter(n != 1)

nem_genes <- filter(nem_genes, !(Tspiralis %in% Tspiralis_remove$Tspiralis),
                      !(Celegans %in% Celegans_remove$Celegans),
                      !(Tspiralis %in% Sratti_remove$Sratti))




out <- nem_genes$Tspiralis
write.table(out, file = "nematodes/get_orthologs/nem_orthologs.csv",sep = ",", row.names = F)


