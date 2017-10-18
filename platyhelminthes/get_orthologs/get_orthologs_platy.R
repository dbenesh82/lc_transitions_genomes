library(biomaRt)
library(dplyr)
library(ggplot2)

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "parasite.wormbase.org")

bm_filters <- listFilters(mart)
bm_attr <- listAttributes(mart)


# get non-exhaustive list of genes for platyhelminthes by querying biomart
# we'll choose 3 distantly related flatworms: S. mansoni, E. multilocularis, S.mediterranea

# S.mansoni, E.multilocularis orthologs
genes1 <- getBM(mart = mart, 
               filters = c("species_id_1010", "only_multiprjeb122_homologue", "with_paralog"),
               value = list("mansoprjea36577", TRUE, FALSE),
               attributes = c("wbps_gene_id", "multiprjeb122_gene"))

# S.mansoni, S.mediterranea orthologs
genes2 <- getBM(mart = mart, 
                filters = c("species_id_1010", "only_meditprjna12585_homologue", "with_paralog"),
                value = list("mansoprjea36577", TRUE, FALSE),
                attributes = c("wbps_gene_id", "meditprjna12585_gene"))

# E.multilocularis, S.mediteranea orthologs
genes3 <- getBM(mart = mart, 
                filters = c("species_id_1010", "only_meditprjna12585_homologue", "with_paralog"),
                value = list("multiprjeb122", TRUE, FALSE),
                attributes = c("wbps_gene_id", "meditprjna12585_gene"))


# make column names more user-friendly
names(genes1) <- c("Smansoni", "Emulti")
names(genes2) <- c("Smansoni", "Smedit")
names(genes3) <- c("Emulti", "Smedit")

# join tables; creates NAs
platy_genes <- full_join(genes1, genes2)
platy_genes <- full_join(platy_genes, genes3)

platy_genes <- na.omit(platy_genes) # remove cases where orthologs in one species pair but not another

# still returns duplicates
Smansoni_remove <- platy_genes%>%
  group_by(Smansoni)%>%
  summarize(n = n())%>%
  filter(n != 1)

Emulti_remove <- platy_genes%>%
  group_by(Emulti)%>%
  summarize(n = n())%>%
  filter(n != 1)

Smedit_remove <- platy_genes%>%
  group_by(Smedit)%>%
  summarize(n = n())%>%
  filter(n != 1)

platy_genes <- filter(platy_genes, !(Smansoni %in% Smansoni_remove$Smansoni),
                      !(Emulti %in% Emulti_remove$Emulti),
                      !(Smansoni %in% Smedit_remove$Smedit))




out <- platy_genes$Smansoni
write.table(out, file = "get_orthologs/platy_orthologs.csv",sep = ",", row.names = F)


