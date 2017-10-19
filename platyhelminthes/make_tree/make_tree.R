library(ggplot2)
library(dplyr)
library(ape)
library(phangorn)

platy <- read.phyDat(file = "platyhelminthes/align_sequences/platy_aa_alignment.fasta", 
                     format = 'fasta', type = 'AA')

d <- dist.ml(platy)
nj_tree <- NJ(d)
plot(nj_tree)

up_tree <- upgma(d)
plot(up_tree)


# max likelihood
njfit <- pml(nj_tree, data=platy)

# mt <- modelTest(platy, model=c("JTT", "LG", "WAG"))
# best_mod = "JTT+G+I"

njfit2 <- pml(nj_tree, platy, model = "JTT", k=4, inv=0.2)

# takes a while
ml_fit <- optim.pml(njfit2, rearrangement = "NNI", # if rearrangement = 'stochastic', very long!
                    optInv=TRUE, optGamma=TRUE)

ml_tree <- ml_fit$tree
ml_tree_root <- root(ml_tree, outgroup = c("schmidtea.mediterranea", "macrostomum.lignano", 
                           "protopolystoma.xenopodis", "gyrodactylus.salaris"))
plot(ml_tree_root)

write.nexus(ml_tree_root, file = "platyhelminthes/make_tree/platy_tree.nexus")
