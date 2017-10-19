library(ggplot2)
library(dplyr)
library(ape)
library(phangorn)

nem <- read.phyDat(file = "nematodes/align_sequences/nematode_aa_alignment.fasta", 
                     format = 'fasta', type = 'AA')

d <- dist.ml(nem)
nj_tree <- NJ(d)
plot(nj_tree)

up_tree <- upgma(d)
plot(up_tree)


# max likelihood
# mt <- modelTest(nem, model=c("JTT", "LG", "WAG"))
# best_mod = "JTT+G+I"

njfit2 <- pml(nj_tree, nem, model = "JTT", k=4, inv=0.2)

# takes a while
ml_fit <- optim.pml(njfit2, rearrangement = "NNI", # if rearrangement = 'stochastic', very long!
                    optInv=TRUE, optGamma=TRUE)

ml_tree <- ml_fit$tree
ml_tree_root <- root(ml_tree, outgroup = c("romanomermis.culicivorax"))
plot(ml_tree_root)

write.nexus(ml_tree, file = "nematodes/make_tree/nem_tree.nexus")