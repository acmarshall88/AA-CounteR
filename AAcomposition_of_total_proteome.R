library(bio3d)
library(tidyverse)

# Read in protein sequence file
fasta <- read.fasta("\\\\uniwa.uwa.edu.au\\userhome\\staff7\\00101127\\My Documents\\Sequence data\\Protein\\human_proteome_1protpergene\\2ndQuartOF_UP000005640_9606.fasta")

# all.aas <- length(fasta[["ali"]]!="-")

# extract aa sequences (as character matrix)
aa.seqs <- as.data.frame(fasta$ali)

all.aas <- unlist(aa.seqs, use.names = FALSE) 
all.aas <- all.aas[all.aas!="-"]

# List 20 standard amino acids (A through Y)
aa.list <- LETTERS[c(-2,-10,-15,-21,-24,-26)]

proteome.composition <- 
  as.data.frame(
    sapply(aa.list, function(x) length(all.aas[all.aas==x]))
  )


# prot2 <- proteome.composition
