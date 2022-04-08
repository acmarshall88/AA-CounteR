library(bio3d)
library(tidyverse)

# Read in protein sequence file
fasta <- read.fasta("sequences.fasta")

# extract aa sequence (as character matrix)
aa.seqs <- as.data.frame(fasta$ali)

# List 20 standard amino acids (A through Y)
aa.list <- LETTERS[c(-2,-10,-15,-21,-24,-26)]

# seq1 <- aa.seqs[1,]

# seq <- aa.seqs[1,][aa.seqs[1,]!="-"]


aa.compositions <-
  lapply(c(1:length(fasta$id)), function(x) {
    no.of.aaIDs <- c(1:length(aa.list))
    aa.seq <- aa.seqs[x,][aa.seqs[x,]!="-"]
    sapply(no.of.aaIDs, function(x) {
      sum(aa.seq==aa.list[x])*100/length(aa.seq)
    }
    )
  }
  )

df <- as.data.frame(aa.compositions)
rownames(df)<-(aa.list)

df_longformat <- df %>%
  # best to avoid rownames, and instead work explicitly with a column with the 
  # same ID
  rownames_to_column() %>% 
  rename(aa = rowname) %>% 
  # this takes the data from wide format to long format
  pivot_longer(
    # cols = c("P"),
    cols = colnames(df),
    names_to = "protein",
    values_to = "proportion"
  )


av.aa.composition <- 
  sapply(aa.list,
         function(x) mean(df_longformat$proportion[df_longformat$aa==x]))
