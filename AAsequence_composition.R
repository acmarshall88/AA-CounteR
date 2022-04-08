library(bio3d)
library(ggplot2)
library(tidyverse)

# Read in protein sequence file
fasta <- read.fasta("hNONO_Q15233.fasta.txt")

# extract aa sequence (as character vector)
aa.seq <- as.vector(fasta$ali)

# List 20 standard amino acids (A through Y)
aa.list <- LETTERS[c(-2,-10,-15,-21,-24,-26)]

# #### AA composition of entire protein ####
# 
# aa.pc <- sapply(
#                 c(1:length(aa.list)),
#                 function(x) sum(aa.seq==aa.list[x])*100/length(aa.seq)
#                 )
# 
# overall_composition <- data.frame(aa.list, aa.pc)
# 
# 
# ggplot(data = overall_composition, mapping = aes(x=aa.list, y=aa.pc)) +
#   geom_bar(stat = "identity")
# 
# 
# 
#### AA composition in sliding window #### 

window.width <- 50

composition <- data.frame(aa.list)

for (i in 1:(length(aa.seq)-(window.width-1))) {
  
  aa.seq.fragment <- aa.seq[seq(i, (i+window.width-1))]
  
  aa.pc <- sapply(
    c(1:length(aa.list)), 
    function(x) sum(aa.seq.fragment==aa.list[x])*100/window.width
  )
  
  # append to dataframe:
  composition <- cbind(composition, aa.pc)
}


# Set column headers as aa number at the centre of the window: 
colnames(composition) <- c("AA", round(seq(window.width/2, length(aa.seq)-(window.width/2))))


#### To plot composition for all 20 AAs on different panels: ####
composition %>% 
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(composition))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>% 
  ggplot(mapping = aes(x = as.numeric(position), y = proportion)) +
  geom_line() +
  xlab("aa position (centre of window)") +
  ylab("frequency (%)") +
  labs(title = fasta$id, caption = paste("window width = ", window.width)) +
  facet_wrap(~AA)


#### To plot composition for all 20 AAs on same graph: ####
composition %>%
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(composition))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>%
  ggplot(aes(as.numeric(position), proportion, color = AA)) + 
  # geom_point() +
  # geom_line(size=2) +
  geom_smooth(size=2, se=FALSE, method="gam", formula = y~poly(x,20)) +
  xlab("aa position (centre of window)") +
  ylab("frequency (%)") +
  labs(title = fasta$id, caption = paste("window width = ", window.width))

###########################################################################
#### Normalise composition to aa composition of entire human proteome: ####

proteome_composition <- read.csv("human_proteome_AAcomposition.csv")

enrichment <- cbind(composition$AA,
                    composition[-1]/proteome_composition$pc)

colnames(enrichment)[1] <- "AA"


#### To plot enrichment for all 20 AAs on different panels: ####
enrichment %>% 
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(enrichment))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>% 
  ggplot(mapping = aes(x = as.numeric(position), y = proportion)) +
  geom_line() +
  xlab("aa position (centre of window)") +
  ylab("fold enrichment") +
  labs(title = fasta$id, caption = paste("window width = ", window.width)) +
  facet_wrap(~AA)


#### To plot enrichment for all 20 AAs on same graph: ####
enrichment %>%
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(enrichment))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>%
  ggplot(aes(as.numeric(position), proportion, color = AA)) + 
  # geom_point() +
  # geom_line(size=2) +
  geom_smooth(size=2, se=FALSE, method="gam", formula = y~poly(x,20)) +
  xlab("aa position (centre of window)") +
  ylab("fold enrichment") +
  labs(title = fasta$id, caption = paste("window width = ", window.width))
