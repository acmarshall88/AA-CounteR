library(bio3d)
library(ggplot2)
library(tidyverse)

# Read in human proteome composition table
proteome_composition <- read.csv("human_proteome_AAcomposition.csv")

# Read in protein sequence file
fasta <- read.fasta("hSFPQ_P23246.fasta.txt")

# extract aa sequence (as character vector)
aa.seq <- as.vector(fasta$ali)

# List 20 standard amino acids (A through Y)
aa.list <- LETTERS[c(-2,-10,-15,-21,-24,-26)]

# Create a "named vector" lookup table for different aa values
comp_lookup <- pull(proteome_composition, pc, name = aa)

# helper funciton to classify a AA position into one of the 3 regions
classify_region <- function(x) {
  case_when(
    x %in% 1:275 ~ "nlcr", 
    x %in% 276:598 ~ "core", 
    x %in% 599:707 ~ "clcr"
  ) %>% 
    factor(levels = c("nlcr", "core", "clcr"))
}

# Set up tibble with all AA, their number and their categorised region of the protein
dat <- tibble(
  seq = aa.seq
) %>% 
  mutate(
    resnum = row_number(),
    region = classify_region(resnum)
  )

# Do the region-wise summary for each AA in each region as a relative proportion
# to the overall genomic enrichment
summary_data <- dat  %>% 
  summarise(
    count = n(),
    .by = c(region, seq)
  ) %>% 
  complete(region, seq = aa.list, fill = list(count = 0)) %>% 
  mutate(
    perc = ifelse(count == 0, 0, count / sum(count) * 100 / comp_lookup[seq]), 
    .by = region
  )

# Create a separate dataframe for the grey bars
grey_bar_data <- summary_data %>% 
  filter(perc == 0) %>%
  mutate(ymin = 0, ymax = 1.0)

# Convert seq to factor to ensure correct positioning
summary_data$seq <- factor(summary_data$seq, levels = aa.list)
grey_bar_data$seq <- factor(grey_bar_data$seq, levels = aa.list)

# Plot with grey bars and enrichment bars
plt1 <- ggplot(summary_data, aes(seq, perc, fill = region)) + 
  geom_col(alpha = 0.5) + 
  geom_rect(data = grey_bar_data, 
            aes(xmin = as.numeric(seq) - 0.45, xmax = as.numeric(seq) + 0.45, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2) +
  facet_wrap(~region)  + 
  xlab("amino acid") +
  ylab("enrichment") +
  scale_y_log10(limits = c(0.09, 10), breaks = c(0.1, 0.3, 1.0, 3, 10)) +
  theme_bw()

#### AA composition of entire protein ####
aa.pc <- sapply(
                c(1:length(aa.list)),
                function(x) sum(aa.seq==aa.list[x])*100/length(aa.seq)
                )

# Normalise composition to aa composition of entire human proteome and include
# in final column:
overall_composition <- data.frame(aa.list, aa.pc, aa.pc/proteome_composition$pc)

# Plot of absolute total composition:
ggplot(data = overall_composition, mapping = aes(x=aa.list, y=aa.pc)) +
  geom_bar(stat = "identity") +
  theme_bw()

# Plot of enrichment (or depletion) relative to proteome composition:
ggplot(data = overall_composition, 
       mapping = aes(
         x=aa.list, 
         y=aa.pc.proteome_composition.pc,
         # fill = aa.pc.proteome_composition.pc
         )) +
  geom_bar(stat = "identity") +
  xlab("amino acid") +
  ylab("enrichment") +
  scale_y_log10(limits = c(0.1,10), breaks = c(0.1,0.3,1.0,3,10)) +
  theme_bw()




#### AA composition in sliding window #### 

window.width <- 30

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


# Create parameters for coloring background of plot(s) by region:
tibble(
  num = 1:length(aa.seq)
) %>% 
  mutate(
    region = classify_region(num)
  ) %>% 
  summarise(start = min(num), 
            end = max(num), 
            .by = region) %>% 
  mutate(
    bottom = 0, top = 100
  ) -> box_params


#### Normalise composition to aa composition of entire human proteome: ####

enrichment <- cbind(composition$AA,
                    composition[-1]/proteome_composition$pc)

colnames(enrichment)[1] <- "AA"



# #### To plot composition for all 20 AAs on different panels: ####
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

# ***OPTION TO ONLY INCLUDE AAs OF INTEREST***:

# aas_to_include <- c(1:20)

proportion_cutoff_for_colouring <- 30
enrichment_cutoff_for_colouring <- 5

aas_to_include <- c()
for (i in 1:nrow(enrichment)) {
  if (any(enrichment[i,2:ncol(enrichment)] > enrichment_cutoff_for_colouring)) {
    aas_to_include <- append(aas_to_include, i)
  }
}

# ####To apply colours to significantly enriched AAs only: ####
# 
# 
# if (any(composition[1, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Acolour <- "orange"} else {Acolour <- "grey"}
# if (any(composition[2, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Ccolour <- "yellow"} else {Ccolour <- "grey"}
# if (any(composition[3, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Dcolour <- "red"} else {Dcolour <- "grey"}
# if (any(composition[4, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Ecolour <- "red"} else {Ecolour <- "grey"}
# if (any(composition[5, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Fcolour <- "green"} else {Fcolour <- "grey"}
# if (any(composition[6, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Gcolour <- "orange"} else {Gcolour <- "grey"}
# if (any(composition[7, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Hcolour <- "blue"} else {Hcolour <- "grey"}
# if (any(composition[8, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Icolour <- "orange"} else {Icolour <- "grey"}
# if (any(composition[9, 2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Kcolour <- "blue"} else {Kcolour <- "grey"}
# if (any(composition[10,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Lcolour <- "orange"} else {Lcolour <- "grey"}
# if (any(composition[11,2:ncol(composition)] > proportion_cutoff_for_colouring)
#     || any(enrichment[11,2:ncol(enrichment)] > enrichment_cutoff_for_colouring)) {
#   Mcolour <- "yellow"} else {Mcolour <- "grey"}
# if (any(composition[12,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Ncolour <- "purple"} else {Ncolour <- "grey"}
# if (any(composition[13,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Pcolour <- "orange"} else {Pcolour <- "grey"}
# if (any(composition[14,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Qcolour <- "purple"} else {Qcolour <- "grey"}
# if (any(composition[15,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Rcolour <- "blue"} else {Rcolour <- "grey"}
# if (any(composition[16,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Scolour <- "pink"} else {Scolour <- "grey"}
# if (any(composition[17,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Tcolour <- "pink"} else {Tcolour <- "grey"}
# if (any(composition[18,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Vcolour <- "orange"} else {Vcolour <- "grey"}
# if (any(composition[19,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Wcolour <- "green"} else {Wcolour <- "grey"}
# if (any(composition[20,2:ncol(composition)] > proportion_cutoff_for_colouring)) {
#   Ycolour <- "green"} else {Ycolour <- "grey"}
# #######################################

plt2 <- composition[aas_to_include,] %>%
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(composition))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>%
  mutate(
    proportion = slider::slide_mean(proportion, before = 4, after = 4)
  ) %>% 

  ggplot() + 
  geom_rect(
    data = box_params,
    aes(xmin = start, xmax = end + 1, ymin = bottom, ymax = top, fill = region, group = ""), 
    color = "transparent",
    alpha = 0.15
  ) +
  # geom_point() +
  geom_line(aes(as.numeric(position), 
                proportion,
                color = AA), 
            size=2) +
  # geom_smooth(size=2, se=FALSE, method="gam", formula = y~poly(x,20)) +
  
  scale_color_manual(
    values = c(
      "A" = "orange",
      "C" = "yellow",
      "D" = "red",
      "E" = "red",
      "F" = "green",
      "G" = "orange",
      "H" = "cyan",
      "I" = "orange",
      "K" = "blue",
      "L" = "orange",
      "M" = "yellow",
      "N" = "purple",
      "P" = "hotpink",
      "Q" = "purple",
      "R" = "blue",
      "S" = "pink",
      "T" = "pink",
      "V" = "orange",
      "W" = "green",
      "Y" = "green"
      # "A" = Acolour,
      # "C" = Ccolour,
      # "D" = Dcolour,
      # "E" = Ecolour,
      # "F" = Fcolour,
      # "G" = Gcolour,
      # "H" = Hcolour,
      # "I" = Icolour,
      # "K" = Kcolour,
      # "L" = Lcolour,
      # "M" = Mcolour,
      # "N" = Ncolour,
      # "P" = Pcolour,
      # "Q" = Qcolour,
      # "R" = Rcolour,
      # "S" = Scolour,
      # "T" = Tcolour,
      # "V" = Vcolour,
      # "W" = Wcolour,
      # "Y" = Ycolour
    )
  ) +
  scale_x_continuous(
    expand = expansion()
  ) + 
  scale_y_continuous(
    expand = expansion()
  ) + 
  coord_cartesian(
    ylim = c(0, 55)
  ) +
  xlab("aa position (centre of window)") +
  ylab("frequency (%)") +
  labs(title = fasta$id, caption = paste("window width = ", window.width)) +
  theme_classic()

###########################################################################
# #### Normalise composition to aa composition of entire human proteome: ####
# 
# enrichment <- cbind(composition$AA,
#                     composition[-1]/proteome_composition$pc)
# 
# colnames(enrichment)[1] <- "AA"


# #### To plot enrichment for all 20 AAs on different panels: ####
# enrichment %>%
#   # this takes the data from wide format to long format
#   #   (best to avoid rownames, and instead work explicitly with a column with the
#   #    same ID)
#   pivot_longer(
#     cols = c(colnames(enrichment))[-1],
#     names_to = "position",
#     values_to = "proportion"
#   ) %>%
#   ggplot(mapping = aes(x = as.numeric(position), y = proportion)) +
#   geom_line() +
#   xlab("aa position (centre of window)") +
#   ylab("fold enrichment") +
#   labs(title = fasta$id, caption = paste("window width = ", window.width)) +
#   facet_wrap(~AA)


#### To plot enrichment for all 20 AAs on same graph: ####

plt3 <- enrichment[aas_to_include,] %>%
  # this takes the data from wide format to long format
  #   (best to avoid rownames, and instead work explicitly with a column with the 
  #    same ID)
  pivot_longer(
    cols = c(colnames(enrichment))[-1],
    names_to = "position",
    values_to = "proportion"
  ) %>%
  mutate(
    proportion = slider::slide_mean(proportion, before = 4, after = 4)
  ) %>% 


  ggplot() + 
  geom_rect(
    data = box_params,
    aes(xmin = start, xmax = end, ymin = bottom, ymax = top, fill = region, group = ""), 
    color = "transparent", 
    alpha = 0.15
  ) +
  # geom_point() +
  geom_line(aes(as.numeric(position), 
                proportion,
                color = AA), 
            size=2) +
  # geom_smooth(size=2, se=FALSE, method="gam", formula = y~poly(x,20)) +

  scale_color_manual(
    values = c(
      "P" = "brown",
      "Q" = "orange"
    )
  ) +
  scale_x_continuous(
    expand = expansion()
  ) + 
  scale_y_continuous(
    expand = expansion()
  ) + 
  
  xlab("aa position (centre of window)") +  
  ylab("fold enrichment") +
  labs(title = fasta$id, caption = paste("window width = ", window.width)) + 
  theme_classic()


#### patch desired plots together ####

patchwork::wrap_plots(
  plt1, plt2, ncol = 1
  
)
