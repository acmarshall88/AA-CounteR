library(bio3d)
library(ggplot2)
library(tidyverse)

# Read in human proteome composition table
proteome_composition <- read.csv("human_proteome_AAcomposition.csv")

# Read in protein sequence file
fasta <- read.fasta("hNONO_Q15233.fasta.txt")

# extract aa sequence (as character vector)
aa.seq <- as.vector(fasta$ali)

# List 20 standard amino acids (A through Y)
aa.list <- LETTERS[c(-2,-10,-15,-21,-24,-26)]


# this creates a "named vector", which we can use as a lookup table for different 
# aa values
comp_lookup <- pull(proteome_composition, pc, name = aa)

# helper funciton to classify a AA position into one of the 3 regions
classify_region <- function(x) {
  case_when(
    x %in% 1:52 ~ "nlcr", 
    x %in% 53:371 ~ "core", 
    x %in% 372:471 ~ "clcr"
  ) %>% 
    factor(levels = c("nlcr", "core", "clcr"))
}



# set up tibble / dataframe with all ~700 AA, their number and their cetegorised
# region of the protien
dat <- tibble(
  seq = aa.seq
) %>% 
  mutate(
    resnum = row_number(),
    region = classify_region(resnum)
  )


# do the region-wise summary for each AA in each region as a relative proportion
# to the overall genomic enrichment and plot it

dat  %>% 
  summarise(
    count = n(),
    .by = c(region, seq)
  ) %>% 
  mutate(
    perc = count / sum(count) * 100 / comp_lookup[seq], 
    .by = region
  ) %>% 


  ggplot(aes(seq, perc, fill = region)) + 
  geom_col(alpha = 0.5) + 
  facet_wrap(~region)  + 
  xlab("amino acid") +
  ylab("enrichment") +
  scale_y_log10(limits = c(0.1,10), breaks = c(0.1,0.3,1.0,3,10)) +
  theme_bw() -> plt1

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


# Create parameters for coloring background of plot(s) by region:
tibble(
  num = 1:707
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


# #### To plot composition for all 20 AAs on different panels: ####
# composition %>% 
#   # this takes the data from wide format to long format
#   #   (best to avoid rownames, and instead work explicitly with a column with the 
#   #    same ID)
#   pivot_longer(
#     cols = c(colnames(composition))[-1],
#     names_to = "position",
#     values_to = "proportion"
#   ) %>% 
#   ggplot(mapping = aes(x = as.numeric(position), y = proportion)) +
#   geom_line() +
#   xlab("aa position (centre of window)") +
#   ylab("frequency (%)") +
#   labs(title = fasta$id, caption = paste("window width = ", window.width)) +
#   facet_wrap(~AA)
# 

#### To plot composition for all 20 AAs on same graph: ####

# ***OPTION TO ONLY INCLUDE AAs OF INTEREST***:
aas_to_include <- c(13:14)

composition[aas_to_include,] %>%
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
      "P" = "darkorange",
      "Q" = "darkblue"
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
  # ggplot(aes(as.numeric(position), proportion, color = AA)) + 
  # # geom_point() +
  # geom_line(linewidth=2) +
  # # geom_smooth(linewidth=2, se=FALSE, method="gam", formula = y~poly(x,20)) +
  # 
  xlab("aa position (centre of window)") +
  ylab("frequency (%)") +
  labs(title = fasta$id, caption = paste("window width = ", window.width)) +
  theme_classic() -> plt2

###########################################################################
#### Normalise composition to aa composition of entire human proteome: ####

enrichment <- cbind(composition$AA,
                    composition[-1]/proteome_composition$pc)

colnames(enrichment)[1] <- "AA"
# 
# # Create parameters for coloring background of plot(s) by region:
# tibble(
#   num = 1:707
# ) %>% 
#   mutate(
#     region = classify_region(num)
#   ) %>% 
#   summarise(start = min(num), 
#             end = max(num), 
#             .by = region) %>% 
#   mutate(
#     bottom = 0, top = 100
#   ) -> box_params


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
# 

#### To plot enrichment for all 20 AAs on same graph: ####

enrichment[aas_to_include,] %>%
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
  theme_classic() -> plt3


#### patch desired plots together ####

patchwork::wrap_plots(
  plt2, plt1, ncol = 1
  
)
