library(ggplot2)
library(plyr)
library(dplyr)
library(ggsci)
library(scales)
library(Rmisc)
library(ggbeeswarm)


input.dir <- #input directory
files.list <- list.files(path = input.dir, pattern = "*.txt", full.names = TRUE)

#set threshold for core unitigs
upper.threshold <- 0.95


#read in table and manipulate
ready_table <- function(file.name) {
  #read file
  full.table <- read.table(file = file.name, sep = '\t', header = TRUE)

  #parse species and tool names
  file.base <- basename(file.name)
  split.name <- strsplit(file.base, "_")[[1]]
  species <- split.name[[1]]
  type <- split.name[[2]]

  #edit species names
  if(species  == 'SP') {
    species <- 'S. pneumoniae'
  } else if (species == 'List') {
    species <- 'L. monocytogenes'
  } else {
    species <- 'S. enterica'
  }

  #remove node id field
  full.table <- full.table[,-1]

  #add in species, tool and graph type
  full.table$Species = species
  full.table$Type = type

  #caculate if genes in core or accessory genome
  full.table$gene_type = ifelse(full.table$Allele_freq >= upper.threshold, "core", "accessory")


  return(full.table)
}

#create empty list to add to
data.list = list()

#iteratively add all data to full.table
for (i in 1:length(files.list)) {
  data.list[[i]] <- ready_table(files.list[[i]])
}


#bind all tables together, generate tally of allele
dataset <- dplyr::bind_rows(data.list)

#add factors and levels to dataset
dataset$Species <- factor(dataset$Species, levels = c("L. monocytogenes", "S. enterica", "S. pneumoniae"))
dataset$Type <- factor(dataset$Type, levels = c("full", "rep"))
dataset$gene_type <- factor(dataset$gene_type, levels = c("accessory", "core"))

#create bins for plot
dataset <- dataset %>% mutate(bin=cut_width(Allele_freq, width = 0.05, boundary = 0))
bin_labels = c("[0,0.05]" = "0", "(0.05,0.1]" = "0.05", "(0.1,0.15]" = "0.1", "(0.15,0.2]" = "0.15", "(0.2,0.25]" = "0.2", "(0.25,0.3]" = "0.25", "(0.3,0.35]" = "0.3", "(0.35,0.4]" = "0.35","(0.4,0.45]" = "0.4", "(0.45,0.5]" = "0.45", "(0.5,0.55]" = "0.5", "(0.55,0.6]" = "0.55", "(0.6,0.65]" = "0.6", "(0.65,0.7]" = "0.65", "(0.7,0.75]" = "0.7", "(0.75,0.8]" = "0.75", "(0.8,0.85]" = "0.8", "(0.85,0.9]" = "0.85", "(0.9,0.95]" = "0.9", "(0.95,1]" = "0.95")


#unitig frequency histogram

overall.hist.plot <- ggplot(dataset, aes(x=Allele_freq)) + geom_histogram(binwidth = 0.025, aes(color = 'darkred', fill = 'darkred')) + facet_grid(Species ~ Type, scales = "free_y") + theme_light() + xlab("Unitig frequency") + ylab("Number of unitigs") + scale_fill_manual(values = c("gray70")) + scale_color_manual(values = c("Black")) + theme(text = element_text(size=20), legend.position = "none") # + scale_x_continuous(minor_breaks = seq(0, 1, 0.05), breaks = seq(0, 1, by = 0.2), limits=c(0,1)) + scale_y_continuous(minor_breaks = seq(0, 12.5, 1), breaks = seq(0, 12.5, by = 2), limits=c(0,12.5))

overall.hist.plot
