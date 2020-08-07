library(ggplot2)
library(plyr)
library(dplyr)
library(ggsci)
library(scales)
library(Rmisc)
library(ggbeeswarm)


input.dir <- #input directory
files.list <- list.files(path = input.dir, pattern = "*.Rtab", full.names = TRUE)

#read in table and manipulate
ready_table <- function(file.name) {
  #read file
  rtab <- read.table(file = file.name, sep = '\t', header = TRUE)

  #count number of input files
  total.genomes <- ncol(rtab) - 1

  #parse species and tool names
  file.base <- basename(file.name)
  split.name <- strsplit(file.base, "_")[[1]]
  species <- split.name[[1]]
  type <- split.name[[2]]

  #edit species names
  if(species  == 'Pneumo') {
    species <- 'S. pneumoniae'
  } else if (species == 'List') {
    species <- 'L. monocytogenes'
  } else {
    species <- 'S. enterica'
  }


  #calculate gene frequency and append to dataframe
  rtab$Gene_freq <- ((rowSums(rtab[, -1])) / total.genomes)

  #add in species, tool and graph type
  rtab$Species = species
  rtab$Type = type

  new.rtab <- data.frame(rtab$Gene, rtab$Species, rtab$Type, rtab$Gene_freq)
  names(new.rtab) <- c("Gene", "Species", "Type", "Gene_freq")

  return(new.rtab)
}

#create empty list to add to
data.list = list()

#iteratively add all data to full.table
for (i in 1:length(files.list)) {
  data.list[[i]] <- ready_table(files.list[[i]])
}

#bind all tables together, generate tally of allele
dataset <- dplyr::bind_rows(data.list)

dataset$Species <- factor(dataset$Species, levels = c("L. monocytogenes", "S. enterica", "S. pneumoniae"))
dataset$Type <- factor(dataset$Type, levels = c("full", "rep"))

#bin allele frequencies for histogram plot
dataset <- dataset %>% mutate(bin=cut_width(Gene_freq, width = 0.05, boundary = 0))
bin_labels = c("[0,0.05]" = "0", "(0.05,0.1]" = "0.05", "(0.1,0.15]" = "0.1", "(0.15,0.2]" = "0.15", "(0.2,0.25]" = "0.2", "(0.25,0.3]" = "0.25", "(0.3,0.35]" = "0.3", "(0.35,0.4]" = "0.35","(0.4,0.45]" = "0.4", "(0.45,0.5]" = "0.45", "(0.5,0.55]" = "0.5", "(0.55,0.6]" = "0.55", "(0.6,0.65]" = "0.6", "(0.65,0.7]" = "0.65", "(0.7,0.75]" = "0.7", "(0.75,0.8]" = "0.75", "(0.8,0.85]" = "0.8", "(0.85,0.9]" = "0.85", "(0.9,0.95]" = "0.9", "(0.95,1]" = "0.95")

#gene frequency histogram

overall.hist.plot <- ggplot(dataset, aes(x=Gene_freq)) + geom_histogram(binwidth = 0.025, aes(color = 'darkred', fill = 'darkred')) + facet_grid(Species ~ Type, scales = 'free_y') + theme_light() + xlab("Gene frequency") + ylab("Number of genes") + scale_fill_manual(values = c("gray70")) + scale_color_manual(values = c("Black")) + theme(text = element_text(size=20), legend.position = "none") #+ scale_x_continuous(minor_breaks = seq(0, 1, 0.05), breaks = seq(0, 1, by = 0.2), limits=c(0,1)) + scale_y_continuous(minor_breaks = seq(0, 12.5, 1), breaks = seq(0, 12.5, by = 2), limits=c(0,12.5))

overall.hist.plot
