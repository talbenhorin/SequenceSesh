rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("dada2", version = "3.11")

library(dada2); packageVersion("dada2")
library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(scales)
theme_set(theme_bw())

seqtab <- readRDS("output/seqtab_final.rds")
taxa <- readRDS("output/tax_final.rds")

samples.out <- rownames(seqtab)

samps <- read.csv("samps.csv", fill = FALSE, header = TRUE) 
samdf <- data.frame(date=samps$date,ID=samps$sample) 
rownames(samdf) <- samples.out

## "Phyloseq" object from OTU table
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

## Filtering
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.1 * nsamples(ps)

# Filter out prevalence
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps1 = prune_taxa(keepTaxa, ps)
ps2 = tax_glom(ps1, "Genus", NArm = TRUE) #glom the pruned taxa 

top10 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:10]
ps2.top10 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top10 <- prune_taxa(top10, ps2.top10)

names <- taxa_names(ps2.top10)

plot_bar(ps, fill="Genus")