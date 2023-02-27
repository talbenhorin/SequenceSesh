rm(list=ls(all=TRUE))

## Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")

## Kalle's dada2 fix
## install.packages("devtools")
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
#library("devtools")
#devtools::install_github("benjjneb/dada2")
#library("dada2")
#library(dada2); packageVersion("dada2")
#path <- "DADA2 tutorial//data//MiSeq_SOP"
#list.files(path)

## Load add-on packages
library(dada2); packageVersion("dada2")
library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(scales)
theme_set(theme_bw())

## Filtering sequences
pathF <- "FNB" 
pathR <- "RNB" 

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") 

fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: consider modifying these parameters if issues arise
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# Infer sequence variants
# File parsing
filtpathF <- "FNB/filtered" 
filtpathR <- "RNB/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

## Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF)

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "E:/tbenhor/libraries/Documents/SequenceSesh/Output/seqtab.rds") # CHANGE ME to where you want sequence table saved

## Remove chimeras
st.all <- readRDS("E:/tbenhor/libraries/Documents/SequenceSesh/Output/seqtab.rds")
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "E:/tbenhor/libraries/Documents/SequenceSesh/silva_nr_v132_train_set.fa.gz", multithread=TRUE)#get tax
tax <- addSpecies(tax, "E:/tbenhor/libraries/Documents/SequenceSesh/silva_species_assignment_v132.fa.gz")

taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## Save files
saveRDS(seqtab, "E:/tbenhor/libraries/Documents/SequenceSesh/Output/seqtab_final.rds") 
saveRDS(tax, "E:/tbenhor/libraries/Documents/SequenceSesh/Output/tax_final.rds") 

## Open Files if not using dada2 workflow
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
#  Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * nsamples(ps)

# Filter out prevalence
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps = prune_taxa(keepTaxa, ps)
ps = tax_glom(ps, "Phylum", NArm = TRUE) #glom the pruned taxa 
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

plot_bar(ps, fill="Phylum")

ordu <- ordinate(ps, "PCoA", "bray", weighted=FALSE)
p = plot_ordination(ps, ordu) + geom_point(size=7, alpha=0.75) + 
  geom_text(fontface = "bold", mapping = aes(label = samdf$date), size = 4, vjust = 2) +
  xlim(-0.5, 1.0) + ylim(-0.25, 0.25)


