rm(list=ls(all=TRUE))

## Install dada2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")

library(dada2); packageVersion("dada2")

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

# Sample inference and merger of paired-end reads
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
rm(derepF); rm(derepR)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "E:/tbenhor/libraries/Documents/SequenceSesh/Output/seqtab.rds") # CHANGE ME to where you want sequence table saved

# Remove chimeras
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
