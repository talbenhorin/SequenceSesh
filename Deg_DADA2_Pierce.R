rm(list=ls(all=TRUE))
sessionInfo()
## Install dada2 package
install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")
install.packages("dada2")
library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(scales)
theme_set(theme_bw())

## Filtering sequences
pathF <- "Deg_Forward" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "Deg_Reverse" # CHANGE ME to the directory containing your demultiplexed reverse-read fastqs
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqFs
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
fastqRs
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: consider modifying these parameters if issues arise
out<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(240,200),trimLeft = 19, trimRight = 20, maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# Infer sequence variants
# File parsing
filtpathF <- "Deg_Forward/filtered" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "Deg_Reverse/filtered"
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtFs
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.names
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Inspect Quality
plotQualityProfile(filtFs)
plotQualityProfile(filtRs)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

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
dim(seqtab)
table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "Deg_Seqtab.rds") # CHANGE ME to where you want sequence table saved

# Remove chimeras
st.all <- readRDS("Deg_seqtab.rds")
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)#get tax
tax <- addSpecies(tax, "silva_species_assignment_v132.fa.gz")
taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Track Reads Through Pipeline
getN <- function(x) sum(getUniques(x))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(out, sapply(filtFs, getN), sapply(filtRs, getN), sapply(mergers, getN), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track

## Save files
saveRDS(seqtab, "Deg_Seqtab_Final.rds") 
saveRDS(tax, "Deg_Tax_Final.rds") 
seqtab <- readRDS("Deg_Seqtab_Final.rds")
taxa <- readRDS("Deg_Tax_Final.rds")
samples.out <- rownames(seqtab)

##Adding Metadata
mdata <- read.csv("Deg_Metadata.csv", fileEncoding= "UTF-8-BOM", fill = FALSE, header = TRUE) 
samdf <- data.frame(Date=mdata$Date,Site=mdata$Site,ID=mdata$Sample) 
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
ps1 = prune_taxa(keepTaxa, ps)
ps1
ps2 = tax_glom(ps1, "Genus", NArm = TRUE) #glom the pruned taxa, glomming condenses taxonomic data at specific level

##Plot top 10 genera bar plot
top10 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:10]
ps2.top10 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top10 <- prune_taxa(top10, ps2.top10)
ps2.top10
names <- taxa_names(ps2.top10)
names
plot_bar(ps2.top10, fill="Genus")+ 
  scale_fill_brewer(palette='Paired')+
  theme_classic(base_size=24)+
  scale_x_discrete(limits=c("AH1T0", "AH2T0", "AH3T0", "CR3T0", "CR8", "CR10", "CR14", "IR1T0", "IR2T0", "IR3T0", "LL1T0", "LL2T0", "LL3T0", "MC1T0", "MC2T0","MC3T0", 'CH1T0','CH2T0','CH3T0'))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')


##Plot top 20 genera bar plot
top20<-names(sort(taxa_sums(ps2), decreasing=TRUE))[1:20]
ps2.top20 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2.top20 <- prune_taxa(top20, ps2.top20)
names20<-taxa_names(ps2.top20)
plot_bar(ps2.top20, fill="Genus")

##Plot top 10 genera heat map
p1 = plot_heatmap(ps2.top10, taxa.label = "Genus",low="white", high="#000033", 
                  na.value = "white",taxa.order = taxa_names(ps2.top10),
                  trans = identity_trans())
p1 = p1  + theme(axis.text.x = element_text(size=8, angle=30, hjust=1, vjust=0.9),
                 legend.title = element_text(size = 10)) +
  labs(fill = "Relative\nabundance")
p1

##Plot NMDS plot
ord.nmds.bray <- ordinate(ps2, method="NMDS", distance="bray")
nmds<-plot_ordination(ps2, ord.nmds.bray, color="Site", shape="Date", title=" ")
nmds = nmds + geom_point(size=7, alpha=0.75)
nmds= nmds+ scale_color_brewer(palette="Set1")
nmds = nmds + geom_text(mapping = aes(label = samdf$ID), size = 4, vjust = 1.5) 
nmds = nmds+theme_bw(base_size = 24)
nmds

ps2.top10@tax_table
##Plot Species Richness Data
plot_richness(ps, measures=c("Shannon", "Simpson"), color="Site")+scale_color_brewer(palette = "Set1")

MCY_Results<-read.csv("Deg_MCY_Abun.csv",header=TRUE,sep=',')
ggplot(MCY_Results)+geom_point(aes(x=Rel_Abun,y=log(Tot_MCY),color=Site,shape=Date))+
  geom_abline(slope=10,intercept=0)
  #geom_smooth(aes(x=Avg_Rel_Abun,y=Avg_Tot_MCY), method='lm', formula=(y~exp(x)))
Exp.Model<-lm((log(Avg_Tot_MCY)~Avg_Rel_Abun), MCY_Results)
summary(Exp.Model)  

x<-seq(0,0.7,0.01)
y=(4.473*exp(9.1453*x))
plot(x,y)
MCY_Res<-rbind(MCY_Results[2:6,],MCY_Results[10:16,],MCY_Results[19,])
TTox_label<-expression(Total ~ Microcystin ~ (µg ~ L^"-1"))
scatplot<-ggplot()+geom_point(data=MCY_Res,aes(x=Rel_Abun,y=Tot_MCY,color=Site,shape=Date), size=5)+
  scale_color_brewer(palette='Set1')+
  geom_line(aes(x=x,y=y))+
  ylim(0,1200)+
  ylab(TTox_label)+
  xlab('Rel. Abundance of Microcystis')+
  theme_classic(base_size = 24)
scatplot
plot(x,y,type="l")
ggarrange(scatplot,nmds,
          nrow = 1,ncol = 2,common.legend = TRUE,
          legend = 'right')
#geom_smooth(aes(x=Avg_Rel_Abun,y=Avg_Tot_MCY), method='lm', formula=(y~exp(x)))
