##SKIP TO LINE 104 AFTER LOADING LIBRARIES!!!!!!!!!!!!!!!!!!!!!!!!!!
#Saving File into GitHub repository
rm(list=ls(all=TRUE))
sessionInfo()
## Install dada2 package
install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2); packageVersion("dada2")
install.packages("dada2")
library(dada2)
library(ape)
library(gridExtra)
library(knitr)
library(phyloseq); packageVersion("phyloseq")
citation('phyloseq')
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(scales)
library(ggpubr)
library(tidyr)
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
tax <- assignTaxonomy(seqtab, "silva_nr99_v138.1_train_set.fa", multithread=TRUE)#get tax
tax <- addSpecies(tax, "silva_species_assignment_v138.1.fa")
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
#When Rerunning from an empty environment START HERE and make PierceV4V5 working dir.
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
nsamples(ps)
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
ps2 = tax_glom(ps1, "Genus", NArm = TRUE) #glom the pruned taxa, glomming condenses taxonomic data at specific level
ps0 = tax_glom(ps, "Genus", NArm = TRUE)
otu_table(ps0)

##Attempting to export phyloseq object as dataframe for ESA
OTU1 = as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf
write.csv(OTUdf,"OTU.csv")
##OTU.csv gives you a count of each ASV/OTU in each sample
TAX = as(tax_table(ps), "matrix")
if(taxa_are_rows(ps)){TAX <- t(TAX)}
TAXdf = as.data.frame(TAX)
TAXdf
write.csv(TAXdf,"Taxa.csv")
##Taxa.csv gives the taxonomic assignment for each unique ASV/OTU

#transform read counts to relative abundance data
ps0.rel <- transform_sample_counts(ps0, function(OTU) OTU/sum(OTU))
OTU2 = as(otu_table(ps0.rel), "matrix")
if(taxa_are_rows(ps0.rel)){OTU2 <- t(OTU2)}
OTUdf2 = as.data.frame(OTU2)
write.csv(OTUdf2,"ps0.rel.csv")
##ps0.rel.csv gives the relative abundance of each genus in each sample

#Gather data into long format for easier reading and plotting
OTU.long<-read.csv("ps0.rel.csv",header=TRUE,sep=",")
OTU.long<-gather(OTU.long,Genus,Rel.Abun,Microcystis.PCC.7914:Sphingopyxis,factor_key=TRUE)
write.csv(OTU.long, "ps0.rel.long.csv")
#I manually edited the sheet to include a column for site
OTU.long.2<-read.csv("ps0.rel.long.csv",header=TRUE,sep=",")

##Load in color palette and establish colors for plotting
library(paletteer) 
paletteer_d("ggthemes::Miller_Stone")
new_p<-paletteer_d("ggthemes::Miller_Stone")
new_p3<-c("#4F6980FF","#B66353FF","#7E756DFF")

##Plot top 10 genera bar plot
top10 <- names(sort(taxa_sums(ps0), decreasing=TRUE))[1:10]
ps.top10 <- transform_sample_counts(ps0, function(OTU) OTU/sum(OTU))
ps.top10 <- prune_taxa(top10, ps.top10)
names <- taxa_names(ps.top10)
plot_bar(ps.top10, fill="Genus")+ 
  theme_classic()+
  scale_x_discrete(limits=c("AH1T0", "AH2T0", "AH3T0", "CR3T0", "CR8", "CR10", "CR14", "IR1T0", "IR2T0", "IR3T0", "LL1T0", "LL2T0", "LL3T0", "MC1T0", "MC2T0","MC3T0", 'CH1T0','CH2T0','CH3T0'))+
  ylab('Relative Abundance (%)')+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

OTUtop10 = as(otu_table(ps.top10), "matrix")
if(taxa_are_rows(ps.top10)){OTUtop10 <- t(OTUtop10)}
OTUtop10df = as.data.frame(OTUtop10)
write.csv(OTUtop10df,"OTUtop10.csv")
TAXtop10 = as(tax_table(ps.top10), "matrix")
if(taxa_are_rows(ps.top10)){TAXtop10 <- t(TAXtop10)}
TAXtop10df = as.data.frame(TAXtop10)
write.csv(TAXtop10df,"Taxatop10.csv")
##At this step I manually transposed the names of the top 10 taxa as the column names for the top 10 OTUs and added a column for Other 

#Gather top 10 asv/otu data into long format
OTUtop10.long<-read.csv("OTUtop10.csv",header=TRUE,sep=",")
OTUtop10.long<-gather(OTUtop10.long,Genus,Rel.Abun,Microcystis.PCC.7914:Other,factor_key=TRUE)

#Make genus a factor with specified levels inluding 'Other'
levels(as.factor(OTUtop10.long$Genus))
Genus.Fill<-factor(OTUtop10.long$Genus, levels = c("Other",
                                     "Ellin6067","Flavobacterium","Flectobacillus", "Microcystis.PCC.7914",
                                     "Paucibacter","Phenylobacterium","Pseudanabaena.PCC.7429","Roseomonas","Sediminibacterium","Vogesella"))

##plot top 10 by sample
fig.5<-ggplot(OTUtop10.long)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus.Fill))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr(base_size = 12)+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(limits=c("AH1T0", "AH2T0", "AH3T0", "CR3T0", "CR8", "CR10", "CR14", "IR1T0", "IR2T0", "IR3T0", "LL1T0", "LL2T0", "LL3T0", "MC1T0", "MC2T0","MC3T0", 'CH1T0','CH2T0','CH3T0'))+
  theme(axis.text.x = element_text(angle=60, hjust=1))
fig.5
ggsave("fig.5.prospectus.pdf",plot=fig.5,width=8,height=5,dpi=300)

#Plot individual sample triplicates
ggplot(OTU.long.2)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
AH<-subset(OTU.long.2,Site%in%"AH")
AH<-subset(AH, Rel.Abun>0)
ggplot(AH)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
AH.5<-subset(AH, Rel.Abun>0.05)
ggplot(AH.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("AH1T0", "AH2T0", "AH3T0"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

CH<-subset(OTU.long.2,Site%in%"CH")
CH<-subset(CH, Rel.Abun>0)
ggplot(CH)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
CH.5<-subset(CH, Rel.Abun>0.05)
ggplot(CH.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("CH1T0", "CH2T0", "CH3T0"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

CR<-subset(OTU.long.2,Site%in%"CR")
CR<-subset(CR, Rel.Abun>0)
ggplot(CR)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
CR.5<-subset(CR, Rel.Abun>0.05)
ggplot(CR.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("CR3T0", "CR8", "CR10", "CR14"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

LL<-subset(OTU.long.2,Site%in%"LL")
LL<-subset(LL, Rel.Abun>0)
ggplot(LL)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
LL.5<-subset(LL, Rel.Abun>0.05)
ggplot(LL.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("LL1T0", "LL2T0", "LL3T0"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

IR<-subset(OTU.long.2,Site%in%"IR")
IR<-subset(IR, Rel.Abun>0)
ggplot(IR)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
IR.5<-subset(IR, Rel.Abun>0.05)
ggplot(IR.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("IR1T0", "IR2T0", "IR3T0"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

MC<-subset(OTU.long.2,Site%in%"MC")
MC<-subset(MC, Rel.Abun>0)
ggplot(MC)+geom_tile(aes(x=Sample,y=Genus,fill=Rel.Abun))
MC.5<-subset(MC, Rel.Abun>0.05)
ggplot(MC.5)+geom_col(aes(x=Sample,y=Rel.Abun,fill=Genus))+
  scale_x_discrete(limits=c("MC1T0", "MC2T0", "MC3T0"))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab('Relative Abundance (%)')+
  labs(fill="Genus")+
  scale_fill_manual(values=new_p)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,.25,.50,.75,1), labels=c(0,25,50,75,100))+
  scale_x_discrete(expand=c(0,0))

##Plot top 10 genera heat map
p1 = plot_heatmap(ps.top10, taxa.label = "Genus",low="white", high="#000033", 
                  na.value = "white",taxa.order = taxa_names(ps.top10),
                  trans = identity_trans())
p1 = p1  + theme(axis.text.x = element_text(size=8, angle=30, hjust=1, vjust=0.9),
                 legend.title = element_text(size = 10)) +
  labs(fill = "Relative\nabundance")
p1

##Plot NMDS plot
ord.nmds.bray <- ordinate(ps0, method="NMDS", distance="bray")
ord.nmds.bray
nmds<-plot_ordination(ps0, ord.nmds.bray, color="Site", shape="Date", title=" ")
nmds = nmds + geom_point(size=7, alpha=0.75)
nmds= nmds+ scale_color_manual(values=new_p)
nmds = nmds + annotate('text',label='Stress = 0.07', x=1.5,y=0.6,size = 4, vjust = 1.5) 
nmds = nmds+theme_bw(base_size = 12)
fig.6<-nmds
fig.6
ggsave("fig.6.prospectus.pdf",plot=fig.6,width=7,height=5,dpi=300)

##Plot Species Richness Data
plot_richness(ps0, measures="Shannon", color="Site")+scale_color_manual(values=new_p)+
  theme_pubr()

##Load in Toxin Data to compare MC conc and microcystis rel. abund.()
MCY_Results<-read.csv("Deg_MCY_Abun.csv",header=TRUE,sep=',')
MCY_Res<-rbind(MCY_Results[3,],MCY_Results[6,],MCY_Results[10,],MCY_Results[13,],MCY_Results[16,],MCY_Results[19,])
TTox_label<-expression(Total ~ Microcystin ~ (µg ~ L^"-1"))
MCY_Results_na<-subset(MCY_Results, Part_MCY>0)

#Build log-transformed model
Exp.Model<-lm(log(Tot_MCY)~Avg_Rel_Abun,data=MCY_Results_na)
summary(Exp.Model)

#Plot log-transformed model
ggplot(data=MCY_Results_na)+geom_point(aes(x=Rel_Abun,y=log(Tot_MCY),color=Site,shape=Date), size=5)+
  geom_smooth(aes(x=Rel_Abun,y=Tot_MCY), method='lm', formula=(log(y)~x))+
  scale_color_manual(values=new_p)

#input exponential model
x<-seq(0,0.7,0.01)
y=(2.3845*exp(10.147*x))
plot(x,y)

#plot exponential model
fig.7<-ggplot()+geom_point(data=MCY_Results_na,aes(x=Rel_Abun,y=Tot_MCY,color=Site,shape=Date), size=5)+
  geom_line(aes(x=x,y=y))+
  ylab(TTox_label)+
  xlab('Average Relative Abundance of Microcystis')+
  scale_color_manual(values=new_p)+
  scale_y_continuous(expand =c(0,0),limits=c(0,1200))+
  theme_pubr(base_size = 12)+
  theme(axis.ticks.y = element_line(size=1), axis.ticks.length=unit(5,"pt"))
fig.7
ggsave("fig.7.prospectus.pdf",plot=fig.7,width=7,height=5,dpi=300)
