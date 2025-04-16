setwd("/Volumes/eswanner-lab/Michelle/April2023_Sequences_DADA2/DL2022S")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.19")
BiocManager::install("Rcpp")
BiocManager::install("microbiome")
BiocManager::install("phyloseq", force = TRUE)
install.packages("vegan")
library("vegan")
install.packages("viridis")
library(viridis)
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 
#change the ref argument to get other versions
library(dada2)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse)
library(patchwork)

#set path to working directory
path <- ("/Volumes/eswanner-lab/Michelle/April2023_Sequences_DADA2/DL2022S")
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample name based on filename format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###Check the quality of your profile in samples###
#Used visual inspection of the quality score heatmaps to determine where to truncate.

#Quality of forward reads
plotQualityProfile(fnFs[1:2])
#Quality of reverse reads
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##Filtering and trimming the reads##
#Truncated based on quality scored from ines 60-63
#MaxN: Default 0.
#MaxEE: After truncation, reads with higher than maxEE "expected errors" will be discarded. Standard is (2,5)
#TrunQ:  Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.
#rm.phix Default is True. Removes sequences that match with PhiX genome
#compress: Default TRUE. Whether the output fastq file should be gzip compressed.
#Multithread: Set to TRUE, multithreading is enabled and the number of available threads is automatically determined. 
#On Windows set multithread=FALSE
#trimLeft=c(,) is set equal to the length of the forward and reverse primers


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(190,225),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = c(19, 19)) # On Windows set multithread=FALSE

head(out)


#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot the error rates. Red line represents expected error rates.
plotErrors(errF, nominalQ=TRUE)

#Sample inference. See doi.org/10.1038/nmeth.3869 for info. Returns dada-class object. see help("dada-class") for some info,
#including multiple diagnostics about the quality of each de-noised sequence variant

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Merge paired reads
#verbose: If TRUE, print status to standard output.

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample. #Inspect the merger data.frame from the first sample. Most of your reads should successfully merge. 
#If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?
head(mergers[[1]])

#construct ASV table.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#frequency of chimeras 
sum(seqtab.nochim)/sum(seqtab)

#Sanity Check: Look at # of reads that made it through filtering (head(track)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy using CyanoSeq

taxa <- assignTaxonomy(seqtab.nochim, "/Volumes/eswanner-lab/Michelle/April2023_Sequences_DADA2/CyanoSeq_1.2_SILVA138.1_dada2.fastq",
                       multithread=TRUE)
# Inspect taxonomic assignments
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)



#### Creating phyloseq object from DADA2 ####

theme_set(theme_bw())

#Make data frame based on information in filenames. 
samples.out <- rownames(seqtab.nochim)

#Read in metadata file
MetaData_SequencingDL2022 <- data.frame(read.csv("MetaData_Sequencing_DL2022.csv"))
MetaData_SequencingDL2022

#Assigns rownames in metadata data.frame to the corresponding sample names in the phyloseq object
rownames(MetaData_SequencingDL2022) <-MetaData_SequencingDL2022$SampleID

#Create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(MetaData_SequencingDL2022), 
               tax_table(taxa))

#Filter physloseq object
ps2 <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
ps2


dna <- Biostrings::DNAStringSet(taxa_names(ps2))
names(dna) <- taxa_names(ps2)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps2) <- paste0("ASV", seq(ntaxa(ps2)))

#transform into long-format dataframe
gg_ps2 <- psmelt(ps2) 

# Filtering out Cyanobacteria ASVs using dplyr::filter
Cyanobacteria <- dplyr::filter(gg_ps2, Phylum %in% c("Cyanobacteriota"))
Cyanobacteria$Type <- c("Cyanobacteria")

# Filtering out purple non-sulfur bacteria
PnSB <- dplyr::filter(gg_ps2, Family %in% c("Rhodospirillaceae", "Rhodopila",
                                            "Rhodovastum", "Hyphomicrobiaceae", "Rhodobacteraceae", "Rhodoferax", 
                                            "Rhodopseudomonas", "Rhodoblastus"))
PnSB$Type <- c("PnSB")

#Filtering out purple sulfur bacteria
PSB <- dplyr::filter(gg_ps2, Family %in% c("Chromatiaceae", "Ectothiorhodospiraceae"))
PSB$Type <- c("PSB")

#Filtering out green sulfur bacteria
GSB <- dplyr::filter(gg_ps2, Class %in% c("Chlorobia"))
GSB$Type <- c("GSB")

#stitch dataframes together into new dataframe Phototrophs
Phototrophs <- rbind(Cyanobacteria, PnSB, PSB, GSB)

# function to normalize abundance to 100%

test1 <- function(x) {(x/sum(x)*100)}

#Filter out phototroph sequences from 3.5 meter depth
Depth3.5meter<- filter(Phototrophs, Depth == 3.5)
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth3.5normalize<-Depth3.5meter %>%
  mutate(Normabund = test1(Abundance))
summarize(Depth3.5meter, Abundance)  

#Filter out phototroph sequences from 4 meter depth
Depth4meter <-filter(Phototrophs, Depth == 4| Depth == 4)
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth4normalize<-Depth4meter %>%
  mutate(Normabund = test1(Abundance))

#Filter out phototroph sequences from 5.25 meter depth
Depth5.25meter <-filter(Phototrophs, Depth == 5.25)
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth5.25normalize<-Depth5.25meter %>%
  mutate(Normabund = test1(Abundance))

#Filter out sequences from 5.75 meter depth
Depth5.75meter <-filter(Phototrophs, Depth == 5.75)
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth5.75normalize<-Depth5.75meter %>%
  mutate(Normabund = test1(Abundance))

#Filter out sequences from 6 meter depth
Depth6meter <-filter(Phototrophs, Depth == 6.0) 
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth6normalize<-Depth6meter %>%
  mutate(Normabund = test1(Abundance))

#Filter out sequences form 6.5 meter depth
Depth6.5meter <-filter(Phototrophs, Depth == 6.5) 
#Normalize to 100% and add normalized abundance to dataframe as "Normabund"
Depth6.5normalize<-Depth6.5meter %>%
  mutate(Normabund = test1(Abundance))

#Stich together dataframes from each depth that have normalized abundances
normalizedPhototrophs <- rbind(Depth3.5normalize, Depth4normalize, Depth5.25normalize,
                               Depth5.75normalize, Depth6normalize, Depth6.5normalize)

# MAKING PLOT#
cbbPalette <- c('#1E88E5','#D81B60',"#FFC107", "#004D40")
Phototroph_plot<- ggplot(normalizedPhototrophs, aes(x = Depth, y = Normabund, fill = Type)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cbbPalette ) + #this is an alternative color palette tool
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis"
  labs(x = "Depth (m)", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  coord_flip() +
  scale_x_reverse(breaks =c(1,2,3,4,5,6))

pdf("2022DL_AugPhototroph_plot.pdf", width = 5, height = 5)
Phototroph_plot
dev.off()


###Filter out Cyanobacteria sequences from each depth and normalize to 100%
CyDepth3meter<- filter(Phototrophs, Depth == 3.5) %>%
  filter(Type == "Cyanobacteria")
CyDepth3normalize<-CyDepth3meter %>%
  mutate(Normabund = test1(Abundance))
summarize(CyDepth3meter, Abundance)  

CyDepth4meter <-filter(Phototrophs, Depth == 4) %>%
  filter(Type == "Cyanobacteria")
CyDepth4normalize<-CyDepth4meter %>%
  mutate(Normabund = test1(Abundance))

CyDepth5.25meter <-filter(Phototrophs, Depth == 5.25) %>%
  filter(Type == "Cyanobacteria")
CyDepth5.25normalize<-CyDepth5.25meter %>%
  mutate(Normabund = test1(Abundance))

CyDepth5.75meter <-filter(Phototrophs, Depth == 5.75) %>%
  filter(Type == "Cyanobacteria")
CyDepth5.75normalize<-CyDepth5.75meter %>%
  mutate(Normabund = test1(Abundance))

CyDepth6meter <-filter(Phototrophs, Depth == 6) %>%
  filter(Type == "Cyanobacteria")
CyDepth6normalize<-CyDepth6meter %>%
  mutate(Normabund = test1(Abundance))

CyDepth6.5meter <-filter(Phototrophs, Depth == 6.5) %>%
  filter(Type == "Cyanobacteria")
CyDepth6.5normalize<-CyDepth6.5meter %>%
  mutate(Normabund = test1(Abundance))

#Stich together Cyanobacteria dataframes from each depth that are normalized 
CynormalizedPhototrophs <- rbind(CyDepth3normalize, CyDepth4normalize, CyDepth5.25normalize,
                                 CyDepth5.75normalize, CyDepth6normalize, CyDepth6.5normalize)

#Make plot of Cyanobacteria orders
Cyorder_plot<- ggplot(CynormalizedPhototrophs, aes(x = Depth, y = Normabund, fill = Order)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) + #this is an alternative color palette tool
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis"
  labs(x = "Depth (m)", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  coord_flip() +
  scale_x_reverse(breaks =c(1,2,3,4,5,6))

pdf("2022_DLAugCyorder_plot.pdf", width = 5, height = 5)
Cyorder_plot
dev.off()

#Stitch together the two plots for one figure in manuscript
Sequencing_plot <- Phototroph_plot + Cyorder_plot + 
  plot_layout(ncol = 2)  # Stacks them in one column. Change to ncol = 3 for side by side

# Display the combined plot
Sequencing_plot
ggsave("Sequencing_plot.pdf", width = 9, height = 5, device = "pdf")
Sequencing_plot
dev.off()
dev.off()
