#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2); packageVersion("dada2")
# I have version 1.14.0
library(ShortRead)
packageVersion("ShortRead")
# 1.44.3
library(Biostrings)
packageVersion("Biostrings")
# 2.54.0
path <- "~/Desktop/TVE/16sonly/16s_seqs" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))

#                     Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       0
# REV.ReverseReads       0          0       0       0

#### Visualizing raw data ####

#Quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1,2,3,4)])
plotQualityProfile(fnFs.filtN[c(51,52,53,54)])
#looks mostly good up to 200

#Quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1,2,3,4)])
plotQualityProfile(fnRs.filtN[c(51,52,53,54)])
#180

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 200 bp for both forward & reverse [leaves ~50bp overlap], added "trimleft" to cut off primers [18 for forward, 20 for reverse]
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leaves ~50bp overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(20,21), #N nucleotides to remove from the start of each read
                     minLen = 50,
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
#               reads.in reads.out
# I2A_R1.fastq    10692      4248
# I2B_R1.fastq    18374     14180
# I2C_R1.fastq    12345      9598
# I2D_R1.fastq    28595     26537
# I2E_R1.fastq    38250     33551
# I2F_R1.fastq    24383     21328
tail(out)
#               reads.in reads.out
# O4D_R1.fastq    31528     27031
# O4E_R1.fastq    31604     27057
# O4F_R1.fastq    42131     37694
# O4G_R1.fastq    16493     14698
# O4H_R1.fastq    31566     28558
# O4I_R1.fastq    13211     12071

#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
#100279800 total bases in 501399 reads from 24 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#100693080 total bases in 559406 reads from 26 samples will be used for learning the error rates.

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

#~############################~#
##### Dereplicate reads ########
#~############################~#
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 148 sequence variants were inferred from 1936 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaRs[[1]]
# dada-class: object describing DADA2 denoising results
# 140 sequence variants were inferred from 1844 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 244-264 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(244,264)] #again, being fairly conservative wrt length

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab2.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 968 bimeras out of 9911 input sequences.
dim(seqtab2.nochim)
#[1]   54 8943

sum(seqtab2.nochim)/sum(seqtab)
#0.9596137
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 


#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="ss16s_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

taxa <- assignTaxonomy(seqtab2.nochim, "~/Desktop/TVE/RefData/silva_nr_v138_train_set.fa.gz",tryRC=TRUE)
unname(head(taxa))
taxa.plus <- addSpecies(taxa, "~/Downloads/silva_species_assignment_v138.fa.gz",tryRC=TRUE,verbose=TRUE)
# 135 out of 8943 were assigned to the species level.
# Of which 108 had genera consistent with the input table.

saveRDS(taxa.plus, file="ss16s_taxaplus2.rds")
saveRDS(taxa, file="ss16s_taxa2.rds")
write.csv(taxa.plus, file="ss16s_taxaplus2.csv")
write.csv(taxa, file="ss16s_taxa2.csv")

saveRDS(seqtab2.nochim, file="ss16s_seqtab.nochim2.rds")
write.csv(seqtab2.nochim, file="ss16s_seqtab.nochim2.csv")
write.csv(seqtab2.nochim, file="ss16s_seqtab.nochim_renamed2.csv")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

#BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library(cowplot)

#import dataframe holding sample information
setwd("~/Desktop/TVE/RefData")
samdf<-read.csv("samples.list_Laura1.csv")
head(samdf)
rownames(samdf) <- samdf$id

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab2.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8943 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 8943 taxa by 7 taxonomic ranks ]

#first look at data
quartz()
ps_glom <- tax_glom(ps, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, x="site", fill="Family")+
  theme(legend.position="top right")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab2.nochim))))
#making output fasta file for lulu step & maybe other things, before giving new ids to sequences

colnames(seqtab2.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab2.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))

ps

#### remove mitochondria, chloroplasts, non-bacteria #### 
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #170 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #231 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #253 taxa to remove

ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #8773 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #8542 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #8289 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #155 taxa

seqtab.clean <- data.frame(otu_table(ps.clean))
write.csv(seqtab.clean,file="seqtab.cleaned2.csv")
seqtab.clean <- read.csv("seqtab.cleaned2.csv",row.names=1)

#### rarefy #####
library(vegan)

quartz()
rarecurve(seqtab.clean,step=100,label=FALSE) #after removing contaminats

total <- rowSums(seqtab.clean)
subset(total, total <6000) #6000 for me (kinda arbitrary)
# I2A  I2G  I3A  I3I 
# 2889 3581 1381 5094 
#4 samples
#I2A I2G I3A I3I identified by MCMC.OTU below as being too low

row.names.remove <- c("I2A","I2G","I3A","I3I")
seqtab.less <- seqtab.clean[!(row.names(seqtab.clean) %in% row.names.remove),]
samdf.rare <- samdf[!(row.names(samdf) %in% row.names.remove), ]
dim(samdf.rare)
# 50 samples left

seqtab.rare <- rrarefy(seqtab.less,sample=6000)
rarecurve(seqtab.rare,step=100,label=FALSE)

#phyloseq object but rarefied
ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2))
ps.rare
#8289 taxa
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8289 taxa and 50 samples ]
# sample_data() Sample Data:       [ 50 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 8289 taxa by 8 taxonomic ranks ]

# #removing missing taxa - lost after rarefying
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare
# 7894 taxa
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7894 taxa and 50 samples ]
# sample_data() Sample Data:       [ 50 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 7894 taxa by 8 taxonomic ranks ]

seqtab.rare <- data.frame(otu_table(ps.rare))

#saving
write.csv(seqtab.rare, file="ss16s_seqtab.rare_6k.csv")
write.csv(samdf.rare, file="ss16s_samdf.rare_6k.csv")
#re-reading
samdf.rare <- read.csv("ss16s_samdf.rare_6k.csv",row.names=1)
seqtab.rare <- read.csv("ss16s_seqtab.rare_6k.csv",row.names=1)

#### alpha diversity ####
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).
quartz()
plot_richness(ps.clean, x="site", measures=c("Shannon", "Simpson"), color="reef") + theme_bw()

df <- data.frame(estimate_richness(ps.rare, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))
df <- data.frame(estimate_richness(ps.clean, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

#diversity for core microbiome
df.core <- data.frame(estimate_richness(pseq.core, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df.core$id <- rownames(df.core)
df.core <- merge(df.core,samdf,by="id") #add sample data

write.csv(df.core,file="16s_diversity_core.csv") #saving

df.core$site <- factor(df.core$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))

#diversity for accessory ASVs
df.acc <- data.frame(estimate_richness(ps.acc, split=TRUE, measures=c("Shannon","InvSimpson","Observed")))

df.acc$id <- rownames(df.acc)
df.acc <- merge(df.acc,samdf,by="id") #add sample data

write.csv(df.acc,file="16s_diversity_acc.csv") #saving

df.acc <- read.csv("16s_diversity_acc.csv",row.names=1,header=TRUE) #reading back in

df.acc$site <- factor(df.acc$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))

library(ggpubr)

diver.core.sh <- summarySE(data=df.core,measurevar=c("Shannon"),groupvars=c("reef","site"))
diver.core.sh$site <- factor(diver.core.sh$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))
diver.core.si <- summarySE(data=df.core,measurevar=c("InvSimpson"),groupvars=c("reef","site"))
diver.core.si$site <- factor(diver.core.si$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))

quartz()
gg.sh.core <- ggplot(diver.core.sh, aes(x=site, y=Shannon,color=reef,shape=reef))+
  geom_errorbar(aes(ymin=Shannon-se,ymax=Shannon+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=reef, shape=reef),size=4,position=position_dodge(0.5))+
  #geom_boxplot(outlier.shape = NA)+
  #geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("PD","SP","CI","BN","BS","CA"))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Inshore","Offshore"))+
  theme(text=element_text(family="Times"))+
  scale_colour_manual(values=c("#11F18F","#6B34F5"),labels=c("Inshore","Offshore"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))
gg.sh.core

gg.si <- ggplot(diver.si, aes(x=site, y=InvSimpson,color=reef,shape=reef))+
  geom_errorbar(aes(ymin=InvSimpson-se,ymax=InvSimpson+se),position=position_dodge(0.5),lwd=0.4,width=0.4)+
  geom_point(aes(colour=reef, shape=reef),size=4,position=position_dodge(0.5))+
  #geom_boxplot(outlier.shape = NA)+
  #geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("PD","SP","CI","BN","BS","CA"))+
  xlab("Site")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Inshore","Offshore"))+
  scale_colour_manual(values=c("#11F18F","#6B34F5"),labels=c("Inshore","Offshore"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  theme(text=element_text(family="Times"))
gg.si

quartz()
annotate_figure(ggarrange(gg.sh,gg.si,common.legend=TRUE,legend="right"),
                top = text_grob("16S Diversity", color = "black", face = "bold", size = 14))

quartz()
gg.site.sha <- ggplot(df.div,aes(x=site,y=Shannon,color=site))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("Shannon Diversity")+
  xlab("Site")+
  scale_x_discrete(labels=c("CI","PD","SP","BN","BS","CA"))+
  theme_cowplot()

gg.site.sim <- ggplot(df.div,aes(x=site,y=InvSimpson,color=site))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("Inv. Simpson diversity")+
  xlab("Site")+
  scale_x_discrete(labels=c("CI","PD","SP","BN","BS","CA"))+
  theme_cowplot()

gg.site.obs <- ggplot(df.div,aes(x=site,y=Observed,color=site))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(alpha=0.5)+
  ylab("OTU richness")+
  xlab("Site")+
  scale_x_discrete(labels=c("CI","PD","SP","BN","BS","CA"))+
  theme_cowplot()

quartz()
ggarrange(gg.site.sha, gg.site.sim, gg.site.obs, nrow=1, labels = c("A", "B","C"),
          common.legend = TRUE, legend = "right")

##for poster
quartz()
gg.sh <- ggplot(df.div, aes(x=site, y=Shannon,color=reef,shape=reef))+
  geom_boxplot(width = c(0.3),outlier.shape = NA)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("Bastimentos N"="BN","Bastimentos S" = "BS","Cayo de Agua" = "CA","Cristobal" = "CI","Punta Donato" = "PD","STRI Point" = "SP"))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(~site, values=c(16,17,15,0,1,2))+
  scale_colour_manual(values=c("Inshore"= "#6157A2","Offshore"="#68B08A"))+
  #guides(color=guide_legend(title="Reef zone"))+
  annotate("text", label = "by site, ANOVA p=",size=1, x = 5,y=0.5)+
  #stat_compare_means(method = "anova")
  theme(legend.position="none",
        axis.text=element_text(size=5),
        axis.title = element_text(size=5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size=5)
  )
gg.sh

gg.si <- ggplot(df.div, aes(x=site, y=si.log,color=reef,shape=site))+
  geom_boxplot(width = c(0.3),outlier.shape = NA)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("Bastimentos N"="BN","Bastimentos S" = "BS","Cayo de Agua" = "CA","Cristobal" = "CI","Punta Donato" = "PD","STRI Point" = "SP"))+
  xlab("Site")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(~site, values=c(16,17,15,0,1,2))+
  scale_colour_manual(values=c("Inshore"= "#6157A2","Offshore"="#68B08A"))+
  annotate("text", label = "by site, ANOVA p<0.05",size=1, x = 2,y=7)+
  #stat_compare_means(method = "anova")+
  theme(legend.position="none",
        axis.text=element_text(size=5),
        axis.title = element_text(size=5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size=5)
  )

gg.total <- ggplot(df.div, aes(x=reef, y=si.log,color=reef,shape=reef))+
  geom_boxplot(width = c(0.3),outlier.shape=NA)+
  scale_x_discrete(labels=c("Offshore"="Off","Inshore" = "In"))+
  xlab("Reef")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,1),labels=c("Inshore","Offshore"))+
  scale_colour_manual(values=c("Inshore"= "#6157A2","Offshore"="#68B08A"))+
  guides(color=guide_legend(title="Reef",
                            override.aes = list(size=1),
                            title.theme = element_text(
                              size = 5)),
         shape="none",
         title.theme = element_text(size = 5))+
  geom_jitter(alpha=0.5)+
  annotate("text", label = "by reef,
           ANOVA p<0.05",size=1, x = 1,y=7)+
  theme(legend.position="none",
        axis.text=element_text(size=5),
        axis.title = element_text(size=5),
        legend.text = element_text(size = 5)
  )

annotate_figure(ggarrange(gg.total,gg.si,
                          labels="AUTO",
                          font.label = list(size = 9),
                          common.legend=TRUE,
                          legend="none",
                          widths=c(0.5,1)),
                top = text_grob("16S Alpha Diversity", size=8,color = "black", face = "bold"),
                theme(axis.text=element_text(size=5),
                      axis.title = element_text(size=5),
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(0.4,"line")))

ggsave(filename="16s_diversity_InvSimp.log.pdf", plot = last_plot(), width=3, height=3, units="in",useDingbats=FALSE)

#stats
library(car)

shapiro.test(df.div$Shannon) #p-value = 3.549e-06
df.div$sh.log <- log(df.div$Shannon)
shapiro.test(df.div$sh.log) #p-value = 1.986e-08
leveneTest(df.div$sh.log~reef,data=df.div) #p=0.01946 *
a.div <- aov(sh.log~reef,data=df.div)
summary(a.div) #reef p=0.0162 * ; #site p=0.061 .
TukeyHSD(a.div) #p adj = 0.0162347 between reefs ; #no site differences

shapiro.test(df.div$InvSimpson) #p-value = 0.0674
leveneTest(df.div$InvSimpson~reef,data=df.div) #p-value=0.7184
df.div$si.log <- log(df.div$InvSimpson)
shapiro.test(df.div$si.log) #p-value = 2.695e-05
leveneTest(df.div$si.log~reef,data=df.div) #p=0.03365 *
a.div <- aov(si.log~site,data=df.div)
summary(a.div) #site log transformed p=0.0498 *; #reef p=0.0221 *
TukeyHSD(a.div) #reefs p adj = 0.0438467; no site differences

#non log-transformed
a.div <- aov(df.div$InvSimpson~reef,data=df.div)
summary(a.div) #p=0.0596 .
TukeyHSD(a.div)
a.div <- aov(df.div$InvSimpson~site,data=df.div)
summary(a.div) #p=0.149
TukeyHSD(a.div)

##stats for core microbiome
shapiro.test(df.core$Shannon) #p-value = 0.0003026
df.core$sh.log <- log(df.core$Shannon)
shapiro.test(df.core$sh.log) #p-value = 1.424e-07
leveneTest(df.core$sh.log~reef,data=df.core) #p=0.05306 .
a.core <- aov(sh.log~reef,data=df.core)
summary(a.core) #p=0.0482 * #reef
TukeyHSD(a.core) #p adj = 0.0481762
summary(a.core) #p=0.277 #site
TukeyHSD(a.core) #nothing

shapiro.test(df.core$InvSimpson) #p-value = 0.577
leveneTest(df.core$InvSimpson~reef,data=df.core) #p-value=0.3409
df.core$si.log <- log(df.core$InvSimpson)
shapiro.test(df.core$si.log) #p-value = 0.01881
leveneTest(df.core$si.log~reef,data=df.core) #p=0.03573 *
a.core <- aov(si.log~reef,data=df.core)
summary(a.core) #site p=0.224; #reef p=0.0467 *
TukeyHSD(a.core) #nothing

a.core <- aov(df.core$InvSimpson~reef,data=df.core)
summary(a.core) #p=0.0821 .
TukeyHSD(a.core)
a.core <- aov(df.core$InvSimpson~site,data=df.core)
summary(a.core) #p=0.389
TukeyHSD(a.core)

##stats for accessosry microbiome
shapiro.test(df.acc$Shannon) #p-value = 5.519e-07
df.acc$sh.log <- log(df.acc$Shannon)
shapiro.test(df.acc$sh.log) #p-value = 1.446e-09

leveneTest(df.acc$sh.log~reef,data=df.acc) #p=0.02248 *
a.acc.reef <- aov(sh.log~reef,data=df.acc)
summary(a.acc.reef) #p=0.0452 * #reef
TukeyHSD(a.acc.reef) #p adj = 0.0451554
leveneTest(df.acc$sh.log~site,data=df.acc) #p=0.08894 .
a.acc.site <- aov(sh.log~site,data=df.acc)
summary(a.acc.site) #p=0.0686 . #site
TukeyHSD(a.acc.site) #nothing 

shapiro.test(df.acc$InvSimpson) #p-value = 0.1106
leveneTest(df.acc$InvSimpson~reef,data=df.acc) #p-value=0.04402 *

df.acc$si.log <- log(df.acc$InvSimpson)
shapiro.test(df.acc$si.log) #p-value = 1.894e-06
leveneTest(df.acc$si.log~reef,data=df.acc) #p=0.02327 *
a.acc.reef <- aov(si.log~reef,data=df.acc)
summary(a.acc.reef) #reef p=0.0953 .
TukeyHSD(a.acc.reef) #nothing

a.acc.site <- aov(si.log~site,data=df.acc)
summary(a.acc.site) #site p=0.0665 .
TukeyHSD(a.acc.site) #nothing

#saving diversity data frame - end of alpha diversity stuff
write.csv(df.div, "16s_diversity.csv")
df.div <- read.csv("16s_diversity.csv", header = TRUE)


#### trim underrepresented ASVs ####
library(MCMC.OTU)

#regular
seq.formcmc <- read.csv("~/Desktop/TVE/Rstudio/16S/16Sv2/seqtab.cleaned2_formcmc.csv")
seq.trim <- purgeOutliers(seq.formcmc,count.columns=3:8291,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#2 bad samples - "I2A" "I3A"
#z-scores: -2.627992 -4.036692
#1333 ASVs passing cutoff for reads
#ASVs show up in 2% of samples:
# FALSE  TRUE 
# 364  1138 
#remove bad samples from sample data frame
row.names.remove <- c("I3A","I2A")
samdf.trim <- samdf[!(row.names(samdf) %in% row.names.remove), ]
rownames(seq.trim) <- seq.trim$sample
seq.trim <- seq.trim[,3:1140] #normal

write.csv(seq.trim,file="seq.trim.csv")
seq.trim <- read.csv("seq.trim.csv",row.names=1)

#remake phyloseq objects
ps.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                    sample_data(samdf.trim), 
                    tax_table(taxa2))
ps.trim #1138 taxa and 52 samples

#rarefied
seq.formcmc <- read.csv("~/Desktop/TVE/Rstudio/16S/16Sv2/ss16s_seqtab.rare_6k.csv")
seq.rare <- purgeOutliers(seq.formcmc,count.columns=3:7895,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#no bad samples
#1565 OTUs passing cutoffs
#1565 show up in 2% of samples

#rename rows
rownames(seq.rare) <- seq.rare$sample.names

#remove sample info
seq.rare <- seq.rare[,3:1567] #rarefied

write.csv(seq.rare,file="seq.rare6k.trim.csv")
seq.rare <- read.csv("seq.rare6k.trim.csv",row.names=1)

#ps object - rarefied
ps.rare.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                         sample_data(samdf.rare), 
                         tax_table(taxa2))
ps.rare.trim #1565 taxa and 50 samples

#### pcoa plots ####
library(stats)
library(MCMC.OTU)

#transform to relative abundance
quartz()
ps.trim.rel <- transform_sample_counts(ps.trim, function(x) x / sum(x))
ord <- ordinate(ps.trim.rel, "PCoA", "bray")
plot_ordination(ps.trim.rel, ord,color="reef", shape="site")+
  stat_ellipse()

#saving
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

ord <- ordinate(ps.rare.trim, "PCoA", "bray")
gg.pcoa.site.rare <- plot_ordination(ps.rare.trim, ord, color=NULL,shape=NULL)+
  geom_point(aes(color=site,shape=site))+
  stat_ellipse(aes(color=reef))+
  theme_cowplot()+
  #scale_color_manual(name="Site",values=c("darkslategray3","darkslategray4","#000004"))+
  #scale_shape_manual(name="Site",values=c(8,4,9))+
  xlab("Axis 1")+
  ylab("Axis 2")+
  #annotate(geom="text", x=0.35, y=0.65, label="p < 0.01**",size=4)+
  ggtitle("Rarefied")

ord.rel <- ordinate(ps.trim.rel, "PCoA", "bray")
gg.pcoa.site <- plot_ordination(ps.trim.rel, ord.rel)+
  geom_point(aes(color=site,shape=site))+
  stat_ellipse(aes(color=reef))+
  theme_cowplot()+
  #scale_color_manual(name="Site",values=c("darkslategray3","darkslategray4","#000004"))+
  #scale_shape_manual(name="Site",values=c(8,4,9))+
  xlab("Axis 1")+
  ylab("Axis 2")+
  #annotate(geom="text", x=0.35, y=0.65, label="p < 0.01**",size=4)+
  ggtitle("Relative abundance")
gg.pcoa.site

quartz()
ggarrange(gg.pcoa.site.rare,gg.pcoa.site,labels="AUTO",common.legend=TRUE,legend="right")

#other color options
#scale_color_manual(name="Site",values=c("darkturquoise","darkslategray4","#000004"))+

#by reef zone
ps.inshore <- subset_samples(ps.trim.rel,reef=="Inshore")
#ps.inrare <- subset_samples(ps.rare.trim,reef=="Inshore")
ord.inshore <- ordinate(ps.inshore, "PCoA", "bray")
gg.inshore <- plot_ordination(ps.inshore, ord.inshore,color="site", shape="site")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  #scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("turquoise","violet","tomato"))+
  guides(color=guide_legend(title="Site"),shape=guide_legend(title="Site"))+
  ggtitle("Inshore - Relative")+
  #xlab("Axis 1 (11.2%)")+#rarefied
  #ylab("Axis 2 (9.8%)")+#rarefied
  #xlab("Axis 1 (11.1%)")+#non-rarefied
  #ylab("Axis 2 (9.6%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
gg.inrare

ps.offshore <- subset_samples(ps.trim.rel,reef=="Offshore")
#ps.offrare <- subset_samples(ps.rare.trim,reef=="Offshore")
ord.offshore <- ordinate(ps.offshore, "PCoA", "bray")
gg.offshore <- plot_ordination(ps.offshore, ord.offshore,color="site", shape="site")+
  geom_point(size=2)+
  stat_ellipse()+
  theme_cowplot()+
  #scale_shape_manual(values=c(16,15),labels=c("Back reef","Fore reef"))+
  scale_color_manual(values=c("tan3","springgreen3","slateblue3"))+
  guides(color=guide_legend(title="Site"),shape=guide_legend(title="Site"))+
  ggtitle("Offshore - Relative")+
  #xlab("Axis 1 (10.4%)")+#rarefied
  #ylab("Axis 2 (8.7%)")+#rarefied
  #xlab("Axis 1 (11.5%)")+#non-rarefied
  #ylab("Axis 2 (9.9%)")+#non-rarefied
  theme(axis.text=element_text(size=10))
quartz()
gg.offshore

quartz()
ggarrange(gg.inshore,gg.offshore,nrow=1,legend="right",labels="AUTO")

#### pcoa plot ####
library(stats)
library(MCMC.OTU)

#transform to relative abundance
ps.trim.rel <- transform_sample_counts(ps.trim, function(x) x / sum(x))
ord <- ordinate(ps.trim.rel, "PCoA", "bray")
plot_ordination(ps.trim.rel, ord,color="reef", shape="site")+
  stat_ellipse()
#Axis.1 [8.4%], Axis.2 [7.2%]

library(ggpubr)

### carly's method
df.seq <- as.data.frame(seq.trim)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="bray")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
scorez$id <- rownames(scorez)
pcoa.all <- merge(scorez,samdf.rare,by='id')
pcoa.all <- merge(scorez,samdf.trim,by='id')

## for printing
ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=reef,fill=reef,shape=site))+
  geom_point(size=2)+
  xlab('Axis 1 (8.4%)')+
  ylab('Axis 2 (7.2%)')+
  #stat_ellipse(aes(color=reef, group=reef))+
  stat_ellipse(geom = "polygon", alpha = 0.15, aes(color=reef,group=reef,fill = reef))+
  #facet_wrap(~site)+
  theme_cowplot()+
  scale_shape_manual(~site, values=c(0,1,2,15,16,17))+
  scale_colour_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  scale_fill_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  annotate("text", label = "by reef, Adonis p=0.01", size=1, x = -0.05, y = 0.05) +
  annotate("text", label = "by site, Adonis p<0.001", size=1, x = -0.05, y = 0.04) +
  guides(color=FALSE,
         shape=FALSE,
         # color=guide_legend(title="Reef",
         #                    override.aes = list(size=1),
         #                    title.theme = element_text(
         #                      size = 5)),
         # shape=guide_legend(title="Site",
         #                    override.aes = list(size=1),
         #                    title.theme = element_text(
         #                      size = 5),
         #                    fill="none"),
         fill=FALSE)+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size=5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4,"line")
  )

leg.pcoa <- get_legend(pcoa)
as_ggplot(leg.pcoa)

ggsave(filename="16s_pcoa_noleg.pdf", plot = last_plot(), width=3, height=3, units="in", useDingbats=FALSE)

### pcoa for core microbiome only
#transform to relative abundance
pseq.core.rel <- transform_sample_counts(pseq.core, function(x) x / sum(x))
ord <- ordinate(pseq.core.rel, "PCoA", "bray")
plot_ordination(pseq.core.rel, ord, color="reef", shape="site")
#Axis.1 [8.4%], Axis.2 [7.2%]

#### Stats ####
library(vegan)
#help on adonis here:
#https://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html#more

ps.trim.rel <- transform_sample_counts(ps.trim, function(x) x / sum(x))
seq.trim.rel <- data.frame(otu_table(ps.trim.rel))

dist.seqtab <- vegdist(seq.trim.rel)
dist.seqtab <- vegdist(seq.trim)
bet.all <- betadisper(dist.seqtab,samdf.trim$reef) #p=0.0165 * reef,rel // p=0.006001 ** not rel
bet.all <- betadisper(dist.seqtab,samdf.trim$site) #p=0.09665 . site,rel // p=0.1006 not rel
anova(bet.all)

#bet.all as not rel, $reef
plot(bet.all)
adonis(seq.trim ~ site, data=samdf.trim, permutations=999) #p=0.001
adonis(seq.trim ~ site*reef, data=samdf.trim, permutations=999) #p=0.001 ***
adonis(seq.trim ~ reef, data=samdf.trim, permutations=999) #p=0.001 ***
adonis(seq.trim ~ reef, strata=samdf.trim$site,data=samdf.trim, permutations=999) #p=1
adonis(seq.trim.rel ~ reef, strata=samdf.trim$site, data=samdf.trim, permutations=999) #p=1
adonis(seq.trim.rel ~ reef*reef, data=samdf.trim, permutations=999) #p=0.001 ***
adonis(seq.trim.rel ~ reef, data=samdf.trim, permutations=999) #p=0.001 ***

#rarefied
ps.rare.rel <- transform_sample_counts(ps.rare, function(x) x / sum(x))
seq.rare.rel <- data.frame(otu_table(ps.rare.rel))

dist.seqtab <- vegdist(seq.rare.rel)
dist.seqtab <- vegdist(seq.rare)
bet.all <- betadisper(dist.seqtab,samdf.rare$reef) #p=0.1133 reef,rel // p=0.09012 . not rel
bet.all <- betadisper(dist.seqtab,samdf.rare$site) #p=0.338 site,rel // p=0.2462 not rel
anova(bet.all)

plot(bet.all)
adonis(seq.rare ~ site, data=samdf.rare, permutations=999) #p=0.001 ***
adonis(seq.rare ~ reef,data=samdf.rare, permutations=999) #p=0.001 ***
adonis(seq.rare ~ site*reef, data=samdf.rare, permutations=999) #p=0.001 ***
adonis(seq.rare ~ reef, strata=samdf.rare$site,data=samdf.rare, permutations=999) #p=1
adonis(seq.rare.rel ~ reef, strata=samdf.rare$site, data=samdf.rare, permutations=999) #p=1
adonis(seq.rare.rel ~ reef*reef, data=samdf.rare, permutations=999) #p=0.001 ***
adonis(seq.rare.rel ~ reef, data=samdf.rare, permutations=999) #p=0.001 ***

#by reef
samdf.inshore <- subset(samdf.trim,reef=="Inshore")
sam.inshore <- rownames(samdf.inshore)
seq.inshore <- seq.trim.rel[(row.names(seq.trim.rel) %in% sam.inshore),]

dist.inshore <- vegdist(seq.inshore)
bet.inshore <- betadisper(dist.inshore,samdf.inshore$site,bias.adjust = TRUE,type="median")
anova(bet.inshore)
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     2 0.008770 0.0043849  1.0347  0.372
# Residuals 22 0.093232 0.0042378           
# Sites not significantly different within inshore (STRI, Punta, Cris)
permutest(bet.inshore, pairwise = FALSE, permutations = 99)
plot(bet.inshore)
#not sig

samdf.offshore <- subset(samdf.trim,reef=="Offshore")
sam.offshore <- rownames(samdf.offshore)
seq.offshore <- seq.trim.rel[(row.names(seq.trim.rel) %in% sam.offshore),]

dist.offshore <- vegdist(seq.offshore)
bet.offshore <- betadisper(dist.offshore,samdf.offshore$site,bias.adjust = TRUE,type="median")
anova(bet.offshore)

#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     2 0.011575 0.0057877  0.6458 0.5331
# Residuals 24 0.215086 0.0089619        
# Sites not significantly different within offshore (Drago, Cayo, Basti)
permutest(bet.offshore, pairwise = FALSE, permutations = 99)
plot(bet.offshore)
#not sig

dist.seqtab <- vegdist(seq.trim)
bet.all <- betadisper(dist.seqtab,samdf.trim$reef)
anova(bet.all)
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     1 0.013987 0.0139874  3.4838 0.06822 .
# Residuals 47 0.188702 0.0040149     
plot(bet.all)
adonis(seq.trim ~ reef, strata=samdf.trim$site, data=samdf.trim, permutations=999)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# reef       1     0.584 0.58401  1.7425 0.03575      1
# Residuals 47    15.752 0.33516         0.96425       
# Total     48    16.337                 1.00000       
adonis(seq.trim ~ reef, data=samdf.trim, permutations=999)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# reef       1     0.584 0.58401  1.7425 0.03575  0.006 **
# Residuals 47    15.752 0.33516         0.96425          
# Total     48    16.337                 1.00000      

ss16sSimper = simper(sqrt(seq.trim[, c(1:ncol(seq.trim))]), samdf.trim$reef)
summary(ss16sSimper, ordered = TRUE)

ps.trim

##pcoa and stats for core microbiome


#### rename ASVs ####
library(rlang)
library(stringr)

tax <- as.data.frame(ps.rare.trim@tax_table@.Data)

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

tax_table(ps.rare.trim) <- as.matrix(tax.clean)

#### bar plot summed by reef zone & site ####
library(dplyr)
ps.all <- transform_sample_counts(ps.rare.trim, function(OTU) OTU/sum(OTU))
pa <- psmelt(ps.all)
tb <- psmelt(ps.all)%>%
  filter(!is.na(Abundance))%>%
  group_by(site,reef,Class,OTU)%>%
  summarize_at("Abundance",mean)

quartz()
ggplot(tb,aes(x=site,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  theme_cowplot()+
  #theme(legend.position="none")+
  xlab('Site')+  
  scale_x_discrete(labels=c("Bas","Cay","Cri","Dra","Pun","STRI"))

#### core v accessory microbiome ####
#BiocManager::install("microbiome")
#remotes::install_github("r-lib/rlang")
library(microbiome)

pseq.core <- core(ps.rare.trim, detection = 0, prevalence = .7)
taxa(pseq.core)
core.taxa <- as.data.frame(tax_table(pseq.core))
core.rel.otu <- as.data.frame(ps.core.rel@otu_table)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 10 taxa and 50 samples ]
# sample_data() Sample Data:       [ 50 samples by 4 sample variables ]
# tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]

#write.csv(core.rel.otu,"/Users/lauratsang/Desktop/TVE/Rstudio/16S/16Sv2/core_rel_otu_0.5prev.csv", row.names = TRUE)

#this shows the relative abundance of the core taxa within each sample
core_abundance(
  ps.rare,
  detection = 0.1/100,
  prevalence = 70/100,
  include.lowest = FALSE
)

ps_glom <- tax_glom(pseq.core, "Genus")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "site")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))

quartz()
plot_bar(pseq.core.rel, fill="Family")+
  geom_bar(stat="identity")+
  xlab('Site')+
  theme_cowplot()+
  scale_x_discrete(labels=c("BS","CA","CI","BN","PD","SP"))

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

# calculating core abundances #
core.sqs <- tax_table(pseq.core)
core.sqs.ids <- row.names(core.sqs)
core.sqs.ids

ps.rare.trim.rel <- transform_sample_counts(ps.rare.trim, function(x) x / sum(x))
seq.rare.rel <- data.frame(otu_table(ps.rare.trim.rel))
tax.core <- tax_table(ps.rare.trim.rel)

library('tidyverse')
#colMeans tells you the average abundance of each ASV across all samples - even the most abundant ASV (sq4) is not very abundant
seq.core <- seq.rare.rel %>% select(all_of(core.sqs.ids))
colMeans(seq.core)

#checking what they were in the non-core table:
colMeans(seq.rare.rel)

ps.core.rel <- phyloseq(otu_table(seq.core, taxa_are_rows=FALSE), 
                        sample_data(samdf.rare))

tax_table(ps.core.rel) <- as.matrix(tax.clean)
ps.core.rel #10 taxa

plot_bar(ps.core.rel, x="site",fill="Genus")+
  geom_bar(stat="identity")+
  theme_cowplot()

#DESEQ to find differentially abundant cores 
#BiocManager::install("DESeq2")
library(DESeq2)

ds.all = phyloseq_to_deseq2(pseq.core, ~ reef)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")
res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all = cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.all), ], "matrix"))
sigtab.all
dim(sigtab.all) #sq 53

#stats
pseq.core.otu <- data.frame(pseq.core@otu_table)
pseq.core.sam <- data.frame(pseq.core@sam_data)

dist.all <- vegdist(pseq.core.otu)
row.names(pseq.core.otu) == row.names(pseq.core.sam)
bet.all <- betadisper(dist.all,pseq.core.sam$reef,bias.adjust = TRUE,type="median")
anova(bet.all) #p=0.01627 *
permutest(bet.all, pairwise = FALSE, permutations = 99) #p=0.01 **
plot(bet.all)
adonis(pseq.core.otu ~ reef, data=pseq.core.sam, permutations=999) #p=0.043 *

#plotting the diff abundant sq53 - doesn't really do anything cuz there's just one ASV
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#stats
library(vegan)
dist.all <- vegdist(seq.core)
row.names(seq.core) == row.names(samdf.rare)
bet.all <- betadisper(dist.all,samdf.rare$reef,bias.adjust = TRUE,type="median")
anova(bet.all) #Pr(>F) = 0.04774 *
permutest(bet.all, pairwise = FALSE, permutations = 99)
plot(bet.all) #Pr(>F) = 0.05 *
adonis(seq.core ~ reef, data=samdf.rare, permutations=999)
#Pr(>F) = 0.069 .

plot_ordination(pseq.core, ordinate(pseq.core, "PCoA"), type="biplot", color="Genus", shape="site", title="biplot")

#### accessory ####
ps.rare.trim.otu <- data.frame(ps.rare.trim@otu_table)
core.tax <- data.frame(pseq.core@tax_table)
core.ids <- c(rownames(core.tax))
ps.rare.trim.acc.otu <- ps.rare.trim.otu[,!colnames(ps.rare.trim.otu) %in% core.ids ]

#remake phyloseq object
ps.acc <- phyloseq(otu_table(ps.rare.trim.acc.otu, taxa_are_rows=FALSE), 
                   sample_data(samdf.rare), 
                   tax_table(taxa2))
ps.acc #1128 taxa accessory

#stats
ps.acc.otu <- data.frame(ps.acc@otu_table)
ps.acc.sam <- data.frame(ps.acc@sam_data)

dist.acc <- vegdist(ps.acc.otu)
row.names(ps.acc.otu) == row.names(ps.acc.sam)
bet.acc <- betadisper(dist.all,ps.acc.sam$reef,bias.adjust = TRUE,type="median")
anova(bet.acc) #p=0.03837 *
permutest(bet.acc, pairwise = FALSE, permutations = 99) #p=0.04 *
plot(bet.acc)
adonis(ps.acc.otu ~ reef, data=ps.acc.sam, permutations=999) #p=0.001 ***

#differential abundance for accessory microbiome
ds.acc = phyloseq_to_deseq2(ps.acc, ~ reef)
dds.acc <- estimateSizeFactors(ds.acc,type="poscounts")
stat.acc = DESeq(dds.acc, test="Wald", fitType="parametric")
res = results(stat.acc, cooksCutoff = FALSE)
alpha = 0.05
sigtab.acc = res[which(res$padj < alpha), ]
sigtab.acc = cbind(as(sigtab.acc, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.acc), ], "matrix"))
sigtab.acc
dim(sigtab.acc)

#plotting differential accessory ASVs
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab.acc$log2FoldChange, sigtab.acc$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab.acc$Phylum = factor(as.character(sigtab.acc$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab.acc$log2FoldChange, sigtab.acc$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.acc$Genus = factor(as.character(sigtab.acc$Genus), levels=names(x))
ggplot(sigtab.acc, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#making differentially abundant accessory ASVs into phyloseq object
acc.da.ids <- c(rownames(sigtab.acc))
acc.da.otu <- data.frame(ps.rare.trim.otu[,colnames(ps.rare.trim.otu) %in% acc.da.ids ])
acc.da.tax <- data.frame(tax[rownames(tax) %in% acc.da.ids, ])
#failed 2/15/21 - taxa and otu don't match?

#remake phyloseq object
ps.da.acc <- phyloseq(otu_table(acc.da.otu, taxa_are_rows=TRUE), 
                   sample_data(samdf.rare), 
                   tax_table(acc.da.tax))
ps.da.acc


#### Bar-plots ####
ps_glom <- tax_glom(ps.rare.trim, "Family")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "id")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
quartz()
plot_bar(ps2, fill="Family")+
  geom_bar(stat="identity")+
  theme_cowplot()
  #theme(legend.position="none",axis.text.x=element_text(size=rel(0.5)))

totalsums <- colSums(seqtab.rare)
summary(totalsums)

#### indicspecies by all sites ####
#BiocManager::install("indicspecies")
library(indicspecies)

#normal
indval = multipatt(seq.trim, samdf.trim$reef, control = how(nperm=999))
summary(indval)
#39 go to inshore, 59 go to offshore (1138 total, 98 selected)

#rarefied - done 3 times
rownames(seq.rare)
indval_rare_1 = multipatt(seq.rare, samdf.rare$reef, control = how(nperm=999))
summary(indval_rare_1) #40 to in, 51 to off (1565 total, 91 selected)

indval_rare_2 = multipatt(seq.rare, samdf.rare$reef, control = how(nperm=999))
summary(indval_rare_2) #43 to in, 46 to off (89 selected)

indval_rare_3 = multipatt(seq.rare, samdf.rare$reef, control = how(nperm=999))
summary(indval_rare_3) #42 to in, 50 to off (92 selected)

#plotting abundances
sqs1 <- data.frame(indval_rare_1[["sign"]])
sqs_sig1 <- subset(sqs1,p.value < 0.05)
sqs_sig1$sqs <- rownames(sqs_sig1)

sqs2 <- data.frame(indval_rare_2[["sign"]])
sqs_sig2 <- subset(sqs2,p.value < 0.05)
sqs_sig2$sqs <- rownames(sqs_sig2)

sqs3 <- data.frame(indval_rare_3[["sign"]])
sqs_sig3 <- subset(sqs3,p.value < 0.05)
sqs_sig3$sqs <- rownames(sqs_sig3)

sqs_12 <- merge(sqs_sig1,sqs_sig2,by="sqs")
sqs <- merge(sqs_sig3,sqs_12,by="sqs")

#benjamini hochberg correction
sqs.adj <- p.adjust(sqs$p.value,method="BH")
sqs.adj < 0.1

#saving
write.csv(sqs,"sqs.csv")

#multitest correction
sqs.all <- data.frame(indval[["sign"]])
sqs.all.nona <- sqs.all[complete.cases(sqs.all),]
sqs.all.nona$padj <- p.adjust(sqs.all.nona$p.value,method="BH")

#plotting abundances
sqs.all.corr <- subset(sqs.all.nona,padj <= 0.1)
sqs.all.corr$sqs <- rownames(sqs.all.corr)
#saving
write.csv(sqs.all.corr,"sqs_all_padj.csv")
sqs.all.corr <- read.csv("sqs_all_padj.csv",row.names=1)

goodtaxa <- sqs.all.corr$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.core <- prune_taxa(allTaxa, pseq.core)
plot_bar(indic.all,x="reef",fill="Genus")+
  facet_wrap(~Genus,scales="free")
plot_bar(indic.all,x="reef",fill="Genus")

#now by more specific results
sqs_in <- subset(sqs.all.corr,index==1)
sqs_out <- subset(sqs.all.corr,index==2)

taxa_in <- sqs_in$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% taxa_in)]
indic.in <- prune_taxa(allTaxa, ps.rare.trim)
indic.in.rel <- transform_sample_counts(indic.in, function(x) x / sum(x))
plot_bar(indic.in,x="reef",fill="Family")
  #facet_wrap(~Family,scales="free")

taxa_out <- sqs_out$sqs
allTaxa <- taxa_names(ps.rare.trim)
allTaxa <- allTaxa[(allTaxa %in% taxa_out)]
indic.out <- prune_taxa(allTaxa, ps.rare.trim)
plot_bar(indic.out,x="reef",fill="Family")
  #facet_wrap(~Family,scales="free")

quartz()
plot_bar(ex1, x="reef", fill="Class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef zone")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  facet_wrap(~Class)

#### bubble plot ####
library(dplyr)
library(tidyverse)
library(grid)
#install.packages("ggplotify")
library("ggplotify")
library(ggpubr)

melt.rare <- psmelt(indic.rare)
melt.rare.sum <- summarySE(melt.rare,measurevar="Abundance",groupvars = c("OTU","Genus","reef","Species","Family","Order"))

ggplot(melt.rare.sum, aes(x = reef, y = Genus)) + 
  geom_point(aes(size = Abundance,fill=Genus), alpha = 0.75, shape = 21)+
  theme_cowplot()+
  ylab("Genus")+
  theme(axis.text.y=element_text(size=7))+
  xlab('Reef zone')+
  guides(fill=FALSE)+
  scale_x_discrete(limits = rev(levels(melt.rare.sum$reef)))

#creating counts by genus
ps.glom.genus.temp <- tax_glom(ps.rare.trim, "Genus")
ps.glom.genus.temp
seq.glom.genus <- data.frame(ps.glom.genus.temp@otu_table)
tax.glom.genus <- as.data.frame(ps.glom.genus.temp@tax_table@.Data)
sqs <- rownames(tax.glom.genus)

rownames(tax.glom.genus) == sqs #check if names match up!
colnames(seq.glom.genus) <- tax.glom.genus$Genus
rownames(tax.glom.genus) <- make.names(tax.glom.genus$Genus,unique = TRUE)
colnames(seq.glom.genus) <- rownames(tax.glom.genus)
rownames(tax.glom.genus) == colnames(seq.glom.genus)
tax.glom.genus <- as.matrix(tax.glom.genus)

#renamed ps object
ps.glom.genus <- phyloseq(otu_table(seq.glom.genus, taxa_are_rows=FALSE), 
                          sample_data(samdf.rare), 
                          tax_table(tax.glom.genus))
indval.genus <- multipatt(seq.glom.genus, samdf.rare$reef, control = how(nperm=999))
summary(indval.genus)
#inshore - 12, offshore - 19, total - 31

genus.all <- data.frame(indval.genus[["sign"]])
genus.all <- subset(genus.all,p.value <= 0.05)
genus.all$genus <- rownames(genus.all)
#saving
write.csv(genus.all,"genus.csv")
genus.all <- read.csv("genus.csv",row.names=1)

goodtaxa <- genus.all$genus
allTaxa <- row.names(tax.glom.genus)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
indic.all <- prune_taxa(allTaxa, ps.glom.genus)
indic.all.rel <- transform_sample_counts(indic.all, function(x) x / sum(x))

quartz()
plot_bar(indic.all,x="reef",fill="Genus")
  #facet_wrap(~Genus,scales="free")

##
ps.glom.genus.0 <- prune_taxa(taxa_sums(ps.glom.genus)>0,ps.glom.genus)

ds.genus = phyloseq_to_deseq2(ps.glom.genus.0, ~ reef)
dds.genus <- estimateSizeFactors(ds.genus,type="poscounts")
stat.genus = DESeq(dds.genus, test="Wald", fitType="parametric")

res = results(stat.genus, cooksCutoff = FALSE)
alpha = 0.05
sigtab.genus = res[which(res$padj < alpha), ]
sigtab.genus = cbind(as(sigtab.genus, "data.frame"), as(tax_table(ps.glom.genus)[rownames(sigtab.genus), ], "matrix"))
sigtab.genus
dim(sigtab.genus)
write.csv(sigtab.genus,"sigtab.genus.csv")


#### bubble plot - genus ####
library(dplyr)

melt.genus.0 <- psmelt(ps.glom.genus.0)
melt.genus.sum.0 <- summarySE(melt.genus.0,measurevar="Abundance",groupvars = c("OTU","Genus","reef","Family","Order"))

quartz()
ggplot(melt.genus.sum.0, aes(x = reef, y = Genus)) + 
  geom_point(aes(size = Abundance,fill=Genus), alpha = 0.75, shape = 21)+
  theme_cowplot()+
  ylab("Genus")+
  theme(axis.text.y=element_text(size=7))+
  xlab('Reef zone')+
  guides(fill=FALSE)+
  scale_x_discrete(limits = rev(levels(melt.genus.sum.0$reef)))

#### DESEQ to find differentially abundant ASVs from ALL ASVs (rarefied set) ####
library(DESeq2)

ds.all = phyloseq_to_deseq2(ps.rare.trim, ~ reef)
dds.all <- estimateSizeFactors(ds.all,type="poscounts")
stat.all = DESeq(dds.all, test="Wald", fitType="parametric")
res = results(stat.all, cooksCutoff = FALSE)
alpha = 0.05
sigtab.all = res[which(res$padj < alpha), ]
sigtab.all = cbind(as(sigtab.all, "data.frame"), as(tax_table(ps.rare.trim)[rownames(sigtab.all), ], "matrix"))
sigtab.all
dim(sigtab.all)
write.csv(sigtab.all,"sigtab.all_DESEQtotalASV.csv")

goodtaxa <- c(row.names(sigtab.all))
allTaxa = taxa_names(ps.rare)
allTaxa <- allTaxa[(allTaxa %in% goodtaxa)]
ex1 = prune_taxa(allTaxa, ps.rare)

quartz()
plot_bar(ex1, x="reef", fill="Class")+
  #  scale_x_discrete(labels=c("Backreef","Forereef"))+
  xlab("Reef")+
  #  scale_fill_manual(values="red",name="Class")+
  theme(axis.text.x = element_text(angle=-45),legend.position="none")+
  facet_wrap(~Genus)

######trying selbal
#devtools::install_github(repo = "UVic-omics/selbal")
library("selbal")











