---
  title: "ITS2"
output: html_document
---
  
#https://github.com/Nicfall/moorea_holobiont/blob/master/mr_ITS2/mr_ITS2.R

  
#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.14")

library(dada2)
#packageVersion("dada2")
#[1] ‘1.14.0’
library(ShortRead)
#packageVersion("ShortRead")
#[1] ‘1.44.1’
library(Biostrings)
#packageVersion("Biostrings")
#[1] ‘2.54.0’

path <- "/Users/lauratsang/Desktop/TVE/itsonly/fasta"

####Sort#####

fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)
sample.names

##### Check for primers ####
FWD <- "GTGAATTGCAGAACTCCGTG"  ## CHANGE ME to your forward primer sequence
REV <- "CCTCCGCTTACTTATATGCTT"  ## CHANGE ME...

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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

###### Visualizing raw data

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,2,3,4)])
plotQualityProfile(fnFs[c(51,52,53,54)])

#cutoff at 210bp

#Then look at quality profile of R2 reads

plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(51,52,53,54)])

#cutoff at 180bp

#Starts to drop off around 210 bp for forwards and 180 for reverse

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

#changing a bit from default settings - maxEE=1 (1 max expected error, more conservative), truncating length at 210 bp for forward & 180 bp reverse

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(210,180), #leaves overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=c(20,21), #N nucleotides to remove from the start of each read
                     minLen = 50,
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce mtching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
#              reads.in reads.out
# I2A_R1.fastq    11060      9852
# I2B_R1.fastq    17016     16010
# I2C_R1.fastq     8043      6431
# I2D_R1.fastq    14758     13853
# I2E_R1.fastq    15602     14639
# I2F_R1.fastq     9674      9250

tail(out)
# reads.in reads.out
# O4D_R1.fastq    10982     10226
# O4E_R1.fastq    11459     10012
# O4F_R1.fastq    12836     12282
# O4G_R1.fastq    18604     17240
# O4H_R1.fastq     9822      8911
# O4I_R1.fastq     6133      5858
length(out)
#~############################~#
##### Learn Error Rates ########
#~############################~#

#setDadaOpt(MAX_CONSIST=30) #if necessary, increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
#102084150 total bases in 537285 reads from 43 samples will be used for learning the error rates.

errR <- learnErrors(filtRs, multithread=TRUE)
#101425941 total bases in 637899 reads from 51 samples will be used for learning the error rates.

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
dadaRs[[1]]

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
rowSums(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#mostly at 300 bp, which is just what was expected [190-330 bp]

plot(table(nchar(getSequences(seqtab))))

#removing some lengths with low frequency
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(297,310)] 

table(nchar(getSequences(seqtab2)))
dim(seqtab2)
plot(table(nchar(getSequences(seqtab2))))

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab2.nochim)/sum(seqtab2)
#0.9487432

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2.nochim), rowSums(seqtab2.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="its2_reads.csv",row.names=TRUE,quote=FALSE)


#~############################~#
##### Assign Taxonomy ##########
#~############################~#
taxa <- assignTaxonomy(seqtab2.nochim, "~/Desktop/TVE/RefData/GeoSymbio_ITS2_LocalDatabase_verForPhyloseq.fasta",tryRC=TRUE,minBoot=80,verbose=TRUE)
unname(head(taxa))

#to come back to later
saveRDS(seqtab2.nochim, file="~/Desktop/TVE/itsonly/its2_seqtab2.nochim.rds")
saveRDS(taxa, file="~/Desktop/TVE/itsonly/its2_taxa.rds")

write.csv(seqtab2.nochim, file="its2_seqtab2.nochim.csv")
write.csv(taxa, file="~/Desktop/TVE/itsonly/its2_taxa.csv")

##### Reading in prior data files ####
setwd("~/Desktop/TVE/itsonly")
seqtab2.nochim <- readRDS("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/its2_seqtab2.nochim.rds")
taxa <- readRDS("its2_taxa.rds")

##### Phyloseq ####

#BiocManager::install("") # <- missing dependency

library('phyloseq')
library('ggplot2')
library('Rmisc')

#import dataframe holding sample information
setwd("~/Desktop/TVE/RefData")
samdf<-read.csv("samples.list_Laura1.csv")
head(samdf)
tail(samdf)

rownames(samdf) <- samdf$id

setwd("~/Desktop/TVE/itsonly")

#phyloseq object with shorter names
ids <- paste0("sq", seq(1, length(colnames(seqtab2.nochim))))
ids

taxa2 <- cbind(taxa, rownames(taxa)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

colnames(seqtab2.nochim) <- ids
seqtab2.nochim

#removed sequences that were host
seqtab2.nc.nohost<-read.csv("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/its2_seqtab2.nc.nohost.csv")
as.data.frame(seqtab2.nc.nohost)
rownames(seqtab2.nc.nohost)<-seqtab2.nc.nohost$X
seqtab2.nc.nohost$X<-NULL
seqtab2.nc.nohost<-data.matrix(seqtab2.nc.nohost, rownames.force = NA)

taxa_nohost<-as.data.frame(read.csv("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/taxtable_nohost.csv"))
rownames(taxa_nohost)<-taxa_nohost$Seq_ID
taxa_nohost$Seq_ID<-NULL
as.matrix(taxa_nohost)

rownames(seqtab2.nc.nohost) == samdf$id #check if names match up!
rownames(taxa_nohost) == colnames(seqtab2.nc.nohost)

ps_nohost <- phyloseq(otu_table(seqtab2.nc.nohost, taxa_are_rows = FALSE), 
               sample_data(samdf),
               tax_table(taxa_nohost))
ps_nohost


#### Alpha diversity #####
library(ggplot2)
#install.packages("extrafont")
library(extrafontdb)
library(extrafont)
library(Rmisc)
library(cowplot)
library(ggpubr)
#font_import(paths = "/Library/Fonts/")
#loadfonts()

#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)

#Shannon: Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).

#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

plot_richness(ps_nohost, x="site", measures=c("Shannon", "Simpson"))

df <- data.frame(estimate_richness(ps_nohost, split=TRUE, measures =c("Observed","Shannon","InvSimpson")))
df

df$id <- rownames(df)
df.div <- merge(df,samdf,by="id") #add sample data

write.csv(df.div,file="~/Desktop/TVE/Rstudio/ITS/ITS_dada2/its_diversity_nohost.csv") #saving
df.div <- read.csv("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/its_diversity_nohost.csv") #reading back in 
df.div$X <- NULL

diver.sh <- summarySE(data=df.div,measurevar=c("Shannon"),groupvars=c("reef","site"))
diver.sh$site <- factor(diver.sh$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))
diver.si <- summarySE(data=df.div,measurevar=c("InvSimpson"),groupvars=c("reef","site"))
diver.si$site <- factor(diver.si$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))

df.div$site <- factor(df.div$site,levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))

quartz()
gg.sh <- ggplot(df.div, aes(x=site, y=Shannon,color=reef,shape=reef))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("Bastimentos N"="BN","Bastimentos S" = "BS","Cayo de Agua" = "CA","Cristobal" = "CI","Punta Donato" = "PD","STRI Point" = "SP"))+
  xlab("Site")+
  ylab("Shannon diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,15),labels=c("Inshore","Offshore"))+
  scale_colour_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))+
  stat_compare_means(method = "anova")
gg.sh

gg.si <- ggplot(df.div, aes(x=site, y=si.log,color=reef,shape=site))+
  geom_boxplot(width = c(0.3),outlier.shape = NA)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(labels=c("Bastimentos N"="BN","Bastimentos S" = "BS","Cayo de Agua" = "CA","Cristobal" = "CI","Punta Donato" = "PD","STRI Point" = "SP"))+
  xlab("Site")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(~site, values=c(16,17,15,0,1,2))+
  scale_colour_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  #guides(color=guide_legend(title="Reef zone"))+
  annotate("text", label = "by site, ANOVA p=0.00261",size=1, x = 4,y=4.5)+
  theme(legend.position="none",
        axis.text=element_text(size=5),
        axis.title = element_text(size=5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size=5)
  )
  
gg.total <- ggplot(df.div, aes(x=reef, y=InvSimpson,color=reef,shape=reef))+
  geom_boxplot(width = c(0.3),outlier.shape=NA)+
  scale_x_discrete(labels=c("Offshore"="Off","Inshore" = "In"))+
  xlab("Reef")+
  ylab("Inv. Simpson diversity")+
  theme_cowplot()+
  scale_shape_manual(values=c(16,1),labels=c("Inshore","Offshore"))+
  scale_colour_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  guides(color=guide_legend(title="Reef",
                            override.aes = list(size=1),
                            title.theme = element_text(
                              size = 5)),
         shape="none",
         title.theme = element_text(size = 5))+
  geom_jitter(alpha=0.5)+
  annotate("text", label = "by reef,
           ANOVA p=0.0438",size=1, x = 1.5,y=4.5)+
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
                top = text_grob("ITS Alpha Diversity", size=8,color = "black", face = "bold"),
                theme(axis.text=element_text(size=5),
                      axis.title = element_text(size=5),
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(0.4,"line")))

ggsave(filename="its_diversity_InvSimp1.pdf", plot = last_plot(), width=3, height=3, units="in", useDingbats=FALSE)

quartz()
gg.si

pdf(file="its_diversity_InvSimp2.pdf",w=5,h=5)
print(plot)
dev.off()

library(car)

shapiro.test(df.div$Shannon) #p-value = 0.0429
df.div$sh.log <- log(df.div$Shannon)
shapiro.test(df.div$sh.log) #p-value = 0.01977
leveneTest(df.div$sh.log~reef,data=df.div) #p=0.2492
a.div <- aov(sh.log~reef,data=df.div)
summary(a.div) #p=0.363
TukeyHSD(a.div) #p adj = 0.3625943
a.div <- aov(sh.log~site,data=df.div)
summary(a.div) #p=0.0963 .
TukeyHSD(a.div)
a.div <- aov(df.div$Shannon~reef,data=df.div)
summary(a.div) #p=0.464
a.div <- aov(df.div$Shannon~site,data=df.div)
summary(a.div) #p=0.108

shapiro.test(df.div$InvSimpson) #p-value = 2.472e-05
leveneTest(df.div$InvSimpson~reef,data=df.div) #p-value=0.8949
df.div$si.log <- log(df.div$InvSimpson)
shapiro.test(df.div$si.log) #p-value = 0.005914
leveneTest(df.div$si.log~reef,data=df.div) #p=0.4415 reef
a.div <- aov(si.log~site,data=df.div)
summary(a.div) #site p=0.00261 **
TukeyHSD(a.div) #p adj = 0.0438467
a.div <- aov(si.log~reef,data=df.div)
summary(a.div) #reef p=0.0438 *

a.div <- aov(df.div$InvSimpson~reef,data=df.div)
summary(a.div) #p=0.0839 .
TukeyHSD(a.div)
a.div <- aov(df.div$InvSimpson~site,data=df.div)
summary(a.div) #p=0.00609 **
TukeyHSD(a.div)


#### rarefy ####
library(vegan)

rarecurve(seqtab2.nc.nohost,step=100,label=FALSE) #after clustering

total <- rowSums(seqtab2.nc.nohost)
total
subset(total, total <1994) #I3G
summary(total)

row.names.remove <- c("I3G")
seqtab2.nc.nohost.rare <- seqtab2.nc.nohost[!(row.names(seqtab2.nc.nohost) %in% row.names.remove),]
samdf.rare <- samdf[!(row.names(samdf) %in% row.names.remove), ]
#53 samples left

seq.rare <- rrarefy(seqtab2.nc.nohost.rare,sample=1994)
rarecurve(seq.rare,step=100,label=FALSE)

#save
write.csv(seqtab2.nc.nohost.rare,"~/Desktop/TVE/Rstudio/ITS/ITS_dada2/seqtab2.nc.nohost.rare.csv")
#read back in
seq.rare <- read.csv("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/seqtab2.nc.nohost.rare.csv",row.names=1,header=TRUE)

#phyloseq object
ps.rare <- phyloseq(otu_table(seq.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf), 
                    tax_table(taxa_nohost))
ps.rare


##############################
#### Bar plot - raw table ####
##############################

quartz()
top <- names(sort(taxa_sums(ps_nohost), decreasing=TRUE))[1:30]
ps.top <- transform_sample_counts(ps_nohost, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(top, ps.top)
plot_bar(ps.top, x="Sample",fill="Phylum") + facet_wrap(~reef, scales="free_x")
write.table(taxa2, "taxtable.txt", sep="\t")

quartz()
bot <- names(sort(taxa_sums(ps_nohost), decreasing=TRUE))[74:104]
ps.bot <- transform_sample_counts(ps_nohost, function(OTU) OTU/sum(OTU))
ps.bot <- prune_taxa(bot, ps.bot)
plot_bar(ps.bot, x="Sample",fill="Phylum") #+ facet_wrap(~ColonyID+Timepoint, scales="free_x")

#bar plot but without the lines between OTUs & no "NAs"
quartz()
ps_glom <- tax_glom(ps_nohost, "Phylum")
ps0 <- transform_sample_counts(ps_glom, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "id")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")+
  theme_dark()+
  facet()


####getting raw average relative abundances
quartz()
ps.rel <- transform_sample_counts(ps_nohost, function(x) x / sum(x))
plot_bar(ps.rel,fill="Phylum")+
  facet_wrap(~reef,scales="free_x")
ps.glom <- tax_glom(ps.rel, "Phylum")


#### Pcoa - raw data ####

library(vegan)
#install.packages("ggforce")
library(ggforce)
library(ggpubr)
library(cowplot)
library(stats)
library(MCMC.OTU)

#transform to relative abundance
ps.rare.rel <- transform_sample_counts(ps.rare, function(x) x / sum(x))
ord <- ordinate(ps.rare.rel, "PCoA", "bray")
plot_ordination(ps.rare.rel, ord,color="reef", shape="reef")+
  stat_ellipse()
#Axis.1 [78.5%], Axis.2 [20%]

# creating a log-transfromed normalized dataset for PCoA:
df.seq <- as.data.frame(seqtab2.nc.nohost)
all.log=logLin(data=df.seq)

# computing Manhattan distances (sum of all log-fold-changes) and performing PCoA:
all.dist=vegdist(all.log,method="bray")
all.pcoa=pcoa(all.dist)

# plotting:
scores=all.pcoa$vectors[,1:2]
scorez <- as.data.frame(scores)
summary(scorez)
scorez$id <- rownames(scorez)
pcoa.all <- merge(scorez,samdf)

quartz()
pcoa <- ggplot(pcoa.all,aes(x=Axis.1,y=Axis.2,color=reef,fill=reef,shape=site))+
  geom_point(size=2)+
  xlab('Axis 1 (78.5%)')+
  ylab('Axis 2 (20%)')+
  #stat_ellipse(aes(color=reef, group=reef))+
  stat_ellipse(geom = "polygon", alpha = 0.15, aes(color=reef,group=reef,fill = reef))+
  #facet_wrap(~site)+
  theme_cowplot()+
  scale_shape_manual(~site, values=c(0,1,2,15,16,17))+
  scale_colour_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  scale_fill_manual(values=c("Inshore"= "#6B34F5","Offshore"="#11F18F"))+
  annotate("text", label = "by reef, Adonis p=0.001***", size=1, x = -0.2, y = 0.2) +
  annotate("text", label = "by site, Adonis p=0.001***", size=1, x = -0.2, y = 0.185) +
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
pcoa

leg.pcoa <- get_legend(pcoa)
as_ggplot(leg.pcoa)

ggsave(filename="its_pcoa_noleg.pdf", plot = last_plot(), width=3, height=3, units="in", useDingbats=FALSE)

adonis(seqtab2.nc.nohost ~ reef, data=samdf, permutations=999)
#0.001 ***
adonis(seqtab2.nc.nohost ~ site, data=samdf, permutations=999)
#0.001 ***

#### MCMC.OTU to remove underrepresented ASVs ####
library(MCMC.OTU)

# #added a column with a blank name in the beginning, with numbers in the column, mcmc.otu likes this
# #also removed the X from the beginning of sample names
seqtab2.rare.mcmc <- read.csv("~/Desktop/TVE/Rstudio/ITS/ITS_dada2/seqtab2.nc.nohost.rare_mcmc.csv",header=TRUE)
row.names(seqtab2.rare.mcmc) <- seqtab2.rare.mcmc[,2]
# #& reading back in things
# 
goods <- purgeOutliers(seqtab2.rare.mcmc, count.columns=c(3:length(seqtab2.rare.mcmc[1,])))
goods <- purgeOutliers(seqtab2.rare.mcmc$sq1,count.columns=3:82,otu.cut=0.001,zero.cut=0.02) #rare
# #otu.cut = 0.1% of reads represented by ASV 
# #zero.cut = present in more than 1 sample (2% of samples)
# colnames(goods)
# #sq 1, 2, 3, 6, 7, 12, 18, 24, 32 with min 99% matching in lulu

#not rare
# lulu.mcmc <- read.csv("~/moorea_holobiont/mr_ITS2/lulu_output_mcmc.csv",header=TRUE)
# goods <- purgeOutliers(lulu.mcmc,count.columns=3:21,otu.cut=0.001,zero.cut=0.02) #not rare
# #otu.cut = 0.1% of reads represented by ASV 
# #zero.cut = present in more than 1 sample (2% of samples)
# colnames(goods)
# #sq 1, 2, 3, 6, 7, 12, 18, 24, 32 with min 99% matching in lulu
# 
# rownames(goods) <- goods$sample
# counts <- goods[,3:11]
# 
# #mcmc.otu removed 3 undersequenced samples: "513" "530" "76", "87" bad to begin with
# remove <- c("513","530","76","87")
# samdf.mcmc <- samdf[!row.names(samdf)%in%remove,]
# 
# #write.csv(samdf.mcmc,"~/Desktop/mr_samples.csv")
# write.csv(counts,file="seqtab_lulu.trimmed.csv")
# write.csv(counts,file="seqtab_lulu.rare.trimmed.csv")
# counts <- read.csv("seqtab_lulu.trimmed.csv",row.names=1,header=TRUE)
# counts <- read.csv("seqtab_lulu.rare.trimmed.csv",row.names=1,header=TRUE)
# 
# ps.mcmc <- phyloseq(otu_table(counts, taxa_are_rows=FALSE), 
#                     sample_data(samdf.mcmc), 
#                     tax_table(taxa2))
# ps.mcmc
# 
# ps.rare.mcmc <- phyloseq(otu_table(counts, taxa_are_rows=FALSE), 
#                          sample_data(samdf), 
#                          tax_table(taxa2))
# ps.rare.mcmc #9 taxa

#### Stats ####
library(vegan)
#help on adonis here:
#https://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html#more

#all
dist.seqtab <- vegdist(seqtab2.nc.nohost.rare)
anova(betadisper(dist.seqtab,samdf.rare$reef))
#Pr(>F) 0.1087
adonis(seqtab2.nc.nohost.rare ~ site, data=samdf.rare, permutations=999)
#0.001 ***
adonis(seqtab2.nc.nohost.rare ~ reef, data=samdf.rare, permutations=999)
#0.001 ***
samdf.rare <- data.frame(sample_data(ps.rare))

#stats - sites
samdf.rare.mcmc <- data.frame(ps.rare@sam_data)
dist.rare <- vegdist(seqtab2.nc.nohost.rare)
bet <- betadisper(dist.rare,samdf.rare.mcmc$site)
anova(bet) #Pr(>F) 0.08359 .
permutest(bet, pairwise = FALSE, permutations = 99) #Pr(>F) = 0.09 .
plot(bet)
adonis(seqtab2.nc.nohost.rare ~ site, data=samdf.rare.mcmc, permutations=999)
#Pr(>F) 0.001 ***

#stats - reef
bet.reef <- betadisper(dist.rare,samdf.rare.mcmc$reef)
anova(bet.reef) #Pr(>F) 0.1087
permutest(bet.reef, pairwise = FALSE, permutations = 99) #Pr(>F) = 0.13
plot(bet.reef)
adonis(seqtab2.nc.nohost.rare ~ reef, data=samdf.rare.mcmc, permutations=999)
#Pr(>F) 0.001 ***

#install.packages("remotes")
#remotes::install_github("Jtrachsel/funfuns")
library("funfuns")
pairwise.adonis(seqtab2.nc.nohost.rare, factors = samdf.rare.mcmc$site, permutations = 999)


```



