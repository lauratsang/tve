setwd("/Users/lauratsang/Desktop/TVE/Rstudio/ITS/Symportal")

library('ggplot2')
library('Rmisc')

#install.packages("MCMC.OTU")
#install.packages(c("Iso", "xtable"))

library(MCMC.OTU)
library(vegan)

dat <- read.csv(file="absolute_counts.csv", sep=",", header=TRUE)
#getting rid of outliers
goods=purgeOutliers(dat, count.columns=c(5:length(dat[1,])))
#What is the proportion of samples with data for these OTUs?
apply(goods[,5:length(goods[1,])],2,function(x){sum(x>0)/length(x)})
#What percentage of global total counts each OTU represents?
apply(goods[,5:length(goods[1,])],2,function(x){sum(x)/sum(goods[,5:length(goods[1,])])})

gss=otuStack(data=goods,count.columns=c(5:length(goods[1,])),condition.columns=c(1:4))
head(gss)

mm=mcmc.otu(
  fixed="reef",
  data=gss,
  nitt=3000,thin=50,burnin=2000 # a long MCMC chain to improve modeling of rare OTUs
)

acpass=otuByAutocorr(mm,gss)
smm0=OTUsummary(mm,gss,otus=acpass,summ.plot=FALSE)
smmA=padjustOTU(smm0)
sigs=signifOTU(smmA)
smm1=OTUsummary(mm,gss,otus=sigs)
gss$count=gss$count+1
smmA$otuWise[sigs]

goods.no0 <- goods[rowSums(goods[,c(5:39)])!=0, ]
goods.no0<-cbind(goods.no0[0:4],"sitename_type"=paste(goods.no0$site,goods.no0$reef, sep="_"),goods.no0[14:ncol(goods.no0)])  

nl=startedLog(data=goods.no0,count.columns=6:length(names(goods.no0)), logstart=1)
head(nl)
goods.dist=vegdist(nl, method="bray")
goods.pcoa=pcoa(goods.dist)

pcp=prcomp(nl, retx=TRUE, center=TRUE)
Sym_scores=goods.pcoa$vectors
# as.data.frame(Sym_scores,row.names = TRUE)
loadings=goods.pcoa$rotation
summary(goods.pcoa)

conditions=goods.no0[,1:5]
summary(pcp)
head(pcp)

write.csv(Sym_scores, file="Sym_scores.csv")
# manually added a column called samples.1 in EXCEL to merge with conditions
Sym_scores<-read.csv("Sym_scores.csv")

symportal.pcoa.all <- merge(conditions,Sym_scores)
Sym_scores$id<-NULL

quartz()
margin=0.01
plot(Sym_scores[,1], Sym_scores[,2],type="n",
     # xlim=c(min(Sym_scores[,1])-margin,max(Sym_scores[,1])+margin),
     # xlim=c(-0.04-margin, 0.25+margin),
     # ylim=c(min(Sym_scores[,2])-margin,max(Sym_scores[,2])+margin),
     xlim = c(min(-0.2),max(0.8+margin)),
     ylim = c(min(-0.2),max(0.8+margin)),
     mgp=c(2.3,1,0),
     xlab=paste("PC1 (",round(summary(pcp)$importance[2,1]*100,1),"%)",sep=""),
     ylab=paste("PC2 (",round(summary(pcp)$importance[2,2]*100,1),"%)",sep=""))
# main="PCA colored by Environment")
points(Sym_scores[conditions$reef=="inshore",1],Sym_scores[conditions$reef=="inshore",2], pch = 16, col="#11F18F", cex=1)
points(Sym_scores[conditions$reef=="offshore",1],Sym_scores[conditions$reef=="offshore",2], pch=16, col="#6B34F5", cex=1)
text(0.4,0.4,labels="**by reef, Adonis p=0.001")
legend("topright",c("Inshore","Offshore"), col=c("#11F18F", "#6B34F5"), pch=c(16, 16), cex=1)

adonis(Sym_scores~reef,data=conditions,method="euclidean")

# ggplot(symportal.pcoa.all,aes(x=Axis.1,y=Axis.2,color=reef,shape=reef))+
#   geom_point()+
#   xlab('Axis 1')+
#   ylab('Axis 2')+
#   stat_ellipse()+
#   #facet_wrap(~site)+
#   #theme_cowplot()+
#   scale_shape_manual(values=c(16,15),labels=c("Inshore","Offshore"))+
#   scale_color_manual(values=c("#11F18F","#6B34F5"),labels=c("Inshore","Offshore"))+
#   guides(color=guide_legend(title="Reef zone"),shape=guide_legend(title="Reef zone"))

#Another method to make barplots from https://ryaneckert.github.io/
if (!require("pacman")) install.packages("pacman")

pacman::p_load("dplyr", "edgeR", "ggplot2", "MCMC.OTU", "pairwiseAdonis", "rcartocolor", "RColorBrewer", "Redmonder", "reshape2", "vegan")

pacman::p_load_gh("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

if (!require("edgeR")){BiocManager::install("edgeR", update = FALSE)
  library(edgeR)}

options("scipen" = 10)

setwd("~/Desktop/TVE/Rstudio/ITS/Symportal")

its2Seq = read.delim("absolute_clean.txt", header = TRUE, check.names = FALSE)
head(its2Seq)

its2MetaData = read.csv("samples.list_Laura1.csv", header = TRUE, check.names = FALSE)
head(its2MetaData)

its2Seq = cbind(its2Seq[1], its2MetaData[,2:4], its2Seq[,c(2:length(its2Seq))])
head(its2Seq)

#purge sequences and transpose the data to work with edgeR - remove low abundance (< 0.01%) sequences
goods = purgeOutliers(its2Seq, count.columns = 5:length(its2Seq), otu.cut = 0.0001, sampleZcut = -5)
its2SeqTransposed = t(goods[,5:length(goods[1, ])])
its2SeqList = DGEList(counts = its2SeqTransposed)
head(its2SeqList$samples)

#TMM normalization in edgeR - normalize sequence counts with weighted trimmed mean of M-values (TMM; Robinson and Oshlack 2010)
its2SeqNorm =  calcNormFactors(its2SeqList, method = "TMM")
head(its2SeqNorm$samples)
its2TMM = t(cpm(its2SeqNorm, normalized.lib.sizes = TRUE))
its2SeqNorm = cbind(its2Seq[,c(2:4)], its2TMM)

#calculate relative abundances
colOrder = order(colSums(its2SeqNorm[4:length(its2SeqNorm[1,])]), decreasing = FALSE) + 3

its2SeqPerc = cbind(its2SeqNorm[,c(1:3)], its2SeqNorm[,c(colOrder)])

its2SeqPerc$sum = apply(its2SeqPerc[, c(4:length(its2SeqPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2SeqPerc = cbind(its2SeqPerc[, c(1:3)], (its2SeqPerc[, c(5:(ncol(its2SeqPerc)-1))] / its2SeqPerc$sum))
#sanity check - should all equal 1
apply(its2SeqPerc[, c(4:(ncol(its2SeqPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

head(its2SeqPerc)

gssSeq = otuStack(its2SeqPerc, count.columns = c(4:length(its2SeqPerc[1, ])),
                  condition.columns = c(1:3))
#remove summ
gssSeq<-gssSeq[-grep("summ", gssSeq$otu),]

levels(gssSeq$otu)

library(cowplot)

#plotting - barplot
colorCount = length(c(4:length(its2SeqPerc[1,])))
getPalette = colorRampPalette(redmonder.pal(8, "qPBI"), bias = 1.7)

its2SeqPlotA = ggplot(gssSeq, aes(x = id, y = count, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black",
           size = 0.25) +
  ylab("Proportion")+
  scale_fill_manual(values=rev(getPalette(colorCount)))+
  labs(fill = expression(paste(italic("ITS2"), " sequence"))) +
  guides(fill = guide_legend(ncol = 9, reverse = TRUE)) +
  theme_bw()

its2SeqPlot = its2SeqPlotA +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        legend.key.size = unit(0.4,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )
quartz()
its2SeqPlot

##now profile data
profiledat <- read.csv(file="profiles_counts.csv", sep=",", header=TRUE)
#getting rid of outliers
profilegoods=purgeOutliers(profiledat, count.columns=c(4:length(profiledat[1,])))
#What is the proportion of samples with data for these OTUs?
apply(profilegoods[,4:length(profilegoods[1,])],2,function(x){sum(x>0)/length(x)})
#What percentage of global total counts each profile represents?
apply(profilegoods[,4:length(profilegoods[1,])],2,function(x){sum(x)/sum(profilegoods[,4:length(profilegoods[1,])])})

#calculate relative abundances
colOrder = order(colSums(profilegoods[4:length(profilegoods[1,])]), decreasing = FALSE) + 3

profilegoodsPerc = cbind(profilegoods[,c(1:3)], profilegoods[,c(colOrder)])

profilegoodsPerc$sum = apply(profilegoodsPerc[, c(4:length(profilegoodsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

profilegoodsPerc = cbind(profilegoodsPerc[, c(1:3)], (profilegoodsPerc[, c(4:(ncol(profilegoodsPerc)-1))] / profilegoodsPerc$sum))
site.order = factor(profilegoodsPerc$site, levels = c("Punta Donato", "STRI Point", "Cristobal","Bastimentos N","Bastimentos S","Cayo de Agua"))
profilegoodsPerc<-arrange(profilegoodsPerc, site.order)

#sanity check - should all equal 1
apply(profilegoodsPerc[, c(4:(ncol(profilegoodsPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

gss.rel=otuStack(data=profilegoodsPerc,count.columns=c(4:length(profilegoodsPerc[1,])),condition.columns=c(1:3))
head(gss.rel)

gss.rel<-gss.rel[-grep("summ", gss.rel$otu),]

# # color palette
# colors <- c(
#   `purple`     = "#6B34F5",
#   `green`      = "#11F18F",
#   `light purple`  = "#B69EF3",
#   `light green`     = "#A6EACD",
#   `pink` = "#FF6F6F",
#   `orange`= "#EC9B46")
# 
# ccols <- function(...) {
#   cols <- c(...)
#   
#   if (is.null(cols))
#     return (colors)
#   
#   colors[cols]
# }
# 
# color_palettes <- list(
#   `main`  = ccols("light purple", "light green","purple","green"),
#   
#   `inshore`  = ccols("purple", "light purple"),
#   
#   `offhsore`   = ccols("green", "light green"),
#   
#   `mixed` = ccols("purple", "green", "light purple", "light green", "blue", "orange"),
#   
#   `grey`  = ccols("light grey", "dark grey")
# )
# 
# color_pal <- function(palette = "mixed", reverse = FALSE, ...) {
#   pal <- color_palettes[[palette]]
#   
#   if (reverse) pal <- rev(pal)
#   
#   colorRampPalette(pal, ...)
# }
# 
# scale_color <- function(palette = "mixed", discrete = TRUE, reverse = FALSE, ...) {
#   pal <- color_pal(palette = palette, reverse = reverse)
#   
#   if (discrete) {
#     discrete_scale("colour", paste0("color_", palette), palette = pal, ...)
#   } else {
#     scale_color_gradientn(colours = pal(256), ...)
#   }
# }
# 
# scale_fill <- function(palette = "mixed", discrete = TRUE, reverse = FALSE, ...) {
#   pal <- color_pal(palette = palette, reverse = reverse)
#   
#   if (discrete) {
#     discrete_scale("fill", paste0("color_", palette), palette = pal, ...)
#   } else {
#     scale_fill_gradientn(colours = pal(256), ...)
#   }
# }


## DIV barplot good
gss.rel$site_f = factor(gss.rel$site, levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))
reef.labs <- c("Inshore", "Offshore")
names(reef.labs) <- c("inshore", "offshore")

its2SeqPlotA = ggplot(gss.rel, aes(x = sample, y = count, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black",
           size = 0.25) +
  scale_fill_manual("DIVs", values = c("D1.D4.D4c.D1c.D2" = "#487DBC", "C1.C1c.C1b" = "#E666A5", "C3" = "#E69ED5", "A4" = "#E8E544", "B5" = "#CB4E2F")) +
  ylab("Proportion")+
  #scale_fill(discrete = TRUE, palette = "main")+
  facet_wrap(reef~site_f,scales="free_x", drop=TRUE, nrow=1, labeller = labeller(reef = reef.labs))+
  #labs(fill = expression(paste(italic("ITS2"), " sequence"))) +
  guides(fill = guide_legend(ncol = 9, reverse = TRUE)) +
  theme_bw()

its2SeqPlot = its2SeqPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 8),
        legend.key = element_blank(),
        legend.key.size = unit(0.4,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 7),
        strip.text.y = element_text(size = 7),
        strip.background = element_rect(fill = "white", size = 0.9)
  )
quartz()
its2SeqPlot

#install.packages('svglite')
ggsave("symportal_DIV.barplot.pdf", plot = last_plot(), width = 6, height = 3, units = c("in"), useDingbats=FALSE)


# ##re-doing profile data with Eckert methods (produces the same graph)
# its2Profs = read.csv("profiles_counts.csv", header = TRUE, check.names = FALSE)
# head(its2Profs)
# 
# its2Profs$site = factor(its2Profs$site, levels = c("Punta Donato", "STRI Point", "Cristobal", "Bastimentos N","Bastimentos S","Cayo de Agua"))
# levels(its2Profs$site)
# 
# ##normalization
# its2ProfsTransposed = t(its2Profs[,4:length(its2Profs[1, ])])
# its2ProfsList = DGEList(counts = its2ProfsTransposed)
# head(its2ProfsList$samples)
# 
# its2ProfsNorm =  calcNormFactors(its2ProfsList, method = "TMM")
# head(its2ProfsNorm$samples)
# its2TMM = t(cpm(its2ProfsNorm, normalized.lib.sizes = TRUE))
# its2ProfsNorm = cbind(its2Profs[,c(1:3)], its2TMM)
# 
# ##plotting
# colOrder2 = order(colSums(its2ProfsNorm[4:length(its2ProfsNorm[1,])]), decreasing = TRUE) + 3
# 
# its2ProfsPerc = cbind(its2ProfsNorm[,c(1:3)],its2ProfsNorm[,c(colOrder2)])
# its2ProfsPerc$sum = apply(its2ProfsPerc[, c(4:length(its2ProfsPerc[1,]))], 1, function(x) {
#   sum(x, na.rm = T)
# })
# 
# its2ProfsPerc = cbind(its2ProfsPerc[, c(1:3)], (its2ProfsPerc[, c(4:(ncol(its2ProfsPerc)-1))]
#                                                 / its2ProfsPerc$sum))
# head(its2ProfsPerc)
# 
# apply(its2ProfsPerc[, c(4:(ncol(its2ProfsPerc)))], 1, function(x) {
#   sum(x, na.rm = T)
# }) #should equal 1
# 
# gssProf = otuStack(its2ProfsPerc, count.columns = c(4:length(its2ProfsPerc[1, ])),
#                    condition.columns = c(1:3))
# 
# gssProf<-gssProf[-grep("summ", gssProf$otu),] # remove summ rows
# 
# levels(gssProf$otu)
# 
# colorCount2 = length(c(4:length(its2ProfsPerc[1,]))) +1
# 
# its2ProfsPlotA = ggplot(gssProf, aes(x = sample_name, y = count, fill = factor(otu))) +
#   geom_bar(position = "stack", stat = "identity", color = "black", size = 0.25) +
#   ylab("Proportion") +
#   scale_fill_manual(values = c(getPalette2(colorCount2)[2:12],
#                                carto_pal(n = 7, "Tropic")[c(6,7)]))+
#   labs(fill = expression(paste(italic("ITS2"), " type profile"))) +
#   guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
#   facet_wrap(gssProf$site, scales = "free_x",nrow = 1) +
#   theme_bw()
# 
# its2ProfsPlot = its2ProfsPlotA +
#   theme(axis.title.x = element_blank(),
#         #axis.text.x = element_blank(),
#         #axis.ticks.x = element_blank(),
#         axis.title.y = element_text(color = "black", size = 12),
#         axis.text.y = element_text(color = "black", size = 12),
#         legend.position = "bottom",
#         legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
#         legend.text = element_text(color = "black", size = 10),
#         legend.key = element_blank(),
#         legend.key.size = unit(0.75,"line"),
#         legend.background = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         plot.background = element_blank(),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12),
#         strip.background = element_rect(fill = "white", size = 0.9)
#   )
# quartz()
# its2ProfsPlot

#stats
set.seed(694) #setting seed allows repetition of randomized processes
#per site
its2dispS = betadisper(vegdist(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$site)

anova(its2dispS) #Pr(>F) = 0.01902 * #shows effect of site on beta diversity

set.seed(694)
#per reef
its2dispR = betadisper(vegdist(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$reef)

anova(its2dispR) #Pr(>F) = 0.01637 * #shows effect of reef on beta diversity

its2PermTestS = permutest(its2dispS, permutations = 9999, pairwise = T, model = "full")
its2PermTestS
# Pr(>F) 0.0191 *
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#                Punta Donato STRI Point  Cristobal Bastimentos N Bastimentos S Cayo de Agua
# Punta Donato               0.10840000 0.96450000    0.00060000    0.96490000       0.0164
# STRI Point      0.10977381            0.10580000    0.25040000    0.22680000       0.8150
# Cristobal       0.96794130 0.10079853               0.00050000    0.98670000       0.0122
# Bastimentos N   0.00187554 0.24816959 0.00084348                  0.02690000       0.0139
# Bastimentos S   0.96562616 0.22439945 0.98843742    0.03271102                     0.1196
# Cayo de Agua    0.02067671 0.81280970 0.01346982    0.02000048    0.11757357             

its2PermTestR = permutest(its2dispR, permutations = 9999, pairwise = T, model = "full")
its2PermTestR
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# inshore offshore
# inshore             0.0188
# offshore 0.016371         

its2PermTestS$statistic
its2PermTestR$statistic
#site
its2Adonis = adonis(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))] ~ site,
                    data = its2ProfsNorm, permutations = 9999, method = "bray")
its2Adonis #Pr(>F) = 0.0001 ***
#reef
its2Adonis = adonis(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))] ~ reef,
                    data = its2ProfsNorm, permutations = 9999, method = "bray")

its2Adonis #Pr(>F) = 0.0001 ***

#Pairwise PERMANOVA
#site
its2PWAdonis = pairwise.adonis(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))],
                               factors = its2ProfsNorm$site,
                               sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonis

#SIMPER test to see which ITS profiles contribute most difference between reef zones
its2SimperR = simper(sqrt(its2ProfsNorm[, c(4:ncol(its2ProfsNorm))]), its2ProfsNorm$reef)
summary(its2SimperR)

#### Alpha diversity #####
library(ggplot2)
library(phyloseq)
#install.packages("extrafontdb")
library(extrafontdb)
library(extrafont)
library(Rmisc)
library(cowplot)
library(ggpubr)




