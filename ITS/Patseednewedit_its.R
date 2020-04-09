#################################### newseedtest_patdata #####################################
# Date: April 6th 2020
# By : AF. Bintarti
# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("sjmisc") 
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("userfriendlyscience")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
install.packages("cowplot")
install.packages("dplyr")
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(cowplot)
library(ggplot2)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(userfriendlyscience)
library(ggpmisc)
library(tibble)
library(fitdistrplus)
# SET THE WORKING DIRECTORY
rm(list=ls())
set.seed(3)
setwd('/Users/arifinabintarti/Documents/PatSeedData/newseedtest/newseedtest_its')
wd <- print(getwd())
otu.its.pat <- read.table('04022020_pat.OpenRef_OTU_table_rare.txt', sep='\t', header=T, row.names = 1)
colnames(otu.its.pat)
#write.csv(taxonomy.rare.pat, file = "taxonomy.rare.pat.csv")
dim(otu.its.pat)
plant.its.data = read.csv("planthealthpat.its.csv", header=T)
# checking the otu table
sort(colSums(otu.its.pat, na.rm = FALSE, dims = 1), decreasing = TRUE)
head(sort(rowSums(otu.its.pat, na.rm = FALSE, dims = 1), decreasing = FALSE))
rarecurve(t(otu.its.pat), step = 20, col = "blue", cex = 0.6) #produce the rarefaction curve
# load the map
map.pat.its <- plant.its.data

# FUNGAL COMPOSITION
BiocManager::install("phyloseq")
library(phyloseq)
taxonomy.its.pat = read.csv("consensus_taxonomy_its_filt.csv", header=T)
View(taxonomy.its.pat)
dim(taxonomy.its.pat)
otu.its.pat <- rownames_to_column(otu.its.pat, var = "OTU_ID")
# merge otu table with taxonomy
otu.tax.its <- merge(otu.its.pat, taxonomy.its.pat, by = "OTU_ID")
dim(otu.tax.its)
# separate both otu table and taxonomy
otu.its.filt <- otu.tax.its[,1:23]
dim(otu.its.filt)
tax.its.filt <- otu.tax.its[,c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
dim(tax.its.filt)
# make OTU_ID as the rowname instead of the first colname
otu.its.filt <- column_to_rownames(otu.its.filt, var = "OTU_ID")
dim(otu.its.filt)

tax.its.filt <- column_to_rownames(tax.its.filt, var = "OTU_ID")
tax.its.filt

rownames(otu.its.filt) <- rownames(tax.its.filt)
# make phyloseq otu table and taxonomy
otu.its.phyl = otu_table(otu.its.filt, taxa_are_rows = TRUE)
tax.its.phyl = tax_table(as.matrix(tax.its.filt))
# make phyloseq map
rownames(map.pat.its) <- map.pat.its$Plant.Name
map.pat.phyl <- sample_data(map.pat.its)
# make phyloseq object
phyl.its <- merge_phyloseq(otu.its.phyl,tax.its.phyl,map.pat.phyl)
phyl.its

# merge taxa by genus
# 1. genus - Fungi
fg.ra <- transform_sample_counts(phyl.its, function(x) x/sum(x))
fg.genus <- tax_glom(fg.ra, taxrank = "Genus", NArm = F)
fg.genus
fg.cum.abun=taxa_sums(fg.genus)
fg.cum.abun # abundance per OTU
df.fg.genus.taxasum=as.data.frame(fg.cum.abun)
df.fg.genus.taxasum
df.fg.genus.taxasum$fg.relabun=df.fg.genus.taxasum$fg.cum.abun/22 # mean relative abundance per OTU
sort.df.fg.genus.taxasum=df.fg.genus.taxasum[order(df.fg.genus.taxasum$fg.relabun, decreasing = T),]
view(sort.df.fg.genus.taxasum)
#write.table(sort.df.bac.phylum.taxasum, file = 'relabun_phylum_bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.fg.genus <- psmelt(fg.genus)
summary(df.fg.genus)
df.fg.genus$Genus <- as.character(df.fg.genus$Genus)
df.fg.genus$Genus[df.fg.genus$Abundance < 0.01] <- "Other"
#df.bac.phylum$Phylum[is.na(df.bac.phylum$Phylum)] <- "Other"


#write.table(df.fg.genus, file = 'df.fg.genus.txt', sep = '\t', row.names= F, quote = FALSE)
df.fg.genus.edit <- read.csv("df.fg.genus.csv", header=T)
df.fg.genus.edit$Treatment = factor(df.fg.genus.edit$Treatment, levels=c('Control','Water withholding','Nutrient addition'))
fg.genus <- ggplot(data=df.fg.genus.edit, aes(x=Plant.Name, y=Abundance, fill=Genus))
barplot.fg.genus <- fg.genus + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(title="B. Fungi",y= "Mean Relative Abundance")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           #axis.line.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           #axis.text.x = element_text(hjust = 1),
                           axis.title.y =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA, 
          size = 0.2))+
  labs(fill="Genus (or lowest taxonomic level)")+
  facet_grid(~Treatment, switch = "x", scales = "free_x")
barplot.fg.genus
# other: either the taxa cannot be classified after Domain or unaasigned taxa.
ggsave("040620_ITS_barplot_patdata.eps",
      barplot.fg.genus, device = "eps",
       width = 11, height =5, 
       units= "in", dpi = 600)


# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR Fungi
# dissimilarity indices for community ecologist to make a distance structure (Jaccard distance between samples)
otu.its.pa <- 1*(otu.its.filt>0)
otu_its.dist <- vegdist(t(otu.its.pa), method='jaccard',binary = T)
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_its_pcoa <- cmdscale(otu_its.dist, eig=T)
#env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores.its=otu_its_pcoa$points[,1]
ax2.scores.its=otu_its_pcoa$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1.its <- otu_its_pcoa$eig[1]/sum(otu_its_pcoa$eig)
ax2.its <- otu_its_pcoa$eig[2]/sum(otu_its_pcoa$eig)
map.its=cbind(map.pat.its,ax1.scores.its,ax2.scores.its)
library(ggrepel)
mult <-.25
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# PCoA Plot by site
its.pcoa <- ggplot(data = map.its, aes(x=ax1.scores.its, y=ax2.scores.its))+
  theme_bw()+
  geom_point(data = map.its, aes(x = ax1.scores.its, y = ax2.scores.its, color=Treatment),size=5,shape=20, alpha=0.5)+
  scale_colour_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#e6194b', '#3cb44b', '#ffe119'))+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1.its,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2.its,3)*100,"% var. explained", sep=""))+
  coord_fixed() + 
  labs(title = "B. Fungi")+
  theme(legend.position="bottom",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=12),
       legend.title = element_blank(),
       legend.spacing.x = unit(0.05, 'cm'))
its.pcoa 
set.seed(3)
adonis(otu_its.dist~map.its$Treatment)
ggsave("its.pcoa.tiff",
       its.pcoa, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)

# 2. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR Fungi --Bray Curtis
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist.its.bc <- vegdist(t(otu.its.filt), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa.its.bc <- cmdscale(otu_dist.its.bc, eig=T)
#env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores.its.bc=otu_pcoa.its.bc$points[,1]
ax2.scores.its.bc=otu_pcoa.its.bc$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1.its.bc <- otu_pcoa.its.bc$eig[1]/sum(otu_pcoa.its.bc$eig)
ax2.its.bc <- otu_pcoa.its.bc$eig[2]/sum(otu_pcoa.its.bc$eig)
mapPat.its.bc=cbind(map.pat.its,ax1.scores.its.bc,ax2.scores.its.bc)
# PCoA Plot by site
pat.pcoa.its.bc <- ggplot(data = mapPat.its.bc, aes(x=ax1.scores.its.bc, y=ax2.scores.its.bc))+
  theme_bw()+
  geom_point(data = mapPat.its.bc, aes(x = ax1.scores.its.bc, y = ax2.scores.its.bc, color=Treatment),size=5,shape=20, alpha=0.5)+
  scale_colour_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#e6194b', '#3cb44b', '#ffe119'))+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  coord_fixed() + 
  labs(title = "A. Bacteria/archaea")+
  theme(legend.position="bottom",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=12),
       legend.title = element_blank(),
       legend.spacing.x = unit(0.05, 'cm'))
pat.pcoa.its.bc
set.seed(3)
adonis(otu_dist.its.bc~mapPat.its.bc$Treatment)











