#################################### newseedtest_patdata #####################################
# Date: November 22nd 2019
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
setwd('/Users/arifinabintarti/Documents/PatSeedData/newseedtest/')
wd <- print(getwd())
otu.rare.pat <- read.table('single_rare_filt.txt', sep='\t', header=T, row.names = 1)
colnames(otu.rare.pat)
taxonomy.rare.pat <- otu.rare.pat[,'taxonomy']
str(taxonomy.rare.pat)
#write.csv(taxonomy.rare.pat, file = "taxonomy.rare.pat.csv")
dim(otu.rare.pat)
otu.rare.pat <- otu.rare.pat[,-25]
dim(otu.rare.pat)
plant.data = read.csv("planthealthpat.csv", header=T)



sort(colSums(otu.rare.pat, na.rm = FALSE, dims = 1), decreasing = TRUE)
head(sort(rowSums(otu.rare.pat, na.rm = FALSE, dims = 1), decreasing = FALSE))

rarecurve(t(otu.rare.pat), step = 20, col = "blue", cex = 0.6) #produce the rarefaction curve

#calculate alpha diversity
otu.rare.pat <- as.data.frame(otu.rare.pat)
s <- specnumber(otu.rare.pat, MARGIN = 2) # richness
rich <- as.data.frame(s)
h <- diversity(t(otu.rare.pat), index = 'shannon') # Shannon index
shannon <- as.data.frame(h)
j <- h/log(s) # Pielou's evenness
pie <- as.data.frame(j)
map.pat <- plant.data
map.pat$Richness <- s
map.pat$Shannon <- h
map.pat$Pielou <- j
map.pat

#1. Compare bacterial and archaeal richness among different treatments
rich.pat <- aov(map.pat$Richness ~ Treatment, data = map.pat)
summary(rich.pat) #F=2.32, p-val=0.123
# testing assumptions
# Generate residual and predicted values
rich.pat.resids <- residuals(rich.pat)
rich.pat.preds <- predict(rich.pat)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(rich.pat.resids) # p-value = 0.047, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.0055 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Richness ~ Treatment, data=map.pat) # no signif differences

#2. Compare bacterial and archaeal shannon among different treatments 
sha.pat <- aov(map.pat$Shannon ~ Treatment, data = map.pat)
summary(sha.pat) #F=6.009, p-val=0.0086
# testing assumptions
# Generate residual and predicted values
sha.pat.resids <- residuals(sha.pat)
sha.pat.preds <- predict(sha.pat)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(sha.pat.resids) # p-value = 0.02, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.043 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Shannon ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 10.64, df = 2, p-value = 0.004893
# there are signif differences
# Do Post Hoc Dunn's Test
sha.dt <- dunnTest(Shannon~Treatment, map.pat, method = "bh", kw=TRUE)
print(sha.dt,dunn.test.results=TRUE)
sha.dt$res
sha.dt.df <- as.data.frame(sha.dt$res)
sha.dt_letter = cldList(P.adj ~ Comparison,
        data = sha.dt$res,
        threshold = 0.05)
colnames(sha.dt_letter)[colnames(sha.dt_letter)=="Group"] <- "Treatment"
sha_plant.summarized <- map.pat %>% group_by(Treatment) %>% summarize(max.sha=max(Shannon))
sha_plant.summ=left_join(sha.dt_letter,sha_plant.summarized, by='Treatment') 
sha_plant.summ
# create new label
label <- c("Control", "Water\nwithholding", "Nutrient\naddition")
# plot
sha<-ggplot(map.pat, aes(x=Treatment, y=Shannon, fill=Treatment)) +
  geom_boxplot()+
  geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
  theme_bw()+
 expand_limits(x = 0, y = 0)+
  geom_text(data=sha_plant.summ ,aes(x=Treatment,y=0.1+max.sha,label=sha_plant.summ$Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  #labs(title = "D. Pod Number", y="Pod count")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text=element_text(size=14), 
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
sha

rich<-ggplot(map.pat, aes(x=Treatment, y=Richness, fill=Treatment)) +
  geom_boxplot()+
  geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
  theme_bw()+
 expand_limits(x = 0, y = 0)+
  #geom_text(data=sha_plant.summ ,aes(x=Treatment,y=0.1+max.sha,label=sha_plant.summ$Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  #labs(title = "D. Pod Number", y="Pod count")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text=element_text(size=14), 
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
rich

plot <- ggarrange(rich, sha,                                  
             ncol = 2, nrow = 1)
plot <- annotate_figure(plot,
                bottom = text_grob("Treatment", face = "bold", size = 20))
ggsave("alpha.tiff",
       plot, device = "tiff",
       width = 8, height= 4, 
       units= "in", dpi = 600)


# BACTERIA COMPOSITION
BiocManager::install("phyloseq")
library(phyloseq)
taxonomy.rare.pat = read.csv("taxa.rare.edit.csv", header=T)
View(taxonomy.rare.pat)
dim(taxonomy.rare.pat)
rownames(taxonomy.rare.pat) <- rownames(otu.rare.pat)
# make phyloseq otu table and taxonomy

otu.pat.phyl = otu_table(otu.rare.pat, taxa_are_rows = TRUE)
tax.pat.phyl = tax_table(as.matrix(taxonomy.rare.pat))
# make phyloseq map
rownames(map.pat) <- map.pat$Plant.Name
map.pat.phyl <- sample_data(map.pat)
# make phyloseq object
phyl.pat.obj <- merge_phyloseq(otu.pat.phyl,tax.pat.phyl,map.pat.phyl)
phyl.pat.obj

# merge taxa by genus
# 1. genus - Bacteria
bac.ra <- transform_sample_counts(phyl.pat.obj, function(x) x/sum(x))
bac.genus <- tax_glom(bac.ra, taxrank = "Genus", NArm = F)
bac.genus
cum.abun=taxa_sums(bac.genus)
cum.abun # abundance per OTU
df.bac.genus.taxasum=as.data.frame(cum.abun)
df.bac.genus.taxasum
df.bac.genus.taxasum$relabun=df.bac.genus.taxasum$cum.abun/24 # mean relative abundance per OTU
sort.df.bac.genus.taxasum=df.bac.genus.taxasum[order(df.bac.genus.taxasum$relabun, decreasing = T),]
view(sort.df.bac.genus.taxasum)
#write.table(sort.df.bac.phylum.taxasum, file = 'relabun_phylum_bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.bac.genus <- psmelt(bac.genus)
summary(df.bac.genus)
df.bac.genus$Genus <- as.character(df.bac.genus$Genus)
#df.bac.phylum$Phylum[df.bac.phylum$Abundance < 0.01] <- "Other"
#df.bac.phylum$Phylum[is.na(df.bac.phylum$Phylum)] <- "Other"
df.bac.genus$Treatment = factor(df.bac.genus$Treatment, levels=c('Control','Water withholding','Nutrient addition'))
genus <- ggplot(data=df.bac.genus, aes(x=Plant.Name, y=Abundance, fill=Genus))
barplot.genus <- genus + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(title="A. Bacteria/archaea",y= "Mean Relative Abundance")+
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
  facet_grid(~Treatment, switch = "x", scales = "free_x")
barplot.genus
# other: either the taxa cannot be classified after Domain or unaasigned taxa.
ggsave("040220_16S_barplot_patdata2.eps",
      barplot.genus, device = "eps",
       width = 11, height =5, 
       units= "in", dpi = 600)

# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA --Jaccard
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu.pat.pa <- 1*(otu.rare.pat>0)
otu_dist <- vegdist(t(otu.pat.pa), method='jaccard',binary = T)
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
#env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
mapPat=cbind(map.pat,ax1.scores,ax2.scores)
library(ggrepel)
mult <-.25
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# PCoA Plot by site
pat.pcoa <- ggplot(data = mapPat, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
  geom_point(data = mapPat, aes(x = ax1.scores, y = ax2.scores, color=Treatment),size=5,shape=20, alpha=0.5)+
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
pat.pcoa 
set.seed(3)
adonis(otu_dist~mapPat$Treatment)
ggsave("pat.pcoa.tiff",
       pat.pcoa, device = "tiff",
       width = 5, height =4, 
       units= "in", dpi = 600)


# 2. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA --Bray Curtis
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist.bc <- vegdist(t(otu.rare.pat), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa.bc <- cmdscale(otu_dist.bc, eig=T)
#env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores.bc=otu_pcoa.bc$points[,1]
ax2.scores.bc=otu_pcoa.bc$points[,2] 
#env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1.bc <- otu_pcoa.bc$eig[1]/sum(otu_pcoa.bc$eig)
ax2.bc <- otu_pcoa.bc$eig[2]/sum(otu_pcoa.bc$eig)
mapPat.bc=cbind(map.pat,ax1.scores.bc,ax2.scores.bc)
# PCoA Plot by site
pat.pcoa.bc <- ggplot(data = mapPat.bc, aes(x=ax1.scores.bc, y=ax2.scores.bc))+
  theme_bw()+
  geom_point(data = mapPat.bc, aes(x = ax1.scores.bc, y = ax2.scores.bc, color=Treatment),size=5,shape=20, alpha=0.5)+
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
pat.pcoa.bc 
set.seed(3)
adonis(otu_dist.bc~mapPat.bc$Treatment)





# Figure 1. Make Bar plot for bacteria and fungi in the same panel
barplot <- ggarrange(barplot.genus,barplot.fg.genus, nrow=2,ncol = 1, align = "hv")
ggsave("Fig.1 Bacteria and fungi_barplot.eps",
       barplot, device = "eps",
       width = 12, height = 11, 
       units= "in", dpi = 600)

# Figure 2. Make PCoA plot for bacteria and fungi in the same panel
library(gridExtra)
plot <- ggarrange(pat.pcoa, its.pcoa,common.legend=T, legend="bottom",nrow=1, ncol=2, align = "h")
plot
ggsave("Fig.2 Bacteria and fungi_PCoAplot.pat.tiff",
       plot, device = "tiff",
       width = 6.5, height = 3.5, 
       units= "in", dpi = 600)











