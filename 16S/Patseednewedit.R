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
install.packages("data.table")
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
library(data.table)
# SET THE WORKING DIRECTORY
rm(list=ls())
set.seed(3)
setwd('/Users/arifinabintarti/Documents/BioRxiv_Seed_Microbiome_2020/16S/')
wd <- print(getwd())
otu.rare.pat <- read.table('single_rare.txt', sep='\t', header=T, row.names = 1)
colnames(otu.rare.pat)
taxonomy.rare.pat <- otu.rare.pat[,'taxonomy']
str(taxonomy.rare.pat)
#write.csv(taxonomy.rare.pat, file = "taxonomy.rare.pat.csv")
dim(otu.rare.pat)
otu.rare.pat <- otu.rare.pat[,-25]
dim(otu.rare.pat)
plant.data = read.csv("planthealthpat.csv", header=T)
map.pat <- plant.data

#checking the colsum and rowsum and produce rarefaction curve
sort(colSums(otu.rare.pat, na.rm = FALSE, dims = 1), decreasing = TRUE)
head(sort(rowSums(otu.rare.pat, na.rm = FALSE, dims = 1), decreasing = FALSE))
rarecurve(t(otu.rare.pat), step = 20, col = "blue", cex = 0.6) #produce the rarefaction curve

################################################################################################################

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
#add Faith's PD index
pd.index <- read.table('PD_whole_tree.txt', sep='\t', header=T, row.names = 1)
pd.index <- rownames_to_column(pd.index, "Plant.Name")
#join aith's PD index to the map 
map.pat <- merge(map.pat, pd.index, by="Plant.Name", all = T)

#1. Compare bacterial and archaeal richness among different treatments
rich.pat <- aov(map.pat$Richness ~ Treatment, data = map.pat)
summary(rich.pat) #F=3.074, p-val=0.0675 
# testing assumptions
# Generate residual and predicted values
rich.pat.resids <- residuals(rich.pat)
rich.pat.preds <- predict(rich.pat)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(rich.pat.resids) # p-value = 0.048, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.006 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Richness ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 3.8317, df = 2, p-value = 0.1472 no signif differences
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
#do plot
rich<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Richness, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  labs(title = "A. Bacteria/archaea")+
  expand_limits(x = 0, y = 0)+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14),
       axis.text.x=element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
rich

#2. Compare bacterial and archaeal Faith's PD index among different treatments 
pd.pat <- aov(PD_whole_tree ~ Treatment, data = map.pat)
summary(pd.pat) #F=2.578 p-val=0.0997
# testing assumptions
# Generate residual and predicted values
pd.pat.resids <- residuals(pd.pat)
pd.pat.preds <- predict(pd.pat)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(pd.pat.resids) # p-value = 0.036, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(PD_whole_tree ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.021 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(PD_whole_tree ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 3.42, df = 2, p-value = 0.1809
# there are no signif differences
# do plot
pd<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=PD_whole_tree, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  labs(y = "Faith's Phylogenetic Diversity")+
  expand_limits(x = 0, y = 0)+
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
pd

#3. Compare bacterial and archaeal shannon among different treatments 
sha.pat <- aov(map.pat$Shannon ~ Treatment, data = map.pat)
summary(sha.pat) #F=6.781 p-val=0.00535 **
# testing assumptions
# Generate residual and predicted values
sha.pat.resids <- residuals(sha.pat)
sha.pat.preds <- predict(sha.pat)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(sha.pat.resids) # p-value = 0.029, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shannon ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.038 variances among group are not homogenous
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
sha.dt_letter = column_to_rownames(sha.dt_letter, "Group" )
rownames(sha.dt_letter)
rownames(sha.dt_letter) <- c("Control","Nutrient addition","Water withholding")
sha.dt_letter = rownames_to_column(sha.dt_letter, "Treatment" )
sha_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.sha=max(Shannon))
sha_plant.summ=left_join(sha.dt_letter,sha_plant.summarized, by='Treatment') 
sha_plant.summ
# plot
#scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+

sha<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Shannon, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"), values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  #labs(title = "B. Shannon Index")+
  geom_text(data=sha_plant.summ ,aes(x=Treatment,y=0.1+max.sha,label=Letter),vjust=0)+
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


alpha <- ggarrange(rich, pd,                                  
             ncol = 2, nrow = 1)
alpha <- annotate_figure(alpha,
                bottom = text_grob("Treatment", face = "bold", size = 20))
ggsave("alpha.tiff",
       alpha, device = "tiff",
       width = 8, height= 4, 
       units= "in", dpi = 600)
################################################################################################################

# BACTERIA COMPOSITION
BiocManager::install("phyloseq")
library(phyloseq)
taxonomy.rare.pat = read.csv("taxa.rare.edit2.csv", header=T)
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
ntaxa(phyl.pat.obj)



#######################################################################################################################
# calculating relative abundance of each OTU in each sample
otu_rel <- decostand(otu.rare.pat, method="total", MARGIN=2)
otu_rel
# calculating mean of relative abundance of each OTU in all samples # "sum(otu_rel[1,])/24"
Mean_abund <- apply(otu_rel, 1, mean)
Mean_abund # similar to "taxa_sum(bac.ra)/24"
#######################################################################################################################

# merge taxa by genus without merge samples into treatment
# 1. genus - Bacteria
bac.genus <- tax_glom(phyl.pat.obj, taxrank = "Genus", NArm = F)
bac.genus
bac.ra <- transform_sample_counts(bac.genus, function(x) x/sum(x))
bac.ra

###create dataframe from phyloseq object
df.bac.genus <- data.table(psmelt(bac.ra))
dim(df.bac.genus)
###convert genus to character vector from a factor
df.bac.genus$Genus <- as.character(df.bac.genus$Genus)
###group dataframe by Phylum, calculate median rel. abundance
mean.relabund.genus <- df.bac.genus[, mean:=mean(Abundance), by ="Genus"]
###change name of genus whose mean. rel. abundance is less than 5%
other <- mean.relabund.genus[(mean <= 0.05), Genus := "other < 5%"]
##checking the total genus abundance per sample
tot.genus.abundance <- df.bac.genus %>%
  group_by(Sample) %>%
  summarise(tot.abund=sum(Abundance)) 

## OR USE DPLYR

df.bac <- psmelt(bac.ra) %>%
  group_by(Sample, Treatment, Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)

df.bac$Genus <- as.character(df.bac$Genus)
df.bac$Genus[df.bac$Mean < 0.01] <- "Other (mean relative abundance < 1%)"
df.bac$Treatment = factor(df.bac$Treatment, levels=c('Control','Water withholding','Nutrient addition'))
# 1. barplot of bacterial/archaeal composition across samples
genus <- ggplot(data=df.bac, aes(x=Sample, y=Mean, fill=Genus))
barplot.genus <- genus + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     #scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c','#f58231', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', 'lightslateblue', '#000000', 'tomato','hotpink2'))+
                     scale_fill_manual(values=c('#911eb4','#ffd8b1','#008080','tomato','#4363d8'))+
                     theme(legend.position="bottom") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(title="A. Bacteria/archaea",y= "Mean Relative Abundance")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           axis.line.x = element_blank(),
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
                           panel.border = element_rect(colour = "black", fill = NA,size = 0.2))+
                           facet_grid(~Treatment, switch = "x", scales = "free_x")+
                           guides(fill=guide_legend(nrow=2,byrow=TRUE))
                           

ggsave("041620_16S_barplot.eps",
      barplot.genus, device = "eps",
       width = 11, height =5, 
       units= "in", dpi = 600)
# other: either the taxa cannot be classified after Domain or unaasigned taxa.
#ggsave("040220_16S_barplot_patdata2.eps",
      #barplot.genus, device = "eps",
       #width = 11, height =5, 
       #units= "in", dpi = 600)


# 2. barplot of bacterial/archaeal composition across treatments
bac.genus <- tax_glom(phyl.pat.obj, taxrank = "Genus", NArm = F)
bac.genus
bac.ra <- transform_sample_counts(bac.genus, function(x) x/sum(x))
bac.ra
bac.trt <- merge_samples(bac.ra, "Treatment")
bac.ra.trt <- transform_sample_counts(bac.trt, function(x) x / sum(x))

df.trt <- psmelt(bac.ra.trt) %>%
  group_by(Sample, Treatment, Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
df.trt$Genus <- as.character(df.trt$Genus)
df.trt$Genus[df.trt$Mean < 0.1] <- "Other"
df.trt$Sample = factor(df.trt$Sample, levels=c('Control','Water withholding','Nutrient addition'))
##checking the total genus abundance per treatment
df.bac %>%
  group_by(Sample) %>%
  summarise(tot.abund=sum(Mean))

df.trt %>%
  group_by(Sample) %>%
  summarise(tot.abund=sum(Mean))
#cum.abun=taxa_sums(bac.ra)
#cum.abun # relatif abundance per genus
#df.bac.genus.taxasum=as.data.frame(cum.abun)
#df.bac.genus.taxasum
#df.bac.genus.taxasum$relabun=df.bac.genus.taxasum$cum.abun/24 # mean relative abundance per genus
#sort.df.bac.genus.taxasum=df.bac.genus.taxasum[order(df.bac.genus.taxasum$relabun, decreasing = T),]
#view(sort.df.bac.genus.taxasum)
#write.table(sort.df.bac.genus.taxasum, file = 'relabun_genus_bac.txt', sep = '\t', col.names = TRUE, row.names = F, quote = FALSE)
#relabun_genus_bac_taxon <- read.csv("relabun_genus_bac_taxon.csv", header=T)

genus.trt <- ggplot(data=df.trt, aes(x=Sample, y=Mean, fill=Genus))
barplot.genus.treat <- genus.trt + 
                     geom_bar(aes(), stat="identity", position="fill") + 
                     scale_fill_manual(values=c( '#911eb4','#008080', '#fffac8', 'tomato'))+
                     theme(legend.position="right") + 
                     guides(fill=guide_legend(nrow=5))+
                     labs(title="A. Bacteria/archaea",y= "Mean Relative Abundance")+
                     theme(plot.title = element_text(size = rel(1.5), face="bold"),
                           #axis.line.y = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text=element_text(size=12, face = "bold"),
                           axis.line.x = element_blank(),
                           #axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           axis.text.x = element_text(hjust = 1),
                           axis.title.y =element_text(size=15,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.title = element_text(size=11, face = "bold"),
                           panel.grid = element_blank(), 
                           panel.background = element_blank(),
                           strip.text.x = element_text(size = 12, face = "bold"),
                           panel.border = element_rect(colour = "black", fill = NA, size = 0.2))
                           #guides(fill=guide_legend(nrow=2,byrow=TRUE))
#  facet_grid(~Treatment, switch = "x", scales = "free_x")
barplot.genus.treat

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
  scale_colour_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8','tomato','#008080'))+
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


################### plant maternal health pat's data #############################

### 1. shoot

shoot.aov <- aov(Shoot_mass_g ~ Treatment, data = map.pat)
summary(shoot.aov) #61.29 1.71e-09 ***
# testing assumptions
# Generate residual and predicted values
shoot.resids <- residuals(shoot.aov)
shoot.preds <- predict(shoot.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(shoot.resids) # p-value = 0.00433, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Shoot_mass_g ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.0014 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Shoot_mass_g ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 20.165, df = 2, p-value = 4.18e-05
# there are signif differences
# Do Post Hoc Dunn's Test
# Do Post Hoc Dunn's Test
shoot.dt <- dunnTest(Shoot_mass_g~Treatment, map.pat, method = "bh", kw=TRUE)
print(shoot.dt,dunn.test.results=TRUE)
shoot.dt$res
shoot.dt.df <- as.data.frame(shoot.dt$res)
shoot.dt_letter = cldList(P.adj ~ Comparison,
        data = shoot.dt$res,
        threshold = 0.05)
shoot.dt_letter = column_to_rownames(shoot.dt_letter, "Group" )
rownames(shoot.dt_letter)
rownames(shoot.dt_letter) <- c("Control","Nutrient addition","Water withholding")
shoot.dt_letter = rownames_to_column(shoot.dt_letter, "Treatment" )
shoot_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.shoot=max(Shoot_mass_g))
shoot_plant.summ=left_join(shoot.dt_letter,shoot_plant.summarized, by='Treatment') 
shoot_plant.summ
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
shoot<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Shoot_mass_g, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=shoot_plant.summ ,aes(x=Treatment,y=2+max.shoot,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="A. Shoot Mass", y="Mass (gram)")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
shoot

### 2. root

root.aov <- aov(Root_mass_g ~ Treatment, data = map.pat)
summary(root.aov) # 10.43 0.000717 ***

# testing assumptions
# Generate residual and predicted values
root.resids <- residuals(root.aov)
root.preds <- predict(root.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(root.resids) # p-value = 5.284e-05, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Root_mass_g ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.146 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Root_mass_g ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 15.365, df = 2, p-value = 0.0004608
# there are signif differences
# Do Post Hoc Dunn's Test
root.dt <- dunnTest(Root_mass_g~Treatment, map.pat, method = "bh", kw=TRUE)
print(root.dt,dunn.test.results=TRUE)
root.dt$res
root.dt.df <- as.data.frame(root.dt$res)
root.dt_letter = cldList(P.adj ~ Comparison,
        data = root.dt$res,
        threshold = 0.05)
root.dt_letter = column_to_rownames(root.dt_letter, "Group" )
rownames(root.dt_letter)
rownames(root.dt_letter) <- c("Control","Nutrient addition","Water withholding")
root.dt_letter = rownames_to_column(root.dt_letter, "Treatment" )
root_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.root=max(Root_mass_g))
root_plant.summ=left_join(root.dt_letter,root_plant.summarized, by='Treatment') 
root_plant.summ
# plot
root<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Root_mass_g, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=root_plant.summ ,aes(x=Treatment,y=2+max.root,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="B. Root Mass", y="Mass (gram)")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
root

## 3. pod mass

podmass.aov <- aov(Pod_mass_g ~ Treatment, data = map.pat)
summary(podmass.aov) # 24.31 3.42e-06 ***

# testing assumptions
# Generate residual and predicted values
podmass.resids <- residuals(podmass.aov)
podmass.preds <- predict(podmass.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(podmass.resids) # p-value = 0.01098, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Pod_mass_g ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.01924 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Pod_mass_g ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 18.305, df = 2, p-value = 0.000106
# there are signif differences
# Do Post Hoc Dunn's Test
podmass.dt <- dunnTest(Pod_mass_g~Treatment, map.pat, method = "bh", kw=TRUE)
print(podmass.dt,dunn.test.results=TRUE)
podmass.dt$res
podmass.dt.df <- as.data.frame(podmass.dt$res)
podmass.dt_letter = cldList(P.adj ~ Comparison,
        data = podmass.dt$res,
        threshold = 0.05)
podmass.dt_letter = column_to_rownames(podmass.dt_letter, "Group" )
rownames(podmass.dt_letter)
rownames(podmass.dt_letter) <- c("Control","Nutrient addition","Water withholding")
podmass.dt_letter = rownames_to_column(podmass.dt_letter, "Treatment" )
podmass_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.podmass=max(Pod_mass_g))
podmass_plant.summ=left_join(podmass.dt_letter,podmass_plant.summarized, by='Treatment') 
podmass_plant.summ
# plot
podmass<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Pod_mass_g, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=podmass_plant.summ ,aes(x=Treatment,y=1+max.podmass,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="C. Pod Mass", y="Mass (gram)")+
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
podmass

## 4. pod number

podnum.aov <- aov(N_pods ~ Treatment, data = map.pat)
summary(podnum.aov) #18.89 2.02e-05 ***
# testing assumptions
# Generate residual and predicted values
podnum.resids <- residuals(podnum.aov)
podnum.preds <- predict(podnum.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(podnum.resids) # p-value = 0.0001188, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(N_pods ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.05903 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(N_pods ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 17.973, df = 2, p-value = 0.0001251
# there are signif differences
# Do Post Hoc Dunn's Test
podnum.dt <- dunnTest(N_pods~Treatment, map.pat, method = "bh", kw=TRUE)
print(podnum.dt,dunn.test.results=TRUE)
podnum.dt$res
podnum.dt.df <- as.data.frame(podnum.dt$res)
podnum.dt_letter = cldList(P.adj ~ Comparison,
        data = podnum.dt$res,
        threshold = 0.05)
podnum.dt_letter = column_to_rownames(podnum.dt_letter, "Group" )
rownames(podnum.dt_letter)
rownames(podnum.dt_letter) <- c("Control","Nutrient addition","Water withholding")
podnum.dt_letter = rownames_to_column(podnum.dt_letter, "Treatment" )
podnum_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.podnum=max(N_pods))
podnum_plant.summ=left_join(podnum.dt_letter,podnum_plant.summarized, by='Treatment') 
podnum_plant.summ
# plot
podnum<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=N_pods, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=podnum_plant.summ ,aes(x=Treatment,y=2+max.podnum,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="D. Pod Number",y="Count")+
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
podnum




## Betadisper 
groups <- factor(c(rep("Control",8),rep("Water withholding",8), rep("Nutrient addition",8)))
otu_dist <- vegdist(t(otu.pat.pa), method='jaccard',binary = T)
mod <- betadisper(otu_dist, groups)
boxplot(mod)
# Null hypothesis of no difference in dispersion between groups
anova(mod) # there is significant differences in dispersion between groups
# the variances among groups are not homogenous,
hsd=TukeyHSD(mod) #which groups differ in relation to their variances
plot(hsd)

## Faith's Phylogenetic Diversity
install.packages("picante")
library(picante)
# load phylogenetic tree
library(ape)
tree <- read.tree("rep_set.tre")
tree$tip.label
# load phyloseq object
phyl.pat.obj
phy_tree(phyl.pat.obj)

## Rhizosphere soil chemistry statistical test

### 1. ph

ph.aov <- aov(pH ~ Treatment, data = map.pat)
summary(ph.aov) #15.6 7.03e-05 ***
# testing assumptions
# Generate residual and predicted values
ph.resids <- residuals(ph.aov)
ph.preds <- predict(ph.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(ph.resids) # p-value = 0.04001, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(pH ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.56 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(pH ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 14.634, df = 2, p-value = 0.0006642
# there are signif differences
# Do Post Hoc Dunn's Test
# Do Post Hoc Dunn's Test
ph.dt <- dunnTest(pH~Treatment, map.pat, method = "bh", kw=TRUE)
print(ph.dt,dunn.test.results=TRUE)
ph.dt$res
ph.dt.df <- as.data.frame(ph.dt$res)
ph.dt_letter = cldList(P.adj ~ Comparison,
        data = ph.dt$res,
        threshold = 0.05)
ph.dt_letter = column_to_rownames(ph.dt_letter, "Group" )
rownames(ph.dt_letter)
rownames(ph.dt_letter) <- c("Control","Nutrient addition","Water withholding")
ph.dt_letter = rownames_to_column(ph.dt_letter, "Treatment" )
ph_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.ph=max(pH))
ph_plant.summ=left_join(ph.dt_letter,ph_plant.summarized, by='Treatment') 
ph_plant.summ
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
ph<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=pH, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=ph_plant.summ ,aes(x=Treatment,y=0.3+max.ph,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="A. pH", y="pH")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
ph

### 2. P

p.aov <- aov(P ~ Treatment, data = map.pat)
summary(p.aov) #22.97 5.17e-06 ***
# testing assumptions
# Generate residual and predicted values
p.resids <- residuals(p.aov)
p.preds <- predict(p.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(p.resids) # p-value = 0.02, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(P ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.001 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(P ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 15.613, df = 2, p-value = 0.0004071
# there are signif differences
# Do Post Hoc Dunn's Test
# Do Post Hoc Dunn's Test
p.dt <- dunnTest(P~Treatment, map.pat, method = "bh", kw=TRUE)
print(p.dt,dunn.test.results=TRUE)
p.dt$res
p.dt.df <- as.data.frame(p.dt$res)
p.dt_letter = cldList(P.adj ~ Comparison,
        data = p.dt$res,
        threshold = 0.05)
p.dt_letter = column_to_rownames(p.dt_letter, "Group" )
rownames(p.dt_letter)
rownames(p.dt_letter) <- c("Control","Nutrient addition","Water withholding")
p.dt_letter = rownames_to_column(p.dt_letter, "Treatment" )
p_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.p=max(P))
p_plant.summ=left_join(p.dt_letter,p_plant.summarized, by='Treatment') 
p_plant.summ
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
p<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=P, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=p_plant.summ ,aes(x=Treatment,y=5+max.p,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="B. Phosphorus (P)", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
p

### 3. K

k.aov <- aov(K ~ Treatment, data = map.pat)
summary(k.aov) #28.74 9.74e-07 ***
# testing assumptions
# Generate residual and predicted values
k.resids <- residuals(p.aov)
k.preds <- predict(p.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(k.resids) # p-value = 0.02, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(K ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.001 variances among group are not homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(K ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 19.172, df = 2, p-value = 6.868e-05
# there are signif differences
# Do Post Hoc Dunn's Test
# Do Post Hoc Dunn's Test
k.dt <- dunnTest(K~Treatment, map.pat, method = "bh", kw=TRUE)
print(k.dt,dunn.test.results=TRUE)
k.dt$res
k.dt.df <- as.data.frame(k.dt$res)
k.dt_letter = cldList(P.adj ~ Comparison,
        data = k.dt$res,
        threshold = 0.05)
k.dt_letter = column_to_rownames(k.dt_letter, "Group" )
rownames(k.dt_letter)
rownames(k.dt_letter) <- c("Control","Nutrient addition","Water withholding")
k.dt_letter = rownames_to_column(k.dt_letter, "Treatment" )
k_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.k=max(K))
k_plant.summ=left_join(k.dt_letter,k_plant.summarized, by='Treatment') 
k_plant.summ
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
k<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=K, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=k_plant.summ ,aes(x=Treatment,y=8+max.k,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="C. Potassium (K)", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
k

### 4 Calcium (Ca)

ca.aov <- aov(Ca ~ Treatment, data = map.pat)
summary(ca.aov) #2.035  0.156
# testing assumptions
# Generate residual and predicted values
ca.resids <- residuals(ca.aov)
ca.preds <- predict(ca.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(ca.resids) # p-value = 0.17, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Ca ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.59 variances among group are homogenous
### RESULT: There are no differences of Ca content among treatments

# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
ca<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Ca, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  #geom_text(data=ca_plant.summ ,aes(x=Treatment,y=1+max.ca,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="D. Calcium (Ca)", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
ca

### 5. Mg
mg.aov <- aov(Mg ~ Treatment, data = map.pat)
summary(mg.aov) #1.692  0.208
# testing assumptions
# Generate residual and predicted values
mg.resids <- residuals(mg.aov)
mg.preds <- predict(mg.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(mg.resids) # p-value = 0.02, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Mg ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.40 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(Mg ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 4.8514, df = 2, p-value = 0.08842
#There are no differences of Mg content among treatments
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
mg<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=Mg, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="E. Magnesium (Mg)", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text.y=element_text(size=14), 
       axis.text.x = element_text(size=14),
       axis.title.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
mg

### 6. OM

om.aov <- aov(OM ~ Treatment, data = map.pat)
summary(om.aov) #0.177  0.839
# testing assumptions
# Generate residual and predicted values
om.resids <- residuals(om.aov)
om.preds <- predict(om.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(om.resids) # p-value = 0.002, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(OM ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.89 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(OM ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 0.53985, df = 2, p-value = 0.7634
# there are no signif differences
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
om<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=OM, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="F. Organic Matter (OM)", y="Percent (%)")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text=element_text(size=14), 
       axis.title.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
om

### 7. NO3

no3.aov <- aov(NO3 ~ Treatment, data = map.pat)
summary(no3.aov) #8.472 0.00201 **
# testing assumptions
# Generate residual and predicted values
no3.resids <- residuals(no3.aov)
no3.preds <- predict(no3.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(no3.resids) # p-value = 0.003, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(NO3 ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.11 variances among group are homogenous
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.
kruskal.test(NO3 ~ Treatment, data=map.pat) # Kruskal-Wallis chi-squared = 14.261, df = 2, p-value = 0.0008002
# there are signif differences
# Do Post Hoc Dunn's Test
no3.dt <- dunnTest(NO3~Treatment, map.pat, method = "bh", kw=TRUE)
print(no3.dt,dunn.test.results=TRUE)
no3.dt$res
no3.dt.df <- as.data.frame(no3.dt$res)
no3.dt_letter = cldList(P.adj ~ Comparison,
        data = no3.dt$res,
        threshold = 0.05)
no3.dt_letter = column_to_rownames(no3.dt_letter, "Group" )
rownames(no3.dt_letter)
rownames(no3.dt_letter) <- c("Control","Nutrient addition","Water withholding")
no3.dt_letter = rownames_to_column(no3.dt_letter, "Treatment" )
no3_plant.summarized <- map.pat %>% 
  group_by(Treatment) %>% 
  summarize(max.no3=max(NO3))
no3_plant.summ=left_join(no3.dt_letter,no3_plant.summarized, by='Treatment') 
no3_plant.summ
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
no3<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=NO3, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  geom_text(data=no3_plant.summ ,aes(x=Treatment,y=3+max.no3,label=Letter),vjust=0)+
  #geom_errorbar(aes(ymin = mean.podnum - sd.podnum, ymax = mean.podnum + sd.podnum), width=0.2)+
  labs(title="G. Nitrate (NO3(-))", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text=element_text(size=14), 
       axis.title.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
no3

### 7. NH4
nh4.aov <- aov(NH4 ~ Treatment, data = map.pat)
summary(nh4.aov) #0.343  0.713
# testing assumptions
# Generate residual and predicted values
nh4.resids <- residuals(nh4.aov)
nh4.preds <- predict(nh4.aov)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(nh4.resids) # p-value = 0.99, data errors are normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(NH4 ~ Treatment, data=map.pat, na.action=na.exclude) # p-val=0.108 variances among group are homogenous
# there are no signif differences
# create new level order and label
level_order <- c("Control","Water withholding","Nutrient addition")
label <- c("Control","Water\nwithholding","Nutrient\naddition")
# plot
nh4<-ggplot(map.pat, aes(x=factor(Treatment, level=level_order), y=NH4, fill=Treatment)) +
  geom_boxplot()+
  scale_fill_manual(labels = c("Control", "Water withholding", "Nutrient addition"),values=c('#4363d8', '#008080', 'tomato'))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0, seed=13), alpha=0.5)+
  theme_bw()+
  expand_limits(x = 0, y = 0)+
  labs(title="H. Ammonium (NH4(+))", y="ppm")+
  scale_x_discrete(labels= label)+
  theme(legend.position="none",
       plot.background = element_blank(),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = 15, face="bold"),
       axis.text=element_text(size=14), 
       axis.title.x = element_blank(),
       axis.title.y=element_text(size=14,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
nh4

### FIGURE

# Figure 1. Plant Biomass
biomass <- ggarrange(shoot, root, podmass, podnum,                                 
             ncol = 2, nrow = 2, align = "hv")
ggsave("Plant biomass.tiff",
       biomass, device = "tiff",
       width = 8.3, height= 7, 
       units= "in", dpi = 600)

# Figure 2. Make Bar plot for bacteria and fungi in the same panel
barplot <- ggarrange(barplot.genus,barplot.fg.genus, nrow=2,ncol = 1, align = "hv")
ggsave("Fig.2 Bacteria and fungi_barplot.edit.eps",
       barplot, device = "eps",
       width = 12, height = 11, 
       units= "in", dpi = 600)

# Figure 3. Make PCoA plot for bacteria and fungi in the same panel
library(gridExtra)
plot <- ggarrange(pat.pcoa, its.pcoa,common.legend=T, legend="bottom",nrow=1, ncol=2, align = "h")
plot
ggsave("Fig.3 Bacteria and fungi_PCoAplot.pat.tiff",
       plot, device = "tiff",
       width = 6.5, height = 3.5, 
       units= "in", dpi = 600)

# Figure 4. Alpha diversity
alpha <- ggarrange(rich, rich.its, sha, sha.its,                                 
             ncol = 2, nrow = 2, align = "hv")
ggsave("Alpha diversity.tiff",
       alpha, device = "tiff",
       width = 8.3, height= 7, 
       units= "in", dpi = 600)

# Figure 5. Soil Chemistry
soilchem <- ggarrange(ph, p, k, ca, mg, om, no3, nh4,                                 
             ncol = 4, nrow = 2, align = "hv")
ggsave("Fig.2 Soil Chemistry.tiff",
       soilchem, device = "tiff",
       width = 16.5, height= 7.5, 
       units= "in", dpi = 600)




