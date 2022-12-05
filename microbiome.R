#################################
#load data
#################################

setwd("//melbngs/Bioinfo/bioinfo-proj-zakrzem/contracts/divpro/scabies")

library (plyr)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(ampvis2)

#import count table
asv	<-	read.table(file	=  "CiP_count_filtered.tsv",	header	=	TRUE, row.names=1)
head(asv)

#import taxonomy
asv.phylo	<-	otu_table(as.matrix(asv),	taxa_are_rows	=	TRUE)
tax <-	read.table(file	="CiP_taxonomy_filtered.tsv",header=T , sep="\t", fill=TRUE, row.names=1)
tax.phylo <- tax_table(as.matrix(tax))
head(tax)

#import metadata
meta <- read.csv(file="metadata.csv",row.names=1)
meta <- meta[match( colnames(asv),rownames(meta)), ]
meta.phylo	<-	sample_data(meta)

#convert into phyloseq object
physeq	<-	phyloseq::phyloseq(asv.phylo,	tax.phylo,	meta.phylo)
physeq

########################################
# data transformation and normalization
########################################
physeq.SyDN <- subset_samples(physeq, Infection %in% c( "SyDn"))
physeq.study <- subset_samples(physeq, Infection %in% c("control", "SyDn", "SyDy"))
relAb.SyDN  <- transform_sample_counts(physeq.SyDN, function(x) x / sum(x) )
relAb.study  <- transform_sample_counts(physeq.study, function(x) x / sum(x) )


#remove taxa that have a mean relative Abundance < 1e-5 , which is 0.001%
relAb.meanFilt.SyDN  <- filter_taxa(relAb.SyDN, function(x) mean(x) > 1e-5, TRUE)
relAb.meanFilt.study  <- filter_taxa(relAb.study, function(x) mean(x) > 1e-5, TRUE)

#sqrt transform
relAb.meanFilt.sqrt.SyDN   <- transform_sample_counts(relAb.meanFilt.SyDN , function(x) sqrt(x) )
relAb.meanFilt.sqrt.study   <- transform_sample_counts(relAb.meanFilt.study , function(x) sqrt(x) )

# get only taxa from the absolute counts physeq with mean abundance 0.005%
keep <- taxa_names(filter_taxa(relAb.SyDN , function(x) mean(x) > 1e-5, TRUE))
physeq.filtered.SyDN  = prune_taxa(keep, physeq.SyDN )

keep <- taxa_names(filter_taxa(relAb.study , function(x) mean(x) > 1e-5, TRUE))
physeq.filtered.study  = prune_taxa(keep, physeq.study )

#rarefy to the min sample size, so all sample have the same sample size
minsampleSize.SyDN <- min(sample_sums(physeq.filtered.SyDN ))
minsampleSize.study <- min(sample_sums(physeq.filtered.study ))
physeq.rarefy.SyDN <- rarefy_even_depth(physeq.filtered.SyDN, sample.size = minsampleSize.SyDN, rngseed = 123, replace = FALSE)
sample_sums(physeq.rarefy.SyDN)
physeq.rarefy.study <- rarefy_even_depth(physeq.filtered.study, sample.size = minsampleSize.study, rngseed = 123, replace = FALSE)
sample_sums(physeq.rarefy.study)




##########################################
# visualizations
##########################################

#piechart
phylum<- tax_glom(physeq.filtered.SyDN, "phylum")
df.phylum <- psmelt(phylum)

df.phylum.converted <- df.phylum %>% 
  group_by(phylum) %>% 
  summarize(total_counts_per_phylum = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(total_counts_percentage = total_counts_per_phylum / sum(total_counts_per_phylum) * 100) 


df.phylum.converted <- df.phylum %>% 
  group_by(phylum) %>% 
  summarize(total_counts_per_phylum = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(total_counts_percentage = total_counts_per_phylum / sum(total_counts_per_phylum) * 100)

ggplot(df.phylum.converted, aes(x="", y=total_counts_percentage, fill=phylum)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank(), legend.position="bottom", legend.text = element_text(size = 8))


#######################################
#bubble plot
########################################
genus<- tax_glom(relAb.meanFilt.SyDN, "genus")
df.genus <- psmelt(genus)
taxa_sel <- df.genus %>% 
  group_by(genus) %>% 
  summarize(mean = mean(Abundance))%>% arrange(-mean) %>% 
  pull(genus) %>% 
  head(n=15)
df.genus.top <- df.genus[ df.genus$genus %in% taxa_sel, ] 
df.genus.top$label <- factor(df.genus.top$label, levels=c("M1.1", "M2.1","M3.1","M4.1","M5.1", "M1.7", "M2.7","M3.7", "M4.7", "M5.7", "M1.21", "M2.21",   "M3.21",  "M4.21",  "M5.21"))

ggplot(df.genus.top, aes(x = label, y = genus)) + 
  geom_point(aes(size = Abundance, fill = genus), alpha = 0.9, shape = 22, show.legend = F) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,24), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom") +  
  #scale_fill_manual(values = black, guide = FALSE) + 
  scale_x_discrete(limits = rev(levels(df.genus$label)))+
  scale_y_discrete(limits = rev(levels(df.genus$genus))) 

#########################################################
#barplot
#######################################################
ggplot(df.genus.top, aes(x = label, y = Abundance, fill = genus)) +
  geom_bar(position = "stack", stat = "identity", width = 0.9, color = "white") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title= 'Bacterial community composition at genus level', x= '', y= 'Proportion of the total counts')



#####################################################
# Heatmap
############################################################
ampvisMeta <- cbind(samples=rownames(sample_data(physeq.filtered.study)), sample_data(physeq.filtered.study)) 
d <- amp_load(
  otutable = as.matrix(otu_table(physeq.filtered.study)),
  metadata = ampvisMeta,
  taxonomy = tax_table(physeq.filtered.study)
  
)

amp_heatmap(
  data = d, 
  tax_aggregate = "Genus",
  facet_by="Timepoint",
  plot_values=FALSE
)
############################################################
#boxplot
############################################################
genus<- tax_glom(relAb.meanFilt.study, "genus")
df.genus <- psmelt(genus)
df.genus.Staph <- df.genus[df.genus$genus=="g__Staphylococcus",]
model.data <-data.frame(
  pig=factor(df.genus.Staph$Pair),
  infected=factor(df.genus.Staph$Infection),
  timepoint=factor(df.genus.Staph$Timepoint),
  taxon= df.genus.Staph$Abundance,
  groupTP=factor(paste0(df.genus.Staph$Infection,"_",df.genus.Staph$Timepoint)))

ggplot(model.data, aes(x=groupTP, y=taxon)) +   geom_boxplot()  + 
  theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Timepoint") + ylab(paste0("Relative abundance Staphylococcus"))

################################################################
#stats: paired t-test
#################################################################
genus<- tax_glom(relAb.meanFilt.sqrt.SyDN, "genus")
df.genus <- psmelt(genus)
df.genus

genera <- unique(df.genus$genus)
taxon  <- c()
pval <- c()
for (i in genera){
  df.genus.taxon <- df.genus[df.genus$genus==i,]
  model.data <-data.frame(
      pig=factor(df.genus.taxon$Pair),
      infected=factor(df.genus.taxon$Infection),
      timepoint=factor(df.genus.taxon$Timepoint),
      taxon= df.genus.taxon$Abundance)
  
  model.data
  model.data.split <-split(model.data, model.data$timepoint)

  model.data.merge <- merge(model.data.split$`TP-1`,model.data.split$`TP-2`,by="pig")
  model.data.merge
  

  if(median(model.data$taxon) > 0.01 & mean(c(model.data.merge$taxon.x, model.data.merge$taxon.y))){
    ttest.results <- t.test(model.data.merge$taxon.x,model.data.merge$taxon.y, paired = TRUE)
    pval <- c(pval, ttest.results$p.value)
    taxon <- c(taxon, i)
  }
}
pval.corr <- p.adjust(pval,method="fdr")
summary.df <- data.frame(pvalue=pval, taxon=taxon, fdr=pval.corr)
summary.df <- summary.df[!is.na(summary.df$pvalue),]

head(summary.df[order(summary.df$fdr),])



##########################################
# beta diversity dendrogram
###########################################

distance <- phyloseq::distance(physeq.rarefy.SyDN, method="bray", type="samples")
distance.matrix <- as.matrix(distance)
rownames(distance.matrix) <- paste0(labels(distance),"_",sample_data(physeq.rarefy.SyDN)$Timepoint)
colnames(distance.matrix) <- paste0(labels(distance),"_",sample_data(physeq.rarefy.SyDN)$Timepoint)
distance <- dist(distance.matrix)
dendro     <- hclust(distance, method="average")
plot(dendro, main="Hierarchical clustering of bray curtis distance", xlab="Sample", ylab="distance")


###############################################################
# pcoa
##################################################################

amp_ordinate(
  d, 
  type = "pcoa",
  distmeasure = "bray",
  sample_color_by = "Timepoint",
  sample_shape_by = "Infection",
  sample_colorframe = TRUE,
  sample_colorframe_label = "Timepoint",
) 


############################################################################
#alpha diversity - all time points
##############################################################################
comps <- list( c("control", "SyDn"), c("control", "SyDy"), c("SyDn", "SyDy") )
p <- plot_richness(physeq.rarefy.study, x = "Infection", color = "Infection", measures = c("Observed", "Shannon"))+ geom_boxplot() + theme_bw()
p + stat_compare_means(
  comparisons = comps,
  label = "p.signif",
  tip.length = 0.05,
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
    symbols = c("xxxx", "***", "**", "*", "ns")
  ),
  method = "wilcox.test")

############################################################################
#alpha diversity - facet
##############################################################################

alpha.metrics <- estimate_richness(physeq.rarefy.study, measures=c("Observed", "Shannon"))

model.data <-data.frame(
  infected=factor(sample_data(physeq.rarefy.study)$Infection),
  timepoint=factor(sample_data(physeq.rarefy.study)$Timepoint),
  diversity=alpha.metrics$Shannon,
  pig=factor(sample_data(physeq.rarefy.study)$Pair),
  it=factor(paste0(sample_data(physeq.rarefy.study)$Infection,"_",sample_data(physeq.rarefy.study)$Timepoint)))
ggplot(model.data, aes(x=timepoint, y=diversity)) +   geom_boxplot() +  facet_wrap(~ infected)


