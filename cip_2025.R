# workshop_analysis.R
# Author: ASP
# Requirements: phyloseq, tidyverse, ampvis2, ggpubr, patchwork

# -----------------------------
# 1) Libraries
# -----------------------------
library(phyloseq)
library(tidyverse)
library(ampvis2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(plotly)
# -----------------------------
# 2) Setup output folders
# -----------------------------
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# -----------------------------
# 3) Read data
# -----------------------------
asv <- read.table(url("https://raw.githubusercontent.com/marthazak/CiP/main/CiP_count_filtered.tsv"), header = TRUE, row.names = 1, check.names = FALSE)
# taxa in rows, samples in columns
asv_phylo <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)

taxonomy <- read.table(url("https://raw.githubusercontent.com/marthazak/CiP/main/CiP_taxonomy_filtered.tsv"), header = TRUE, sep = "\t", row.names = 1, fill = TRUE, check.names = FALSE)
tax_phylo <- tax_table(as.matrix(taxonomy))

metadata <- read.csv(url("https://raw.githubusercontent.com/marthazak/CiP/main/metadata.csv"), row.names = 1, check.names = FALSE)
# match sample order
metadata <- metadata[match(colnames(asv), rownames(metadata)), , drop = FALSE]
meta_phylo <- sample_data(metadata)

# phyloseq object
physeq <- phyloseq(asv_phylo, tax_phylo, meta_phylo)
physeq <- subset_samples(physeq, Infection %in% c("control", "SyDn"))

# -----------------------------
# 4) Filter & Normalize
# -----------------------------
# relative abundance
relAb <- transform_sample_counts(physeq, function(x) x/sum(x))

# filter taxa with mean relative abundance > 0.001%
meanThreshold <- 1e-5
relAbFilt <- filter_taxa(relAb, function(x) mean(x) > meanThreshold, TRUE)
physeqFilt <- prune_taxa(taxa_names(relAbFilt), physeq)

# rarefy to minimum sample size
minSample <- min(sample_sums(physeqFilt))
physeqRare <- rarefy_even_depth(physeqFilt, sample.size = minSample, rngseed = 123, replace = FALSE)

# sqrt transform (for stats)
relAbSqrt <- transform_sample_counts(relAbFilt, sqrt)

# -----------------------------
# 5) Phylum pie/donut plot
# -----------------------------
phylum <- tax_glom(physeqFilt, "phylum")
df_phylum <- psmelt(phylum) %>% 
  group_by(phylum) %>% 
  summarise(total = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(percent = total/sum(total)*100)

p_pie <- ggplot(df_phylum, aes(x = 1, y = percent, fill = phylum)) +
  geom_col(color = "white") + coord_polar(theta = "y") +
  theme_void() + labs(title = "Phylum composition") +
  theme(legend.title = element_blank())
p_pie
#ggsave("plots/phylum_pie.png", p_pie, width = 6, height = 6, dpi = 300)

# -----------------------------
# 6) Bubble plot top 15 genera
# -----------------------------
genus <- tax_glom(relAbFilt, "genus")
df_genus <- psmelt(genus)
top15 <- df_genus %>% group_by(genus) %>% summarise(mean_abund = mean(Abundance)) %>% slice_max(mean_abund, n=15) %>% pull(genus)
df_top <- df_genus %>% filter(genus %in% top15)
if(!"label" %in% colnames(df_top)) df_top$label <- sample_names(genus)

p_bubble <- ggplot(df_top, aes(x = label, y = genus, size = Abundance, fill = genus)) +
  geom_point(shape=21, alpha=0.9) + scale_size(range=c(2,18)) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), legend.position = "none") +
  labs(title = "Top 15 Genera — Bubble Plot", x="Sample", y="Genus")
p_bubble
#ggsave("plots/genus_bubble.png", p_bubble, width=11, height=6, dpi=300)

# -----------------------------
# 7) Stacked bar plot top genera
# -----------------------------
df_bar <- df_top %>% group_by(label, genus) %>% summarise(Abundance=sum(Abundance))

p_bar <- ggplot(df_bar, aes(x=label, y=Abundance, fill=genus)) +
  geom_col(color="white") +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(title="Genus composition", x="Sample", y="Relative abundance")
p_bar
#ggsave("plots/genus_barplot.png", p_bar, width=11, height=6, dpi=300)

# -----------------------------
# 8) Heatmap (ggplot)
# -----------------------------
df_heat <- df_top %>% group_by(label, genus) %>% summarise(Abundance=sum(Abundance))


p_heat <- ggplot(df_heat, aes(x=label, y=genus, fill=Abundance)) +
  geom_tile(color="white") +
  scale_fill_gradient(low="white", high="steelblue") +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.text.y = element_text(size=8)) +
  labs(title="Genus Heatmap", x="Sample", y="Genus")
p_heat
#ggsave("plots/genus_heatmap.png", p_heat, width=12, height=8, dpi=300)
# -----------------------------
# 9) Boxplot example for a specific taxon (e.g., Staphylococcus)
# -----------------------------
df_genus_staph <- df_genus %>% filter(str_detect(genus, "Staphylococcus"))
if(nrow(df_genus_staph)>0){
  df_box <- df_genus_staph %>% mutate(group = paste(Infection, Timepoint, sep="_"))
  p_box <- ggplot(df_box, aes(x=group, y=Abundance, fill=Infection)) + geom_boxplot() +
    theme(axis.text.x=element_text(angle=90, hjust=1)) + labs(title="Staphylococcus Relative Abundance", x="Group", y="Rel. abundance")
#  ggsave("plots/staph_boxplot.png", p_box, width=9, height=6, dpi=300)
}
p_box
# -----------------------------
# 10) Alpha diversity (Observed, Shannon)
# -----------------------------
alpha_metrics <- estimate_richness(physeqRare, measures=c("Observed", "Shannon"))
df_alpha <- data.frame(alpha_metrics,sample_data(physeqRare))

p_alpha <- ggplot(df_alpha, aes(x=Group1, y=Shannon, color=Infection)) + geom_boxplot(alpha=0.7) + labs(title="Alpha diversity")
p_alpha
#ggsave("plots/alpha_diversity.png", p_alpha, width=10, height=6, dpi=300)

# -----------------------------
# 11) Beta diversity dendrogram (Bray-Curtis)
# -----------------------------
bc_dist <- phyloseq::distance(physeqRare, method="bray")
dist_mat <- as.matrix(bc_dist)
rownames(dist_mat) <- sample_data(physeqRare)$Group1
hc <- hclust(as.dist(dist_mat), method="average")
#png("plots/bray_dendrogram.png", width=1200, height=800)
plot(hc, main="Bray-Curtis Dendrogram", xlab="Sample")
#dev.off()

# -----------------------------
# 12) PCoA 
# -----------------------------
bc_dist <- phyloseq::distance(physeqRare, method="bray")
dist_mat <- as.matrix(bc_dist)
rownames(dist_mat) <-sample_data(physeqRare)$Group1
hc <- hclust(as.dist(dist_mat), method="average")
plot(hc, main="Bray-Curtis Dendrogram", xlab="Sample")

pcoa_res <- ordinate(physeqRare, method="PCoA", distance="bray")
df_pcoa <- plot_ordination(physeqRare, pcoa_res, type="samples") +
  geom_point(aes(color=Infection, shape=Timepoint), size=4) +
  theme_minimal() + labs(title="PCoA 2D")
df_pcoa
#ggsave("plots/pcoa_2D.png", df_pcoa, width=8, height=6, dpi=300)


# 3D interactive PCoA using plotly
coords <- as.data.frame(pcoa_res$vectors[,1:3])
coords$Sample <- rownames(coords)
coords <- cbind(coords, metadata[coords$Sample, ])

coords$TP_Inf <- paste(coords$Timepoint, coords$Infection, sep = "_")
p_pcoa3D <- plot_ly(coords, x=~Axis.1, y=~Axis.2, z=~Axis.3, color=~TP_Inf, symbol=~TP_Inf, text=~Sample, type='scatter3d', mode='markers', marker=list(size=5))
p_pcoa3D
#htmlwidgets::saveWidget(p_pcoa3D, "plots/pcoa_3D.html")

# -----------------------------
# 13) Paired t-tests per genus (SyDn or subset)
# -----------------------------
genera <- unique(df_genus$genus)
results <- tibble::tibble(taxon=character(), pvalue=double(), fdr=double())
for(g in genera){
  dsub <- df_genus %>% filter(genus==g)
  if(!all(c("TP-1","TP-2") %in% unique(dsub$Timepoint))) next
  d1 <- dsub %>% filter(Timepoint=="TP-1") %>% select(Pair, Abundance)
  d2 <- dsub %>% filter(Timepoint=="TP-2") %>% select(Pair, Abundance)
  merged <- merge(d1,d2,by="Pair",suffixes=c('.x','.y'))
  if(nrow(merged)<3) next
  t_res <- try(t.test(merged$Abundance.x, merged$Abundance.y, paired=TRUE), silent=TRUE)
  if(inherits(t_res,"htest")) results <- bind_rows(results, tibble::tibble(taxon=g, pvalue=t_res$p.value, fdr=NA_real_))
}
if(nrow(results)>0){ results$fdr <- p.adjust(results$pvalue, method="fdr"); write_csv(results,"results/paired_ttests.csv") }

message("Workshop analysis completed. All plots in 'plots/', stats in 'results/'.")
