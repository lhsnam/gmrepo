
library(phyloseq)
library(dplyr)
library(ggplot2)

# import data
asv_table <- read.csv('asv_table.csv')
tax_table <- read.csv('tax_table.csv', sep = '\t')
met_table <- read.csv('GMrepoMetaData.csv', sep = '\t')

# ASV table
asv_table2 <- asv_table[,-1]
rownames(asv_table2) <- asv_table[,1]
asv_mat <- as.matrix(asv_table2)
ASV <- otu_table(t(asv_mat), taxa_are_rows = T)

# TAX table
tax_table2 <- tax_table[,-1]
tax_table2 <- tax_table2[,-1]

rownames(tax_table2) <- tax_table[,2]
colnames(tax_table2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_mat <- as.matrix(tax_table2)
tax_mat <- tax_mat[rownames(tax_mat) %in% rownames(t(asv_mat)), ]
TAX <- tax_table(tax_mat)

# Sample data
met_table2 <- met_table[,-1]
rownames(met_table2) <- met_table[,1]
META <- sample_data(met_table2)

# phyloseq
ps <- phyloseq(ASV, TAX, META)

plot_bar(ps, fill = "Phylum")

#PCoA
ord <- ordinate(ps, 'PCoA', 'bray')
plot_ordination(ps, ord, type = 'samples', color = 'Country')


# Alpha diversity
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(ps, measures = "Simpson"),
  "Chao1" = phyloseq::estimate_richness(ps, measures = "Chao1"),
  "Age" = sample_data(ps)$Host.age,
  "Disease" = sample_data(ps)$Disease.name)

adiv$Disease <- as.factor(adiv$Disease)

ggplot(adiv, aes(Age, Shannon, color = Disease)) +
  geom_point() +
  geom_smooth()

#Relative abundance
library(microbiome)
library(microbiomeutilities)

ps_phylum <- microbiome::aggregate_rare(ps, level = "Phylum", detection = 0.1, prevalence = 0)

abund_df <- microbiomeutilities::phy_to_ldf(ps_phylum, 
                                            transform.counts = "compositional")

ggplot(abund_df, aes(x = Sam_rep, y = Abundance, fill = Phylum)) +
  geom_bar(position = 'fill', stat = 'identity')

#Pick top phylum
top_phylum <- abund_df %>% 
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

top_phylum <- top_phylum %>% arrange(desc(mean_abundance))
top_phylum <- as.data.frame(top_phylum)
top.phyla <- c("p__Bacteroidota", "p__Firmicutes_A", "p__Proteobacteria", "p__Firmicutes_C", "p__Actinobacteriota")

abund_df <- abund_df %>% dplyr::mutate(pick_Phylum = ifelse(Phylum %in% top.phyla, yes = Phylum, "Others"))

# Plot
abund_df$pick_Phylum <- as.factor(abund_df$pick_Phylum)
abund_df$pick_Phylum <- factor(abund_df$pick_Phylum, levels = c(top.phyla, "Others")) 

abund_df <- abund_df[order(abund_df$Host.age, decreasing = F),]

abund_plot <- 
ggplot(abund_df, aes(x = Sam_rep, y = Abundance, fill = pick_Phylum)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual(values = c("p__Bacteroidota" = "#6388B4",
                               "p__Firmicutes_A" = "#FFAE34",
                               "p__Proteobacteria" = "#EF6F6A",
                               "p__Firmicutes_C" = "#8CC2CA",
                               "p__Actinobacteriota" = "#F1FAEE",
                               "Others" = "grey75"),
                    breaks = c(top.phyla, "Others"),
                    ) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

library(plotly)
ggplotly(abund_plot)

#genus level
ps_genus <- microbiome::aggregate_rare(ps, level = "Genus", detection = 0.1, prevalence = 0)

abund_df <- microbiomeutilities::phy_to_ldf(ps_genus, 
                                            transform.counts = "compositional")

ggplot(abund_df, aes(x = Sam_rep, y = Abundance, fill = Genus)) +
  geom_bar(position = 'fill', stat = 'identity')

#Pick top genus
top_genus <- abund_df %>% 
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance), se = se(Abundance), median = median(Abundance), q1 = quantile(Abundance)[2], q3 = quantile(Abundance)[4]) %>% 
  ungroup()

top_genus <- top_genus %>% arrange(desc(mean_abundance))
top_genus <- as.data.frame(top_genus)
top.genera <- c("g__Bacteroides", "g__Phocaeicola", "g__Prevotella", "g__Agathobacter", "g__Alistipes")

abund_df <- abund_df %>% dplyr::mutate(pick_Genus = ifelse(Genus %in% top.genera, yes = Genus, "Others"))

# Plot
abund_df$pick_Genus <- as.factor(abund_df$pick_Genus)
abund_df$pick_Genus <- factor(abund_df$pick_Genus, levels = c(top.genera, "Others")) 

abund_df <- abund_df[order(abund_df$Host.age, decreasing = F),]

abund_plot <- 
  ggplot(abund_df, aes(x = Sam_rep, y = Abundance, fill = pick_Genus)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual(values = c("g__Bacteroides" = "#6388B4",
                               "g__Phocaeicola" = "#FFAE34",
                               "g__Prevotella" = "#EF6F6A",
                               "g__Agathobacter" = "#8CC2CA",
                               "g__Alistipes" = "#F1FAEE",
                               "Others" = "grey75"),
                    breaks = c(top.genera, "Others"),
  ) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggplotly(abund_plot)


