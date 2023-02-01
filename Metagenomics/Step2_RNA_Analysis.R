# Libraries
library("ranacapa")
library("phyloseq")
library("ggplot2")
library("stringr")
library("plyr")
library("reshape2")
library("reshape")
library("dplyr")
library("tidyr")
library("doBy")
library("plyr")
library("microbiome")
library("ggpubr")
library("vegan")
library("tidyverse")
library("magrittr")
library("cowplot")
library("dendextend")
library("WGCNA")
library("metagenomeSeq")
library("decontam")
library("RColorBrewer")
library("ampvis2")


# Load data
raw <- import_biom("/Users/shashankgupta/Desktop/ImprovAFish/exported-feature-table/feature-table_taxonomy.biom")
tree <- read_tree("/Users/shashankgupta/Desktop/ImprovAFish/exported-feature-table/tree.nwk")
refseq <- Biostrings::readDNAStringSet("/Users/shashankgupta/Desktop/ImprovAFish/exported-feature-table/dna-sequences.fasta", use.names = TRUE)
dat <- read.table("/Users/shashankgupta/Desktop/ImprovAFish/metadata.txt", header = TRUE,row.names = 1, sep = "\t")
# Merge into one complete phyloseq object
all <- merge_phyloseq(raw, sample_data(dat), tree, refseq)
tax <- data.frame(tax_table(all), stringsAsFactors = FALSE)
tax <- tax[,1:7] # No info in col 8-15
# Set informative colnames
colnames(tax) <- c("Kingdom", "Phylum","Class","Order","Family","Genus", "Species")
library(stringr)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""), 
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""), 
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
# - Clean rank by rank
# Kingdom - Remove the unassigned completely
# Phylum
table(tax.clean$Phylum)
# Class
table(tax.clean$Class)
# Remove extra info about origin from some bacteria
# Remove all fields that contain "uncultured", "Unknown" or "Ambigious"
bad <- c("Ambiguous_taxa","uncultured", "Subgroup_21")
tax.clean[tax.clean$Class %in% bad,3:7] <- ""

# Order
table(tax.clean$Order)

bad <- c("0319-6G20","1-20","11-24", "ADurb.Bin180","D8A-2", "Group_1.1c", "JGI_0000069-P22","Marine_Group_II",
         "Pla3_lineage","Run-SP154",
         "Ambiguous_taxa", "uncultured", "UBA10353_marine_group",
         "Subgroup_17", "SAR86_clade", "SAR11_clade", "SAR202_clade", "Chloroplast")
tax.clean[tax.clean$Order %in% bad,4:7] <- ""

# Family
table(tax.clean$Family)
bad <- c("Ambiguous_taxa","11-24","67-14", "uncultured", "SAR116_clade", "Run-SP154", 
         "Marine_Group_II", "env.OPS_17", "SAR116_clade", "S085", "S-70", "NS9_marine_group", "Mitochondria")
tax.clean[tax.clean$Family %in% bad,5:7] <- ""

# Genus
table(tax.clean$Genus)
bad <- c("Ambiguous_taxa","Unknown_Family","uncultured","Subgroup_10", "1174-901-12", "67-14")
tax.clean[tax.clean$Genus %in% bad,6:7] <- ""

# Species
table(tax.clean$Species)
bad <- c("Ambiguous_taxa","marine_metagenome","low_GC","wastewater_metagenome","unidentified", 
         "uncultured_synthetic", "uncultured_organism")
tax.clean[tax.clean$Species %in% bad,6:7] <- ""

#tax.clean[grepl("uncultured", tax.clean$Species),"Species"] <- ""
#tax.clean[grepl("unidentified", tax.clean$Species),"Species"] <- ""

# Remove remove ".", change "-" and " " to "_"
for (i in 1:ncol(tax.clean)){
  tax.clean[,i] <- str_replace_all(tax.clean[,i], "[.]","")
  tax.clean[,i] <- str_replace_all(tax.clean[,i], "[(]","")
  tax.clean[,i] <- str_replace_all(tax.clean[,i], "[)]","")
  tax.clean[,i] <- str_replace_all(tax.clean[,i], "-","_")
  tax.clean[,i] <- str_replace_all(tax.clean[,i], " ","_")
}

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
# File holes in the tax table
for (i in 1:nrow(tax.clean)){
  #  Fill in missing taxonomy
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
    tax.clean$Species[i] <- paste("Genus_",tax.clean$Genus[i], sep = "_")
  }
}

rm(bad, class, family, i, kingdom,new,order,phylum,uncul)

tax_table(all) <- as.matrix(tax.clean)
all
#Remove Unnassigned
#all.clean <- subset_taxa(all, Kingdom != "Archaea")
#all.clean <- subset_taxa(all.clean, Kingdom != "Eukaryota")
all.clean <- subset_taxa(all, Kingdom != "Unassigned")
all.clean <- prune_taxa(taxa_sums(all.clean) > 0, all.clean)
all.clean

# Rarefaction plot
p <- ggrare(all.clean, step = 1000, color = "New_Diet", label = "samplingTime", se = FALSE)
p + facet_wrap(~New_Diet)+ theme_bw()

# Supplementary figures
# Enter variables for calculating rarefaction curves
psdata <- all.clean  # The phyloseq object containing the data
measures <- c("Observed", "Shannon") #Which measures should the rarefaction curves be calculated for (phyloseq naming)
depths <- rep(c(1, 100, c(1, 2, 3, 4, 5, 7.5, 10, 12.5, 15)*1000), each = 2)  #Which depths should be included
variables <- c("samplingTime", "New_Diet") #First variable for colour, second for faceting

# this enables automatic addition of the Depth to the output by ldply
names(depths) <- depths 

#This step runs the estimate rarefied richness function on the samples at all values in depth
rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none')) 

# convert Depth from factor to numeric
rarefaction_curve_data$Depth <- as.numeric(as.character(rarefaction_curve_data$Depth))
rarefaction_curve_data$Sample <- as.character(rarefaction_curve_data$Sample)
rarefaction_curve_data$Sample <- gsub('\\.', '-', rarefaction_curve_data$Sample)

#Add sample data
rarefaction_curve_verbose <- merge(rarefaction_curve_data, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')

#Summarize alpha diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_verbose, c("Depth", variables, "Measure"), 
                                        summarise, Alpha_diversity_mean = mean(Alpha_diversity), 
                                        Alpha_diversity_sd = sd(Alpha_diversity),
                                        Alpha_diversity_sem = sd(Alpha_diversity)/sqrt(length(Alpha_diversity)))

#Create plot
if (length(variables) == 2){
  p = ggplot( data = rarefaction_curve_data_summary, 
              mapping = aes( x = Depth, y = Alpha_diversity_mean, 
                             ymin = Alpha_diversity_mean - Alpha_diversity_sd, 
                             ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = get(variables[1]), group = get(variables[1]))) + 
    geom_line() + geom_errorbar(width = 750) + 
    facet_grid(facets = Measure ~ get(variables[2]), scales = 'free_y') +
    ylab("Alpha diversity (mean)") + 
    labs(colour = variables[1]) +
    xlab("Number of reads") 
} else if (length(variables == 1)){
  p = ggplot( data = rarefaction_curve_data_summary, 
              mapping = aes( x = Depth, y = Alpha_diversity_mean, 
                             ymin = Alpha_diversity_mean - Alpha_diversity_sd, 
                             ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = get(variables[1]), group = get(variables[1]))) + 
    geom_line() + geom_errorbar(width = 750) + 
    ylab("Alpha diversity (mean)") + 
    labs(colour = variables[1]) +
    xlab("Number of reads") 
}

library(RColorBrewer)
cols  <- c(brewer.pal(8,"Set1"), brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(7,"Accent"),brewer.pal(12,"Paired"),"gray")

# Plot and remember to save
p <- p + theme_bw() +   scale_fill_manual(values =cols) + scale_colour_manual( values = cols)
p


# Load Required Libraries
library(ggpubr)

# Estimate richness of all.clean
shannon.div <- estimate_richness(all.clean, measures = c("Shannon", "Simpson", "Observed","Chao1"))

# Get sample data
sampledata1<- data.frame(sample_data(all.clean))

# Rename row names
row.names(shannon.div) <- gsub("[.]","-", row.names(shannon.div))

# Merge data
sampleData <- merge(sampledata1, shannon.div, by = 0 , all = TRUE)

# Factorize New_Diet
sampleData$New_Diet <- factor(sampleData$New_Diet, levels=c('ext-ctrl', 'CTR', 'MC1', 'MC2', 'MN3'))

# List of comparisons
my_comparisons <- list( c("ext-ctrl", "MC1"), c("ext-ctrl", "MC2"), c("ext-ctrl", "MN3"),
                        c("CTR", "MC1"), c("CTR", "MC2"),
                        c("CTR", "MN3"))

# Create Observed Richness Plot
p1 <- ggboxplot(sampleData, x = "New_Diet", y = "Observed",
                color = "New_Diet", palette = "jco", legend = "none") + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 400) +
  geom_jitter(aes(colour = New_Diet), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank()) +  
  scale_fill_manual(values = cols) + 
  scale_colour_manual( values = cols)

# Create Shannon Plot
ggboxplot(sampleData, x = "New_Diet", y = "Shannon",
          color = "New_Diet", palette = "jco", legend = "none") + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 8) +
  geom_jitter(aes(colour = New_Diet), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())


# The code appears to be an analysis of microbial community data using the 
# R packages phyloseq and ggplot2. The script performs a Principal Coordinate Analysis (PCoA) 
# on the Bray-Curtis and weighted UniFrac distances between samples and visualizes the 
# results as scatter plots colored by a categorical variable "New_Diet". The last line 
# performs an analysis of variance (Adonis) on the Bray-Curtis distance matrix, with the
# response variable being the "New_Diet". The results are shown as a scatter plot with 
# ellipses representing 95% confidence intervals.

set.seed(1)

# Calculate PCoA on Bray distance
PCoA_bray <- ordinate(physeq = all.clean, method = "PCoA", distance = "bray")

# Plot PCoA on Bray distance
PCoA_bray_plot <- plot_ordination(
  physeq = all.clean, 
  ordination = PCoA_bray, 
  color = "New_Diet"
) + 
  geom_point(shape = 19, alpha=0.7) + 
  theme_bw() + ggtitle("PCoA Plot - Bray") + 
  xlab("PCoA 1 [17.4 %]") + ylab("PCoA 2 [8.8 %]") + 
  stat_ellipse() + scale_fill_manual(values =cols) + 
  scale_colour_manual( values = cols)

# Calculate PCoA on weighted unifrac distance
PCoA_wunifrac <- ordinate(physeq = all.clean, method = "PCoA", distance = "wunifrac")

# Plot PCoA on weighted unifrac distance
PCoA_wunifrac_plot <- plot_ordination(
  physeq = all.clean, 
  ordination = PCoA_wunifrac, 
  color = "New_Diet"
) + 
  geom_point(shape = 19, alpha=0.7) + 
  theme_bw() + ggtitle("PCoA Plot - weighted unifrac") + 
  xlab("PCoA 1 [68.2 %]") + ylab("PCoA 2 [8.3 %]") + 
  stat_ellipse() + scale_fill_manual(values =cols) + 
  scale_colour_manual( values = cols)

# Grid plot of PCoA plots
bottom_row <- plot_grid(p1, PCoA_bray_plot, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1.3))
plot_grid(p, bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1.2))

# Run adonis test
sampledf <- data.frame(sample_data(all.clean))
bcdist <- phyloseq::distance(all.clean, method="bray",normalized=TRUE) 
adonis2(bcdist ~ New_Diet, data = sampledf, permutations = 9999)


#Contamination removal
df <- as.data.frame(sample_data(all.clean)) # Put sample_data into a ggplot-friendly data.frame
df$Sample_or_Control <- ifelse( df$New_Diet  %in% c("ext-ctrl"), "Control_Sample", "True_Sample")
sample_data(all.clean) <- df
df$LibrarySize <- sample_sums(all.clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() + theme_bw()

# Identify Contaminants - Prevalence
sample_data(all.clean)$is.neg <- sample_data(all.clean)$Sample_or_Control == "Control_Sample"
contamdf.prev <- isContaminant(all.clean, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

#more aggressive classification threshold rather than the default. i.e., 0.05
contamdf.prev05 <- isContaminant(all.clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
all.noncontam <- prune_taxa(!contamdf.prev05$contaminant, all.clean)
all.noncontam

# Bray-Curtis PCoA after contamination removal
PCoA_bray <- ordinate(physeq = all.noncontam, method = "PCoA", distance = "bray")
PCoA_bray_plot <- plot_ordination(
  physeq = all.noncontam, 
  ordination = PCoA_bray, 
  color = "New_Diet"
) + 
  geom_point(shape = 19, alpha = 0.7) + 
  theme_bw() + 
  ggtitle("PCoA Plot - Bray") + 
  xlab("PCoA 1 [17.6 %]") + 
  ylab("PCoA 2 [9 %]") + 
  stat_ellipse()


# Weighted UniFrac PCoA
PCoA_wunifrac <- ordinate(physeq = all.noncontam, method = "PCoA", distance = "unifrac")
PCoA_wunifrac_plot <- plot_ordination(
  physeq = all.noncontam, 
  ordination = PCoA_wunifrac, 
  color = "New_Diet"
) + 
  geom_point(shape = 19, alpha = 0.7) + 
  theme_bw() + 
  ggtitle("PCoA Plot - Weighted UniFrac") + 
  xlab("PCoA 1 [3.3 %]") + 
  ylab("PCoA 2 [2.3 %]") + 
  stat_ellipse()



psdata <- subset_samples(all.noncontam, Sample_or_Control=="True_Sample")
psdata <- prune_taxa(taxa_sums(psdata) > 0, psdata)
psdata


# Inspect the number of reads per sample and compare to rarefaction curves
sample_data(psdata)$reads <- unlist(sample_sums(psdata))
dat1 <- data.frame(sample_data(psdata))
table(dat1$reads)

# Set a cutoff where Shannon diversity dont increase, and observed increases markedly slower.
# At the same time don't throw too many samples out. Set the cutoff value accordingly
cutoff <- 0

# See which samples have been removed
dat1[dat1$reads < cutoff,]

# Remove the samples with fewer reads than the cutoff
psdata.p <- prune_samples(sample_sums(psdata) > cutoff, psdata)
psdata.p <- prune_taxa(taxa_sums(psdata.p) >0, psdata.p)
psdata.p

# Rarefy the samples using the function multiple_rarefy
psdata.r <- multiple_rarefy(psdata.p)
#psdata.r = rarefy_even_depth(psdata.p, rngseed=1, sample.size=0.9*min(sample_sums(psdata.p)), replace=F)
#psdata.r = rarefy_even_depth(psdata.p)
psdata.r<- transform_sample_counts(psdata.p, function(x) x / sum(x) )
psdata.r <- prune_taxa(taxa_sums(psdata.r) > 0, psdata.r)
psdata.r
table(sample_sums(psdata.r))
rm(all)
rm(all.clean)
rm(all.noncontam)

ggrare(psdata.p, step = 1000, color = "New_Diet", label = "samplingTime", se = FALSE)+ facet_wrap(~New_Diet)+ theme_bw()

#### 
# Remove Time point T0
psdata.r <- subset_samples(psdata.r, samplingTime != "T0")
psdata.r <- prune_taxa(taxa_sums(psdata.r) > 0, psdata.r)
psdata.r

################# RNA
#Bar plot
#phylum plot
# Final.RNA_phylum <- psdata.r %>%
#   aggregate_taxa(level = "Phylum") %>%
#   transform(transform = "compositional")
# Final.RNA_phylum <- phyloseq_filter_taxa_tot_fraction(Final.RNA, frac = 0.01)
# Final.RNA_phylum


#Final.RNA <- transform(psdata.r, "compositional")
Final.RNA <- aggregate_rare(psdata.r, level = "Phylum", detection = 1/100, prevalence = 20/100)

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(10, "Dark2")) 
PhylaPalette = getPalette(10)



Final.RNA_phylum_plot<- plot_composition(Final.RNA, sample.sort = "Proteobacteria",otu.sort = "abundance", verbose = TRUE)
Final.RNA_phylum_plot <- Final.RNA_phylum_plot + theme_bw() + theme(axis.text.x=element_blank(),
                                                                    axis.ticks.x=element_blank())+ scale_fill_manual(values = PhylaPalette)
Final.RNA_phylum_plot


#Bacterial Community Composition for Manuscript
Final.seq.melt.RNA <- psmelt(tax_glom(psdata.r, "Species"))
tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

for (rank in tax_ranks) {
  n_unique <- length(unique(Final.seq.melt.RNA[[rank]]))
  message(paste(rank, ": ", n_unique, sep = ""))
}


table(grepl("Phylum", unique(Final.seq.melt.RNA$Class)))
table(grepl("Class|Phylum", unique(Final.seq.melt.RNA$Order)))
table(grepl("Order|Class|Phylum", unique(Final.seq.melt.RNA$Family)))
table(grepl("Family|Order|Class|Phylum", unique(Final.seq.melt.RNA$Genus)))
table(grepl("Family|Order|Class|Phylum|Genus", unique(Final.seq.melt.RNA$Species)))


library(doBy)
# Final.psdata.r.RNA.phylum <- tax_glom(Final.RNA, "Phylum")
# Final.psdata.r.RNA.phylum.psmelt <- psmelt(Final.psdata.r.RNA.phylum)
Phylum_df <- summaryBy(Abundance~Phylum, data=Final.seq.melt.RNA, FUN=sum)
Phylum_df$Percent <- round(Phylum_df$Abundance.sum/sum(Phylum_df$Abundance.sum)*100, 4)
Phylum_df <- plyr::arrange(Phylum_df, plyr::desc(Percent))
Phylum_df$PercentageRound <- round(Phylum_df$Percent, digits = 2)
head(Phylum_df)


# Final.psdata.r.RNA.class <- tax_glom(Final.RNA, "Class")
# Final.psdata.r.RNA.class.psmelt <- psmelt(Final.psdata.r.RNA.class)
class_df <- summaryBy(Abundance~Phylum+Class, data=Final.seq.melt.RNA, FUN=sum)
class_df$Percent <- round(class_df$Abundance.sum/sum(class_df$Abundance.sum)*100, 4)
class_df <- plyr::arrange(class_df, plyr::desc(Percent))
class_df$Round <- round(class_df$Percent, digits = 2)
head(class_df)


# Final.psdata.r.RNA.order <- tax_glom(Final.RNA, "Order")
# Final.psdata.r.RNA.order.psmelt <- psmelt(Final.psdata.r.RNA.order)
order_df <- summaryBy(Abundance~Order, data=Final.seq.melt.RNA, FUN=sum)
order_df$Percent <- round(order_df$Abundance.sum/sum(order_df$Abundance.sum)*100, 4)
order_df <- plyr::arrange(order_df, plyr::desc(Percent))
order_df

# Final.psdata.r.RNA.genus <- tax_glom(Final.RNA, "Genus")
# Final.psdata.r.RNA.genus.psmelt <- psmelt(Final.psdata.r.RNA.genus)
genus_df <- summaryBy(Abundance~Genus, data=Final.seq.melt.RNA, FUN=sum)
genus_df$Percent <- round(genus_df$Abundance.sum/sum(genus_df$Abundance.sum)*100, 4)
genus_df <- plyr::arrange(genus_df, plyr::desc(Percent))
genus_df$Round <- round(genus_df$Percent, digits = 2)

View(genus_df)


library(ggpubr)
psdata1 <- subset_samples(psdata, samplingTime %in% c("T1", "T2", "T3"))
psdata1 <- prune_taxa(taxa_sums(psdata1) > 0, psdata1)
shannon.div <- estimate_richness(psdata1, measures = c("Shannon", "Observed"))
sampledata1<- data.frame(sample_data(psdata1))
row.names(shannon.div) <- gsub("[.]","-", row.names(shannon.div))
sampleData <- merge(sampledata1, shannon.div, by = 0 , all = TRUE)


sampleData$New_Diet <- factor(sampleData$New_Diet, levels=c( 'CTR', 'MC1', 'MC2', 'MN3'))
levels(sampleData$New_Diet)

sampleData$samplingTime <- factor(sampleData$samplingTime, levels=c('T1', 'T2', 'T3'))
levels(sampleData$samplingTime)


library(cowplot)

my_comparisons <- list( c("CTR", "MC1"), c("CTR", "MC2"), c("CTR", "MN3"), 
                        c("MC1", "MC2"), c("MC1", "MN3"), c("MC2", "MN3") )

p1 <- ggboxplot(sampleData, x = "New_Diet", y = "Shannon",
                color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6.5) +
  #geom_jitter(aes(colour = New_Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) + 
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank()) + facet_wrap("samplingTime")


#my_comparisons <- list( c("T0", "T1"), c("T0", "T2"), c("T0", "T3"), 
#                        c("T1", "T2"), c("T1", "T3"), c("T2", "T3") )

my_comparisons <- list(c("T1", "T2"), c("T1", "T3"), c("T2", "T3") )


p4 <- ggboxplot(sampleData, x = "samplingTime", y = "Shannon",
                color = "samplingTime", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 5.5) +
  #geom_jitter(aes(colour = samplingTime), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = samplingTime), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank()) + facet_wrap(. ~ New_Diet, ncol = 4)


plot_grid(p4, p1, labels = c('A', 'B'), label_size = 12, ncol = 2)




#Beta Diversity Unifrac PCoA
set.seed(1)
PCoA_bray <- ordinate(physeq = psdata.r, method = "PCoA", distance = "unifrac")
p2 <- plot_ordination(
  physeq = psdata.r,
  ordination = PCoA_bray, 
  color = "samplingTime"
) + 
  geom_point(shape = 19, alpha=0.7) + theme_bw() + ggtitle("PCoA Plot - Unifrac") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("PCoA 1 [4.5 %]") + ylab("PCoA 2 [2.6 %]") + 
  scale_color_manual(values = c("#0073c2cc", "#efc000cc", "#868686cc", "#cd534ccc")) + stat_ellipse()


p2 



set.seed(1)
PCoA_bray <- ordinate(physeq = psdata.r, method = "PCoA", distance = "Unifrac")
p3 <- plot_ordination(
  physeq = psdata.r, 
  ordination = PCoA_bray, 
  color = "New_Diet"
) + 
  geom_point(shape = 19, alpha=0.7) + theme_bw() + ggtitle("PCoA Plot - Unifrac")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("PCoA 1 [4.5%]") + ylab("PCoA 2 [2.6 %]")+ 
  scale_color_manual(values = c("#0073c2cc", "#efc000cc", "#868686cc", "#cd534ccc")) + stat_ellipse()

p3  


p5 <- plot_grid(p1, p4, labels = c('A', 'B'), label_size = 12, ncol = 1, align = "h")

p6 <- plot_grid(p3, p2, labels = c('C', 'D'), label_size = 12, ncol = 1, align = "hv")

plot_grid(p5, p6, ncol = 2, rel_widths = c(1.5, 1))

#Adonis
library(vegan)
sampledf <- data.frame(sample_data(psdata.r))
bcdist <- phyloseq::distance(psdata.r, method="bray",normalized=TRUE) 
adonis2(bcdist ~ samplingTime, 
        data = sampledf, permutations = 9999)

wudist <- phyloseq::distance(psdata.r, method="Unifrac",normalized=TRUE) 
adonis2(wudist ~ samplingTime, 
        data = sampledf, permutations = 9999)

adonis2(wudist ~ New_Diet, 
        data = sampledf, permutations = 9999)

adonis2(wudist ~ New_Diet + samplingTime + samplingTank, 
        data = sampledf, permutations = 9999)




my_comparisons <- list( c("T0", "T1"), c("T0", "T2"), c("T0", "T3"), 
                        c("T1", "T2"), c("T1", "T3"), c("T2", "T3") )
ggboxplot(sampleData, x = "samplingTime", y = "Shannon",
          color = "samplingTime", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6) +
  geom_boxplot(aes(fill = samplingTime), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  + facet_grid(~New_Diet) 



my_comparisons <- list( c("CTR", "MC1"), c("CTR", "MC2"), c("CTR", "MN3"), 
                        c("MC1", "MC2"), c("MC1", "MN3"), c("MC2", "MN3") )
ggboxplot(sampleData, x = "New_Diet", y = "Observed",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 350) + #geom_jitter(aes(colour = New_Diet), size = 3, alpha = 0.8) + 
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  + facet_grid(~samplingTime) 

new_df=sampleData[!grepl("CTR",sampleData$sampleDiet),]
my_comparisons <- list(c("MC1", "MC2"), c("MC1", "MN3"), c("MC2", "MN3") )

ggboxplot(new_df, x = "New_Diet", y = "Shannon",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 5.5) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  + facet_grid(~samplingTime) 


my_comparisons <- list(  c("T1", "T2"), c("T1", "T3"), c("T2", "T3") )
ggboxplot(new_df, x = "samplingTime", y = "Shannon",
          color = "samplingTime", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 5.5) +
  geom_boxplot(aes(fill = samplingTime), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  + facet_grid(~New_Diet) 


library(microbiomeutilities)
library(viridis)
library(RColorBrewer)


heat.sample <- plot_taxa_heatmap(psdata.r,
                                 subset.top = 50,
                                 VariableA = "New_Diet",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")




ampvis2_obj <- phyloseq_to_ampvis2(psdata.r)

amp_heatmap(ampvis2_obj,
            group_by = "New_Diet",
            facet_by = "samplingTime",
            tax_aggregate = "Genus",
            tax_add = c("Kingdom", "Phylum"),
            tax_show = 50,
            color_vector = c("white", "dark red"),
            plot_colorscale = "sqrt",
            plot_values = T) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")



# Differential abundance analysis
rabuplot(phylo_ob = psdata.r, predictor= "New_Diet", type = "Phylum", facet_wrap   ="samplingTime")
sample_data(psdata.r)$merged <- paste(sample_data(psdata.r)$New_Diet,"-",sample_data(psdata.r)$samplingTime)
rabuplot(phylo_ob = psdata.r, predictor= "merged", type = "Phylum", N_taxa= 20 )

ps1 <- subset_samples(psdata.r, samplingTime %in% c("T1"))
ps1 <- prune_taxa(taxa_sums(ps1) >0, ps1)
A1 <- rabuplot(phylo_ob = ps1, predictor= "New_Diet", type = "Phylum",facet_wrap   ="samplingTime") + theme(legend.position = "none")

ps1 <- subset_samples(psdata.r, samplingTime %in% c("T2"))
ps1 <- prune_taxa(taxa_sums(ps1) >0, ps1)
A2 <-rabuplot(phylo_ob = ps1, predictor= "New_Diet", type = "Phylum",facet_wrap   ="samplingTime") + theme(legend.position = "none")

ps1 <- subset_samples(psdata.r, samplingTime %in% c("T3"))
ps1 <- prune_taxa(taxa_sums(ps1) >0, ps1)
A3 <-rabuplot(phylo_ob = ps1, predictor= "New_Diet", type = "Phylum",facet_wrap   ="samplingTime")

cowplot::plot_grid(A1, A2, A3, ncol = 3, rel_widths = c(1,1,  1.3))


###################################################################################################################################################################
#host-microbiome interactions
# Set seed as a precaution for reproducibility as some methods are non-deterministic.
set.seed(13118)
Run_analysis <- TRUE # if FALSE it tries to load data instead of running
save_plots <- TRUE
plot_labeling_size = 15
prefix <- "m"
save_TOM <- TRUE
pam_stage <- FALSE    # Partitioning around medoids, tries to put more OTUs into modules where it is not directly clear from the dendrogram cutting.
set_verbose = 1 # How much detail from the WGCNA functions? Higher means more detail.
omics_type = "otu" 
take_average <- F
max_block_size = 10000      #  OTU data is small enough that this is not necessary.
applied_filter = 0.005      #  filter lowest x percent of total sum,  100 is 100 percent.
parameter_sets <- list(set_1 = list(applied_norm = "TSS", applied_transf = "CLR", assoc_measure = "bicor"),
                       set_2 = list(applied_norm = "CSS", applied_transf = "log2", assoc_measure = "bicor"))
chosen_parameter_set <- parameter_sets$set_2



#host-microbiome interactions
psdata.p <- prune_samples(sample_sums(psdata) > cutoff, psdata)
psdata.p <- filter_taxa(psdata.p, function(x) sum(x > 3) > (0.01*length(x)), TRUE)# Remove taxa not seen more than 3 times in at least 5% of samples
psdata.p

psdata.p <- transform_sample_counts(psdata.p, function(x) x / sum(x) )
#psdata.p = filter_taxa(psdata.p, function(x) mean(x) > 1e-4, TRUE)
psdata.p



# Rarefy the samples using the function multiple_rarefy
#psdata.r <- multiple_rarefy(psdata.p)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(face="bold"))

gut_16SrRNA_mapping_file <- sample_data(psdata.p)
omics_data <- otu_table(psdata.p)
# host_gut_mapping_file <- read.table("../host_data/host_gut_mapping.txt", sep = "\t", row.names = 1, header = T)
# host_gut_rawCounts_total1 <- read.table("../host_data/host_gut_count.txt")

host_gut_mapping_file <- sample_info
host_gut_rawCounts_total1 <- omics_data_host

gut_16SrRNA_mapping_file$New_Diet <- tolower(gut_16SrRNA_mapping_file$New_Diet)
gut_16SrRNA_mapping_file$ID_New <- paste(gut_16SrRNA_mapping_file$samplingTime, "_",
                                         gut_16SrRNA_mapping_file$New_Diet, "_", 
                                         gut_16SrRNA_mapping_file$samplingTank, "_", "0",
                                         gut_16SrRNA_mapping_file$sampleNumber, 
                                         sep = "")


#gut_16SrRNA_mapping_file$ID_New <- gsub("[-]","_", gut_16SrRNA_mapping_file$sampleName)
#gut_16SrRNA_mapping_file$ID_New <- gsub("_HGm","_", gut_16SrRNA_mapping_file$ID_New)


# host_gut_mapping_file$sampleName <- row.names(host_gut_mapping_file)
# host_gut_mapping_file$common <- host_gut_mapping_file$sampleName
# host_gut_mapping_file$common <- gsub("_HGh","_", host_gut_mapping_file$common)


table(gut_16SrRNA_mapping_file$ID_New %in% host_gut_mapping_file$ID_New)
table(host_gut_mapping_file$ID_New %in% gut_16SrRNA_mapping_file$ID_New)

host_gut_mapping_file <- host_gut_mapping_file[host_gut_mapping_file$ID_New %in% gut_16SrRNA_mapping_file$ID_New, ]
gut_16SrRNA_mapping_file <- gut_16SrRNA_mapping_file[gut_16SrRNA_mapping_file$ID_New %in% host_gut_mapping_file$ID_New, ]


gut_16SrRNA_mapping_file$rownames <- row.names(gut_16SrRNA_mapping_file)
row.names(gut_16SrRNA_mapping_file) <- gut_16SrRNA_mapping_file$ID_New



ID_New_name_list <-  gut_16SrRNA_mapping_file$rownames
gut_16SrRNA_OTU_table <- omics_data[, colnames(omics_data) %in% ID_New_name_list]

table(colnames(gut_16SrRNA_OTU_table) ==  gut_16SrRNA_mapping_file$rownames)
colnames(gut_16SrRNA_OTU_table) <- gut_16SrRNA_mapping_file$ID_New
table(colnames(gut_16SrRNA_OTU_table) ==  gut_16SrRNA_mapping_file$ID_New)

omics_data <- gut_16SrRNA_OTU_table

####################################################################################################
####################################################################################################
#Step 1

dim(omics_data)
omics_data <- t(omics_data)
dim(omics_data) %>% paste(c( "Samples", "OTUs"))

HellingerData<-decostand(omics_data,method = "hellinger")
omics_data <- HellingerData

dim(omics_data) %>% paste(c("Samples", "OTUs"))

powers <- c(1:10, seq(12,20,2))
suppressWarnings(sft <- pickSoftThreshold(omics_data, 
                                          powerVector = powers, 
                                          verbose = set_verbose, 
                                          networkType = "signed",
                                          corFn= chosen_parameter_set$assoc_measure))


# Find the soft thresholding power beta to which co-expression similarity is raised to calculate adjacency.
# based on the criterion of approximate scale-free topology.

idx <- min(which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) > 0.90))
if(is.infinite(idx)){
  idx <- min(which((-sign(sft$fitIndices[,3])*sft$fitIndices[,2]) > 0.80))
  if(!is.infinite(idx)){
    st <- sft$fitIndices[idx,1]
  } else{
    idx <- which.max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
    st <- sft$fitIndices[idx,1]
  }
} else{
  st <- sft$fitIndices[idx,1]
}


# Plot Scale independence measure and Mean connectivity measure

# Scale-free topology fit index as a function of the soft-thresholding power
data.frame(Indices = sft$fitIndices[,1],
           sfApprox = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]) %>% 
  ggplot() + 
  geom_hline(yintercept = 0.9, color = "red", alpha = 0.6) + # corresponds to R^2 cut-off of 0.9
  geom_hline(yintercept = 0.8, color = "red", alpha = 0.2) + # corresponds to R^2 cut-off of 0.8
  geom_line(aes(x = Indices, y = sfApprox), color = "red", alpha = 0.1, size = 2.5) +
  geom_text(mapping = aes(x = Indices, y = sfApprox, label = Indices), color = "red", size = 4) +
  ggtitle("Scale independence") +
  xlab("Soft Threshold (power)") +
  ylab("SF Model Fit,signed R^2") +
  xlim(1,20) +
  ylim(-1,1) +
  geom_segment(aes(x = st, y = 0.25, xend = st, yend = sfApprox[idx]-0.05), 
               arrow = arrow(length = unit(0.2,"cm")), 
               size = 0.5)-> scale_independence_plot 




# Mean connectivity as a function of the soft-thresholding power

data.frame(Indices = sft$fitIndices[,1],
           meanApprox = sft$fitIndices[,5]) %>% 
  ggplot() + 
  geom_line(aes(x = Indices, y = meanApprox), color = "red", alpha = 0.1, size = 2.5) +
  geom_text(mapping = aes(x = Indices, y = meanApprox, label = Indices), color = "red", size = 4) +
  xlab("Soft Threshold (power)") +
  ylab("Mean Connectivity") +
  geom_segment(aes(x = st-0.4, 
                   y = sft$fitIndices$mean.k.[idx], 
                   xend = 0, 
                   yend = sft$fitIndices$mean.k.[idx]),
               arrow = arrow(length = unit(0.2,"cm")), 
               size = 0.4) +
  ggtitle(paste0("Mean connectivity: ", 
                 round(sft$fitIndices$mean.k.[idx],2))) -> mean_connectivity_plot


cowplot::plot_grid(scale_independence_plot, mean_connectivity_plot, ncol = 2, align = "h", labels = c("A", "B"), label_size = plot_labeling_size) -> si_mc_plot

si_mc_plot





allowWGCNAThreads()

modules.omics.Y <- blockwiseModules(omics_data,
                                    power = st, 
                                    networkType = "signed", 
                                    TOMType = "signed",
                                    corType = chosen_parameter_set$assoc_measure,
                                    #maxPOutliers = 0.05,
                                    #deepSplit = 4, # Default 2
                                    minModuleSize = 12, # 30
                                    #minCoreKME = 0.5,      # Default 0.5
                                    #minCoreKMESize = 2,    # Default minModuleSize/3,
                                    #minKMEtoStay = 0.5,    # Default 0.3
                                    #reassignThreshold = 0, # Default 1e-6
                                    #mergeCutHeight = 0.4,  # Default 0.15
                                    #pamStage = pam_stage, 
                                    #pamRespectsDendro = TRUE,
                                    #replaceMissingAdjacencies = TRUE,
                                    numericLabels = TRUE,
                                    saveTOMs = save_TOM,
                                    saveTOMFileBase = "TOM",
                                    verbose = 3,
                                    maxBlockSize=8000, nThreads = 10)


rownames(modules.omics.Y$MEs) <- rownames(omics_data)
names(modules.omics.Y$colors) <- colnames(omics_data)
names(modules.omics.Y$unmergedColors) <- colnames(omics_data)

hubs <- chooseTopHubInEachModule(omics_data, modules.omics.Y$colors, power = st, omitColors = "0")

stage2results_Y <- list(modules = modules.omics.Y, 
                        hubs = hubs)




# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_Y$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_Y$modules$goodSamples) == length(stage2results_Y$modules$goodSamples)
genes_good <- sum(stage2results_Y$modules$goodGenes) == length(stage2results_Y$modules$goodGenes)

ME_good <- sum(stage2results_Y$modules$MEsOK) == length(stage2results_Y$modules$MEsOK)
table(stage2results_Y$modules$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename(Module = Var1, Size = Freq) %>%
  dplyr::mutate(Module_color = labels2colors(as.numeric(as.character(Module)))) -> module_size



module_size %>% 
  ggplot(aes(x = Module, y = Size, fill = Module)) +
  geom_col(color =  "#000000") +
  ggtitle("Number of genes in each module") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = setNames(module_size$Module_color,module_size$Module)) +
  geom_text(aes(label = Size),vjust = 0.5, hjust = -0.18, size = 3.5) +
  ylim(0, max(module_size$Size)*1.1) +
  theme(plot.margin = margin(2, 2, 2, 2, "pt")) +
  coord_flip()-> module_size_barplot

module_size_barplot + theme_bw()


table(stage2results_Y$modules$colors) %>% as.data.frame() -> res
res$`Module color` <- WGCNA::labels2colors(as.numeric(as.character(res$Var1)))
res <- res[, c(1,3,2)]
colnames(res) <- c("Module", "Module color", "Number of ASVs")

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_Y$modules$dendrograms)){
  plotDendroAndColors(stage2results_Y$modules$dendrograms[[i]], merged_colors[stage2results_Y$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_Y$modules$blockGenes[[i]]),
                                    " OTUs"))
}
MEs <- stage2results_Y$modules$MEs

# Module correlation to other modules
MEs_R <- bicor(MEs, MEs, maxPOutliers = 0.05)

idx.r <- which(rownames(MEs_R) == "ME0")
idx.c <- which(colnames(MEs_R) == "ME0")

MEs_R_noME0 <- MEs_R[-idx.r, -idx.c]


MEs_R_noME0[upper.tri(MEs_R_noME0)] %>% 
  as.data.frame() %>% 
  dplyr::rename("correlation" = ".") %>% 
  ggplot(aes(x=correlation)) + 
  geom_histogram(bins = 20) +
  #geom_density() + 
  xlim(-1, 1) +
  ggtitle(paste0(prefix,"ME correlation\n w/o ",prefix ,"ME0")) -> MEs_R_density

pheatmap::pheatmap(MEs_R_noME0, color = colorRampPalette(c("Blue", "White", "Red"))(100),
                   silent = T, 
                   breaks = seq(-1,1,length.out = 101),
                   treeheight_row = 5, 
                   treeheight_col = 5,
                   main = paste0(prefix,"ME correlation heatmap w/o ",prefix ,"ME0"),
                   labels_row = paste0(prefix, rownames(MEs_R)),
                   labels_col = paste0(prefix, colnames(MEs_R))) -> MEs_R_Corr

cowplot::plot_grid(MEs_R_density, MEs_R_Corr$gtable, labels = c("D", "E"), label_size = plot_labeling_size, rel_widths = c(0.6, 1)) -> density_eigen

density_eigen

all(rownames(omics_data) == rownames(MEs))
dim(omics_data) %>% paste0(c(" samples", " OTUs"))
kME <- bicor(omics_data, MEs, maxPOutliers = 0.05)
dim(kME) %>% paste0(c(" OTUs", " modules"))

intra_cor <- c()
for (i in 1:ncol(omics_data)) {
  m <- stage2results_Y$modules$colors[i]
  intra_cor[i] <- kME[i, paste0("ME", m)]
  if(m != 0){
    intra_cor[i] <- kME[i, paste0("ME", m)]
  } else{
    intra_cor[i] <- NA
  }
  
}

idx <- which(is.na(intra_cor))
intra_cor <- intra_cor[-idx]

plot(density(intra_cor), main = "Correlations with module-eigenASVs (within module correlation)\nNo ME0", xlim = c(-1,1))



# Corr within modules
corr_within_module <- function(omics_data, modules, module_x = 1){
  idx.omics_data <- which(modules$colors == module_x)
  idx.me <- which(colnames(modules$MEs) == paste0("ME",module_x))
  kME_x <- bicor(omics_data[,idx.omics_data], modules$MEs[,idx.me], maxPOutliers = 0.05)
  kME_x
}

ggplot.list <- list()

for(m in colnames(stage2results_Y$modules$MEs))
  {
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(omics_data = omics_data, modules = stage2results_Y$modules, module_x = h))) %>% 
    ggplot() + 
    #geom_density(aes(x = x), fill = labels2colors(h), color = "black", alpha = 0.5) + 
    geom_histogram(aes(x), fill = labels2colors(h), color = "black", alpha = 0.5, bins = 20) + 
    xlim(-1, 1) +
    xlab("ASV correlation")+
    ggtitle(paste0(prefix,m)) -> da_plot
  
  ggplot.list[[m]] <- da_plot
}

ggplot.list <- ggplot.list[ggplot.list %>% names() %>% sub("ME", "", .) %>% as.numeric() %>% order()]

cowplot::plot_grid(plotlist = ggplot.list, ncol = 6) -> density_all_plot
density_all_plot


# comine to one plot
cowplot::plot_grid(si_mc_plot , density_eigen, ncol = 1, rel_heights = c(0.8,1)) -> part_1


cowplot::plot_grid(part_1, module_size_barplot, labels = c("", "C"), label_size = plot_labeling_size, rel_widths = c(1,0.5)) -> part_2


cowplot::plot_grid(part_2, density_all_plot, ncol = 1, rel_heights = c(0.8,1), labels = c("", "F"), label_size = plot_labeling_size)




stage2results_Y$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("OTU_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubOTUs

library(phyloseq)
otu_taxonomy <- data.frame(tax_table(psdata.p))
dplyr::left_join(hubOTUs, 
                 (otu_taxonomy %>%
                    tibble::rownames_to_column(var = "OTU_name")), 
                 by = "OTU_name") -> hubOTUs_Before

hubOTUs_Before$ModuleOTU <- paste0("ME" ,hubOTUs_Before$Module)




degrees <- intramodularConnectivity.fromExpr(omics_data, colors = modules.omics.Y$colors, power = st,
                                             networkType = "signed", distFnc = chosen_parameter_set$assoc_measure)

degrees$OTU_name <- colnames(omics_data)
degrees$Module <- modules.omics.Y$colors
degrees <- degrees[,c(6,5,1:4)]


dplyr::left_join(degrees, 
                 (otu_taxonomy %>%
                    tibble::rownames_to_column(var = "OTU_name")), 
                 by = "OTU_name") -> degrees


degrees %>%
  dplyr::filter(Module == "8")



moduleLabels <- modules.omics.Y$colors
moduleColors <- labels2colors(modules.omics.Y$colors)
table(moduleColors)
table(moduleLabels)


# Further downstream analysis
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html

module_eigengenes <- modules.omics.Y$MEs

all.equal(sample_info$common, rownames(module_eigengenes))
sample_info %>% data.frame() %>%
  mutate(samplingTime = factor(samplingTime)) %>% 
  mutate(New_Diet = factor(New_Diet)) %>%
  mutate(Day = ifelse(samplingTime == "T1", paste0("T1_", New_Diet),
                      ifelse(samplingTime == "T2", paste0("T2_", New_Diet), paste0("T3_", New_Diet)))) -> meta
head(meta)


des_mat <- model.matrix(~ meta$Day)
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)

stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

stats_df


#Let’s make plot of module 7
module_19_df <- module_eigengenes %>%
  tibble::rownames_to_column("common") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(meta %>%
                      dplyr::select(common, Day, samplingTime, New_Diet),
                    by = c("common" = "common"))


ggplot(module_19_df, aes(x = New_Diet, y = ME18, color = New_Diet)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic() + facet_wrap("samplingTime")


ggplot(module_19_df, aes(x = Day, y = ME18, color = samplingTime)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


#What genes are a part of module 7?
gene_module_key <- tibble::enframe(modules.omics.Y$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
table(gene_module_key$module)

gene_module_key %>%
  dplyr::filter(module == "ME18")



my_comparisons <- list( c("CTR", "MC1"), c("CTR", "MC2"), c("CTR", "MN3"), 
                        c("MC1", "MC2"), c("MC1", "MN3"), c("MC2", "MN3") )
ggboxplot(module_19_df, x = "New_Diet", y = "ME24",
          color = "samplingTime", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.5) + #geom_jitter(aes(colour = New_Diet), size = 3, alpha = 0.8) + 
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  + facet_grid(~samplingTime) 







#lefse_analysis_cleaned_prevalence
#https://rpubs.com/mohsen/lefse_analysis_cleaned_prevalence_phyloseq
library(tidyverse)
library(phyloseq)
library(microbiomeMarker)
library(microbial)



#T1
phy_MC1<-subset_samples(psdata, samplingTime %in% c("T1") & New_Diet %in% c("CTR", "MC1"))
phy_MC1 <- prune_taxa(taxa_sums(phy_MC1) >0, phy_MC1)
table(sample_data(phy_MC1)$New_Diet)


lef_out1<-run_lefse(phy_MC1, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                   kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")



#plot_ef_bar(lef_out)

a1 <- plot_abundance(lef_out1, group = "New_Diet")
table(marker_table(lef_out1)$enrich_group)
#rabuplot(phy_MC1, predictor= "New_Diet",N_taxa = 25)

#marker_1 <- microbiomeMarker::abundances(lef_out1)
data.frame(marker_table(lef_out1)) %>%
  filter(enrich_group != "CTR") %>%
  select(feature) %>%
  `rownames<-`( NULL )

phy_MC2<-subset_samples(psdata, samplingTime %in% c("T1") & New_Diet %in% c("CTR", "MC2"))
phy_MC2 <- prune_taxa(taxa_sums(phy_MC2) >0, phy_MC2)
table(sample_data(phy_MC2)$New_Diet)
lef_out2<-run_lefse(phy_MC2, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                   kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")

#plot_ef_bar(lef_out)
a2 <- plot_abundance(lef_out2, group = "New_Diet")
table(marker_table(lef_out2)$enrich_group)


phy_MN3<-subset_samples(psdata, samplingTime %in% c("T1") & New_Diet %in% c("CTR", "MN3"))
phy_MN3 <- prune_taxa(taxa_sums(phy_MN3) >0, phy_MN3)
table(phyloseq::sample_data(phy_MN3)$New_Diet)
lef_out3<-microbiomeMarker::run_lefse(phy_MN3, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                     kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")
#plot_ef_bar(lef_out)
plot_abundance(lef_out, group = "New_Diet")
table(marker_table(lef_out)$enrich_group)

plot_grid(a1, a2, labels = "AUTO", ncol = 1)



#Code for all time points and diets
output_table <- data.frame()

for (i in c("T1","T2","T3")){
  for (j in c("MC1","MC2", "MN3")){
    phy<-subset_samples(psdata, samplingTime %in% c(i) & New_Diet %in% c("CTR", j))
    phy <- prune_taxa(taxa_sums(phy) >0, phy)
    lef_out<-run_lefse(phy, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                       kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")
    assign(paste0("p", i, j), plot_abundance(lef_out, group = "New_Diet"))
    output_table[i,j]<-sum(marker_table(lef_out)$enrich_group == j)
  }
}
colnames(output_table)<-c("MC1", "MC2", "MN3")
rownames(output_table) <- c("T1","T2","T3")
output_table

plot_grid(pT1MC1, pT1MC2, pT1MN3, pT2MC1, pT2MC2, pT2MN3, pT3MC1, pT3MC2, p3MN3, labels = "AUTO", ncol = 1)
plot_grid(pT1MC1, pT1MC2,  pT2MC1, pT2MC2,  pT3MC1, pT3MC2, labels = "AUTO", ncol = 1, align = "hv")





#Code for all time points and diets but ignoring "T1MN3"
for (i in c("T1","T2","T3")){
  for (j in c("MC1","MC2", "MN3")){
    if (i == "T1" & j == "MN3"){
      next
    } else {
      phy<-subset_samples(psdata, samplingTime %in% c(i) & New_Diet %in% c("CTR", j))
      phy <- prune_taxa(taxa_sums(phy) >0, phy)
      lef_out<-run_lefse(phy, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                         kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")
      assign(paste0("metagenomics", i, j), plot_abundance(lef_out, group = "New_Diet"))
      output_table[i,j]<-sum(marker_table(lef_out)$enrich_group == j)
    }
  }
}
colnames(output_table)<-c("MC1", "MC2", "MN3")
rownames(output_table) <- c("T1","T2","T3")
t(output_table)

plot_grid(metagenomicsT1MC1, metagenomicsT1MC2, 
          metagenomicsT2MC1, metagenomicsT2MC2, metagenomicsT2MN3, 
          metagenomicsT3MC1, metagenomicsT3MC2, metagenomicsT3MN3, 
          labels = "AUTO", ncol = 1, align = "hv")




#Code for all time points and diets but ignoring "T1MN3"
# Control
for (i in c("T1","T2","T3")){
  for (j in c("MC1","MC2", "MN3")){
    if (i == "T1" & j == "MN3"){
      next
    } else {
      phy<-subset_samples(psdata, samplingTime %in% c(i) & New_Diet %in% c("CTR", j))
      phy <- prune_taxa(taxa_sums(phy) >0, phy)
      lef_out<-run_lefse(phy, group = "New_Diet", taxa_rank = "Genus", transform = "log10",
                         kw_cutoff = 0.05, lda_cutoff = 2, strict = "2")
      assign(paste0("metagenomics", i, j), plot_abundance(lef_out, group = "New_Diet"))
      output_table[i,j]<-sum(marker_table(lef_out)$enrich_group == "CTR")
    }
  }
}
colnames(output_table)<-c("MC1", "MC2", "MN3")
rownames(output_table) <- c("T1","T2","T3")
t(output_table)



