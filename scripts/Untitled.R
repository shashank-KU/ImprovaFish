

library(gtsummary)
library(dplyr)

df1 <- data.frame(sample_info1)
df2<- data.frame(sample_info)
rownames(df2) <- df2$ID_New
df1$ID_New <- rownames(df1)



sample_info2<- merge(df1, df2, by = "ID_New")


trial2 <- sample_info2 %>% select(BW, Length, Gut_light, Anus_bleed, Gut_bleed, New_Diet)


table2 <- 
  tbl_summary(
    trial2,
    by = New_Diet, # split table by group
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p() %>% # sample_info for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 

table2
head(sample_info)


library(ggpubr)

sample_info$New_Diet <- factor(sample_info$New_Diet, levels=c( 'ctr', 'mc1', 'mc2', 'mn3'))
levels(sample_info$New_Diet)

sample_info$Time <- factor(sample_info$Time, levels=c('T1', 'T2', 'T3'))
levels(sample_info$Time)

table(sample_info$Group)

sample_info$Group <- factor(sample_info$Group, levels=c('ctr_T1', 'ctr_T2', 'ctr_T3', 
                                          'mc1_T1', 'mc1_T2', 'mc1_T3',
                                          'mc2_T1', 'mc2_T2', 'mc2_T3',
                                          'mn3_T1', 'mn3_T2', 'mn3_T3'))
levels(sample_info$Group)

#compare_means(Shannon ~ race, data = sample_info, method = "kruskal.sample_info")
my_comparisons <- list( c("ctr", "mc1"), c("ctr", "mc2"),
                        c("ctr", "mn3"), c("mc1", "mn2"),
                        c("mc1", "mn3"), c("mc2", "mn3"))

my_comparisons <- list( c("ctr_T1", "mc1_T1"), c("ctr_T1", "mc2_T1"),c("ctr_T1", "mn3_T1"),
                        c("ctr_T2", "mc1_T2"), c("ctr_T2", "mc2_T2"),c("ctr_T2", "mn3_T2"),
                        c("ctr_T3", "mc1_T3"), c("ctr_T3", "mc2_T3"),c("ctr_T3", "mn3_T3"))
library(cowplot)

p9 <- ggboxplot(sample_info, x = "Group", y = "BW",
                color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 600) + 
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value


p10 <- ggboxplot(sample_info, x = "Group", y = "Length",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 45) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value

plot_grid(p9, p10, labels = c('A', 'B' ), label_size = 12, ncol = 1)


p11 <- ggboxplot(sample_info, x = "Group", y = "CSI",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.45) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  


p12 <- ggboxplot(sample_info, x = "Group", y = "HSI",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6.5) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  

plot_grid(p11, p12, labels = c('A', 'B' ), label_size = 12, ncol = 1)


ggboxplot(sample_info, x = "Group", y = "Length",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) +facet_wrap("Sex")+
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value






ggboxplot(sample_info, x = "Group", y = "Shannon",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 7) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value





# a diversity correlation
ITS_sample <- sample_data(Final.ITS)
RNA_sample <- sample_data(Final.RNA)

ITS_alpha <- data.frame(ITS_sample[,which(colnames(ITS_sample) %in% c("abcno","ObservedRichness_mean","SDI_mean"))])
RNA_alpha <- data.frame(RNA_sample[,which(colnames(RNA_sample) %in% c("abcno","ObservedRichness_mean","SDI_mean"))])
colnames(ITS_alpha)[c(2:3)] <- c("ITS_OR_mean","ITS_SDI_mean")
colnames(RNA_alpha)[c(2:3)] <- c("RNA_OR_mean","RNA_SDI_mean")
ITS_RNA_alpha <- merge(ITS_alpha,RNA_alpha,by="abcno")

ggscatter(sample_info, x = "BW", y = "Richness",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "BW", ylab = "RNA_Richness")

ggscatter(sample_info, x = "Length", y = "Richness",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Length", ylab = "RNA_Richness")


p13 <- ggscatter(sample_info, x = "BW", y = "Shannon",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "BW", ylab = "Shannon Diversity Index")

p14 <- ggscatter(sample_info, x = "Length", y = "Shannon",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Length", ylab = "Shannon Diversity Index")



plot_grid(p13, p14, labels = c('A', 'B' ), label_size = 12, ncol = 1)






ah <- AnnotationHub()
SsalOrg <- subset(ah, ah$species == "Salmo salar" & ah$rdataclass=="OrgDb")
SsalOrg <- ah[["AH96017"]]

library(gridExtra)
library(clusterProfiler)
library("stringr")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("tibble")
library("AnnotationHub")
library("AnnotationDbi")
library("readr")
library("dplyr")
library("clusterProfiler")
library("GenomicFeatures")
library("GO.db")
library("BiocParallel")
library("DT")
library("WGCNA")
library("cowplot")
library("gridExtra")
library("EnhancedVolcano")
library("plyr")
library("png")
library("grid")
library("gridExtra")


library(clusterProfiler)

# Read in the module colors and MEs from your WGCNA analysis
moduleLabels <- modules.omics_X_Liver$colors
moduleColors <- labels2colors(modules.omics_X_Liver$colors)
MEs <- modules.omics_X_Liver$MEs

# Read in the gene annotations and create a vector of all gene IDs
annot1 <- read.csv("/Users/shashankgupta/Desktop/ImprovAFish/github/Input/Salmo_salar-GCA_905237065.2_gene_annotations.csv", header = TRUE)
probes <- colnames(omics_data_host) 
probes2annot <- match(probes, annot1$gene_id)
allLLIDs <- annot1$v2.gene_id.NCBI[probes2annot]

# Create a list of data frames, one for each module, containing the gene IDs in that module
intModules <- unique(moduleColors)
geneLists <- lapply(intModules, function(module) {
  modGenes <- (moduleColors == module) 
  modLLIDs <- allLLIDs[modGenes]
  data.frame(gene = modLLIDs, module = module)
})

# Combine all gene lists into a single data frame
geneList <- do.call(rbind, geneLists)

# Perform GO enrichment analysis on each module
resGO <- compareCluster(geneClusters = geneList, fun = "enrichGO", pvalueCutoff = 0.05, qvalueCutoff = 0.05, OrgDb = "SsalOrg", ont = "BP")

# Create a single dotplot for all enriched GO terms across all modules
resGO_dotplot <- dotplot(resGO, showCategory = 20, title = "GO Enrichment Analysis")

# Perform KEGG pathway enrichment analysis on each module
resKEGG <- compareCluster(geneClusters = geneList, fun = "enrichKEGG", pvalueCutoff = 0.05, qvalueCutoff = 0.05, organism = "sasa")

# Create a single dotplot for all enriched KEGG pathways across all modules
resKEGG_dotplot <- dotplot(resKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

# Combine the dotplots for GO and KEGG into a single plot with multiple panels, one for each module
all_res_dotplot <- list()
for (i in intModules) {
  all_res_dotplot[[i]] <- grid.arrange(get(paste0("resGO_dotplot_", i)))
}

# Combine all panels into a single plot with multiple columns, one for each module
grid.arrange(grobs = all_res_dotplot, ncol = 3)



#
library(clusterProfiler)
# Making a list of module names
MEList <- moduleEigengenes(modules.omics_X_Liver, colors = dynamicColors)
MEs
modNames <- substring(names(MEs), 3)

# Correlating each genes expression profile with the module eigengenes in order to create module gene sets
geneModuleMembership <- as.data.frame(cor(omics_data_host, MEs, use = "p"))
# "For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile." - WGCNA tutorial

# Iteratively creating a list of module genesets to test. These are in ensembl ids
moduleGeneSets<-lapply(modNames,function(module){
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  rownames(geneModuleMembership[moduleGenes,])
})
names(moduleGeneSets)<-modNames

# Trimming the module gene sets so that the final two digits after the "." are removed
moduleGeneSets.trimmed<-lapply(moduleGeneSets,function(x){
  str_split_fixed(x,"\\.",2)[,1]
})

# Looking up the ENTREZ id for each gene
moduleGeneSets.Entrez<-lapply(moduleGeneSets.trimmed,function(x){
  bitr(x,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Mm.eg.db")$ENTREZID
})





# Extract module labels and colors
moduleLabels <- modules.omics_X_Liver$colors
moduleColors <- labels2colors(modules.omics_X_Liver$colors)

# Create a table of module colors
table(moduleColors)

# Extract module eigengenes and gene tree
MEs <- modules.omics_X_Liver$MEs
geneTree <- modules.omics_X_Liver$dendrograms[[1]]

# Load gene annotations from file
annot1 <- read.csv("/Users/shashankgupta/Desktop/ImprovAFish/github/Input/Salmo_salar-GCA_905237065.2_gene_annotations.csv", header = TRUE)

# Map probe IDs to gene IDs using annotations
probes <- colnames(omics_data_host) 
probes2annot <- match(probes, annot1$gene_id)
allLLIDs <- annot1$v2.gene_id.NCBI[probes2annot]

# Extract unique module colors
intModules <- unique(moduleColors)

# Write LocusLinkIDs for each module to a file
for (module in intModules) 
{
  # Select module probes 
  modGenes = (moduleColors==module) 
  # Get their LocusLink IDs 
  modLLIDs = allLLIDs[modGenes]
  # Write them into a file 
  fileName = paste("LocusLinkIDs-", module, ".txt", sep=""); 
  write.table(as.character(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE) 
}

# Write all LocusLinkIDs to a single file
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.character(allLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)

# Perform gene ontology and pathway enrichment analysis for each module
list_names <- unique(moduleColors)
all.df <- read.table("LocusLinkIDs-all.txt")
all.df$V1 <- as.character(all.df$V1)
for (i in list_names){
  assign(paste(i,"_df", sep=""), read.table(paste0("LocusLinkIDs-",i,".txt")))
  assign(paste("resGO_",i, sep=""), enrichGO(gene = get(paste0(i,"_df"))$V1, ont = "BP", universe = all.df$V1, OrgDb = "SsalOrg", pvalueCutoff  = 0.05,qvalueCutoff  = 0.05))
  assign(paste("resGO_dotplot_",i, sep=""), dotplot(get(paste0("resGO_",i))))
  assign(paste("resKEGG_",i, sep=""), enrichKEGG(gene = get(paste0(i,"_df"))$V1, universe = all.df$V1,  organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05))
  assign(paste("resKEGG_dotplot_",i, sep=""), dotplot(get(paste0("resKEGG_",i))))
}
















######
merge <- modules.omics_X_Liver
moduleColors <- merge$colors
mergedMEs <- merge$MEs
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs <- mergedMEs
MEs <- MEs[,c(1:76)] 
modNames <- substring(names(MEs), 3)

# Correlating each genes expression profile with the module eigengenes in order to create module gene sets
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
# "For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile." - WGCNA tutorial

# Iteratively creating a list of module genesets to test. These are in ensembl ids
moduleGeneSets<-lapply(modNames,function(module){
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  rownames(geneModuleMembership[moduleGenes,])
})
names(moduleGeneSets)<-modNames

# Trimming the module gene sets so that the final two digits after the "." are removed
moduleGeneSets.trimmed<-lapply(moduleGeneSets,function(x){
  str_split_fixed(x,"\\.",2)[,1]
})

# Looking up the ENTREZ id for each gene
moduleGeneSets.Entrez<-lapply(moduleGeneSets.trimmed,function(x){
  bitr(x,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Mm.eg.db")$ENTREZID
})






































