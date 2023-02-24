#WGCNA Host analysis

library(stringr)
library(dplyr)
library(DESeq2)

# Set seed as a precaution for reproducibility as some methods are non-deterministic.
set.seed(13118)
Run_analysis <- TRUE    # if FALSE it tries to load data instead of running the module creation step.
prefix <- "h"
save_TOM <- FALSE       # Can get rather big
pam_stage <- FALSE      # Partitioning around medoids, tries to put more genes into modules where it is not directly clear from the dendrogram cutting.
set_verbose = 1         # How much detail from the WGCNA functions? A higher number means more detail.
omics_type = "rna"      # Just a name in the form of a string
take_average = F        # Take the average of each sample group
assoc_measure = "bicor"
max_block_size = 10000      #  13000 ~ 2-2.5 h |
applied_norm = "TMM"        #  TPM or (TMM, choose raw counts below)
applied_transf = "log2"     #  log2 or CLR
applied_filter = "sd"       #  sd or mad
pcCorrection <- F
if(pcCorrection){
  estimate_n.pc = T
  if(!estimate_n.pc){
    number_of_pcs_to_remove = 1 # Does not matter when pcCorrection is FALSE
  }
}
theme_set(theme_bw())


#host_gut_mapping_file[host_gut_mapping_file$sampleName == "T1_mn3_tk2_HGh07", "New_Diet"] <- "mc2"
omics_data_host <- host_gut_rawCounts_total1
sample_info<- host_gut_mapping_file

#sample_info$New_Diet <- tolower(sample_info$New_Diet)
#Rename the sample names based on new diet
sample_info$ID_New <- paste(sample_info$Time, "_",
                            sample_info$New_Diet, "_", 
                            sample_info$Tank_number, "_", "0",
                            sample_info$Sample_Number, 
                                      sep = "")
head_change <- subset(sample_info, select=c("common" ,"ID_New"))
row.names(sample_info) <- sample_info$ID_New



#omics_data_host <- t(omics_data_host)
#omics_data_host$names <- row.names(omics_data_host)


omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% sample_info$ID_New]
omics_data_host <- omics_data_host1
all(colnames(omics_data_host) == sample_info$ID_New)


omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% sample_info$ID_New]
omics_data_host <- omics_data_host1

# rename colnames of omics_data_host using new ID
library(dplyr)
library(tibble)


df1 <- subset(sample_info, select= c("common", "ID_New"))
omics_data_host %>% 
  rename_with(~deframe(df1)[.x], .cols = df1$ID_New) %>% 
  select(any_of(df1$ID_New)) -> omics_data_host_new
omics_data_host <- omics_data_host_new

all(colnames(omics_data_host) == sample_info$ID_New)
reorder_idx <- match(sample_info$ID_New, colnames(omics_data_host))
omics_data_host <- omics_data_host[ , reorder_idx]
all(colnames(omics_data_host) %in% sample_info$ID_New)
all(colnames(omics_data_host) == sample_info$ID_New)
table( sample_info$ID_New %in% colnames(omics_data_host))


#test <- subset(sample_info, !(colnames(omics_data_host) %in% sample_info$ID_New ))
#sample_info <- sample_info[sample_info$ID_New %in% colnames(omics_data_host), ]




#rownames(sample_info) <- sample_info$ID_New
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
df <- round(omics_data_host) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50)

all(colnames(df) == rownames(sample_info))


dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = sample_info, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)


dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data
rm(dds_norm)
#normalized_counts[1:2, 1:2]
omics_data_host <- normalized_counts 
rm(normalized_counts)

dim(omics_data_host) %>% paste(c("Samples", "Genes"))


powers <- c(c(1:10), seq(from = 12, to=20, by=2))

library(WGCNA)
allowWGCNAThreads()

sft <- pickSoftThreshold(omics_data_host, 
                         powerVector = powers, 
                         verbose = set_verbose, 
                         networkType = "signed",
                         corFn= chosen_parameter_set$assoc_measure)

library(ggplot2)

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
  ylim(-1,1) + theme_bw() +
  geom_segment(aes(x = st, y = 0.25, xend = st, yend = sfApprox[idx]-0.05), 
               arrow = arrow(length = unit(0.2,"cm")), 
               size = 0.5)-> scale_independence_plot 



data.frame(Indices = sft$fitIndices[,1],
           meanApprox = sft$fitIndices[,5]) %>% 
  ggplot() + 
  geom_line(aes(x = Indices, y = meanApprox), color = "red", alpha = 0.1, size = 2.5) +
  geom_text(mapping = aes(x = Indices, y = meanApprox, label = Indices), color = "red", size = 4) +
  xlab("Soft Threshold (power)") +
  ylab("Mean Connectivity") +theme_bw() +
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
st

modules.omics_X <- blockwiseModules(omics_data_host, 
                                  power = st, 
                                  networkType = "signed", 
                                  TOMType = "signed",
                                  corType = assoc_measure,
                                  #maxPOutliers = 0.05,
                                  #deepSplit = 4, # Default 2
                                  minModuleSize = 20,           # 30
                                  #minCoreKME = 0.5,            # Default 0.5
                                  #minCoreKMESize = 2,          # Default minModuleSize/3,
                                  #minKMEtoStay = 0.5,          # Default 0.3
                                  #reassignThreshold = 0,       # Default 1e-6
                                  #mergeCutHeight = 0.2,        # Default 0.15
                                  #pamStage = pam_stage, 
                                  #pamRespectsDendro = TRUE,
                                  #replaceMissingAdjacencies = TRUE,
                                  nThreads = 11,
                                  numericLabels = TRUE,
                                  saveTOMs = save_TOM,
                                  saveTOMFileBase = "HOST_TOM",
                                  verbose = 3,
                                  maxBlockSize=8000)



rownames(modules.omics_X$MEs) <- rownames(omics_data_host)
names(modules.omics_X$colors) <- colnames(omics_data_host)
names(modules.omics_X$unmergedColors) <- colnames(omics_data_host)

hubs <- chooseTopHubInEachModule(omics_data_host, modules.omics_X$colors, omitColors = 0)
hubs <- chooseTopHubInEachModule(omics_data_host, modules.omics_X$colors, power = st, omitColors = "0")

stage2results_X <- list(modules = modules.omics_X, 
                      hubs = hubs)

stage2results_X
# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_X$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_X$modules$goodSamples) == length(stage2results_X$modules$goodSamples)
genes_good <- sum(stage2results_X$modules$goodGenes) == length(stage2results_X$modules$goodGenes)

ME_good <- sum(stage2results_X$modules$MEsOK) == length(stage2results_X$modules$MEsOK)
table(stage2results_X$modules$colors) %>% 
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

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_X$modules$dendrograms)){
  plotDendroAndColors(stage2results_X$modules$dendrograms[[i]], merged_colors[stage2results_X$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_X$modules$blockGenes[[i]]),
                                    " genes"))
}


MEs <- stage2results_X$modules$MEs


# Module correlation to other modules
MEs_R <- bicor(MEs, MEs, maxPOutliers = 0.05)

idx.r <- which(rownames(MEs_R) == "ME0")
idx.c <- which(colnames(MEs_R) == "ME0")
MEs_R_noME0 <- MEs_R[-idx.r, -idx.c]

MEs_R_noME0[upper.tri(MEs_R_noME0)] %>% 
  as.data.frame() %>% 
  dplyr::rename("correlation" = ".") %>% 
  ggplot(aes(correlation)) + 
  geom_histogram(bins = 20) +
  #geom_density() + 
  xlim(-1, 1) +
  ggtitle(paste0(prefix,"ME correlation\n w/o ",prefix ,"ME0")) -> MEs_R_density

pheatmap::pheatmap(MEs_R_noME0, color = colorRampPalette(c("Blue", "White", "Red"))(100),
                   silent = T, 
                   breaks = seq(-1,1,length.out = 101),
                   treeheight_row = 10, 
                   treeheight_col = 10,
                   main = paste0(prefix,"ME correlation heatmap w/o ", prefix ,"ME0"),
                   labels_row = rep("", length(rownames(MEs_R))), #paste0(prefix, rownames(MEs_R)),
                   labels_col = rep("", length(colnames(MEs_R))) #paste0(prefix, colnames(MEs_R))
) -> MEs_R_Corr

cowplot::plot_grid(MEs_R_density, MEs_R_Corr$gtable, labels = c("D", "E"), label_size = plot_labeling_size, rel_widths = c(0.6, 1)) -> density_eigen

density_eigen

all(rownames(omics_data_host) == rownames(MEs))
dim(omics_data_host) %>% paste0(c(" samples", " genes"))
kME <- bicor(omics_data_host, MEs, maxPOutliers = 0.05)
dim(kME) %>% paste0(c(" genes", " modules"))



# Corr within modules
corr_within_module <- function(omics_data_host, modules, module_x = 1){
  idx.omics_data_host <- which(modules$colors == module_x)
  idx.me <- which(colnames(modules$MEs) == paste0("ME",module_x))
  kME_x <- bicor(omics_data_host[,idx.omics_data_host], modules$MEs[,idx.me], maxPOutliers = 0.05)
  kME_x
}

ggplot.list <- list()

modules <- 
  colnames(stage2results_X$modules$MEs)[colnames(stage2results_X$modules$MEs) %>% sub("ME", "", .) %>% as.numeric() %>% order()]

for(m in colnames(stage2results_X$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(omics_data_host = omics_data_host, modules = stage2results_X$modules, module_x = h))) %>% 
    ggplot() + 
    #geom_density(aes(x = x), fill = labels2colors(h), color = "black", alpha = 0.5) +
    geom_histogram(aes(x), fill = labels2colors(h), color = "black", alpha = 0.5, bins = 20) + 
    xlim(-1, 1) +
    xlab("Gene correlation") +
    ggtitle(paste0(prefix,m)) -> da_plot
  
  ggplot.list[[m]] <- da_plot
}

ggplot.list <- ggplot.list[ggplot.list %>% names() %>% sub("ME", "", .) %>% as.numeric() %>% order()]

cowplot::plot_grid(plotlist = ggplot.list) -> density_all_plot

#density_all_plot



# comine to one plot
cowplot::plot_grid(si_mc_plot , density_eigen, ncol = 1, rel_heights = c(0.8,1)) -> part_1


cowplot::plot_grid(part_1, module_size_barplot, labels = c("", "C"), label_size = plot_labeling_size, rel_widths = c(1,0.5)) -> part_2


cowplot::plot_grid(part_2, density_all_plot, ncol = 1, rel_heights = c(0.8,1), labels = c("", "F"), label_size = plot_labeling_size) -> part_3

part_3 # 1300 * 1300

#Hub genes
stage2results_X$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("gene_id" = ".") %>% 
  tibble::rownames_to_column(var = "Module") -> top_genes

top_genes



#



# Further downstream analysis
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
#modules.omics_X.backup <- modules.omics_X 
module_eigengenes <- modules.omics_X$MEs

all.equal(sample_info$ID_New, rownames(module_eigengenes))
sample_info %>% 
  mutate(Time = factor(Time)) %>% 
  mutate(New_Diet = factor(New_Diet)) %>%
  mutate(Day = ifelse(Time == "T1", paste0("T1_", New_Diet),
                      ifelse(Time == "T2", paste0("T2_", New_Diet), paste0("T3_", New_Diet)))) -> meta
head(meta)


des_mat <- model.matrix(~ meta$Day)
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
fit <- limma::eBayes(fit)

stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

stats_df


#Let’s make plot of module 7
module_19_df <- module_eigengenes %>%
  tibble::rownames_to_column("ID_New") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(meta %>%
                      dplyr::select(ID_New, Day, Time, New_Diet),
                    by = c("ID_New" = "ID_New"))


ggplot(module_19_df, aes(x = Time, y = ME23, color = Day)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic() + facet_wrap("New_Diet")


ggplot(module_19_df, aes(x = Day, y = ME23, color = Time)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


#What genes are a part of module 7?
gene_module_key <- tibble::enframe(modules.omics_X$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME1")





#GO term and KEGG pathway enrichment of Atlantic salmon genes using the clusterProfiler package
library(AnnotationHub)
library(AnnotationDbi)
library(readr)
library(dplyr)
library(clusterProfiler)
library(GenomicFeatures)
library(GO.db)

#https://gitlab.com/sandve-lab/salmon-go-and-kegg-enrichment/-/blob/master/Ssal_GO_and_KEGG_ClusterProfiler.Rmd
# get OrgDb from AnnotationHub
ah <- AnnotationHub()
# list OrgDb's available for Salmo salar
SsalOrg <- subset(ah, ah$species == "Salmo salar" & ah$rdataclass=="OrgDb")
SsalOrg <- ah[["AH96017"]]


# Extract each modules information and NCBI gene ID associated with them 
moduleLabels <- modules.omics_X$colors
moduleColors <- labels2colors(modules.omics_X$colors)
table(moduleColors)
table(moduleLabels)

table(moduleColors, moduleLabels)

MEs <- modules.omics_X$MEs
geneTree <- modules.omics_X$dendrograms[[1]]


#https://gitlab.com/garethgillard/ssalv3annotation
annot1 <- read.csv("/Users/shashankgupta/Desktop/ImprovAFish/ImprovAFish_Final/Salmo_salar-GCA_905237065.2_gene_annotations.csv", header = TRUE)

probes <- colnames(omics_data_host) 
probes2annot <- match(probes, annot1$gene_id)
allLLIDs <- annot1$v2.gene_id.NCBI[probes2annot]

intModules <- unique(moduleColors)

setwd("/Users/shashankgupta/Desktop/ImprovAFish/ImprovAFish/resGO/")
for (module in intModules) 
{
  # Select module probes 
  modGenes = (moduleColors==module) 
  # Get their entrez ID codes 
  modLLIDs = allLLIDs[modGenes]
  # Write them into a file 
  fileName = paste("LocusLinkIDs-", module, ".txt", sep=""); 
  write.table(as.character(modLLIDs), file = fileName, row.names = FALSE, col.names = FALSE) 
}

fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.character(allLLIDs), file = fileName, row.names = FALSE, col.names = FALSE)

#table(moduleColors)
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

resGO_dot_plots <- ls(pattern="resGO_dotplot_*")
resKEGG_dot_plots <- ls(pattern="resKEGG_dotplot_*")












# save all the plots
library(cowplot)
library(gridExtra)

# PNG
plots <- ls(pattern = "^resKEGG_dotplot_.*")
for (i in 1:length(plots)) {
  tryCatch({
    png(paste0(plots[i], ".png"))
    print(get(plots[i]))
    dev.off()
  }, error = function(e) {
    print(paste("Error printing", plots[i], ":", e))
    file.remove(paste0(plots[i], ".png"))
  })
}



# PDF
dir.create("plots")
setwd("plots/")
library(plotflow)
library(pdftools)
# Get all the objects in the environment that start with the string "resKEGG_dotplot_"
plots <- ls(pattern = "^resKEGG_dotplot_.*")
for (i in 1:length(plots)) {
  tryCatch({
    pdf(paste0(plots[i], ".pdf"))
    print(get(plots[i]))
    dev.off()
  }, error = function(e) {
    print(paste("Error printing", plots[i], ":", e))
    file.remove(paste0(plots[i], ".pdf"))
  })
}

pdf_list <- list() 
path <- getwd() 
pdf_list <- list.files(path = path, pattern = "*.pdf", full.names = TRUE)
pdftools::pdf_combine(input = pdf_list,
                      output = "merge_resKEGG.pdf")

# Remove the individual PDF files
pdf_files_to_remove <- setdiff(pdf_list, "merge_resKEGG.pdf")
file.remove(pdf_files_to_remove)














#Deseq2 analysis
library(DESeq2)
library(BiocParallel)

dds$group <- factor(paste0(dds$Time, dds$New_Diet))
design(dds) <- ~ group-1

dds <- DESeq(dds, parallel = T)
resultsNames(dds)

annot_function <- read_tsv("../../../Salmo_salar-GCA_905237065.2_gene_annotations.tsv")
annot_function <- annot_function[, c("gene_id", "v2.gene_id.NCBI","v2.gene_name.ensembl", "v2.product" )]

# T1
res <- results(dds,
               contrast = list("groupT1mc1","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj)

res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements

head(res_tbl, n = 25)



res <- results(dds,
               contrast = list("groupT1mc2","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT1mn3","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements


# T2
res <- results(dds,
               contrast = list("groupT2mc1","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



res <- results(dds,
               contrast = list("groupT2mc2","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT2mn3","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements


# T3
res <- results(dds,
               contrast = list("groupT3mc1","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



res <- results(dds,
               contrast = list("groupT3mc2","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT3mn3","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(v2.product)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



# loop
# Create the dataframe
df <- data.frame(Group=character(), a1=character(), a2=character(), a3=character())

# Loop through each group and store the results in the dataframe
for(treatment in c("T1", "T2", "T3")){
  for(group in c("mc1","mc2","mn3")){
    sample_name <- paste0("group",treatment,group)
    ctr_sample <- paste0("group",treatment,"ctr")
    
    # Calculate results
    res <- results(dds, contrast = list(sample_name,ctr_sample))
    
    # Filter the results and count the number of positive elements
    res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
      filter(padj <0.05)%>%
      arrange(padj)
    a <- length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0])
    
    # Store the results in the dataframe
    df <- rbind(df, c(sample_name,a))
  }
}
names(df) <- c("Group", "a1", "a2", "a3")
df












