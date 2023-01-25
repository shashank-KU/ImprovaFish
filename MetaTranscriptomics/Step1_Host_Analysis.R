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
annot1 <- read.csv("../../ImprovAFish_Final/Salmo_salar-GCA_905237065.2_gene_annotations.csv", header = TRUE)
probes <- colnames(omics_data_host) 
probes2annot <- match(probes, annot1$gene_id)
allLLIDs <- annot1$v2.gene_id.NCBI[probes2annot]

intModules <- unique(moduleColors)

setwd("resGO/")
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


#read a text file in R for all the unique names present in "moduleColors" and read the file with unique names. e.g.,

setwd("../../ImprovAFish/resGO/")
table(moduleColors)

# in R make a loop-
#   
# black.df <- read.table("LocusLinkIDs-black.txt")
# resGO.black <- enrichGO(gene = black.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
# resGO.black.dotplot <-  dotplot(resGO.black)
# r1 <- as.character(resGO.black$Description)
# resKEGG.black <- enrichKEGG(gene = black.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
# resKEGG.black.dotplot <- dotplot(resKEGG.black, showCategory = 10)
# 
# 
# do it for all the names in list_names
# list_names <- unique(moduleColors)
# "grey"          "blue"          "magenta"       "red"           "turquoise"     "midnightblue" 
#  "cyan"                  "purple"        "pink"          "lightyellow"   "darkgreen"    
#  "brown"         "green"         "darkturquoise" "darkred"       "yellow"        "greenyellow"  
# "lightcyan"     "salmon"        "tan"           "darkgrey"      "grey60"        "royalblue"    
# "lightgreen" 
# 
# In R, make a loop for all the names in list_names and perform all the commands below by changing the names
# black_df <- read.table("LocusLinkIDs-black.txt")
# resGO_black <- enrichGO(gene = black_df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
# resGO_black_dotplot <-  dotplot(resGO_black)
# 
# 
# change the names of black.df, resGO.black, resGO.black.dotplot to list_names
# 
# 
# r1 <- as.character(resGO.black$Description)
# resKEGG.black <- enrichKEGG(gene = black.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
# resKEGG.black.dotplot <- dotplot(resKEGG.black, showCategory = 10)






black.df <- read.table("LocusLinkIDs-black.txt")
blue.df <- read.table("LocusLinkIDs-blue.txt")
brown.df <- read.table("LocusLinkIDs-brown.txt")
cyan.df <- read.table("LocusLinkIDs-cyan.txt")
darkgreen.df <- read.table("LocusLinkIDs-darkgreen.txt")
darkgrey.df <- read.table("LocusLinkIDs-darkgrey.txt")
darkmagenta.df <- read.table("LocusLinkIDs-darkmagenta.txt")
darkolivegreen.df <- read.table("LocusLinkIDs-darkolivegreen.txt")
darkorange.df <- read.table("LocusLinkIDs-darkorange.txt")
darkred.df <- read.table("LocusLinkIDs-darkred.txt")
darkturquoise.df <- read.table("LocusLinkIDs-darkturquoise.txt")
green.df <- read.table("LocusLinkIDs-green.txt")
greenyellow.df <- read.table("LocusLinkIDs-greenyellow.txt")
grey60.df <- read.table("LocusLinkIDs-grey60.txt")
lightcyan.df <- read.table("LocusLinkIDs-lightcyan.txt")
lightgreen.df <- read.table("LocusLinkIDs-lightgreen.txt")
lightyellow.df <- read.table("LocusLinkIDs-lightyellow.txt")
magenta.df <- read.table("LocusLinkIDs-magenta.txt")
#mediumpurple3.df <- read.table("LocusLinkIDs-mediumpurple3.txt")
midnightblue.df <- read.table("LocusLinkIDs-midnightblue.txt")
orange.df <- read.table("LocusLinkIDs-orange.txt")
#orangered4.df <- read.table("LocusLinkIDs-orangered4.txt")
paleturquoise.df <- read.table("LocusLinkIDs-paleturquoise.txt")
pink.df <- read.table("LocusLinkIDs-pink.txt")
#plum1.df <- read.table("LocusLinkIDs-plum1.txt")
purple.df <- read.table("LocusLinkIDs-purple.txt")
red.df <- read.table("LocusLinkIDs-red.txt")
royalblue.df <- read.table("LocusLinkIDs-royalblue.txt")
saddlebrown.df <- read.table("LocusLinkIDs-saddlebrown.txt")
salmon.df <- read.table("LocusLinkIDs-salmon.txt")
#sienna3.df <- read.table("LocusLinkIDs-sienna3.txt")
skyblue.df <- read.table("LocusLinkIDs-skyblue.txt")
#skyblue3.df <- read.table("LocusLinkIDs-skyblue3.txt")
steelblue.df <- read.table("LocusLinkIDs-steelblue.txt")
tan.df <- read.table("LocusLinkIDs-tan.txt")
turquoise.df <- read.table("LocusLinkIDs-turquoise.txt")
violet.df <- read.table("LocusLinkIDs-violet.txt")
white.df <- read.table("LocusLinkIDs-white.txt")
yellow.df <- read.table("LocusLinkIDs-yellow.txt")
#yellowgreen.df <- read.table("LocusLinkIDs-yellowgreen.txt")

all.df <- read.table("LocusLinkIDs-all.txt")
all.df$V1 <- as.character(all.df$V1)
 
resGO.black <- enrichGO(gene = black.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.blue <- enrichGO(gene = blue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.brown <- enrichGO(gene = brown.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.cyan <- enrichGO(gene = cyan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkgreen <- enrichGO(gene = darkgreen.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkgrey <- enrichGO(gene = darkgrey.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkmagenta <- enrichGO(gene = darkmagenta.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkolivegreen <- enrichGO(gene = darkolivegreen.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkorange <- enrichGO(gene = darkorange.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkred <- enrichGO(gene = darkred.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.darkturquoise <- enrichGO(gene = darkturquoise.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.green <- enrichGO(gene = green.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.greenyellow <- enrichGO(gene = greenyellow.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.grey60 <- enrichGO(gene = grey60.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.lightcyan <- enrichGO(gene = lightcyan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.lightgreen <- enrichGO(gene = lightgreen.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.lightyellow <- enrichGO(gene = lightyellow.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.magenta <- enrichGO(gene = magenta.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#resGO.mediumpurple3 <- enrichGO(gene = mediumpurple3.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.midnightblue <- enrichGO(gene = midnightblue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.orange <- enrichGO(gene = orange.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.orangered4 <- enrichGO(gene = orangered4.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.paleturquoise <- enrichGO(gene = paleturquoise.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.pink <- enrichGO(gene = pink.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#resGO.plum1 <- enrichGO(gene = plum1.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.purple <- enrichGO(gene = purple.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.red <- enrichGO(gene = red.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.royalblue <- enrichGO(gene = royalblue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.saddlebrown <- enrichGO(gene = saddlebrown.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.salmon <- enrichGO(gene = salmon.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#resGO.sienna3 <- enrichGO(gene = sienna3.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.skyblue <- enrichGO(gene = skyblue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#resGO.skyblue3 <- enrichGO(gene = skyblue3.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.steelblue <- enrichGO(gene = steelblue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.tan <- enrichGO(gene = tan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.turquoise <- enrichGO(gene = turquoise.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.violet <- enrichGO(gene = violet.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.white <- enrichGO(gene = white.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
resGO.yellow <- enrichGO(gene = yellow.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
#resGO.yellowgreen <- enrichGO(gene = yellowgreen.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)

#DT::datatable(dplyr::select(resGO.turquoise@result,-geneID),rownames = F)
a1 <-  dotplot(resGO.black)
a2 <- dotplot(resGO.blue, showCategory = 10 )
a3 <- dotplot(resGO.brown, showCategory = 10 )
 dotplot(resGO.cyan, showCategory = 10 )
 dotplot(resGO.darkgreen, showCategory = 10 )
a23 <- dotplot(resGO.darkgrey, showCategory = 10 )
 dotplot(resGO.darkmagenta, showCategory = 10 )
a4 <- dotplot(resGO.darkolivegreen, showCategory = 10 )
dotplot(resGO.darkorange, showCategory = 10 )
a5 <- dotplot(resGO.darkred, showCategory = 10 )
 dotplot(resGO.darkturquoise, showCategory = 10 )
a6 <- dotplot(resGO.green, showCategory = 10 )
a7 <- dotplot(resGO.greenyellow, showCategory = 10 )
a8 <- dotplot(resGO.grey60, showCategory = 10 )
a9 <- dotplot(resGO.lightcyan, showCategory = 10 )
a10 <- dotplot(resGO.lightgreen, showCategory = 10 )
a11 <- dotplot(resGO.lightyellow, showCategory = 10 )
a12 <- dotplot(resGO.magenta, showCategory = 10 )
 dotplot(resGO.mediumpurple3, showCategory = 10 )
a13 <- dotplot(resGO.midnightblue, showCategory = 10 )
 dotplot(resGO.orange, showCategory = 10 )
a14 <- dotplot(resGO.orangered4, showCategory = 10 )
 dotplot(resGO.paleturquoise, showCategory = 10 )
a15 <- dotplot(resGO.pink, showCategory = 10 )
a16 <- dotplot(resGO.plum1, showCategory = 10 )
a17 <- dotplot(resGO.purple, showCategory = 10 )
 dotplot(resGO.red, showCategory = 10 )
 dotplot(resGO.royalblue, showCategory = 10 )
 dotplot(resGO.saddlebrown, showCategory = 10 )
 dotplot(resGO.salmon, showCategory = 10 )
 dotplot(resGO.sienna3, showCategory = 10 )
 dotplot(resGO.skyblue, showCategory = 10 )
 dotplot(resGO.skyblue3, showCategory = 10 )
a18 <- dotplot(resGO.steelblue, showCategory = 10 )
a19 <- dotplot(resGO.tan, showCategory = 10 )
a20 <- dotplot(resGO.turquoise, showCategory = 10 )
 dotplot(resGO.violet, showCategory = 10 )
 dotplot(resGO.white, showCategory = 10 )
a21 <- dotplot(resGO.yellow, showCategory = 10 )
a22 <- dotplot(resGO.yellowgreen, showCategory = 10 )

cowplot::plot_grid(a1, a2,a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, ncol=4, labels = "auto")
cowplot::plot_grid(a13,a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, ncol=4, labels = "auto")

# Bubble chart
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

r1 <- as.character(resGO.black$Description)
r2 <- as.character(resGO.blue$Description)
r3 <- as.character(resGO.brown$Description)
r4 <- as.character(resGO.cyan$Description)
r5 <- as.character(resGO.darkgreen$Description)
r6 <- as.character(resGO.darkgrey$Description)
r7 <- as.character(resGO.darkmagenta$Description)
r8 <- as.character(resGO.darkolivegreen$Description)
r9 <- as.character(resGO.darkorange$Description)
r10 <- as.character(resGO.darkred$Description)
r11 <- as.character(resGO.darkturquoise$Description)
r12 <- as.character(resGO.green$Description)
r13 <- as.character(resGO.greenyellow$Description)
r15 <- as.character(resGO.grey60$Description)
r16 <- as.character(resGO.lightcyan$Description)
r17 <- as.character(resGO.lightgreen$Description)
r18 <- as.character(resGO.lightyellow$Description)
r19 <- as.character(resGO.magenta$Description)
r20 <- as.character(resGO.mediumpurple3$Description)
r21 <- as.character(resGO.midnightblue$Description)
r22 <- as.character(resGO.orange$Description)
r23 <- as.character(resGO.orangered4$Description)
r24 <- as.character(resGO.paleturquoise$Description)
r25 <- as.character(resGO.pink$Description)
r26 <- as.character(resGO.plum1$Description)
r27 <- as.character(resGO.purple$Description)
r28 <- as.character(resGO.red$Description)
r29 <- as.character(resGO.royalblue$Description)
r30 <- as.character(resGO.saddlebrown$Description)
r31 <- as.character(resGO.salmon$Description)
r32 <- as.character(resGO.sienna3$Description)
r33 <- as.character(resGO.skyblue$Description)
r34 <- as.character(resGO.skyblue3$Description)
r35 <- as.character(resGO.steelblue$Description)
r36 <- as.character(resGO.tan$Description)
r37 <- as.character(resGO.turquoise$Description)
r38 <- as.character(resGO.violet$Description)
r39 <- as.character(resGO.white$Description)
r40 <- as.character(resGO.yellow$Description)
r41 <- as.character(resGO.yellowgreen$Description)
#http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know
full <- c(r1, r2, 	r3, 	r4, 	r5, 	r6, 	r7, 	r8, 	r9, 	r10, 	r11, 	r12, 	r13, 	r15, 	r16, 	r17, 	r18, 	r19, 	r20, 	r21, 	r22, 	r23, 	r24, 	r25, 	r26, 	r27, 	r28, 	r29, 	r30, 	r31, 	r32, 	r33, 	r34, 	r35, 	r36, 	r37, 	r38, 	r39, 	r40, 	r41)
docs <- Corpus(VectorSource(full))
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)
set.seed(12345)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))
#https://www.presentermedia.com/customize?id=25632#


#heatplot(resGO.turquoise, showCategory=30) # analyse overlap of categories
#enrichplot::upsetplot(resGO.turquoise) # analyse overlap of categories
goplot(resGO.turquoise) # GO terms hierarchy
cnetplot(resGO.turquoise) # this shows the genes connected to the enriched KEGG terms
edox <- setReadable(resGO.turquoise, OrgDb = SsalOrg)
a14 <- cnetplot(resGO.turquoise)
a15 <- cnetplot(edox, categorySize="pvalue", foldChange=turquoise.df$V1)
a16 <- cnetplot(edox,  circular = T, colorEdge = TRUE) 
cowplot::plot_grid(a14, a15, a16, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8,  1.2))



resKEGG.black <- enrichKEGG(gene = black.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.blue <- enrichKEGG(gene = blue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.brown <- enrichKEGG(gene = brown.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.cyan <- enrichKEGG(gene = cyan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkgreen <- enrichKEGG(gene = darkgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkgrey <- enrichKEGG(gene = darkgrey.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkmagenta <- enrichKEGG(gene = darkmagenta.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkolivegreen <- enrichKEGG(gene = darkolivegreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkorange <- enrichKEGG(gene = darkorange.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkred <- enrichKEGG(gene = darkred.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.darkturquoise <- enrichKEGG(gene = darkturquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.green <- enrichKEGG(gene = green.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.greenyellow <- enrichKEGG(gene = greenyellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.grey60 <- enrichKEGG(gene = grey60.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.lightcyan <- enrichKEGG(gene = lightcyan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.lightgreen <- enrichKEGG(gene = lightgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.lightyellow <- enrichKEGG(gene = lightyellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.magenta <- enrichKEGG(gene = magenta.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
#resKEGG.mediumpurple3 <- enrichKEGG(gene = mediumpurple3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.midnightblue <- enrichKEGG(gene = midnightblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.orange <- enrichKEGG(gene = orange.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.orangered4 <- enrichKEGG(gene = orangered4.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.paleturquoise <- enrichKEGG(gene = paleturquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.pink <- enrichKEGG(gene = pink.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
#resKEGG.plum1 <- enrichKEGG(gene = plum1.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.purple <- enrichKEGG(gene = purple.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.red <- enrichKEGG(gene = red.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.royalblue <- enrichKEGG(gene = royalblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.saddlebrown <- enrichKEGG(gene = saddlebrown.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.salmon <- enrichKEGG(gene = salmon.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
#resKEGG.sienna3 <- enrichKEGG(gene = sienna3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.skyblue <- enrichKEGG(gene = skyblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
#resKEGG.skyblue3 <- enrichKEGG(gene = skyblue3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.steelblue <- enrichKEGG(gene = steelblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.tan <- enrichKEGG(gene = tan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.turquoise <- enrichKEGG(gene = turquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.violet <- enrichKEGG(gene = violet.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.white <- enrichKEGG(gene = white.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
resKEGG.yellow <- enrichKEGG(gene = yellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)
#resKEGG.yellowgreen <- enrichKEGG(gene = yellowgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.05, qvalueCutoff  = 0.05)


#DT::datatable(dplyr::select(resKEGG.turquoise@result,-geneID),rownames = F)

a24 <- dotplot(resKEGG.black, showCategory = 10)
a25 <- dotplot(resKEGG.blue, showCategory = 10)
a26 <- dotplot(resKEGG.brown, showCategory = 10)
 dotplot(resKEGG.cyan, showCategory = 10)
 dotplot(resKEGG.darkgreen, showCategory = 10)
 dotplot(resKEGG.darkgrey, showCategory = 10)
 dotplot(resKEGG.darkmagenta, showCategory = 10)
 dotplot(resKEGG.darkolivegreen, showCategory = 10)
a27 <-  dotplot(resKEGG.darkorange, showCategory = 10)
 dotplot(resKEGG.darkred, showCategory = 10)
 dotplot(resKEGG.darkturquoise, showCategory = 10)
a28 <- dotplot(resKEGG.green, showCategory = 10)
a29 <- dotplot(resKEGG.greenyellow, showCategory = 10)
a30 <- dotplot(resKEGG.grey60, showCategory = 10)
a31 <- dotplot(resKEGG.lightcyan, showCategory = 10)
a32 <- dotplot(resKEGG.lightgreen, showCategory = 10)
 dotplot(resKEGG.lightyellow, showCategory = 10)
a33 <- dotplot(resKEGG.magenta, showCategory = 10)
 dotplot(resKEGG.mediumpurple3, showCategory = 10)
a43 <- dotplot(resKEGG.midnightblue, showCategory = 10)
 dotplot(resKEGG.orange, showCategory = 10)
 dotplot(resKEGG.orangered4, showCategory = 10)
 dotplot(resKEGG.paleturquoise, showCategory = 10)
 dotplot(resKEGG.pink, showCategory = 10)
 dotplot(resKEGG.plum1, showCategory = 10)
a35 <- dotplot(resKEGG.purple, showCategory = 10)
 dotplot(resKEGG.red, showCategory = 10)
a44 <- dotplot(resKEGG.royalblue, showCategory = 10)
dotplot(resKEGG.saddlebrown, showCategory = 10)
a37 <- dotplot(resKEGG.salmon, showCategory = 10)
 dotplot(resKEGG.sienna3, showCategory = 10)
 dotplot(resKEGG.skyblue, showCategory = 10)
 dotplot(resKEGG.skyblue3, showCategory = 10)
 dotplot(resKEGG.steelblue, showCategory = 10)
a38 <- dotplot(resKEGG.tan, showCategory = 10)
a39 <- dotplot(resKEGG.turquoise, showCategory = 10)
 dotplot(resKEGG.violet, showCategory = 10)
a40 <- dotplot(resKEGG.white, showCategory = 10)
a41 <- dotplot(resKEGG.yellow, showCategory = 10)
a42 <- dotplot(resKEGG.yellowgreen, showCategory = 10)

cowplot::plot_grid(a24, a25,a26, a27, a28, a29, a30, a31, a32, a33, a35,  ncol=4, labels = "auto")
cowplot::plot_grid(a37, a38, a39, a40, a41, a42, a43, a44, ncol=4, labels = "auto")


#Turquoise modules: GO and KEGG
aplot::plot_list(a20, a39, tag_levels='A')


library(DOSE)
library(enrichplot)
resKEGG.turquoise <- pairwise_termsim(resKEGG.turquoise)
emapplot(resKEGG.turquoise)


edox <- setReadable(resKEGG.turquoise, OrgDb = SsalOrg)

a50 <- cnetplot(edox, categorySize="pvalue")
a51 <- cnetplot(edox, foldChange=turquoise.df$V1, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(a50, a51, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.4))

edox2 <- pairwise_termsim(edox)
treeplot(edox2)
treeplot(edox2, hclust_method = "average")
aplot::plot_list(a20, a39, tag_levels='A')




##
#Deseq2 analysis

library(DESeq2)
library(BiocParallel)

dds$group <- factor(paste0(dds$Time, dds$New_Diet))
design(dds) <- ~ group-1

dds <- DESeq(dds, parallel = T)
resultsNames(dds)



# T1
res <- results(dds,
               contrast = list("groupT1mc1","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements

head(res_tbl, n = 25)



res <- results(dds,
               contrast = list("groupT1mc2","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT1mn3","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

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

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



res <- results(dds,
               contrast = list("groupT3mc2","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT3mn3","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

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












