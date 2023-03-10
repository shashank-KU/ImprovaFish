#PC correction

#WGCNA Host analysis
omics_data_host <- host_gut_rawCounts_total1
sample_info<- host_gut_mapping_file

omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% sample_info$common]
omics_data_host <- omics_data_host1
all(colnames(omics_data_host) == sample_info$common)

theme_set(theme_bw())




# Set seed as a precaution for reproducibility as some methods are non-deterministic.
set.seed(13118)
Run_analysis <- TRUE    # if FALSE it tries to load data instead of running the module creation step.

save_plot <- TRUE
plot_labeling_size = 20
prefix <- "h"

save_TOM <- FALSE       # Can get rather big
pam_stage <- FALSE      # Partitioning around medoids, tries to put more genes into modules where it is not directly clear from the dendrogram cutting.
set_verbose = 1         # How much detail from the WGCNA functions? A higher number means more detail.

omics_type = "rna"      # Just a name in the form of a string

tissue = "G"            # G for gut, or L for liver

what_samples = "FW_and_SW"     # FW, SW or FW_and_SW

take_average = F        # Take the average of each sample group

####

assoc_measure = "bicor"

max_block_size = 10000      #  13000 ~ 2-2.5 h |

####

applied_norm = "TMM"        #  TPM or (TMM, choose raw counts below)


applied_transf = "log2"     #  log2 or CLR


applied_filter = "sd"       #  sd or mad

pcCorrection <- T
if(pcCorrection){
  estimate_n.pc = T
  if(!estimate_n.pc){
    number_of_pcs_to_remove = 1 # Does not matter when pcCorrection is FALSE
  }
}





library(DESeq2)
omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% sample_info$common]
omics_data_host <- omics_data_host1
all(colnames(omics_data_host) == sample_info$common)
reorder_idx <- match(sample_info$common, colnames(omics_data_host))
omics_data_host <- omics_data_host[ , reorder_idx]
all(colnames(omics_data_host) %in% sample_info$common)
all(colnames(omics_data_host) == sample_info$common)


rownames(sample_info) <- sample_info$common
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
normalized_counts[1:2, 1:2]
omics_data_host<- normalized_counts 



omics_data_host <- t(omics_data_host)
dim(omics_data_host) %>% paste(c("Genes", "Samples"))

if(pcCorrection){
  library(sva)
  if(estimate_n.pc){
    
    print(t(omics_data_host) %>% dim() %>% paste(c("Samples", "genes"))) 
    mod=matrix(1,nrow=nrow(t(omics_data_host)),ncol=1)
    colnames(mod)="Intercept"
    
    print(omics_data_host %>% dim() %>% paste(c("genes", "Samples"))) 
    n.pc = num.sv(dat = omics_data_host, mod = mod,  method="be")
  } else{
    n.pc = number_of_pcs_to_remove
  }
  
  print(omics_data_host %>% dim() %>% paste(c("genes", "Samples")))
  omics_data_host <- sva_network(omics_data_host, n.pc = n.pc)
}

if(pcCorrection){
  report <- paste0("PC-correction with ", n.pc, " PC(s) removed.")
  
  cowplot::plot_grid(plot_hist(omics_data_host, "sd"),
                     plot_hist(omics_data_host, "mad"),
                     plot_hist(omics_data_host, "no"),
                     align = "h",
                     ncol = 3)
} else{
  report <- "There was no PC-correction performed."
}
report

omics_data_host <- omics_data_host %>% t()
dim(omics_data_host) %>% paste(c("Samples", "Genes"))

powers <- c(c(1:10), seq(from = 12, to=20, by=2))
library(WGCNA)

sft <- pickSoftThreshold(omics_data_host, 
                         powerVector = powers, 
                         verbose = set_verbose, 
                         networkType = "signed",
                         corFn= chosen_parameter_set$assoc_measure)



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
  ylim(-1,1) +
  geom_segment(aes(x = st, y = 0.25, xend = st, yend = sfApprox[idx]-0.05), 
               arrow = arrow(length = unit(0.2,"cm")), 
               size = 0.5)-> scale_independence_plot 



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



st
modules.omics_X_PC_corrected  <- blockwiseModules(omics_data_host, 
                                    power = st, 
                                    networkType = "signed", 
                                    TOMType = "signed",
                                    corType = assoc_measure,
                                    maxPOutliers = 0.05,
                                    deepSplit = 4, # Default 2
                                    minModuleSize = 20,           # 30
                                    minCoreKME = 0.5,            # Default 0.5
                                    minCoreKMESize = 2,          # Default minModuleSize/3,
                                    minKMEtoStay = 0.5,          # Default 0.3
                                    reassignThreshold = 0,       # Default 1e-6
                                    mergeCutHeight = 0.2,        # Default 0.15
                                    pamStage = pam_stage, 
                                    pamRespectsDendro = TRUE,
                                    replaceMissingAdjacencies = TRUE,
                                    nThreads = 10,
                                    numericLabels = TRUE,
                                    saveTOMs = save_TOM,
                                    saveTOMFileBase = "HOST_TOM",
                                    verbose = 3,
                                    maxBlockSize=8000)




rownames(modules.omics_X_PC_corrected$MEs) <- rownames(omics_data_host)
names(modules.omics_X_PC_corrected$colors) <- colnames(omics_data_host)
names(modules.omics_X_PC_corrected$unmergedColors) <- colnames(omics_data_host)

hubs <- chooseTopHubInEachModule(omics_data_host, modules.omics_X_PC_corrected$colors, omitColors = 0)
hubs <- chooseTopHubInEachModule(omics_data_host, modules.omics_X_PC_corrected$colors, power = st, omitColors = "0")

stage2results_X_PC_corrected <- list(modules = modules.omics_X_PC_corrected , 
                        hubs = hubs)


# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_X_PC_corrected$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_X_PC_corrected$modules$goodSamples) == length(stage2results_X_PC_corrected$modules$goodSamples)
genes_good <- sum(stage2results_X_PC_corrected$modules$goodGenes) == length(stage2results_X_PC_corrected$modules$goodGenes)

ME_good <- sum(stage2results_X_PC_corrected$modules$MEsOK) == length(stage2results_X_PC_corrected$modules$MEsOK)
table(stage2results_X_PC_corrected$modules$colors) %>% 
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
for(i in seq_along(stage2results_X_PC_corrected$modules$dendrograms)){
  plotDendroAndColors(stage2results_X_PC_corrected$modules$dendrograms[[i]], merged_colors[stage2results_X_PC_corrected$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_X_PC_corrected$modules$blockGenes[[i]]),
                                    " genes"))
}


MEs <- stage2results_X_PC_corrected$modules$MEs


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

#density_eigen

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
  colnames(stage2results_X_PC_corrected$modules$MEs)[colnames(stage2results_X_PC_corrected$modules$MEs) %>% sub("ME", "", .) %>% as.numeric() %>% order()]

for(m in colnames(stage2results_X_PC_corrected$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(omics_data_host = omics_data_host, modules = stage2results_X_PC_corrected$modules, module_x = h))) %>% 
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


cowplot::plot_grid(part_2, density_all_plot, ncol = 1, rel_heights = c(0.8,1), labels = c("", "F"), label_size = plot_labeling_size) #1000*1000 size of plot for png

#Hub genes
stage2results_X_PC_corrected$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("gene_id" = ".") %>% 
  tibble::rownames_to_column(var = "Module") -> top_genes

top_genes







# Further downstream analysis
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
#modules.omics_X_PC_corrected .backup <- modules.omics_X_PC_corrected  
module_eigengenes <- modules.omics_X_PC_corrected $MEs

all.equal(sample_info$common, rownames(module_eigengenes))
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


#Let???s make plot of module 7
module_19_df <- module_eigengenes %>%
  tibble::rownames_to_column("common") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(meta %>%
                      dplyr::select(common, Day, Time, New_Diet),
                    by = c("common" = "common"))


ggplot(module_19_df, aes(x = Time, y = ME1, color = Day)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic() + facet_wrap("New_Diet")


ggplot(module_19_df, aes(x = Day, y = ME1, color = Time)) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


#What genes are a part of module 7?
gene_module_key <- tibble::enframe(modules.omics_X_PC_corrected $colors, name = "gene", value = "module") %>%
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
moduleLabels <- modules.omics_X_PC_corrected$colors
moduleColors <- labels2colors(modules.omics_X_PC_corrected$colors)
table(moduleColors)
table(moduleLabels)

table(moduleColors, moduleLabels)

MEs <- modules.omics_X_PC_corrected$MEs
geneTree <- modules.omics_X_PC_corrected$dendrograms[[1]]


#https://gitlab.com/garethgillard/ssalv3annotation
annot1 <- read.csv("Salmo_salar-GCA_905237065.2_gene_annotations.csv", header = TRUE)
probes <- colnames(omics_data_host) 
probes2annot <- match(probes, annot1$gene_id)
allLLIDs <- annot1$v2.gene_id.NCBI[probes2annot]

intModules <- unique(moduleColors)

setwd("resGO_after_PC_corrected/")
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



table(moduleColors)
getwd()

black.df <- read.table("LocusLinkIDs-black.txt")
blue.df <- read.table("LocusLinkIDs-blue.txt")
brown.df <- read.table("LocusLinkIDs-brown.txt")
cyan.df <- read.table("LocusLinkIDs-cyan.txt")
green.df <- read.table("LocusLinkIDs-green.txt")
greenyellow.df <- read.table("LocusLinkIDs-greenyellow.txt")
grey60.df <- read.table("LocusLinkIDs-grey60.txt")
lightcyan.df <- read.table("LocusLinkIDs-lightcyan.txt")
lightgreen.df <- read.table("LocusLinkIDs-lightgreen.txt")
magenta.df <- read.table("LocusLinkIDs-magenta.txt")
midnightblue.df <- read.table("LocusLinkIDs-midnightblue.txt")
pink.df <- read.table("LocusLinkIDs-pink.txt")
purple.df <- read.table("LocusLinkIDs-purple.txt")
red.df <- read.table("LocusLinkIDs-red.txt")
salmon.df <- read.table("LocusLinkIDs-salmon.txt")
tan.df <- read.table("LocusLinkIDs-tan.txt")
turquoise.df <- read.table("LocusLinkIDs-turquoise.txt")
yellow.df <- read.table("LocusLinkIDs-yellow.txt")

all.df <- read.table("LocusLinkIDs-all.txt")
all.df$V1 <- as.character(all.df$V1)

resGO.black <- enrichGO(gene = black.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.blue <- enrichGO(gene = blue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.brown <- enrichGO(gene = brown.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.cyan <- enrichGO(gene = cyan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.green <- enrichGO(gene = green.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.greenyellow <- enrichGO(gene = greenyellow.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.grey60 <- enrichGO(gene = grey60.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.lightcyan <- enrichGO(gene = lightcyan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.lightgreen <- enrichGO(gene = lightgreen.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.magenta <- enrichGO(gene = magenta.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.midnightblue <- enrichGO(gene = midnightblue.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.pink <- enrichGO(gene = pink.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.purple <- enrichGO(gene = purple.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.red <- enrichGO(gene = red.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.salmon <- enrichGO(gene = salmon.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.tan <- enrichGO(gene = tan.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.turquoise <- enrichGO(gene = turquoise.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
resGO.yellow <- enrichGO(gene = yellow.df$V1, ont = "BP", universe = all.df$V1, OrgDb = SsalOrg, pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)

#DT::datatable(dplyr::select(resGO.turquoise@result,-geneID),rownames = F)
a1 <- dotplot(resGO.black)
a2 <- dotplot(resGO.blue, showCategory = 10 )
a3 <- dotplot(resGO.brown, showCategory = 10 )
a4 <- dotplot(resGO.cyan, showCategory = 10 )
a12 <- dotplot(resGO.green, showCategory = 10 )
a13 <- dotplot(resGO.greenyellow, showCategory = 10 )
dotplot(resGO.grey60, showCategory = 10 )
dotplot(resGO.lightcyan, showCategory = 10 )
dotplot(resGO.lightgreen, showCategory = 10 )
dotplot(resGO.magenta, showCategory = 10 )
dotplot(resGO.midnightblue, showCategory = 10 )
dotplot(resGO.pink, showCategory = 10 )
a27 <- dotplot(resGO.purple, showCategory = 10 )
a28 <- dotplot(resGO.red, showCategory = 10 )
a29 <- dotplot(resGO.salmon, showCategory = 10 )
dotplot(resGO.tan, showCategory = 10 )
a37 <- dotplot(resGO.turquoise, showCategory = 10 )
a40 <- dotplot(resGO.yellow, showCategory = 10 )
cowplot::plot_grid(a1,a2, a3, a4, a12, a13, a27, a28, a29, a37, a40, ncol=4, labels = "auto")

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



resKEGG.black <- enrichKEGG(gene = black.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.blue <- enrichKEGG(gene = blue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.brown <- enrichKEGG(gene = brown.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.cyan <- enrichKEGG(gene = cyan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkgreen <- enrichKEGG(gene = darkgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkgrey <- enrichKEGG(gene = darkgrey.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkmagenta <- enrichKEGG(gene = darkmagenta.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkolivegreen <- enrichKEGG(gene = darkolivegreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkorange <- enrichKEGG(gene = darkorange.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkred <- enrichKEGG(gene = darkred.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.darkturquoise <- enrichKEGG(gene = darkturquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.green <- enrichKEGG(gene = green.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.greenyellow <- enrichKEGG(gene = greenyellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.grey60 <- enrichKEGG(gene = grey60.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.lightcyan <- enrichKEGG(gene = lightcyan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.lightgreen <- enrichKEGG(gene = lightgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.lightyellow <- enrichKEGG(gene = lightyellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.magenta <- enrichKEGG(gene = magenta.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.mediumpurple3 <- enrichKEGG(gene = mediumpurple3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.midnightblue <- enrichKEGG(gene = midnightblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.orange <- enrichKEGG(gene = orange.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.orangered4 <- enrichKEGG(gene = orangered4.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.paleturquoise <- enrichKEGG(gene = paleturquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.pink <- enrichKEGG(gene = pink.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.plum1 <- enrichKEGG(gene = plum1.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.purple <- enrichKEGG(gene = purple.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.red <- enrichKEGG(gene = red.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.royalblue <- enrichKEGG(gene = royalblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.saddlebrown <- enrichKEGG(gene = saddlebrown.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.salmon <- enrichKEGG(gene = salmon.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.sienna3 <- enrichKEGG(gene = sienna3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.skyblue <- enrichKEGG(gene = skyblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.skyblue3 <- enrichKEGG(gene = skyblue3.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.steelblue <- enrichKEGG(gene = steelblue.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.tan <- enrichKEGG(gene = tan.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.turquoise <- enrichKEGG(gene = turquoise.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.violet <- enrichKEGG(gene = violet.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.white <- enrichKEGG(gene = white.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.yellow <- enrichKEGG(gene = yellow.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)
resKEGG.yellowgreen <- enrichKEGG(gene = yellowgreen.df$V1, universe = all.df$V1, organism = "sasa", pvalueCutoff = 0.01, qvalueCutoff  = 0.05)


#DT::datatable(dplyr::select(resKEGG.turquoise@result,-geneID),rownames = F)

a42 <- dotplot(resKEGG.black, showCategory = 10)
a43 <- dotplot(resKEGG.blue, showCategory = 10)
a44 <- dotplot(resKEGG.brown, showCategory = 10)
a45 <- dotplot(resKEGG.cyan, showCategory = 10)
dotplot(resKEGG.darkgreen, showCategory = 10)
dotplot(resKEGG.darkgrey, showCategory = 10)
dotplot(resKEGG.darkmagenta, showCategory = 10)
dotplot(resKEGG.darkolivegreen, showCategory = 10)
dotplot(resKEGG.darkorange, showCategory = 10)
dotplot(resKEGG.darkred, showCategory = 10)
dotplot(resKEGG.darkturquoise, showCategory = 10)
a53 <- dotplot(resKEGG.green, showCategory = 10)
a54 <- dotplot(resKEGG.greenyellow, showCategory = 10)
a56 <- dotplot(resKEGG.grey60, showCategory = 10)
a57 <- dotplot(resKEGG.lightcyan, showCategory = 10)
a58 <- dotplot(resKEGG.lightgreen, showCategory = 10)
dotplot(resKEGG.lightyellow, showCategory = 10)
a60 <- dotplot(resKEGG.magenta, showCategory = 10)
a61 <- dotplot(resKEGG.mediumpurple3, showCategory = 10)
dotplot(resKEGG.midnightblue, showCategory = 10)
dotplot(resKEGG.orange, showCategory = 10)
dotplot(resKEGG.orangered4, showCategory = 10)
dotplot(resKEGG.paleturquoise, showCategory = 10)
a66 <- dotplot(resKEGG.pink, showCategory = 10)
dotplot(resKEGG.plum1, showCategory = 10)
a68 <- dotplot(resKEGG.purple, showCategory = 10)
a69 <- dotplot(resKEGG.red, showCategory = 10)
dotplot(resKEGG.royalblue, showCategory = 10)
dotplot(resKEGG.saddlebrown, showCategory = 10)
dotplot(resKEGG.salmon, showCategory = 10)
dotplot(resKEGG.sienna3, showCategory = 10)
dotplot(resKEGG.skyblue, showCategory = 10)
dotplot(resKEGG.skyblue3, showCategory = 10)
a76 <- dotplot(resKEGG.steelblue, showCategory = 10)
a77 <- dotplot(resKEGG.tan, showCategory = 10)
a78 <- dotplot(resKEGG.turquoise, showCategory = 10)
dotplot(resKEGG.violet, showCategory = 10)
dotplot(resKEGG.white, showCategory = 10)
a81 <- dotplot(resKEGG.yellow, showCategory = 10)
a82 <- dotplot(resKEGG.yellowgreen, showCategory = 10)


library(DOSE)
library(enrichplot)
resKEGG.yellow <- pairwise_termsim(resKEGG.yellow)
emapplot(resKEGG.yellow)
gseKEGG(turquoise.df$V1, organism = "sasa")

edox <- setReadable(resGO.yellowgreen, OrgDb = SsalOrg)

a29 <- cnetplot(edox, categorySize="pvalue")
a30 <- cnetplot(edox, foldChange=yellow.df$V1, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(a1, a2, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.4))

edox <- setReadable(resGO.yellowgreen, OrgDb = SsalOrg)
edox2 <- pairwise_termsim(edox)
a31 <- treeplot(edox2)
a32 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(a31, a32, tag_levels='A')






#Select the gene modules
power = 14
adjacency = adjacency(omics_data_host, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");
rm(dissTOM)
rm(TOM)
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


table(merged_colors)
plotDendroAndColors(geneTree, merged_colors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")


MEList = moduleEigengenes(omics_data_host, colors = merged_colors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(omics_data_host, merged_colors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
#pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(merged_colors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  


# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(omics_data_host)
  inModule = is.finite(match(merged_colors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = merged_colors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(omics_data_host)
  inModule = is.finite(match(merged_colors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = merged_colors[inModule]);
}

rm(modTOM) # remove large object

#if(!"RCy3" %in% installed.packages()){
#  install.packages("BiocManager")
#  BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-brown.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-brown.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)








































