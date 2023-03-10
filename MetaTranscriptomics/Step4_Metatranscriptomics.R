setwd("/Users/shashankgupta/Desktop/ImprovAFish/metatranscriptomics/")
library(DESeq2)
library(WGCNA)
library(stringr)
library(dplyr)
load("/Users/shashankgupta/Desktop/ImprovAFish/metatranscriptomics/CountsKallistoQuantificationDeNovoSTxDRAM/dds.Stx.RData")
dds <- dds[ , dds$TimePoint != "T0"]
keep <- rowSums(counts(dds) >=5) >= 5
dds <- dds[keep,]
dds
dds_norm <- vst(dds)
metatranscriptomics <- assay(dds_norm) %>%
  t() # Transpose this data
rm(dds_norm)

#keep <- rowSums(counts(dds) >=5) >= 5
#modules.omics.MetaTrans.backup <-  modules.omics.MetaTrans

# ######
# dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
# dds <- dds[idx,]
# dds
# dds_norm <- vst(dds)
# metatranscriptomics <- assay(dds_norm) %>%
#   t() # Transpose this data
# rm(dds_norm)
# # modules.omics.MetaTrans.bk <- modules.omics.MetaTrans


dim(metatranscriptomics) %>% paste(c("Samples", "metatranscriptomics"))


powers <- c(c(1:10), seq(from = 12, to=20, by=2))

library(WGCNA)
allowWGCNAThreads()

sft <- pickSoftThreshold(metatranscriptomics, 
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



modules.omics.MetaTrans <- blockwiseModules(metatranscriptomics, 
                                    power = st, 
                                    networkType = "signed", 
                                    TOMType = "signed",
                                    corType = assoc_measure,
                                    #maxPOutliers = 0.05,
                                    #deepSplit = 4, # Default 2
                                    minModuleSize = 21,           # 30
                                    #minCoreKME = 0.5,            # Default 0.5
                                    #minCoreKMESize = 2,          # Default minModuleSize/3,
                                    #minKMEtoStay = 0.5,          # Default 0.3
                                    #reassignThreshold = 0,       # Default 1e-6
                                    #mergeCutHeight = 0.2,        # Default 0.15
                                    #pamStage = pam_stage, 
                                    #pamRespectsDendro = TRUE,
                                    #replaceMissingAdjacencies = TRUE,
                                    #nThreads = 10,
                                    numericLabels = TRUE,
                                    saveTOMs = save_TOM,
                                    saveTOMFileBase = "HOST_TOM",
                                    verbose = 3,
                                    maxBlockSize=8000)



rownames(modules.omics.MetaTrans$MEs) <- rownames(metatranscriptomics)
names(modules.omics.MetaTrans$colors) <- colnames(metatranscriptomics)
names(modules.omics.MetaTrans$unmergedColors) <- colnames(metatranscriptomics)
table(modules.omics.MetaTrans$colors)
hubs <- chooseTopHubInEachModule(metatranscriptomics, modules.omics.MetaTrans$colors, power = st, omitColors = "0")

stage2results_MetaTrans <- list(modules = modules.omics.MetaTrans, 
                        hubs = hubs)

library(ggplot2)
# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_MetaTrans$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_MetaTrans$modules$goodSamples) == length(stage2results_MetaTrans$modules$goodSamples)
genes_good <- sum(stage2results_MetaTrans$modules$goodGenes) == length(stage2results_MetaTrans$modules$goodGenes)

ME_good <- sum(stage2results_MetaTrans$modules$MEsOK) == length(stage2results_MetaTrans$modules$MEsOK)
table(stage2results_MetaTrans$modules$colors) %>% 
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


table(stage2results_MetaTrans$modules$colors) %>% as.data.frame() -> res
res$`Module color` <- WGCNA::labels2colors(as.numeric(as.character(res$Var1)))
res <- res[, c(1,3,2)]
colnames(res) <- c("Module", "Module color", "Number of metatranscriptomics")

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_MetaTrans$modules$dendrograms)){
  plotDendroAndColors(stage2results_MetaTrans$modules$dendrograms[[i]], merged_colors[stage2results_MetaTrans$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_MetaTrans$modules$blockGenes[[i]]),
                                    " metatranscriptomics"))
}
MEs <- stage2results_MetaTrans$modules$MEs

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

all(rownames(metatranscriptomics) == rownames(MEs))
dim(metatranscriptomics) %>% paste0(c(" samples", " metatranscriptomics"))
kME <- bicor(metatranscriptomics, MEs, maxPOutliers = 0.05)
dim(kME) %>% paste0(c(" metatranscriptomics", " modules"))

intra_cor <- c()
for (i in 1:ncol(metatranscriptomics)) {
  m <- stage2results_MetaTrans$modules$colors[i]
  intra_cor[i] <- kME[i, paste0("ME", m)]
  if(m != 0){
    intra_cor[i] <- kME[i, paste0("ME", m)]
  } else{
    intra_cor[i] <- NA
  }
  
}

idx <- which(is.na(intra_cor))
intra_cor <- intra_cor[-idx]

plot(density(intra_cor), main = "Correlations with module-eigenmetatranscriptomics (within module correlation)\nNo ME0", xlim = c(-1,1))



# Corr within modules
corr_within_module <- function(metatranscriptomics, modules, module_x = 1){
  idx.metatranscriptomics <- which(modules$colors == module_x)
  idx.me <- which(colnames(modules$MEs) == paste0("ME",module_x))
  kME_x <- bicor(metatranscriptomics[,idx.metatranscriptomics], modules$MEs[,idx.me], maxPOutliers = 0.05)
  kME_x
}

ggplot.list <- list()

for(m in colnames(stage2results_MetaTrans$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(metatranscriptomics = metatranscriptomics, modules = stage2results_MetaTrans$modules, module_x = h))) %>% 
    ggplot() + 
    #geom_density(aes(x = x), fill = labels2colors(h), color = "black", alpha = 0.5) + 
    geom_histogram(aes(x), fill = labels2colors(h), color = "black", alpha = 0.5, bins = 20) + 
    xlim(-1, 1) +
    xlab("metatranscriptomics correlation")+
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









#stage2results_MetaTrans_backup <- stage2results_MetaTrans
stage2results_MetaTrans

# Remove ME0 from analysis
idx <- which(colnames(stage2results_MetaTrans$modules$MEs) == "ME0")
stage2results_MetaTrans$modules$MEs <- stage2results_MetaTrans$modules$MEs[,-idx]


if(take_average){
  stage2results_MetaTrans$modules$MEs %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_HGm[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_MetaTrans_eigengenes
  
  stage2results_MetaTrans_eigengenes <- as.data.frame(stage2results_MetaTrans_eigengenes)
} else{
  stage2results_MetaTrans_eigengenes <- stage2results_MetaTrans$modules$MEs
}

rownames(stage2results_MetaTrans_eigengenes) <- stage2results_MetaTrans_eigengenes$common
stage2results_MetaTrans_eigengenes$common <- NULL

dim(stage2results_MetaTrans_eigengenes); dim(stage2results_X_eigengenes)



table(rownames(stage2results_X_eigengenes) %in% rownames(stage2results_MetaTrans_eigengenes))

stage2results_X_eigengenes_New <- stage2results_X_eigengenes[rownames(stage2results_X_eigengenes) %in% rownames(stage2results_MetaTrans_eigengenes), ]


table(rownames(stage2results_X_eigengenes_New) %in% rownames(stage2results_MetaTrans_eigengenes))



#Correlate modules from transcriptomics and metagenomics.

# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(stage2results_MetaTrans_eigengenes) == rownames(stage2results_X_eigengenes_New))
if(!same_order){
  stage2results_MetaTrans_eigengenes <- stage2results_MetaTrans_eigengenes[order(rownames(stage2results_MetaTrans_eigengenes)),]
  stage2results_X_eigengenes_New <- stage2results_X_eigengenes_New[order(rownames(stage2results_X_eigengenes_New)),]
  same_order <- all(rownames(stage2results_MetaTrans_eigengenes) == rownames(stage2results_X_eigengenes_New))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


#####
p.value_matr <- corr.value_matr <- matrix(ncol = ncol(stage2results_MetaTrans_eigengenes), 
                                          nrow = ncol(stage2results_X_eigengenes_New), 
                                          dimnames = list(colnames(stage2results_X_eigengenes_New), 
                                                          colnames(stage2results_MetaTrans_eigengenes)))


for(i in 1:ncol(stage2results_X_eigengenes_New)){
  for(j in 1:ncol(stage2results_MetaTrans_eigengenes)){
    cor.res <- cor.test(stage2results_X_eigengenes_New[,i], stage2results_MetaTrans_eigengenes[,j])
    p.value_matr[i, j] <- cor.res$p.value
    corr.value_matr[i, j] <- cor.res$estimate
  }
}

# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_eigengenes_New), colnames(stage2results_MetaTrans_eigengenes))


# Add significance level.  
# One star means a p-value of less than 0.05; Two stars is less than 0.01, and three, is less than 0.001.

signif_matrix <- rep("", length(p.value_matr))
three_star <- which( p.value_matr <= 0.001)
signif_matrix[three_star] <- "***"
two_star <- which((p.value_matr <= 0.01) & (p.value_matr > 0.001))
signif_matrix[two_star] <- "**"
one_star <- which((p.value_matr <= 0.05) & (p.value_matr > 0.01))
signif_matrix[one_star] <- "*"
dim(signif_matrix) = dim(p.value_matr) # Give textMatrix the correct dimensions 


# Collect all results into a list.
MetaTrans_corr_X <- list(p_value = p.value_matr, 
                 p_value_adj = p.value_matr.adjust,
                 signif_matrix = signif_matrix,
                 correlation = corr.value_matr)
rm(p.value_matr, p.value_matr.adjust, signif_matrix, corr.value_matr)



heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)



#Annotation for Eigen Metatranscriptomics
meta_info <- read_tsv("CountsKallistoQuantificationDeNovoSTxDRAM/Total.DRAM.STx.bac.k2.total.annotations.tsv")
rownames(meta_info) <- meta_info$...1
stage2results_MetaTrans$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("metatranscriptomics_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubmetatranscriptomics


dplyr::left_join(hubmetatranscriptomics, 
                 (meta_info %>%
                    tibble::rownames_to_column(var = "metatranscriptomics_name")), 
                 by = "metatranscriptomics_name") -> hubmetatranscriptomics_Before

hubmetatranscriptomics_Before$Modulemetatranscriptomic <- paste0("ME" ,hubmetatranscriptomics_Before$Module)
hubmetatranscriptomics_Before <- hubmetatranscriptomics_Before[, c(33, 15)]

#Match with Eigen Metatranscriptomics
test <- MetaTrans_corr_X$correlation
match(hubmetatranscriptomics_Before[, "Modulemetatranscriptomic"], colnames(test))
colnames(test)[match(hubmetatranscriptomics_Before[,"Modulemetatranscriptomic"], colnames(test))] = hubmetatranscriptomics_Before[,"uniref_taxonomy"]

MetaTrans_corr_X$correlation <- test

pheatmap::pheatmap(MetaTrans_corr_X$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,  # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = MetaTrans_corr_X$signif_matrix, 
                   fontsize_number = 8, #10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F,
                   show_rownames = F, legend = F,
                   labels_row = paste0(prefix_OTUs, rownames(MetaTrans_corr_X$correlation)),
                   labels_col = paste0(colnames(MetaTrans_corr_X$correlation)),
                   main = paste("Eigen\n", "Metatranscriptomics" )) -> MetaTrans_corr_X_plot





##
#Deseq2 analysis

library(DESeq2)
library(BiocParallel)
load("/Users/shashankgupta/Desktop/ImprovAFish/metatranscriptomics/CountsKallistoQuantificationDeNovoSTxDRAM/dds.Stx.RData")
dds <- dds[ , dds$TimePoint != "T0"]
keep <- rowSums(counts(dds) >=5) >= 5
dds <- dds[keep,]
dds$group <- factor(paste0(dds$TimePoint, dds$Treatment))
design(dds) <- ~ group-1

dds <- DESeq(dds, parallel = T)
resultsNames(dds)

annot_function_metatrans <- read_tsv("../../../metatranscriptomics/CountsKallistoQuantificationDeNovoSTxDRAM/Total.DRAM.STx.bac.k2.total.annotations.tsv")
annot_function_metatrans <- annot_function_metatrans[, c("ID", "kegg_hit", "uniref_taxonomy","pfam_hits", "cazy_hits" )]


# T1
res <- results(dds,
               contrast = list("groupT1mc1","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements

head(res_tbl, n = 25)



res <- results(dds,
               contrast = list("groupT1mc2","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)


length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT1mn3","groupT1ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



# T2
res <- results(dds,
               contrast = list("groupT2mc1","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements

head(res_tbl, n = 25)



res <- results(dds,
               contrast = list("groupT2mc2","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 

res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)

length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT2mn3","groupT2ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)



length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements


# T3
res <- results(dds,
               contrast = list("groupT3mc1","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)


length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements



res <- results(dds,
               contrast = list("groupT3mc2","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj)
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)


length(res_tbl$log2FoldChange[res_tbl$log2FoldChange>0]) # count positive elements
length(res_tbl$log2FoldChange[res_tbl$log2FoldChange<0]) # count negative elements




res <- results(dds,
               contrast = list("groupT3mn3","groupT3ctr"))
table(res$padj < 0.05)

res_tbl <- as_tibble(res, rownames = "ENSEMBL") %>%
  filter(padj <0.05)%>%
  arrange(padj) 
res_tbl <- merge(res_tbl, annot_function_metatrans, by.x = "ENSEMBL", by.y = "ID", all.x = TRUE )
res_tbl %>% 
  filter(log2FoldChange > 0) %>%
  select(pfam_hits, uniref_taxonomy)


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

res <- results(dds,contrast = list("groupT3mn3","groupT3ctr")) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat("DEGs: ", nrow(res))
vsd <- vst(dds)
library(stringi)
p1 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% res$Genes[1:23]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(Treatment = stri_replace_all_regex(Sample, colnames(dds), dds$Treatment, vectorize=FALSE)) %>%
  ggplot(aes(x = Treatment, y = Expression, color = Treatment)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes)) + theme_bw()

p2 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% res$Genes[1:23]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(TimePoint = stri_replace_all_regex(Sample, colnames(dds), dds$TimePoint, vectorize=FALSE)) %>%
  ggplot(aes(x = TimePoint, y = Expression, color = TimePoint)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes)) + theme_bw()

cowplot::plot_grid(plotlist = list(p1,p2), ncol = 2)











