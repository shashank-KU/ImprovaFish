library("MetaboAnalystR")
getwd()
mSet<-InitDataObjects("conc", "mf", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, "../../Metabolomics/Lipidomics/Lipidomics Data/Lipidomics.csv", "colmf")
mSet<-ReadMetaData(mSet, "../../Metabolomics/Lipidomics/Lipidomics Data/Metadata.csv")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-SanityCheckMeta(mSet, 1)
mSet<-SetDataTypeOfMeta(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MeanNorm", "SrNorm", "NULL", ratio=FALSE, ratioNum=20)
Lipidomics <- data.frame(mSet$dataSet$norm)
Lipidomics <- t(Lipidomics)

#Look for outliers by examining tree of Lipidomics  
dim(Lipidomics) %>% paste(c("Lipidomics", "Samples"))
sampleTree = hclust(dist(Lipidomics), method = "average");
plot(sampleTree, main = "Lipidomics clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# outlier<-c("N.Lauroylsarcosine...271.21434.Da.608.71.s")
# Lipidomics<- Lipidomics %>%
#   as.data.frame %>%
#   filter(!row.names(Lipidomics) %in% outlier)
# sampleTree = hclust(dist(Lipidomics), method = "average");
# plot(sampleTree, main = "Lipidomics clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#Look for outlier samples by examining tree of samples
#Transpose such that samples are in rows and Lipidomics are in columns.  
Lipidomics <- t(Lipidomics)
dim(Lipidomics) %>% paste(c("Samples", "Lipidomics"))
sampleTree = hclust(dist(Lipidomics), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# outlier<-c("T3_ctr_tk3_VFA07")
# Lipidomics<- Lipidomics %>%
#   as.data.frame %>%
# filter(!row.names(Lipidomics) %in% outlier)

#Transpose such that samples are in rows and Lipidomics are in columns.  
dim(Lipidomics) %>% paste(c("Samples", "Lipidomics"))
powers <- c(1:10, seq(12,30,2))
sft <- pickSoftThreshold(Lipidomics, 
                         powerVector = powers, 
                         verbose = set_verbose, 
                         networkType = "signed",
                         corFn= chosen_parameter_set$assoc_measure)

allowWGCNAThreads()


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
  xlim(1,30) + theme_bw() +
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
#write.table(Lipidomics, "Lipidomics.txt", sep = "\t")
#Lipidomics<- read.table("Lipidomics.txt", header = T, row.names = 1, stringsAsFactors = F, fill = T)

modules.omics.Lipidomics <- blockwiseModules(Lipidomics,
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


rownames(modules.omics.Lipidomics$MEs) <- rownames(Lipidomics)
names(modules.omics.Lipidomics$colors) <- colnames(Lipidomics)
names(modules.omics.Lipidomics$unmergedColors) <- colnames(Lipidomics)

hubs_M <- chooseTopHubInEachModule(Lipidomics, modules.omics.Lipidomics$colors, power = st, omitColors = "0")

stage2results_Lipidomics <- list(modules = modules.omics.Lipidomics, 
                         hubs = hubs_M)





# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_Lipidomics$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_Lipidomics$modules$goodSamples) == length(stage2results_Lipidomics$modules$goodSamples)
genes_good <- sum(stage2results_Lipidomics$modules$goodGenes) == length(stage2results_Lipidomics$modules$goodGenes)

ME_good <- sum(stage2results_Lipidomics$modules$MEsOK) == length(stage2results_Lipidomics$modules$MEsOK)
table(stage2results_Lipidomics$modules$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename(Module = Var1, Size = Freq) %>%
  dplyr::mutate(Module_color = labels2colors(as.numeric(as.character(Module)))) -> module_size



module_size %>% 
  ggplot(aes(x = Module, y = Size, fill = Module)) +
  geom_col(color =  "#000000") +
  ggtitle("Number of genes in each module") +
  theme(legend.position = "none") + 
  theme_bw()+
  scale_fill_manual(values = setNames(module_size$Module_color,module_size$Module)) +
  geom_text(aes(label = Size),vjust = 0.5, hjust = -0.18, size = 3.5) +
  ylim(0, max(module_size$Size)*1.1) +
  theme(plot.margin = margin(2, 2, 2, 2, "pt")) +
  coord_flip()-> module_size_barplot

module_size_barplot + theme_bw()


table(stage2results_Lipidomics$modules$colors) %>% as.data.frame() -> res
res$`Module color` <- WGCNA::labels2colors(as.numeric(as.character(res$Var1)))
res <- res[, c(1,3,2)]
colnames(res) <- c("Module", "Module color", "Number of Lipidomics")

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_Lipidomics$modules$dendrograms)){
  plotDendroAndColors(stage2results_Lipidomics$modules$dendrograms[[i]], merged_colors[stage2results_Lipidomics$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_Lipidomics$modules$blockGenes[[i]]),
                                    " Lipidomics"))
}
MEs <- stage2results_Lipidomics$modules$MEs

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
  theme_bw() +
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

all(rownames(Lipidomics) == rownames(MEs))
dim(Lipidomics) %>% paste0(c(" samples", " Lipidomics"))
kME <- bicor(Lipidomics, MEs, maxPOutliers = 0.05)
dim(kME) %>% paste0(c(" Lipidomics", " modules"))

intra_cor <- c()
for (i in 1:ncol(Lipidomics)) {
  m <- stage2results_Lipidomics$modules$colors[i]
  intra_cor[i] <- kME[i, paste0("ME", m)]
  if(m != 0){
    intra_cor[i] <- kME[i, paste0("ME", m)]
  } else{
    intra_cor[i] <- NA
  }
  
}

idx <- which(is.na(intra_cor))
intra_cor <- intra_cor[-idx]

plot(density(intra_cor), main = "Correlations with module-eigenLipidomics (within module correlation)\nNo ME0", xlim = c(-1,1))



# Corr within modules
corr_within_module <- function(Lipidomics, modules, module_x = 1){
  idx.Lipidomics <- which(modules$colors == module_x)
  idx.me <- which(colnames(modules$MEs) == paste0("ME",module_x))
  kME_x <- bicor(Lipidomics[,idx.Lipidomics], modules$MEs[,idx.me], maxPOutliers = 0.05)
  kME_x
}

ggplot.list <- list()

for(m in colnames(stage2results_Lipidomics$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(Lipidomics = Lipidomics, modules = stage2results_Lipidomics$modules, module_x = h))) %>% 
    ggplot() + 
    #geom_density(aes(x = x), fill = labels2colors(h), color = "black", alpha = 0.5) + 
    geom_histogram(aes(x), fill = labels2colors(h), color = "black", alpha = 0.5, bins = 20) + 
    xlim(-1, 1) +
    theme_bw()+
    xlab("Lipidomics correlation")+
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









#stage2results_Lipidomics_backup <- stage2results_Lipidomics
stage2results_Lipidomics 
#rownames(stage2results_Lipidomics$modules$MEs) <- gsub("-MC1-TK1-", "_mc1_tk1_", rownames(stage2results_Lipidomics$modules$MEs))

# Remove ME0 from analysis
#idx <- which(colnames(stage2results_Lipidomics$modules$MEs) == "ME0")
#stage2results_Lipidomics$modules$MEs <- stage2results_Lipidomics$modules$MEs[,-idx]


if(take_average){
  stage2results_Lipidomics$modules$MEs %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_VFA.*", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_Lipidomics_eigengenes
  
  stage2results_Lipidomics_eigengenes <- as.data.frame(stage2results_Lipidomics_eigengenes)
} else{
  stage2results_Lipidomics_eigengenes <- stage2results_Lipidomics$modules$MEs
}

rownames(stage2results_Lipidomics_eigengenes) <- stage2results_Lipidomics_eigengenes$common
stage2results_Lipidomics_eigengenes$common <- NULL

dim(stage2results_Lipidomics_eigengenes); dim(stage2results_X_eigengenes)



table(rownames(stage2results_X_eigengenes) %in% rownames(stage2results_Lipidomics_eigengenes))

stage2results_X_eigengenes_New <- stage2results_X_eigengenes[rownames(stage2results_X_eigengenes) %in% rownames(stage2results_Lipidomics_eigengenes), ]


table( rownames(stage2results_X_eigengenes_New) %in% rownames(stage2results_Lipidomics_eigengenes))
table( rownames(stage2results_X_eigengenes_New) == rownames(stage2results_Lipidomics_eigengenes))



#Correlate modules from transcriptomics and metagenomics.

# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(stage2results_Lipidomics_eigengenes) == rownames(stage2results_X_eigengenes_New))
if(!same_order){
  stage2results_Lipidomics_eigengenes <- stage2results_Lipidomics_eigengenes[order(rownames(stage2results_Lipidomics_eigengenes)),]
  stage2results_X_eigengenes_New <- stage2results_X_eigengenes_New[order(rownames(stage2results_X_eigengenes_New)),]
  same_order <- all(rownames(stage2results_Lipidomics_eigengenes) == rownames(stage2results_X_eigengenes_New))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


#####
p.value_matr <- corr.value_matr <- matrix(ncol = ncol(stage2results_Lipidomics_eigengenes), 
                                          nrow = ncol(stage2results_X_eigengenes_New), 
                                          dimnames = list(colnames(stage2results_X_eigengenes_New), 
                                                          colnames(stage2results_Lipidomics_eigengenes)))


for(i in 1:ncol(stage2results_X_eigengenes_New)){
  for(j in 1:ncol(stage2results_Lipidomics_eigengenes)){
    cor.res <- cor.test(stage2results_X_eigengenes_New[,i], stage2results_Lipidomics_eigengenes[,j])
    p.value_matr[i, j] <- cor.res$p.value
    corr.value_matr[i, j] <- cor.res$estimate
  }
}

# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_eigengenes_New), colnames(stage2results_Lipidomics_eigengenes))


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
Lipidomics_corr_X <- list(p_value = p.value_matr, 
                             p_value_adj = p.value_matr.adjust,
                             signif_matrix = signif_matrix,
                             correlation = corr.value_matr)
rm(p.value_matr, p.value_matr.adjust, signif_matrix, corr.value_matr)



heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)


stage2results_Lipidomics$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("Lipidomics_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubLipidomics

hubLipidomics -> hubLipidomics_Before

hubLipidomics_Before$ModuleLipidomics <- paste0("ME" ,hubLipidomics_Before$Module)


test <- Lipidomics_corr_X$correlation
#hubLipidomics_Before$Lipidomics_name <- replace(hubLipidomics_Before$Lipidomics_name, hubLipidomics_Before$Lipidomics_name=="Ile.Leu...244.17825.Da.264.70.s", "Ile-Leu")
#hubLipidomics_Before$Lipidomics_name <- replace(hubLipidomics_Before$Lipidomics_name, hubLipidomics_Before$Lipidomics_name=="Di.N.octyl.phthalate...390.27623.Da.707.55.s", "Di-N-octylphthalate")
#hubLipidomics_Before$Lipidomics_name <- replace(hubLipidomics_Before$Lipidomics_name, hubLipidomics_Before$Lipidomics_name=="X2_.Deoxyadenosine...251.10146.Da.111.43.s", "2-Deoxyadenosine")
#hubLipidomics_Before$Lipidomics_name <- replace(hubLipidomics_Before$Lipidomics_name, hubLipidomics_Before$Lipidomics_name=="Asp.Leu...246.12114.Da.191.61.s", "Asp-Leu")
#colnames(test)[1:2] <- c("D-Arabinonic acid", "Thr-Thr")

match(hubLipidomics_Before[, "ModuleLipidomics"], colnames(test))
colnames(test)[match(hubLipidomics_Before[,"ModuleLipidomics"], colnames(test))] = hubLipidomics_Before[,"Lipidomics_name"]

Lipidomics_corr_X$correlation <- test
pheatmap::pheatmap(Lipidomics_corr_X$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,  # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = Lipidomics_corr_X$signif_matrix, 
                   fontsize_number = 8, #10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F,
                   show_rownames = F,
                   labels_row = paste0(prefix_OTUs, rownames(Lipidomics_corr_X$correlation)),
                   labels_col = paste0(colnames(Lipidomics_corr_X$correlation)),
                   #main = "Eigen Lipidomics"
                   ) -> Lipidomics_corr_X_plot








#########
library(tidyverse)
library(phyloseq)
library(microbiomeMarker)
raw <- as.matrix(Lipidomics)
OTU = otu_table(raw, taxa_are_rows = TRUE)
dat <- read.csv("/Users/shashankgupta/Desktop/ImprovAFish/Metabolomics/Lipidomics/Lipidomics Data/Metadata.csv")
row.names(dat) <- dat$Sample
dat$New_Diet <- toupper(dat$Diet)
# Merge into one complete phyloseq object
ps <- merge_phyloseq(otu_table(OTU), sample_data(dat))
tt <- as.data.frame(row.names(Lipidomics))
row.names(tt) <- tt$`row.names(Lipidomics)`
colnames(tt)[1] <- "Kingdom"
tax_table(ps) <- as.matrix(tt)

#T2
ps_MC1<-subset_samples(ps, Time %in% c("T2") & New_Diet %in% c("CTR", "MC1"))
ps_MC1 <- prune_taxa(taxa_sums(ps_MC1) >0, ps_MC1)
table(sample_data(ps_MC1)$New_Diet)


lef_out<-run_lefse(ps_MC1, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 2, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)

res <- ldamarker(ps_MC1,group="New_Diet")
plotLDA(res,group=c("MC1","CTR"),lda=0, padj = 0.05, fontsize.y = 6)



ps_MC2<-subset_samples(ps, Time %in% c("T2") & New_Diet %in% c("CTR", "MC2"))
ps_MC2 <- prune_taxa(taxa_sums(ps_MC2) >0, ps_MC2)
table(sample_data(ps_MC2)$New_Diet)


lef_out<-run_lefse(ps_MC2, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)

res <- ldamarker(ps_MC2, group="New_Diet")
plotLDA(res,group=c("MC2","CTR"), lda=0, padj = 0.05, fontsize.y = 6)


ps_MN3<-subset_samples(ps, Time %in% c("T2") & New_Diet %in% c("CTR", "MN3"))
ps_MN3 <- prune_taxa(taxa_sums(ps_MN3) >0, ps_MN3)
table(sample_data(ps_MN3)$New_Diet)


lef_out<-run_lefse(ps_MN3, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)

res <- ldamarker(ps_MN3, group="New_Diet")
plotLDA(res,group=c("MN3","CTR"), lda=0, padj = 0.05, fontsize.y = 6)




#T3
ps_MC1<-subset_samples(ps, Time %in% c("T3") & New_Diet %in% c("CTR", "MC1"))
ps_MC1 <- prune_taxa(taxa_sums(ps_MC1) >0, ps_MC1)
table(sample_data(ps_MC1)$New_Diet)


lef_out<-run_lefse(ps_MC1, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 2, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)

res <- ldamarker(ps_MC1,group="New_Diet")
plotLDA(res,group=c("MC1","CTR"),lda=0, padj = 0.05, fontsize.y = 6)



ps_MC2<-subset_samples(ps, Time %in% c("T3") & New_Diet %in% c("CTR", "MC2"))
ps_MC2 <- prune_taxa(taxa_sums(ps_MC2) >0, ps_MC2)
table(sample_data(ps_MC2)$New_Diet)


lef_out<-run_lefse(ps_MC2, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 2, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)

res <- ldamarker(ps_MC2, group="New_Diet")
plotLDA(res,group=c("MC2","CTR"), lda=2, padj = 0.05, fontsize.y = 6)

nrow(subset(res, res$direction ==  "MC2" & res$p.adj < 0.05))

ps_MN3<-subset_samples(ps, Time %in% c("T3") & New_Diet %in% c("CTR", "MN3"))
ps_MN3 <- prune_taxa(taxa_sums(ps_MN3) >0, ps_MN3)
table(sample_data(ps_MN3)$New_Diet)


lef_out<-run_lefse(ps_MN3, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 2, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)


res <- ldamarker(ps_MN3, group="New_Diet")
plotLDA(res,group=c("MN3","CTR"), lda=2, padj = 0.05, fontsize.y = 6)

nrow(subset(res, res$direction ==  "MN3" & res$p.adj < 0.05))


#loop
output_table <- data.frame()

for (i in c("T2", "T3")){
  for (j in c("MC1","MC2", "MN3")){
    phy<-subset_samples(ps, Time %in% c(i) & New_Diet %in% c("CTR", j))
    phy <- prune_taxa(taxa_sums(phy) >0, phy)
    lef_out<-run_lefse(phy,group = "New_Diet", norm = "CPM", 
                       kw_cutoff = 0.05, lda_cutoff = 2, taxa_rank = "none")
    assign(paste0("Lipid", i, j), plot_abundance(lef_out, group = "New_Diet"))
    output_table[i,j]<-sum(marker_table(lef_out)$enrich_group == j)
  }
}
#colnames(output_table)<-c("MC1", "MC2", "MN3")
#rownames(output_table) <- c("T1","T2","T3")
t(output_table)

plot_grid(LipidT2MC1, LipidT2MC2, LipidT2MN3, LipidT3MC1, LipidT3MC2, LipidT3MN3, labels = "AUTO", ncol = 1, align = "hv")
# portrait 25 *10
#plot_grid(poteomicsT3MC1, poteomicsT3MC2, poteomicsT3MN3, labels = "AUTO", ncol = 1, align = "hv")























MEs = moduleEigengenes(Lipidomics_host, moduleLabels1)$eigengenes

MEs <- orderMEs(MEs)
module_order = names(MEs) %>% gsub("ME","", .)
MEs0 <- MEs
MEs0$Sample.ID = row.names(MEs)

mME = MEs0 %>%
  pivot_longer(-Sample.ID) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )



library(tidyverse)
setwd("/Users/shashankgupta/Desktop/ImprovAFish/ImprovAFish")
sample_info3 <- read.csv("sample_info.csv", row.names = 1)
sample_info3$Sample.ID <- sample_info3$ID_New
mME_meta <- merge(mME, sample_info3, by = "Sample.ID") 
colnames(mME_meta)[2] <-"Module"

mME_meta$New_Diet <- factor(mME_meta$New_Diet, levels = c("ctr", "mc1", "mc2", "mn3"))
mME_meta$Module <- factor(mME_meta$Module, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                      "10", "11", "12", "13", "14", "15", "16", "17",
                                                      "18", "19", "20","21","22","23", "24"))
expression_plots<-mME_meta%>%
  #  group_by(Module) %>%
  ggplot(aes(x=Time, y=value, fill=New_Diet)) +
  facet_wrap(~ Module)+
  ylab("Mean Module Eigenegene Value") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  stat_summary(fun=mean, geom="line", aes(New_Diet=New_Diet, color = New_Diet), position = position_dodge(width = 0.5))  + 
  geom_point(pch = 21, position = position_dodge(width = 0.5)) +
  scale_fill_manual(name="Lifestage", values=c("#8C510A", "#DFC27D","#80CDC1", "#003C30", "#BA55D3")) +
  scale_color_manual(name="Lifestage", values=c("#8C510A", "#DFC27D","#80CDC1", "#003C30", "#BA55D3")) + 
  #xlab("Hours Post-Fertilization") + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 11 , color = "black"),
        axis.title = element_text(size = 16, color = "black"), 
        axis.text.x = element_text(size=11, color="black"), 
        legend.title=element_blank(), 
        legend.text=element_text(color="black", size=12)); expression_plots


lifestage








