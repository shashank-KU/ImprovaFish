library(magrittr) # For pipes %>% 
library(dplyr)    # For data frame manipulations
library(tibble)   # For rownames_to_columns etc
library(tidyr)    # For functions: Separate, ...
library(textshape)

take_average = T
cluster_the_columns = F
col_cut <- 1
row_cut <- 1
prefix_genes <- "h"
prefix_OTUs <- "m"


#Transcriptomics
#stage2results_X_PC_corrected_backup <- stage2results_X_PC_corrected  
################################################################################
stage2results_X_PC_corrected 

# Remove ME0 from analysis
idx <- which(colnames(stage2results_X_PC_corrected$modules$MEs) == "ME0")
stage2results_X_PC_corrected$modules$MEs <- stage2results_X_PC_corrected$modules$MEs[,-idx]
filter.sig.mod <- 2




if(take_average){
  stage2results_X_PC_corrected$modules$MEs %>% 
    rownames_to_column(var = "common") %>%
    mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_X_eigengenes
  
  stage2results_X_eigengenes <- as.data.frame(stage2results_X_eigengenes)
  
} else{
  stage2results_X_eigengenes <- stage2results_X_PC_corrected$modules$MEs
}

rownames(stage2results_X_eigengenes) <- stage2results_X_eigengenes$common
stage2results_X_eigengenes$common <- NULL

# Create a dendrogram of the transcriptomics eigengenes to organise the final plots.
X_ME_dendro <- hclust(as.dist(1 - WGCNA::bicor(stage2results_X_eigengenes, maxPOutliers = 0.05)), method = "ward.D2")

# "..palettes that vary from one hue to another via white may have a more symmetrical appearance in RGB space"
if(take_average){
  heatmap_colors <- colorRampPalette(c("#18b29f","#FFFFFF","#ac6721"), interpolate = "spline", space = "rgb")(51)
} else{
  heatmap_colors <- colorRampPalette(c("#18b29f","#FFFFFF","#ac6721"), interpolate = "spline", space = "rgb")(51)
}

heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)

annotation_col <- 
  rownames(stage2results_X_eigengenes) %>% 
  as.data.frame() %>% 
  dplyr::rename(common = ".") %>% 
  tidyr::separate(col = common,  into = c("Time", "New_Diet"), remove = F, extra = "drop") %>% # If replicates are present, we drop them here
  #dplyr::mutate(Feed = factor(Feed, levels = c("FO", "VOFO", "FOVO", "VO"))) %>% 
  #dplyr::mutate(Day = sub("D","", Day)) %>% 
  #dplyr::mutate(Days_total = ifelse(Water == "FW", paste0("FW", Day), paste0("SW", Day))) %>%
  #dplyr::mutate(Days_total = factor(Days_total, levels= c(paste0("FW", 0:20), paste0("SW", 0:20)))) %>%
  #dplyr::mutate(Body_Weight = as.numeric(BW), Length = as.numeric(Length)) %>% 
  dplyr::select(common,"Time", "New_Diet") %>%  
  # Order here is reflected in the order of annotations, Days_total = new number for each unique day
  tibble::column_to_rownames(var = "common")

colnames(annotation_col) <- c("Time", "Diet")


annotation_colors <- list(
  Time = c(`T1` = "navy", `T2` = "darkgreen", `T3` = "red"),
  Diet = c(`ctr` = "#1D88AF", `mc1` = "darkblue", `mc2` = "darkgreen", `mn3` = "#F08A46")
  #Tank_number = c(`tk1` = "#F08A46", `tk2` = "#8EB470",`tk3` = "#B7CFA4"),
  #Sample_Number =c(`01` = "#1D88AF", `03` = "#12556D", `05` = "#F6B890", `07` = "#F08A46") 
  )


moduleLabels = stage2results_X_PC_corrected$modules$colors
moduleColors = labels2colors(stage2results_X_PC_corrected$modules$colors)

stage2results_X_eigengenes_to_plot <- 
  dplyr::inner_join(annotation_col %>% 
                      rownames_to_column(var = "common"), 
                    stage2results_X_eigengenes %>% 
                      rownames_to_column(var = "common"), 
                    by = "common") %>%
  #dplyr::arrange(Feed, Water, Day) %>%              # The order at which the columns should appear, given that there is no clustering.
  dplyr::select(common, starts_with("ME")) %>% 
  tibble::column_to_rownames(var = "common") %>% 
  t()


pheatmap::pheatmap(stage2results_X_eigengenes_to_plot, 
                   cluster_cols = F,
                   cluster_rows = X_ME_dendro,
                   treeheight_row = 20,
                   cutree_rows = row_cut,
                   cutree_cols = col_cut,
                   color = heatmap_colors,
                   fontsize = 10,
                   fontsize_col = ifelse(take_average, 10, 6),
                   fontsize_row = 8,
                   annotation_colors = annotation_colors,
                   annotation_col = annotation_col, 
                   silent = F, annotation_legend = F,
                   labels_row = paste0(prefix_genes, rownames(stage2results_X_eigengenes_to_plot)),
                   main = paste("Gene Module 'expression'\n", ifelse(take_average, "mean values", ""))) -> X_plot





#stage2results_Y_PC_corrected_backup <- stage2results_Y_PC_corrected
stage2results_Y_PC_corrected

# Remove ME0 from analysis
idx <- which(colnames(stage2results_Y_PC_corrected$modules$MEs) == "ME0")
stage2results_Y_PC_corrected$modules$MEs <- stage2results_Y_PC_corrected$modules$MEs[,-idx]


if(take_average){
  stage2results_Y_PC_corrected$modules$MEs %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_Y_eigengenes
  
  stage2results_Y_eigengenes <- as.data.frame(stage2results_Y_eigengenes)
} else{
  stage2results_Y_eigengenes <- stage2results_Y_PC_corrected$modules$MEs
}

rownames(stage2results_Y_eigengenes) <- stage2results_Y_eigengenes$common
stage2results_Y_eigengenes$common <- NULL

dim(stage2results_Y_eigengenes); dim(stage2results_X_eigengenes)


#Correlate modules from transcriptomics and metagenomics.

# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_eigengenes))
if(!same_order){
  stage2results_Y_eigengenes <- stage2results_Y_eigengenes[order(rownames(stage2results_Y_eigengenes)),]
  stage2results_X_eigengenes <- stage2results_X_eigengenes[order(rownames(stage2results_X_eigengenes)),]
  same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_eigengenes))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


#####
p.value_matr <- corr.value_matr <- matrix(ncol = ncol(stage2results_Y_eigengenes), 
                                          nrow = ncol(stage2results_X_eigengenes), 
                                          dimnames = list(colnames(stage2results_X_eigengenes), 
                                                          colnames(stage2results_Y_eigengenes)))


for(i in 1:ncol(stage2results_X_eigengenes)){
  for(j in 1:ncol(stage2results_Y_eigengenes)){
    cor.res <- cor.test(stage2results_X_eigengenes[,i], stage2results_Y_eigengenes[,j])
    p.value_matr[i, j] <- cor.res$p.value
    corr.value_matr[i, j] <- cor.res$estimate
  }
}

# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_eigengenes), colnames(stage2results_Y_eigengenes))


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
Y_corr_X <- list(p_value = p.value_matr, 
                 p_value_adj = p.value_matr.adjust,
                 signif_matrix = signif_matrix,
                 correlation = corr.value_matr)
rm(p.value_matr, p.value_matr.adjust, signif_matrix, corr.value_matr)



heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)

pheatmap::pheatmap(Y_corr_X$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,  # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = Y_corr_X$signif_matrix, 
                   fontsize_number = 8, #10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F,
                   show_rownames = F,
                   labels_row = paste0(prefix_OTUs, rownames(Y_corr_X$correlation)),
                   labels_col = paste0(prefix_OTUs, colnames(Y_corr_X$correlation)),
                   main = "EigenASVs") -> Y_corr_X_plot


stage2results_Y_PC_corrected$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("OTU_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubOTUs


otu_taxonomy <- data.frame(tax_table(psdata.p))
dplyr::left_join(hubOTUs, 
                 (otu_taxonomy %>%
                    tibble::rownames_to_column(var = "OTU_name")), 
                 by = "OTU_name") -> hubOTUs_Before

hubOTUs_Before$ModuleOTU <- paste0("ME" ,hubOTUs_Before$Module)
head(hubOTUs_Before)

test <- Y_corr_X$correlation
hubOTUs_Before <- hubOTUs_Before[, c(10, 8)]

match(hubOTUs_Before[, "ModuleOTU"], colnames(test))
colnames(test)[match(hubOTUs_Before[,"ModuleOTU"], colnames(test))] = hubOTUs_Before[,"Genus"]

Y_corr_X$correlation <- test

pheatmap::pheatmap(Y_corr_X$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,  # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = Y_corr_X$signif_matrix, 
                   fontsize_number = 8, #10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F,
                   show_rownames = F,
                   labels_row = paste0(prefix_OTUs, rownames(Y_corr_X$correlation)),
                   labels_col = paste0(colnames(Y_corr_X$correlation)),
                   main = "Eigen ASVs") -> Y_corr_X_plot



# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(Z_traits) == rownames(stage2results_X_eigengenes))
if(!same_order){
  stage2results_Y_eigengenes <- stage2results_Y_eigengenes[order(rownames(stage2results_Y_eigengenes)),]
  stage2results_X_eigengenes <- stage2results_X_eigengenes[order(rownames(stage2results_X_eigengenes)),]
  same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_eigengenes))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


p.value_matr <- corr.value_matr <- matrix(ncol = ncol(Z_traits), 
                                          nrow = ncol(stage2results_X_eigengenes), 
                                          dimnames = list(colnames(stage2results_X_eigengenes), 
                                                          colnames(Z_traits)))
for(i in 1:ncol(stage2results_X_eigengenes)){
  for(j in 1:ncol(Z_traits)){
    
    if(colnames(Z_traits)[j] == "Day"){
      cor.res <- cor.test(stage2results_X_eigengenes[,i], Z_traits[,j], method = "spearman", exact = FALSE)
      p.value_matr[i, j] <- cor.res$p.value
      corr.value_matr[i, j] <- cor.res$estimate
    } else{
      cor.res <- cor.test(stage2results_X_eigengenes[,i], Z_traits[,j])
      p.value_matr[i, j] <- cor.res$p.value
      corr.value_matr[i, j] <- cor.res$estimate
    }
  }
}


# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_eigengenes), colnames(Z_traits))


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
Z_corr_X <- list(p_value = p.value_matr, 
                 p_value_adj = p.value_matr.adjust,
                 signif_matrix = signif_matrix,
                 correlation = corr.value_matr)
rm(p.value_matr, p.value_matr.adjust, signif_matrix, corr.value_matr)


heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)

pheatmap::pheatmap(Z_corr_X$correlation, 
                   color = heatmap_colors, 
                   
                   cluster_cols = F,
                   #treeheight_col = 0, 
                   treeheight_row = 0, # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = Z_corr_X$signif_matrix, 
                   fontsize_number = 8, # 10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F, 
                   legend = F,
                   show_rownames = F,
                   main = "Traits and external variables") -> Z_corr_X_plot



########################################################################
########################################################################
# module–trait associations

# Define numbers of genes and samples

nGenes = ncol(omics_data_host);
nSamples = nrow(omics_data_host);
moduleLabels1 = stage2results_X_PC_corrected$modules$colors
moduleColors1 = labels2colors(stage2results_X_PC_corrected$modules$colors)
MEs = stage2results_X_PC_corrected$MEs;

table(moduleColors1)
table(moduleLabels1)

table(moduleLabels1, moduleColors1)


# Recalculate MEs with color labels
MEs0 = moduleEigengenes(omics_data_host, moduleLabels1)$eigengenes
#MEs = orderMEs(MEs0)
head(MEs0)
MEs <- MEs0[, c(2:19)]
head(MEs)
#colnames(MEs)<- paste("h", colnames(MEs), sep = "")
head(MEs)

sample_info<- host_gut_mapping_file

sample_info1 <- sample_info
sample_info1$Time = NULL
sample_info1$Diet = NULL
sample_info1$New_Diet = NULL
sample_info1$Tank_number = NULL
sample_info1$Sample_Type = NULL
sample_info1$New_Diet = NULL
sample_info1$ID = NULL
sample_info1$sampleName = NULL
sample_info1$Sample_Number = NULL
row.names(sample_info1) <- sample_info1$common
sample_info1$common <- NULL

#MEs <- subset(MEs, select = -c(MEgrey))

reorder_idx <- match(row.names(MEs), row.names(sample_info1))
sample_info1 <- sample_info1[ reorder_idx , ]
all(row.names(sample_info1) %in% row.names(MEs))
all(row.names(sample_info1) == row.names(MEs))


moduleTraitCor = cor(MEs, sample_info1, use = "p")
moduleTraitCor <- moduleTraitCor[order(match(rownames(moduleTraitCor), c("ME1",	"ME9",	"ME4",	"ME12",	"ME16",	"ME3",	"ME14",	"ME15",	"ME8",	"ME17",	"ME6",	"ME18",	"ME2",	"ME7",	"ME10",	"ME5",	"ME11",	"ME13"))), ]


moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sample_info1),
               yLabels = rownames(moduleTraitCor),
               ySymbols =rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))






signif_matrix <- rep("", length(moduleTraitPvalue))
three_star <- which( moduleTraitPvalue <= 0.001)
signif_matrix[three_star] <- "***"
two_star <- which((moduleTraitPvalue <= 0.01) & (moduleTraitPvalue > 0.001))
signif_matrix[two_star] <- "**"
one_star <- which((moduleTraitPvalue <= 0.05) & (moduleTraitPvalue > 0.01))
signif_matrix[one_star] <- "*"
dim(signif_matrix) = dim(moduleTraitPvalue) # Give textMatrix the correct dimensions 


pheatmap::pheatmap(moduleTraitCor, 
                   color = heatmap_colors,
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   cluster_cols = F,
                   #treeheight_col = 0, 
                   #treeheight_row = 0, # will be shown on the transcriptomics ME heatmap
                   cluster_rows = F,
                   #cutree_rows = row_cut,
                   display_numbers = signif_matrix, 
                   fontsize_number = 8, # 10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F, 
                   legend = F,
                   show_rownames = F,
                   main = "Traits and external variables") -> MT_plot



cowplot::plot_grid(MT_plot$gtable,
                   Z_corr_X_plot$gtable,
                   X_plot$gtable,
                   Y_corr_X_plot$gtable,
                   ncol = 4,
                   align = 'h',
                   rel_widths = c(1.4, 1.8, 4.3, 1.0) # PNG 1700 (Width) * 1000 (Height)  # pdf 16 * 10 portrait
)+ ggplot2::theme(plot.margin = ggplot2::unit(c(3,0,2.5,1), "cm"))











































