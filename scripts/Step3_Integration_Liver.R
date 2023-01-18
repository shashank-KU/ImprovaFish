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
#stage2results_X_Liver_backup <- stage2results_X_Liver
#stage2results_X_Liver <- stage2results_X_Liver_backup
################################################################################
stage2results_X_Liver 


# Remove ME0 from analysis
idx <- which(colnames(stage2results_X_Liver$modules$MEs) == "ME0")
stage2results_X_Liver$modules$MEs <- stage2results_X_Liver$modules$MEs[,-idx]
filter.sig.mod <- 2


take_average <- TRUE

if(take_average){
  stage2results_X_Liver$modules$MEs %>% 
    rownames_to_column(var = "common") %>%
    mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_X_Liver_eigengenes
  
  stage2results_X_Liver_eigengenes <- as.data.frame(stage2results_X_Liver_eigengenes)
  
} else{
  stage2results_X_Liver_eigengenes <- stage2results_X_Liver$modules$MEs
}

rownames(stage2results_X_Liver_eigengenes) <- stage2results_X_Liver_eigengenes$common
stage2results_X_Liver_eigengenes$common <- NULL

# Create a dendrogram of the transcriptomics eigengenes to organise the final plots.
X_ME_dendro <- hclust(as.dist(1 - WGCNA::bicor(stage2results_X_Liver_eigengenes, maxPOutliers = 0.05)), method = "ward.D2")

# "..palettes that vary from one hue to another via white may have a more symmetrical appearance in RGB space"
if(take_average){
  heatmap_colors <- colorRampPalette(c("#18b29f","#FFFFFF","#ac6721"), interpolate = "spline", space = "rgb")(51)
} else{
  heatmap_colors <- colorRampPalette(c("#18b29f","#FFFFFF","#ac6721"), interpolate = "spline", space = "rgb")(51)
}


annotation_col <- 
  rownames(stage2results_X_Liver_eigengenes) %>% 
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

colnames(annotation_col) <- c("Time" , "Diet")

annotation_colors <- list(
  Time = c(`T1` = "navy", `T2` = "darkgreen", `T3` = "red"),
  Diet = c(`ctr` = "#1D88AF", `mc1` = "darkblue", `mc2` = "darkgreen", `mn3` = "#F08A46")
  #Tank_number = c(`tk1` = "#F08A46", `tk2` = "#8EB470",`tk3` = "#B7CFA4"),
  #Sample_Number =c(`01` = "#1D88AF", `03` = "#12556D", `05` = "#F6B890", `07` = "#F08A46") 
)



stage2results_X_Liver_eigengenes_to_plot <- 
  dplyr::inner_join(annotation_col %>% 
                      rownames_to_column(var = "common"), 
                    stage2results_X_Liver_eigengenes %>% 
                      rownames_to_column(var = "common"), 
                    by = "common") %>%
  #dplyr::arrange(Feed, Water, Day) %>%              # The order at which the columns should appear, given that there is no clustering.
  dplyr::select(common, starts_with("ME")) %>% 
  tibble::column_to_rownames(var = "common") %>% 
  t()
library(WGCNA)
moduleLabels = stage2results_X_Liver$modules$colors
moduleColors = labels2colors(stage2results_X_Liver$modules$colors)
table(moduleLabels)
table(moduleColors)
table(moduleLabels, moduleColors)
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)


pheatmap::pheatmap(stage2results_X_Liver_eigengenes_to_plot, 
                   cluster_cols = F,
                   cluster_rows = X_ME_dendro,
                   treeheight_row = 20,
                   #cluster_rows = T,
                   cutree_rows = row_cut,
                   cutree_cols = col_cut,
                   color = heatmap_colors,
                   fontsize = 10,
                   fontsize_col = ifelse(take_average, 10, 6),
                   fontsize_row = 8,
                   annotation_colors = annotation_colors,
                   annotation_col = annotation_col, 
                   silent = F,annotation_legend = F,
                   labels_row = paste0(prefix_genes, rownames(stage2results_X_Liver_eigengenes_to_plot)),
                   main = paste("Gene Module 'expression'\n", ifelse(take_average, "mean values", ""))) -> X_plot





#stage2results_Y_backup <- stage2results_Y
stage2results_Y

# Remove ME0 from analysis
idx <- which(colnames(stage2results_Y$modules$MEs) == "ME0")
stage2results_Y$modules$MEs <- stage2results_Y$modules$MEs[,-idx]

take_average <- T
if(take_average){
  stage2results_Y$modules$MEs %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_Y_eigengenes
  
  stage2results_Y_eigengenes <- as.data.frame(stage2results_Y_eigengenes)
} else{
  stage2results_Y_eigengenes <- stage2results_Y$modules$MEs
}

rownames(stage2results_Y_eigengenes) <- stage2results_Y_eigengenes$common
stage2results_Y_eigengenes$common <- NULL

dim(stage2results_Y_eigengenes); dim(stage2results_X_Liver_eigengenes)


#Correlate modules from transcriptomics and metagenomics.

# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_Liver_eigengenes))
if(!same_order){
  stage2results_Y_eigengenes <- stage2results_Y_eigengenes[order(rownames(stage2results_Y_eigengenes)),]
  stage2results_X_Liver_eigengenes <- stage2results_X_Liver_eigengenes[order(rownames(stage2results_X_Liver_eigengenes)),]
  same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_Liver_eigengenes))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


#####
p.value_matr <- corr.value_matr <- matrix(ncol = ncol(stage2results_Y_eigengenes), 
                                          nrow = ncol(stage2results_X_Liver_eigengenes), 
                                          dimnames = list(colnames(stage2results_X_Liver_eigengenes), 
                                                          colnames(stage2results_Y_eigengenes)))


for(i in 1:ncol(stage2results_X_Liver_eigengenes)){
  for(j in 1:ncol(stage2results_Y_eigengenes)){
    cor.res <- cor.test(stage2results_X_Liver_eigengenes[,i], stage2results_Y_eigengenes[,j])
    p.value_matr[i, j] <- cor.res$p.value
    corr.value_matr[i, j] <- cor.res$estimate
  }
}

# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_Liver_eigengenes), colnames(stage2results_Y_eigengenes))


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


stage2results_Y$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("OTU_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubOTUs


otu_taxonomy <- data.frame(tax_table(psdata.p))
dplyr::left_join(hubOTUs, 
                 (otu_taxonomy %>%
                    tibble::rownames_to_column(var = "OTU_name")), 
                 by = "OTU_name") -> hubOTUs_Before

hubOTUs_Before$ModuleOTU <- paste0("ME" ,hubOTUs_Before$Module)


test <- Y_corr_X$correlation
hubOTUs_Before <- hubOTUs_Before[, c(10, 8)]
hubOTUs_Before$Genus <- replace(hubOTUs_Before$Genus, hubOTUs_Before$Genus=="Allorhizobium_Neorhizobium_Pararhizobium_Rhizobium", "Rhizobium")
hubOTUs_Before$Genus <- replace(hubOTUs_Before$Genus, hubOTUs_Before$Genus=="Burkholderia_Caballeronia_Paraburkholderia", "Paraburkholderia")


match(hubOTUs_Before[, "ModuleOTU"], colnames(test))
colnames(test)[match(hubOTUs_Before[,"ModuleOTU"], colnames(test))] = hubOTUs_Before[,"Genus"]

Y_corr_X$correlation <- test
#Y_corr_X$correlation <- Y_corr_X$correlation[, order(match(colnames(Y_corr_X$correlation), c("Photobacterium",	"Brachybacterium",	"HOC36",	"Kingdom_Bacteria",	"Lactobacillus",	"Class_Dehalococcoidia",	"Paraburkholderia",	"Teleostei",	"Oceanisphaera"))) ]

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


#####################################################################################
#####################################################################################
#####################################################################################
#Correlate transcriptomics with traits
# simplify by chosing only the columns we need.
common <- sample_info$common
sample_info <- as_tibble(sample_info)
sample_info %>% dplyr::select(common, samplingTime, New_Diet) -> Z_traits
head(Z_traits)


Z_traits %>% 
  mutate(samplingTime = factor(samplingTime, levels = c("T1", "T2", "T3"))) %>% 
  mutate(New_Diet = factor(New_Diet)) %>%
  mutate(Day = ifelse(samplingTime == "T1", paste0("T1_", New_Diet),
                      ifelse(samplingTime == "T2", paste0("T2_", New_Diet), paste0("T3_", New_Diet)))) -> Z_traits
head(Z_traits)

Z_traits$Day<- factor(Z_traits$Day)
Z_traits$samplingTime<- factor(Z_traits$samplingTime)
levels(Z_traits$samplingTime)
levels(Z_traits$Day)
#Z_traits$common <- NULL
#One hot encoding the nominal variables.
rownames(Z_traits) <- Z_traits$common
Z_traits %>%
  dplyr::select(one_of("Day", "New_Diet", "samplingTime" )) %>% 
  #column_to_rownames(var = "common") %>% 
  model.matrix(~.-1, data = .) %>% 
  as.data.frame() %>% 
  #rownames_to_column(var = "common") %>% 
  as_tibble() -> nominal_variables

head(nominal_variables)



Z_traits %>% 
  dplyr::select(one_of("Day", "New_Diet", "samplingTime" , "common")) -> ordinal_and_other_variables
ordinal_and_other_variables


#ordinal_and_other_variables$common <- row.names(ordinal_and_other_variables)

Z_traits <- full_join(nominal_variables, ordinal_and_other_variables, by = "common")
head(Z_traits)
dim(Z_traits)
Z_traits$samplingTimeT1 <- ifelse(Z_traits$samplingTime == "T1", 1, 0)



if(take_average){
  Z_traits %>% 
    mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean), na.rm = TRUE) %>%  # Remove the NA values 
    ungroup() %>% 
    #dplyr::select(common, c(starts_with("New_Diet")))  %>% 
    tibble::column_to_rownames(var = "common")  -> Z_traits
  
  Z_traits <- as.data.frame(Z_traits)
  
} 

Z_traits %>% rmarkdown::paged_table()

Z_traits$Day <- NULL
Z_traits$New_Diet <- NULL
Z_traits$samplingTime <- NULL
Z_traits <- Z_traits[, c(1:15, 18, 16, 17)]
#####################################################################################
#####################################################################################
#####################################################################################


# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(Z_traits) == rownames(stage2results_X_Liver_eigengenes))
if(!same_order){
  stage2results_Y_eigengenes <- stage2results_Y_eigengenes[order(rownames(stage2results_Y_eigengenes)),]
  stage2results_X_Liver_eigengenes <- stage2results_X_Liver_eigengenes[order(rownames(stage2results_X_Liver_eigengenes)),]
  same_order <- all(rownames(stage2results_Y_eigengenes) == rownames(stage2results_X_Liver_eigengenes))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


p.value_matr <- corr.value_matr <- matrix(ncol = ncol(Z_traits), 
                                          nrow = ncol(stage2results_X_Liver_eigengenes), 
                                          dimnames = list(colnames(stage2results_X_Liver_eigengenes), 
                                                          colnames(Z_traits)))
for(i in 1:ncol(stage2results_X_Liver_eigengenes)){
  for(j in 1:ncol(Z_traits)){
    
    if(colnames(Z_traits)[j] == "Day"){
      cor.res <- cor.test(stage2results_X_Liver_eigengenes[,i], Z_traits[,j], method = "spearman", exact = FALSE)
      p.value_matr[i, j] <- cor.res$p.value
      corr.value_matr[i, j] <- cor.res$estimate
    } else{
      cor.res <- cor.test(stage2results_X_Liver_eigengenes[,i], Z_traits[,j])
      p.value_matr[i, j] <- cor.res$p.value
      corr.value_matr[i, j] <- cor.res$estimate
    }
  }
}


# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_Liver_eigengenes), colnames(Z_traits))


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
                   treeheight_col = 0, 
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


# Modules with no significant associations
if (filter.sig.mod == 1) {
  idx1 <- rowSums(#Z_corr_X$signif_matrix == "*" |
    Z_corr_X$signif_matrix == "**" |
      Z_corr_X$signif_matrix == "***") > 0
  idx2 <- rowSums(#Y_corr_X$signif_matrix == "*" |
    Y_corr_X$signif_matrix == "**" | 
      Y_corr_X$signif_matrix == "***") > 0
  sig.modules.h <- rownames(Y_corr_X$p_value)[idx1 & idx2]
  
  sig.modules.m <- colnames(Y_corr_X$p_value)[colSums(#Y_corr_X$signif_matrix[idx1 & idx2,] == "*" |
    Y_corr_X$signif_matrix[idx1 & idx2,] == "**" | 
      Y_corr_X$signif_matrix[idx1 & idx2,] == "***") > 0]
  
  
  save(sig.modules.h, sig.modules.m, file = "../data/module_data/sig_modules.Rdata")
}

cowplot::plot_grid(Z_corr_X_plot$gtable,
                   X_plot$gtable,
                   Y_corr_X_plot$gtable, 
                   ncol = 3, 
                   rel_widths = c(dim(Z_traits)[2]/3, 
                                  ifelse(take_average,dim(stage2results_X_Liver_eigengenes)[1]/2,dim(stage2results_X_Liver_eigengenes)[1]/6),
                                  dim(stage2results_Y_eigengenes)[2]/2.2),
                   align = "h") + ggplot2::theme(plot.margin = ggplot2::unit(c(3,0,2.5,1), "cm")) # Top, Right, Bottom, Left





########################################################################
########################################################################
# module–trait associations

# Define numbers of genes and samples

nGenes = ncol(omics_data_host_Liver);
nSamples = nrow(omics_data_host_Liver);
moduleLabels1 = stage2results_X_Liver$modules$colors
moduleColors1 = labels2colors(stage2results_X_Liver$modules$colors)
MEs = stage2results_X_Liver$MEs;

table(moduleColors1)
table(moduleLabels1)

table(moduleLabels1, moduleColors1)



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(omics_data_host_Liver, moduleLabels1)$eigengenes
#MEs = orderMEs(MEs0)
head(MEs0)
MEs <- MEs0[, c(2:44)]
head(MEs)
#colnames(MEs)<- paste("h", colnames(MEs), sep = "")
head(MEs)

sample_info1 <- host_gut_mapping_file

sample_info1[, c("Time", "Diet","New_Diet","Tank_number", "sampleName", "Sample_Number", "Sample_Type","ID")] <- list(NULL)
row.names(sample_info1) <- sample_info1$common
sample_info1$common <- NULL

reorder_idx <- match(row.names(MEs), row.names(sample_info1))
sample_info1 <- sample_info1[ reorder_idx , ]
all(row.names(sample_info1) %in% row.names(MEs))
all(row.names(sample_info1) == row.names(MEs))

#MEs <- subset(MEs, select = -c(MEgrey))
moduleTraitCor = cor(MEs, sample_info1, use = "p")
moduleTraitCor <- moduleTraitCor[order(match(rownames(moduleTraitCor), c("ME15",	"ME25",	"ME2",	"ME32",	"ME18",	"ME10",	"ME13",	"ME40",	"ME31",	"ME17",	"ME16",	"ME20",	"ME37",	"ME7",	"ME5",	"ME35",	"ME38",	"ME3",	"ME23",	"ME22",	"ME34",	"ME33",	"ME28",	"ME9",	"ME29",	"ME39",	"ME1",	"ME8",	"ME6",	"ME36",	"ME14",	"ME43",	"ME11",	"ME4",	"ME41",	"ME26",	"ME24",	"ME19",	"ME30",	"ME42",	"ME27",	"ME12",	"ME21"))), ]


moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8, 1, 1))
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
                   rel_widths = c(1.2, 1.5, 3.5, 1.1) # PNG 1200 (Width) * 800 (Height)  # pdf 14 * 8 portrait
)+ ggplot2::theme(plot.margin = ggplot2::unit(c(3,0,2.5,1.5), "cm"))



















































