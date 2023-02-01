library("MetaboAnalystR")
getwd()
mSet<-InitDataObjects("conc", "mf", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, "/Users/shashankgupta/Desktop/ImprovAFish/Proteomics-meta/metaP_FpTopN_imputedValues.csv", "colmf")
mSet<-ReadMetaData(mSet, "/Users/shashankgupta/Desktop/ImprovAFish/Proteomics-meta/metadata_proteomics.csv")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-SanityCheckMeta(mSet, 1)
mSet<-SetDataTypeOfMeta(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MeanNorm", "SrNorm", "NULL", ratio=FALSE, ratioNum=20)
Proteomics <- data.frame(mSet$dataSet$norm)
Proteomics <- t(Proteomics)

#Look for outliers by examining tree of Proteomics  
dim(Proteomics) %>% paste(c("Proteomics", "Samples"))
sampleTree = hclust(dist(Proteomics), method = "average");
plot(sampleTree, main = "Proteomics clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# outlier<-c("N.Lauroylsarcosine...271.21434.Da.608.71.s")
# Proteomics<- Proteomics %>%
#   as.data.frame %>%
#   filter(!row.names(Proteomics) %in% outlier)
# sampleTree = hclust(dist(Proteomics), method = "average");
# plot(sampleTree, main = "Proteomics clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#Look for outlier samples by examining tree of samples
#Transpose such that samples are in rows and Proteomics are in columns.  
Proteomics <- t(Proteomics)
dim(Proteomics) %>% paste(c("Samples", "Proteomics"))
sampleTree = hclust(dist(Proteomics), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# outlier<-c("T3_CTR_tk3_VFA07")
# Proteomics<- Proteomics %>%
#   as.data.frame %>%
#   filter(!row.names(Proteomics) %in% outlier)

#Transpose such that samples are in rows and Proteomics are in columns.  
dim(Proteomics) %>% paste(c("Samples", "Proteomics"))
powers <- c(1:10, seq(12,30,2))
sft <- pickSoftThreshold(Proteomics, 
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
#write.table(Proteomics, "Proteomics.txt", sep = "\t")
#Proteomics<- read.table("Proteomics.txt", header = T, row.names = 1, stringsAsFactors = F, fill = T)

modules.omics.Proteomics <- blockwiseModules(Proteomics,
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


rownames(modules.omics.Proteomics$MEs) <- rownames(Proteomics)
names(modules.omics.Proteomics$colors) <- colnames(Proteomics)
names(modules.omics.Proteomics$unmergedColors) <- colnames(Proteomics)

hubs_M <- chooseTopHubInEachModule(Proteomics, modules.omics.Proteomics$colors, power = st, omitColors = "0")

stage2results_Proteomics <- list(modules = modules.omics.Proteomics, 
                                 hubs = hubs_M)





# Convert labels to colors for plotting
merged_colors <- labels2colors(stage2results_Proteomics$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_Proteomics$modules$goodSamples) == length(stage2results_Proteomics$modules$goodSamples)
genes_good <- sum(stage2results_Proteomics$modules$goodGenes) == length(stage2results_Proteomics$modules$goodGenes)

ME_good <- sum(stage2results_Proteomics$modules$MEsOK) == length(stage2results_Proteomics$modules$MEsOK)
table(stage2results_Proteomics$modules$colors) %>% 
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


table(stage2results_Proteomics$modules$colors) %>% as.data.frame() -> res
res$`Module color` <- WGCNA::labels2colors(as.numeric(as.character(res$Var1)))
res <- res[, c(1,3,2)]
colnames(res) <- c("Module", "Module color", "Number of Proteomics")

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_Proteomics$modules$dendrograms)){
  plotDendroAndColors(stage2results_Proteomics$modules$dendrograms[[i]], merged_colors[stage2results_Proteomics$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_Proteomics$modules$blockGenes[[i]]),
                                    " Proteomics"))
}
MEs <- stage2results_Proteomics$modules$MEs

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

all(rownames(Proteomics) == rownames(MEs))
dim(Proteomics) %>% paste0(c(" samples", " Proteomics"))
kME <- bicor(Proteomics, MEs, maxPOutliers = 0.05)
dim(kME) %>% paste0(c(" Proteomics", " modules"))

intra_cor <- c()
for (i in 1:ncol(Proteomics)) {
  m <- stage2results_Proteomics$modules$colors[i]
  intra_cor[i] <- kME[i, paste0("ME", m)]
  if(m != 0){
    intra_cor[i] <- kME[i, paste0("ME", m)]
  } else{
    intra_cor[i] <- NA
  }
  
}

idx <- which(is.na(intra_cor))
intra_cor <- intra_cor[-idx]

plot(density(intra_cor), main = "Correlations with module-eigenProteomics (within module correlation)\nNo ME0", xlim = c(-1,1))



# Corr within modules
corr_within_module <- function(Proteomics, modules, module_x = 1){
  idx.Proteomics <- which(modules$colors == module_x)
  idx.me <- which(colnames(modules$MEs) == paste0("ME",module_x))
  kME_x <- bicor(Proteomics[,idx.Proteomics], modules$MEs[,idx.me], maxPOutliers = 0.05)
  kME_x
}

ggplot.list <- list()

for(m in colnames(stage2results_Proteomics$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(Proteomics = Proteomics, modules = stage2results_Proteomics$modules, module_x = h))) %>% 
    ggplot() + 
    #geom_density(aes(x = x), fill = labels2colors(h), color = "black", alpha = 0.5) + 
    geom_histogram(aes(x), fill = labels2colors(h), color = "black", alpha = 0.5, bins = 20) + 
    xlim(-1, 1) +
    theme_bw()+
    xlab("Proteomics correlation")+
    ggtitle(paste0(prefix,m)) -> da_plot
  
  ggplot.list[[m]] <- da_plot
}

ggplot.list <- ggplot.list[ggplot.list %>% names() %>% sub("ME", "", .) %>% as.numeric() %>% order()]

cowplot::plot_grid(plotlist = ggplot.list, ncol = 6) -> density_all_plot
density_all_plot


# combine to one plot
cowplot::plot_grid(si_mc_plot , density_eigen, ncol = 1, rel_heights = c(0.8,1)) -> part_1


cowplot::plot_grid(part_1, module_size_barplot, labels = c("", "C"), label_size = plot_labeling_size, rel_widths = c(1,0.5)) -> part_2


cowplot::plot_grid(part_2, density_all_plot, ncol = 1, rel_heights = c(0.8,1), labels = c("", "F"), label_size = plot_labeling_size)









#stage2results_Proteomics_backup <- stage2results_Proteomics
stage2results_Proteomics 
#rownames(stage2results_Proteomics$modules$MEs) <- gsub("-MC1-TK1-", "_MC1_tk1_", rownames(stage2results_Proteomics$modules$MEs))

# Remove ME0 from analysis
idx <- which(colnames(stage2results_Proteomics$modules$MEs) == "ME0")
stage2results_Proteomics$modules$MEs <- stage2results_Proteomics$modules$MEs[,-idx]


if(take_average){
  stage2results_Proteomics$modules$MEs %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_HGm.*", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    ungroup() -> stage2results_Proteomics_eigengenes
  
  stage2results_Proteomics_eigengenes <- as.data.frame(stage2results_Proteomics_eigengenes)
} else{
  stage2results_Proteomics_eigengenes <- stage2results_Proteomics$modules$MEs
}

rownames(stage2results_Proteomics_eigengenes) <- stage2results_Proteomics_eigengenes$common
stage2results_Proteomics_eigengenes$common <- NULL

dim(stage2results_Proteomics_eigengenes); dim(stage2results_X_eigengenes)



table(rownames(stage2results_X_eigengenes) %in% rownames(stage2results_Proteomics_eigengenes))

stage2results_X_eigengenes_New <- stage2results_X_eigengenes[rownames(stage2results_X_eigengenes) %in% rownames(stage2results_Proteomics_eigengenes), ]


table( rownames(stage2results_X_eigengenes_New) %in% rownames(stage2results_Proteomics_eigengenes))
table( rownames(stage2results_X_eigengenes_New) == rownames(stage2results_Proteomics_eigengenes))



#Correlate modules from transcriptomics and metagenomics.

# Check that the samples are in the same order. 
# If they are not in order, change their order to match; If they do not match one-to-one, call an error.
same_order <- all(rownames(stage2results_Proteomics_eigengenes) == rownames(stage2results_X_eigengenes_New))
if(!same_order){
  stage2results_Proteomics_eigengenes <- stage2results_Proteomics_eigengenes[order(rownames(stage2results_Proteomics_eigengenes)),]
  stage2results_X_eigengenes_New <- stage2results_X_eigengenes_New[order(rownames(stage2results_X_eigengenes_New)),]
  same_order <- all(rownames(stage2results_Proteomics_eigengenes) == rownames(stage2results_X_eigengenes_New))
  if(!same_order){
    stop("Sample names do not match. Samples should be identical.", call. = F)
  }
} else{cat("Samples match")}


#####
p.value_matr <- corr.value_matr <- matrix(ncol = ncol(stage2results_Proteomics_eigengenes), 
                                          nrow = ncol(stage2results_X_eigengenes_New), 
                                          dimnames = list(colnames(stage2results_X_eigengenes_New), 
                                                          colnames(stage2results_Proteomics_eigengenes)))


for(i in 1:ncol(stage2results_X_eigengenes_New)){
  for(j in 1:ncol(stage2results_Proteomics_eigengenes)){
    cor.res <- cor.test(stage2results_X_eigengenes_New[,i], stage2results_Proteomics_eigengenes[,j])
    p.value_matr[i, j] <- cor.res$p.value
    corr.value_matr[i, j] <- cor.res$estimate
  }
}

# Correct for number of tests
p.value_matr.adjust <- p.adjust(p.value_matr, method = "fdr")
dim(p.value_matr.adjust) <- dim(p.value_matr)
dimnames(p.value_matr.adjust) <- list(colnames(stage2results_X_eigengenes_New), colnames(stage2results_Proteomics_eigengenes))


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
Proteomics_corr_X <- list(p_value = p.value_matr, 
                          p_value_adj = p.value_matr.adjust,
                          signif_matrix = signif_matrix,
                          correlation = corr.value_matr)
rm(p.value_matr, p.value_matr.adjust, signif_matrix, corr.value_matr)



heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 6, name ="RdBu")))(51)


stage2results_Proteomics$hubs %>% 
  as.data.frame() %>% 
  dplyr::rename("Proteomics_name" = ".") %>%
  tibble::rownames_to_column(var = "Module") -> hubProteomics

hubProteomics -> hubProteomics_Before

hubProteomics_Before$ModuleProteomics <- paste0("ME" ,hubProteomics_Before$Module)


test <- Proteomics_corr_X$correlation
hubProteomics_Before$Proteomics_name <- sub("XP_.*", "", hubProteomics_Before$Proteomics_name)
hubProteomics_Before$Proteomics_name <- sub("NP_.*", "", hubProteomics_Before$Proteomics_name)


hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="cytochrome.c.oxidase.subunit.5B..mitochondrial", "Cytochrome c oxidase subunit 5B, mitochondrial")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="myosin.7.like.isoform.X1", "myosin-7 isoform X1")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="creatine.kinase.S.type..mitochondrial", "Creatine kinase S-type, mitochondrial")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="X3.hydroxyisobutyrate.dehydrogenase..mitochondrial", "3-hydroxyisobutyrate dehydrogenase, mitochondrial")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="LOW.QUALITY.PROTEIN..zonadhesin.like", "Zonadhesin like")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="cytosolic.non.specific.dipeptidase", "Cytosol nonspecific dipeptidas")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="protocadherin.Fat.4", "Protocadherin Fat 4")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="cadherin.17.precursor", "Cadherins 17 precursor")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="sushi.domain.containing.protein.2.isoform.X2", "Sushi domain-containing protein 2")
hubProteomics_Before$Proteomics_name <- replace(hubProteomics_Before$Proteomics_name, hubProteomics_Before$Proteomics_name=="intestinal.type.alkaline.phosphatase", "Intestinal-type alkaline phosphatase")


match(hubProteomics_Before[, "ModuleProteomics"], colnames(test))
colnames(test)[match(hubProteomics_Before[,"ModuleProteomics"], colnames(test))] = hubProteomics_Before[,"Proteomics_name"]

Proteomics_corr_X$correlation <- test
pheatmap::pheatmap(Proteomics_corr_X$correlation, 
                   color = heatmap_colors, 
                   treeheight_col = 0, 
                   treeheight_row = 0,  # will be shown on the transcriptomics ME heatmap
                   cluster_rows = X_ME_dendro,
                   cutree_rows = row_cut,
                   display_numbers = Proteomics_corr_X$signif_matrix, 
                   fontsize_number = 8, #10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F,
                   show_rownames = F,
                   legend = F,
                   labels_row = paste0(prefix_OTUs, rownames(Proteomics_corr_X$correlation)),
                   labels_col = paste0(colnames(Proteomics_corr_X$correlation)),
                   main = paste("Eigen\n", "Meta-proteins" )
                   ) -> Proteomics_corr_X_plot














#https://rpubs.com/mohsen/lefse_analysis_cleaned_prevalence_phyloseq
library(tidyverse)
library(phyloseq)
library(microbiomeMarker)
raw <- as.matrix(Proteomics)
OTU = otu_table(raw, taxa_are_rows = TRUE)
dat <- read.csv("/Users/shashankgupta/Desktop/ImprovAFish/Proteomics-meta/metadata_proteomics.csv")
row.names(dat) <- dat$Samples
dat$New_Diet <- toupper(dat$Diet)

# Merge into one complete phyloseq object
ps <- merge_phyloseq(otu_table(OTU), sample_data(dat))
tt <- as.data.frame(row.names(Proteomics))
row.names(tt) <- tt$`row.names(Proteomics)`
colnames(tt)[1] <- "Kingdom"
tax_table(ps) <- as.matrix(tt)

set.seed(12345)
ps_MC1<-subset_samples(ps, New_Diet %in% c("CTR", "MC1"))
ps_MC1 <- prune_taxa(taxa_sums(ps_MC1) >0, ps_MC1)
table(sample_data(ps_MC1)$New_Diet)


lef_out<-run_lefse(ps_MC1, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")


plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)
data.frame(marker_table(lef_out)) %>%
  filter(enrich_group != "CTR") %>%
  select(feature) %>%
  `rownames<-`( NULL )


res <- ldamarker(ps_MC1, group="New_Diet")
plotLDA(res,group=c("MC1","CTR"), lda=1.75, padj =  0.05, fontsize.y = 6)


##Change orientation
# t1<- data.frame(marker_table(lef_out))
# t1 %>%
#   mutate(ef_lda = ifelse(enrich_group == "CTR", ef_lda, -ef_lda)) -> t1
# marker_table(lef_out) <- t1


ps_MC2<-subset_samples(ps, New_Diet %in% c("CTR", "MC2"))
ps_MC2 <- prune_taxa(taxa_sums(ps_MC2) >0, ps_MC2)
table(sample_data(ps_MC2)$Diet)


lef_out<-run_lefse(ps_MC2, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)
data.frame(marker_table(lef_out)) %>%
  filter(enrich_group != "CTR") %>%
  select(feature) %>%
  `rownames<-`( NULL )


res <- ldamarker(ps_MC2, group="Diet")
plotLDA(res,group=c("MC2","CTR"), lda=1.75, padj =  0.05, fontsize.y = 6)


ps_MN3<-subset_samples(ps, New_Diet %in% c("CTR", "MN3"))
ps_MN3 <- prune_taxa(taxa_sums(ps_MN3) >0, ps_MN3)
table(sample_data(ps_MN3)$New_Diet)


lef_out<-run_lefse(ps_MN3, group = "New_Diet", norm = "CPM", 
                   kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")

plot_ef_bar(lef_out)
table(marker_table(lef_out)$enrich_group)
data.frame(marker_table(lef_out)) %>%
  filter(enrich_group != "CTR") %>%
  select(feature) %>%
  `rownames<-`( NULL )



res <- ldamarker(ps_MN3, group="Diet")
plotLDA(res,group=c("MN3","CTR"), lda=1.75, padj =  0.05, fontsize.y = 6)


#loop
output_table <- data.frame()

for (i in c("T3")){
  for (j in c("MC1","MC2", "MN3")){
    phy<-subset_samples(ps, Time %in% c(i) & New_Diet %in% c("CTR", j))
    phy <- prune_taxa(taxa_sums(phy) >0, phy)
    lef_out<-run_lefse(phy,group = "New_Diet", norm = "CPM", 
                       kw_cutoff = 0.05, lda_cutoff = 1.75, taxa_rank = "none")
    assign(paste0("poteomics", i, j), plot_abundance(lef_out, group = "New_Diet"))
    output_table[i,j]<-sum(marker_table(lef_out)$enrich_group == j)
  }
}
#colnames(output_table)<-c("MC1", "MC2", "MN3")
#rownames(output_table) <- c("T1","T2","T3")
output_table

#plot_grid(pT1MC1, pT1MC2, pT1MN3, pT2MC1, pT2MC2, pT2MN3, pT3MC1, pT3MC2, p3MN3, labels = "AUTO", ncol = 1)
plot_grid(poteomicsT3MC1, poteomicsT3MC2, poteomicsT3MN3, labels = "AUTO", ncol = 1, align = "hv")





MEs = moduleEigengenes(Proteomics_host, moduleLabels1)$eigengenes

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

mME_meta$New_Diet <- factor(mME_meta$New_Diet, levels = c("CTR", "MC1", "MC2", "MN3"))
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








