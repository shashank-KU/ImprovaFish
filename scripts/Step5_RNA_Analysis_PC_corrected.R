#host-microbiome interactions
library(vegan)
library(tidyverse)
library(dendextend)
library(magrittr)                  # For pipes %>% 
library(cowplot)                   # For combining plots
library(dendextend)                # Modify dendrograms
library(WGCNA)                     # Weighted gene co-expression network analysis
library(dplyr)                     # Manipulation of data frames
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(genefilter)
psdata.p  = transform_sample_counts(psdata, function(x) sqrt(x / sum(x))) # Hellinger transformation
flist<- filterfun(kOverA(5, 2e-05))

psdata.p = filter_taxa(psdata.p, flist, TRUE)
psdata.p


theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.title = element_text(face="bold"))

dat2 <- sample_data(psdata.p)
omics_data <- otu_table(psdata.p)
host_gut_mapping_file <- read.table("../host_data/host_gut_mapping.txt", sep = "\t", row.names = 1, header = T)
#host_gut_rawCounts_total1 <- read.table("../host_data/host_gut_count.txt")


dat2$common <- gsub("[-]","_", dat2$sampleName)
dat2$common <- gsub("_HGm","_", dat2$common)

host_gut_mapping_file$sampleName <- row.names(host_gut_mapping_file)
host_gut_mapping_file$common <- host_gut_mapping_file$sampleName
host_gut_mapping_file$common <- gsub("_HGh","_", host_gut_mapping_file$common)


table(dat2$common %in% host_gut_mapping_file$common)
table(host_gut_mapping_file$common %in% dat2$common)

host_gut_mapping_file <- host_gut_mapping_file[host_gut_mapping_file$common %in% dat2$common, ]
dat2 <- dat2[dat2$common %in% host_gut_mapping_file$common, ]


dat2$rownames <- row.names(dat2)
row.names(dat2) <- dat2$common



common_name_list <-  dat2$rownames
gut_16SrRNA_OTU_table <- omics_data[, colnames(omics_data) %in% common_name_list]

table(colnames(gut_16SrRNA_OTU_table) ==  dat2$rownames)
colnames(gut_16SrRNA_OTU_table) <- dat2$common
table(colnames(gut_16SrRNA_OTU_table) ==  dat2$common)

# colnames(host_gut_rawCounts_total1) <- gsub("_HGh","_", colnames(host_gut_rawCounts_total1))
# row.names(host_gut_mapping_file) <- host_gut_mapping_file$common


omics_data <- gut_16SrRNA_OTU_table
sample_info <- dat2

all(colnames(omics_data) == sample_info$common)


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


pcCorrection <- T
if(pcCorrection){
  estimate_n.pc = T
  if(!estimate_n.pc){
    number_of_pcs_to_remove = 1 # Does not matter when pcCorrection is FALSE
  }
}



#Step 1
library(WGCNA)
dim(omics_data)
omics_data <- t(omics_data)
dim(omics_data) %>% paste(c( "Samples", "OTUs"))

# omics_data[1:3, 1:3]
# HellingerData<-decostand(omics_data,method = "hellinger")
# omics_data <- HellingerData

omics_data <- t(omics_data)

dim(omics_data) %>% paste(c("OTUs", "Samples"))



if(pcCorrection){
  library(sva)
  if(estimate_n.pc){
    
    print(t(omics_data) %>% dim() %>% paste(c("Samples", "genes"))) 
    mod=matrix(1,nrow=nrow(t(omics_data)),ncol=1)
    colnames(mod)="Intercept"
    
    print(omics_data %>% dim() %>% paste(c("genes", "Samples"))) 
    n.pc = num.sv(dat = omics_data, mod = mod,  method="be")
  } else{
    n.pc = number_of_pcs_to_remove
  }
  
  print(omics_data %>% dim() %>% paste(c("genes", "Samples")))
  omics_data <- sva_network(omics_data, n.pc = n.pc)
}

if(pcCorrection){
  report <- paste0("PC-correction with ", n.pc, " PC(s) removed.")
  
  cowplot::plot_grid(plot_hist(omics_data, "sd"),
                     plot_hist(omics_data, "mad"),
                     plot_hist(omics_data, "no"),
                     align = "h",
                     ncol = 3)
} else{
  report <- "There was no PC-correction performed."
}
report

omics_data <- omics_data %>% t()
dim(omics_data) %>% paste(c("Samples", "OTUs"))

powers <- c(1:10, seq(12,30,2))
suppressWarnings(sft <- pickSoftThreshold(omics_data, 
                                          powerVector = powers, 
                                          verbose = set_verbose, 
                                          networkType = "signed",
                                          corFn= chosen_parameter_set$assoc_measure))



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
  xlim(1,30) +
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
modules.omics_Y_PC_corrected <- blockwiseModules(omics_data,
                                    power = st, 
                                    networkType = "signed", 
                                    TOMType = "signed",
                                    corType = chosen_parameter_set$assoc_measure,
                                    maxPOutliers = 0.05,
                                    deepSplit = 4, # Default 2
                                    minModuleSize = 10, # 30
                                    minCoreKME = 0.5,      # Default 0.5
                                    minCoreKMESize = 2,    # Default minModuleSize/3,
                                    minKMEtoStay = 0.5,    # Default 0.3
                                    reassignThreshold = 0, # Default 1e-6
                                    mergeCutHeight = 0.2,  # Default 0.15
                                    pamStage = pam_stage, 
                                    pamRespectsDendro = TRUE,
                                    replaceMissingAdjacencies = TRUE,
                                    numericLabels = TRUE,
                                    saveTOMs = save_TOM,
                                    saveTOMFileBase = "TOM",
                                    nThreads = 10,
                                    verbose = 3,
                                    maxBlockSize=8000)




rownames(modules.omics_Y_PC_corrected$MEs) <- rownames(omics_data)
names(modules.omics_Y_PC_corrected$colors) <- colnames(omics_data)
names(modules.omics_Y_PC_corrected$unmergedColors) <- colnames(omics_data)

hubs <- chooseTopHubInEachModule(omics_data, modules.omics_Y_PC_corrected$colors, omitColors = 0)

stage2results_Y_PC_corrected <- list(modules = modules.omics_Y_PC_corrected, 
                        hubs = hubs)

merged_colors <- labels2colors(stage2results_Y_PC_corrected$modules$colors)
n_modules <- unique(merged_colors) %>% length()

samples_good <- sum(stage2results_Y_PC_corrected$modules$goodSamples) == length(stage2results_Y_PC_corrected$modules$goodSamples)
OTUs_good <- sum(stage2results_Y_PC_corrected$modules$goodGenes) == length(stage2results_Y_PC_corrected$modules$goodGenes)

ME_good <- sum(stage2results_Y_PC_corrected$modules$MEsOK) == length(stage2results_Y_PC_corrected$modules$MEsOK)

table(stage2results_Y_PC_corrected$modules$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename(Module = Var1, Size = Freq) %>% 
  dplyr::mutate(Module_color = labels2colors(as.numeric(as.character(Module)))) -> module_size

module_size %>% 
  ggplot(aes(x = Module, y = Size, fill = Module)) +
  geom_col(color =  "#000000") +
  ggtitle("Number of AVs in each module") +
  theme(legend.position = "none") + 
  scale_fill_manual(values = setNames(module_size$Module_color,module_size$Module)) +
  geom_text(aes(label = Size),vjust = 0.5, hjust = -0.18, size = 3.5) +
  ylim(0, max(module_size$Size)*1.1) +
  theme(plot.margin = margin(2, 2, 2, 2, "pt")) +
  coord_flip()-> module_size_barplot

module_size_barplot



table(stage2results_Y_PC_corrected$modules$colors) %>% as.data.frame() -> res
res$`Module color` <- WGCNA::labels2colors(as.numeric(as.character(res$Var1)))
res <- res[, c(1,3,2)]
colnames(res) <- c("Module", "Module color", "Number of ASVs")
res %>% stargazer::stargazer(type = "latex", summary = FALSE, rownames = FALSE)

# Plot the dendrogram and the module colors underneath for each block
for(i in seq_along(stage2results_Y_PC_corrected$modules$dendrograms)){
  plotDendroAndColors(stage2results_Y_PC_corrected$modules$dendrograms[[i]], merged_colors[stage2results_Y_PC_corrected$modules$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste0("Cluster Dendrogram\n", 
                                    "for block ", 
                                    i,": ",
                                    length(stage2results_Y_PC_corrected$modules$blockGenes[[i]]),
                                    " OTUs"))
}
MEs <- stage2results_Y_PC_corrected$modules$MEs

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
  m <- stage2results_Y_PC_corrected$modules$colors[i]
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

for(m in colnames(stage2results_Y_PC_corrected$modules$MEs)){
  h <- as.numeric(sub("ME","", m))
  data.frame(x = suppressWarnings(corr_within_module(omics_data = omics_data, modules = stage2results_Y_PC_corrected$modules, module_x = h))) %>% 
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


degrees <- intramodularConnectivity.fromExpr(omics_data, colors = modules.omics_Y_PC_corrected$colors, power = st,
                                             networkType = "signed", distFnc = chosen_parameter_set$assoc_measure)

degrees$OTU_name <- colnames(omics_data)
degrees$Module <- modules.omics_Y_PC_corrected$colors
degrees <- degrees[,c(6,5,1:4)]


dplyr::left_join(degrees, 
                 (otu_taxonomy %>%
                    tibble::rownames_to_column(var = "OTU_name")), 
                 by = "OTU_name") -> degrees


degrees  %>%
  dplyr::filter(Module == "4")









