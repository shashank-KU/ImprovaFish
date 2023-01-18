#Additional figures
set.seed(13118) # Set seed for reproducibility if any of the methods used are non deterministic.
options(stringsAsFactors = F)
library(magrittr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)


stage2results_Y

otus_to_keep <- omics_data %>% rownames()

samples_to_keep <- omics_data %>% colnames()

otu.raw <- omics_data[otus_to_keep, samples_to_keep]

otu.raw %>% dim()

percentage <- function(x){x/sum(x)}

otu.raw <- apply(otu.raw, MARGIN = 2, FUN = percentage)

take_average <- T
if(take_average){
  otu.raw %>% 
    t() %>% 
    data.frame %>% 
    rownames_to_column(var = "common") %>% 
    mutate(common = sub("_[0-9]{1,2}", "", common)) %>% 
    group_by(common) %>% 
    summarise_all(list(mean)) %>% 
    column_to_rownames(var = "common") %>% 
    t() -> otu.raw
}


# 
# if(take_average){
#     otu.raw %>% 
#   
#     data.frame %>% 
#     rownames_to_column(var = "common") %>%
#     mutate(common = sub("_[0-9]{1,2}$", "", common)) %>% 
#     group_by(common) %>% 
#     summarise_all(list(mean))-> otu.raw
# 
# }
# 
# rownames(otu.raw) <- otu.raw$common
# otu.raw$common <- NULL
# otu.raw <- t(otu.raw)
 

stage2results_Y$modules$dendrograms[[1]] %>% as.dendrogram %>% dendextend::hang.dendrogram(., hang = -1) %>% dendextend::as_hclust_fixed(.) -> otu_clusters

otu_annotat_module <- data.frame(otu.names =  names(stage2results_Y$modules$colors), Module = stage2results_Y$modules$colors) %>%
  mutate(Module_color = WGCNA::labels2colors(Module)) %>% 
  mutate(Module = paste0("Module ", Module)) %>% 
  mutate(Module = as.factor(Module))

module_col <- otu_annotat_module[,-1] %>% distinct()
color_vector <- setNames(module_col$Module_color, module_col$Module)


annotation_rows_otu <- data.frame(OTU_name = otu_annotat_module$otu.names, Module = otu_annotat_module$Module)


otu_taxonomy %>% 
  rownames_to_column(var = "OTU_name") -> A

annotation_rows_otu <- left_join(annotation_rows_otu, A, by = "OTU_name") %>% 
  #mutate(`Core OTU` = ifelse(OTU_name %in% c("OTU_1", "OTU_2", "OTU_5", "OTU_10"), "Yes", "No")) %>% 
  dplyr::select(OTU_name, Module, Phylum) %>% column_to_rownames(var = "OTU_name")

tax_color_vector <- annotation_rows_otu %>% 
  dplyr::select("Phylum") %>% 
  mutate_all(.funs = list(as.character)) %>%  # Remove factor
  unique() %>% 
  mutate(tax_color = colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(nrow(.)))

tax_color_vector <- setNames(as.character(tax_color_vector$tax_color), tax_color_vector$Phylum)




annotation_colors <- list(
  Time = c(`T1` = "navy", `T2` = "darkgreen", `T3` = "red"),
  New_Diet = c(`ctr` = "#1D88AF", `mc1` = "darkblue", `mc2` = "darkgreen", `mn3` = "#F08A46"),
  Tank_number = c(`tk1` = "#F08A46", `tk2` = "#8EB470",`tk3` = "#B7CFA4"),
  Sample_Number =c(`01` = "#1D88AF", `03` = "#12556D", `05` = "#F6B890", `07` = "#F08A46"),
  Module = color_vector,
  Phylum = tax_color_vector)

otu.raw %>% t() -> abundata

abundata[1:5,1:5]

abundata %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "common") %>% 
  separate(col = common, into = c("Time", "New_Diet", "Tank_number"), sep = "_", remove = F, fill = "warn") %>% 
  mutate(New_Diet = factor(New_Diet, levels = c("ctr", "mc1", "mc2", "mn3"))) %>% 
  mutate(Time = factor(Time, levels = c("T1", "T2", "T3"))) %>% 
  #mutate(Day = as.integer(sub("D", "", Day))) %>% 
  arrange(New_Diet, Time, Tank_number) %>% 
  tibble::column_to_rownames(var = "common")  -> otu.raw.prop_reordered




otu.raw.prop_reordered %>%
  rownames() %>%
  as.data.frame() %>% 
  dplyr::rename(., "common" = ".") %>%
  separate(col = common, into = c("Time", "New_Diet", "Tank_number"), sep = "_", remove = F, fill = "warn") %>% 
  dplyr::select(common,  Time, New_Diet, Tank_number) %>%
  mutate(New_Diet = factor(New_Diet, levels = c("ctr", "mc1", "mc2", "mn3"))) %>% 
  mutate(Time = factor(Time, levels = c("T1", "T2", "T3"))) %>% 
  #mutate(Day = as.integer(sub("D", "", Day))) %>%
  tibble::column_to_rownames(var = "common") -> annotation_columns


pheatmap::pheatmap(otu.raw.prop_reordered, 
                   #cluster_cols = T,
                   clustering_method = "ward.D2",
                   cluster_rows = otu_clusters,
                   breaks = c(-0.0001, 0, 0.0001, 0.001, 0.005, seq(0.01,1,length.out = 10)),
                   color = c("#F5F5F5", colorRampPalette(colors = c("#ffc7b8","#3f0000"), bias = 1)(12)), 
                   annotation_row = annotation_rows_otu, 
                   annotation_col = annotation_columns,
                   annotation_colors = annotation_colors,
                   show_rownames = F, 
                   show_colnames = F,
                   fontsize = 12,
                   treeheight_row = 100) -> OTU_heatmap
























