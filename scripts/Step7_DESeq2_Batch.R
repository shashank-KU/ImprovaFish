#DESeq2 analysis



omics_data_host <- rawCounts_total2
sample_info<- mapping_file

omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% row.names(sample_info)]
omics_data_host <- omics_data_host1
all(colnames(omics_data_host) == row.names(sample_info))

theme_set(theme_bw())


library(DESeq2)


omics_data_host1 <- omics_data_host[, colnames(omics_data_host) %in% row.names(sample_info)]
omics_data_host <- omics_data_host1
all(colnames(omics_data_host) == row.names(sample_info))
reorder_idx <- match(row.names(sample_info), colnames(omics_data_host))
omics_data_host <- omics_data_host[ , reorder_idx]
all(colnames(omics_data_host) %in% row.names(sample_info))
all(colnames(omics_data_host) == row.names(sample_info))



#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
df <- round(omics_data_host) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 1)

all(colnames(df) == rownames(sample_info))

dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = sample_info, # Data frame with annotation for our samples
  design = ~Time+New_Diet
)

dds$New_Diet <- relevel(dds$New_Diet, ref = "ctr")

dds <- DESeq(dds, parallel = T)
resultsNames(dds)


resCondition <- results(dds, name = "New_Diet_mc1_vs_ctr") %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat("DEGs: ", nrow(resCondition))

vsd <- vst(dds, blind=FALSE)

library(stringi)
library(tidyr)
plot1 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(New_Diet = stri_replace_all_regex(Sample, colnames(dds), dds$New_Diet, vectorize=FALSE)) %>%
  ggplot(aes(x = New_Diet, y = Expression, color = New_Diet)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

plot2 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(Time = stri_replace_all_regex(Sample, colnames(dds), dds$Time, vectorize=FALSE)) %>%
  ggplot(aes(x = Time, y = Expression, color = Time)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)


library(EnhancedVolcano)
resCondition <- results(dds, name = "New_Diet_MC1_vs_CTR")
EnhancedVolcano(resCondition,
                lab = rownames(resCondition),
                x = 'log2FoldChange',
                y = 'padj',
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)



# Diet 2
resCondition <- results(dds, name = "New_Diet_MC2_vs_CTR") %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat("DEGs: ", nrow(resCondition))

vsd <- vst(dds, blind=FALSE)


library(stringi)
plot1 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(New_Diet = stri_replace_all_regex(Sample, colnames(dds), dds$New_Diet, vectorize=FALSE)) %>%
  ggplot(aes(x = New_Diet, y = Expression, color = New_Diet)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

plot2 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(Time = stri_replace_all_regex(Sample, colnames(dds), dds$Time, vectorize=FALSE)) %>%
  ggplot(aes(x = Time, y = Expression, color = Time)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)

resCondition <- results(dds, name = "New_Diet_MC2_vs_CTR")
EnhancedVolcano(resCondition,
                lab = rownames(resCondition),
                x = 'log2FoldChange',
                y = 'padj',
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)






# Diet 3
resCondition <- results(dds, name = "New_Diet_MN3_vs_CTR") %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat("DEGs: ", nrow(resCondition))

vsd <- vst(dds, blind=FALSE)


library(stringi)
plot1 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(New_Diet = stri_replace_all_regex(Sample, colnames(dds), dds$New_Diet, vectorize=FALSE)) %>%
  ggplot(aes(x = New_Diet, y = Expression, color = New_Diet)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

plot2 <- assay(vsd) %>%
  as.data.frame %>%
  add_rownames(var = "Genes") %>%
  filter(Genes %in% resCondition$Genes[1:5]) %>%
  gather(Sample, Expression, -Genes) %>%
  mutate(Time = stri_replace_all_regex(Sample, colnames(dds), dds$Time, vectorize=FALSE)) %>%
  ggplot(aes(x = Time, y = Expression, color = Time)) +
  geom_boxplot() + 
  facet_grid(rows = vars(Genes))

cowplot::plot_grid(plotlist = list(plot1, plot2), ncol = 2)


resCondition <- results(dds, name = "New_Diet_MN3_vs_CTR")
EnhancedVolcano(resCondition,
                lab = rownames(resCondition),
                x = 'log2FoldChange',
                y = 'padj',
                #xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)





#
vsd <- vst(dds)
plotPCA(vsd, "Time") + stat_ellipse()
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Time)
plotPCA(vsd, "Time") + stat_ellipse()






################################################################################
# Liver
library(edgeR)
mapping_file3 <- subset(mapping_file, mapping_file$Sample_Type %in% c("Liver"))
#mapping_file3 <- subset(mapping_file3, mapping_file3$Diet %in% c("mc1", "ctr"))
mapping_file4 <- mapping_file3[ !grepl(c("T0_ctr_") , row.names(mapping_file3)) , ]

rawCounts_total3<- rawCounts_total2
rawCounts_total3 <- rawCounts_total3[ , grepl( "_Lh" , names( rawCounts_total3 ) ) ]
#rawCounts_total3 <- rawCounts_total3[ , grepl( c("_mc1_|_ctr_") , names( rawCounts_total3 ) ) ]
rawCounts_total4 <- rawCounts_total3[ , !grepl(c("T0_ctr_") , names(rawCounts_total3)) ]
all(colnames(rawCounts_total4) == rownames(mapping_file4))

#head(rawCounts_total4)
#head(mapping_file4)
#plotMDS(y)

mapping_file4$group<-paste(mapping_file4$Diet, mapping_file4$Time, sep = ".")
head(mapping_file4)
table(mapping_file4$Time, mapping_file4$Diet)
# Filtering







