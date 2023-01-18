#Time as a covariates




nGenes = ncol(omics_data_host)
nSamples = nrow(omics_data_host)
moduleLabels1 = stage2results_X$modules$colors
moduleColors1 = labels2colors(stage2results_X$modules$colors)
MEs = stage2results_X$modules$MEs

table(moduleColors1)
table(moduleLabels1)
moduleTraitCor = cor(MEs, sample_info1, use = "p")
#moduleTraitCor <- as.data.frame(moduleTraitCor)
moduleTraitCor <- moduleTraitCor[order(match(rownames(moduleTraitCor), c("ME4","ME1","ME24","ME7","ME15","ME9","ME22","ME6","ME19","ME8","ME20","ME10","ME16","ME18","ME13","ME12","ME5","ME11","ME3","ME14","ME21","ME23","ME17","ME2"))), ]
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


# Loop for linear model by keeping TIME as a co-variate
library(purrr)
sample_info1$ID_New <- row.names(sample_info1)
sample_info2 <- left_join(sample_info1, sample_info, by = c("ID_New" = "ID_New"))
row.names(sample_info2) <- sample_info2$ID_New
sample_info2 <- sample_info2[,c(1:16)]
sample_info2$ID_New <- NULL

do_thing <- function(var){
  (results <- map_dbl(1:ncol(MEs),
                      ~{
                        fit <- lm(MEs[[.x]] ~ sample_info2[[var]] * sample_info2$Time)
                        summary(fit)$coefficients[2,4] #To extract the pvalue
                      }) |> set_names(1:ncol(MEs)))
  test1 <- as.data.frame(results)
  
  names(test1) <- var
  test1
}

coldf <- colnames(sample_info2[,c(1:14)])
map_dfc(coldf,
        do_thing) -> FINAL


# fit <- lm(MEs$ME9 ~ sample_info2$HSI * sample_info2$Time)
# summary(fit)$coefficients[2,4]

rownames(FINAL) <- colnames(MEs)
FINAL <- FINAL[order(match(rownames(FINAL),  c("ME4","ME1","ME24","ME7","ME15","ME9","ME22","ME6","ME19","ME8","ME20","ME10","ME16","ME18","ME13","ME12","ME5","ME11","ME3","ME14","ME21","ME23","ME17","ME2"))), ]
FINAL_1 <- as.matrix(FINAL)
signif_matrix1 <- rep("", length(FINAL_1))
three_star <- which( FINAL_1 <= 0.001)
signif_matrix1[three_star] <- "***"
two_star <- which((FINAL_1 <= 0.01) & (FINAL_1 > 0.001))
signif_matrix1[two_star] <- "**"
one_star <- which((FINAL_1 <= 0.05) & (FINAL_1 > 0.01))
signif_matrix1[one_star] <- "*"
dim(signif_matrix1) = dim(FINAL_1) # Give textMatrix the correct dimensions 

pheatmap::pheatmap(moduleTraitCor, 
                   color = heatmap_colors,
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   cluster_cols = F,
                   #treeheight_col = 0, 
                   #treeheight_row = 0, # will be shown on the transcriptomics ME heatmap
                   cluster_rows = F,
                   #cutree_rows = row_cut,
                   display_numbers = signif_matrix1, 
                   fontsize_number = 8, # 10
                   breaks = seq(from = -1, to = 1, length.out = 51), 
                   silent = F, 
                   legend = F,
                   show_rownames = F,
                   main = "Traits and external variables") -> MT_plot



