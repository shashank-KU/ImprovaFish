

library(gtsummary)
library(dplyr)

df1 <- data.frame(sample_info1)
df2<- data.frame(sample_info)
rownames(df2) <- df2$ID_New
df1$ID_New <- rownames(df1)



sample_info2<- merge(df1, df2, by = "ID_New")


trial2 <- sample_info2 %>% select(BW, Length, Gut_light, Anus_bleed, Gut_bleed, New_Diet)


table2 <- 
  tbl_summary(
    trial2,
    by = New_Diet, # split table by group
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p() %>% # sample_info for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 

table2
head(sample_info)


library(ggpubr)

sample_info$New_Diet <- factor(sample_info$New_Diet, levels=c( 'ctr', 'mc1', 'mc2', 'mn3'))
levels(sample_info$New_Diet)

sample_info$Time <- factor(sample_info$Time, levels=c('T1', 'T2', 'T3'))
levels(sample_info$Time)

table(sample_info$Group)

sample_info$Group <- factor(sample_info$Group, levels=c('ctr_T1', 'ctr_T2', 'ctr_T3', 
                                          'mc1_T1', 'mc1_T2', 'mc1_T3',
                                          'mc2_T1', 'mc2_T2', 'mc2_T3',
                                          'mn3_T1', 'mn3_T2', 'mn3_T3'))
levels(sample_info$Group)

#compare_means(Shannon ~ race, data = sample_info, method = "kruskal.sample_info")
my_comparisons <- list( c("ctr", "mc1"), c("ctr", "mc2"),
                        c("ctr", "mn3"), c("mc1", "mn2"),
                        c("mc1", "mn3"), c("mc2", "mn3"))

my_comparisons <- list( c("ctr_T1", "mc1_T1"), c("ctr_T1", "mc2_T1"),c("ctr_T1", "mn3_T1"),
                        c("ctr_T2", "mc1_T2"), c("ctr_T2", "mc2_T2"),c("ctr_T2", "mn3_T2"),
                        c("ctr_T3", "mc1_T3"), c("ctr_T3", "mc2_T3"),c("ctr_T3", "mn3_T3"))
library(cowplot)

p9 <- ggboxplot(sample_info, x = "Group", y = "BW",
                color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 600) + 
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value


p10 <- ggboxplot(sample_info, x = "Group", y = "Length",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 45) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value

plot_grid(p9, p10, labels = c('A', 'B' ), label_size = 12, ncol = 1)


p11 <- ggboxplot(sample_info, x = "Group", y = "CSI",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.45) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  


p12 <- ggboxplot(sample_info, x = "Group", y = "HSI",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6.5) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())  

plot_grid(p11, p12, labels = c('A', 'B' ), label_size = 12, ncol = 1)


ggboxplot(sample_info, x = "Group", y = "Length",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) +facet_wrap("Sex")+
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value






ggboxplot(sample_info, x = "Group", y = "Shannon",
          color = "New_Diet", palette = "jco", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 7) +
  #geom_jitter(aes(colour = Diet), size = 2, alpha = 0.8) +
  geom_boxplot(aes(fill = New_Diet), width=0.7, alpha = 0.5) +
  theme_bw() +  theme(legend.position="none",axis.title.x=element_blank())     # Add global p-value





# a diversity correlation
ITS_sample <- sample_data(Final.ITS)
RNA_sample <- sample_data(Final.RNA)

ITS_alpha <- data.frame(ITS_sample[,which(colnames(ITS_sample) %in% c("abcno","ObservedRichness_mean","SDI_mean"))])
RNA_alpha <- data.frame(RNA_sample[,which(colnames(RNA_sample) %in% c("abcno","ObservedRichness_mean","SDI_mean"))])
colnames(ITS_alpha)[c(2:3)] <- c("ITS_OR_mean","ITS_SDI_mean")
colnames(RNA_alpha)[c(2:3)] <- c("RNA_OR_mean","RNA_SDI_mean")
ITS_RNA_alpha <- merge(ITS_alpha,RNA_alpha,by="abcno")

ggscatter(sample_info, x = "BW", y = "Richness",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "BW", ylab = "RNA_Richness")

ggscatter(sample_info, x = "Length", y = "Richness",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Length", ylab = "RNA_Richness")


p13 <- ggscatter(sample_info, x = "BW", y = "Shannon",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "BW", ylab = "Shannon Diversity Index")

p14 <- ggscatter(sample_info, x = "Length", y = "Shannon",
          add = "reg.line", add.params = list(color = "red"),size = 1, conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Length", ylab = "Shannon Diversity Index")



plot_grid(p13, p14, labels = c('A', 'B' ), label_size = 12, ncol = 1)







