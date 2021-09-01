library(tidyverse)
library(ggplot2)
library(reshape2)
library(factoextra)

# load data ---------------------------------------------------------------
samples.df <- read.csv("sample_names.txt", header = F, col.names = c("Sample"))

rc.df <- read.csv("GBS/POPOOLATION2/F_rc.cov60.MAF0.1.txt", sep = '\t')

rc.df <- rc.df[,colSums(is.na(rc.df))<nrow(rc.df)]
rc.df <- rc.df %>% mutate(SNP = paste0(X..chr, "_", pos))



colnames(rc.df) <- c(colnames(rc.df)[1:9], paste(samples.df[1:6, 1]),  "SNP")

rc.df %>% head(3)


# PCA ---------------------------------------------------------------------

rc.2.df <- rc.df[, 10:length(colnames(rc.df))] 
rownames(rc.2.df) <- rc.2.df$SNP

rc.3.df <- rc.2.df %>% select(-SNP) %>% t()

rc.3.df[, 1:10]



colores.gp <- data.frame( group = c("EB", "LT"), 
                          col = c( "black",  "orange" ))

groups <- data.frame( group = row.names(rc.3.df) ) %>% separate( group, sep = "\\-", into = c("type", "group", "rep") ) %>%
  merge(colores.gp, by = "group")

groups$group <- factor( groups$group, levels = c("EB", "LT") )


fviz_pca_ind(res.pca,
             col.ind = groups$group, # color by groups
             palette = colores.gp$col,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
