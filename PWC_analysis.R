library(tidyverse)
library(ggplot2)
library(reshape2)


# load data ---------------------------------------------------------------
samples.df <- read.csv("sample_names.txt", header = F, col.names = c("Sample"))

rc.df <- read.csv("GBS/POPOOLATION2/F_rc.cov60.txt", sep = '\t')
rc.df <- rc.df %>% mutate(SNP = paste0(X..chr, "_", pos))

rc.2.df <- rc.df[, 10:length(colnames(rc.df))] 

colnames(rc.2.df) <- c(paste(samples.df[1:6, 1], "_M"), paste0(samples.df[1:6, 1], "_m"), "SNP")

rc.3.df <- rc.2.df %>% melt(id.vars = "SNP") %>%
  separate(variable, sep = "_", into = c("Sample", "Allele.Type")) %>%
  filter(Allele.Type == 'm') %>%
  group_by(SNP) %>%
  summarise(Min.MAF = min(value, na.rm = T)) %>%
  filter(Min.MAF > 0.01)

#rc.3.df 

pwc.df <- read.csv("GBS/POPOOLATION2/F_pwc", sep = '\t') %>% 
  mutate(SNP = paste0(X..chr, "_", pos)) %>% 
  filter(SNP %in% rc.3.df$SNP) %>% 
  filter(allele_count == 2)

pwc.df  %>% head(3)




# transform pwc ----------------------------------------------------------

## Creating a data frame with the names of the sample comparisons

sample_iterator_1 = function(x) {
  
  v <- c()
  j = 1
  for (i in sort(c(1:5), decreasing = T )) {
    
    v <- c(v, rep(x[j, 1], i))
 #   print(i)
    j = j + 1
  }
  
  return(v)
}

sample_iterator_2 = function(x) {
  
  v <- c()
  
  for (i in c(2:6)) {
    
    v <- c(v, x[i:6, 1])
 #   print(i)
    
  }
  
  return(v)
}



samples.2.df <- data.frame(
  comparison = colnames(pwc.df)[10:length(colnames(pwc.df))-1],
  S1 = sample_iterator_1(samples.df),
  S2 = sample_iterator_2(samples.df)
  
  
)

sample_iterator_3 = function(x) {
  
  v <- c()
  j = 1
  for (i in sort(c(1:6), decreasing = T )) {
    
    v <- c(v, rep(x[j, 1], i))
    #   print(i)
    j = j + 1
  }
  
  return(v)
}

sample_iterator_4 = function(x) {
  
  v <- c()
  
  for (i in c(1:6)) {
    
    v <- c(v, x[i:6, 1])
    #   print(i)
    
  }
  
  return(v)
}


samples.3.df <- data.frame(
  
  S1 = sample_iterator_3(samples.df),
  S2 = sample_iterator_4(samples.df)
  
  
) %>% filter(S1 == S2) %>%
  mutate(diff.SUM = 0)



## Generating a data frame with the sum of the allele frequency differences for all SNPs grouped by S1 and S2

pwc.2.df <- pwc.df[, c(length(colnames(pwc.df)), 10:length(colnames(pwc.df))-1)] %>%
  melt(id.vars = c("SNP")) %>%
  rename(
    
    "comparison" = "variable",
    "diff" = "value"
    
  ) %>%
  merge(samples.2.df, by = "comparison", all = T) %>%
  group_by(S1, S2) %>%
  summarise(diff.SUM = sum(diff, na.rm = T)) %>%
  rbind(samples.3.df)

# transform to matrix -----------------------------------------------------


pwc.3.df <- pwc.2.df %>% 
  dcast(S1 ~ S2, fun.aggregate = sum, value.var = 'diff.SUM') #%>% select(-S1)

row.names(pwc.3.df) <- pwc.3.df$S1

pwc.3.df <- pwc.3.df %>% select(-S1) %>% as.matrix()
pwc.3.df[lower.tri(pwc.3.df, diag = F)] <- t(pwc.3.df)[lower.tri(pwc.3.df, diag = F)]

# NMDS prep ---------------------------------------------------------------

colores.gp <- data.frame( group = c("EB", "LT"), 
                          col = c( "black",  "orange" ))
groups <- data.frame( group = colnames(pwc.3.df) ) %>% separate( group, sep = "\\-", into = c("type", "group", "rep") ) %>%
  merge(colores.gp, by = "group")

groups$group <- factor( groups$group, levels = c("EB", "LT") )

# NMDS --------------------------------------------------------------------
library(vegan)
fit <- cmdscale(pwc.3.df)
x <- fit[, 1]

x %>% as.data.frame()
y <- fit[, 2]

p1.df <- data.frame(Sample = rownames(as.data.frame(x)), x = as.data.frame(x)$x, y = as.data.frame(y)$y ) %>%
  separate(Sample, sep = '-', into = c("Type", "Site", "Rep"))


hulls.p <- p1.df %>%
  group_by(Site) %>%
  slice(chull(x, y))

# Update the plot with a fill group, and overlay the new hulls



ggplot(data = p1.df, aes(x, y)) + geom_point(size = 4, aes(col = Site)) + 
  aes(fill = factor(Site)) + geom_polygon(data = hulls.p, alpha = 0.5)+ ggp.cust +
  scale_color_manual(values = c("gray75", "orangered")) +
  scale_fill_manual(values = c("gray75", "orangered"))


plot(fit[ , 1], fit[ , 2], xlab = "Dimension 1", ylab = "Dimension 2")



pwc.nms <- metaMDS( pwc.3.df, distance = 'euclidean', k = 3, wascores = F)


orditorp( pwc.nms, "sites" , cex = 0.5)
#dev.off()

ordihull(
  pwc.nms,
  groups$group,
  #island.spp_groups$habitat,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  #border = c("gray0", "gray0", "gray48", "gray48"),
  # lty = c(1, 2, 1, 2),
  lwd = 2.5
)
scrs <-  scores(pwc.nms, display = "sites")
cent <-
  aggregate(scrs ~ groups$group, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c( 8 , 8 , 8, 8, 8),
       col = c("#d1495b", "#edae49", "#8d96a3",  "#00798c","#2e4057"),
       bg = c("black"),
       lwd = 3.0,
       cex = 2.0 # Plots centroids as points on ordination
)


