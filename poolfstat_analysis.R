library(poolfstat)
library(ggplot2)
library(GGally)
library(factoextra)
# sample names ------------------------------------------------------------

samples.df <- read.csv("sample_names.txt", header = F, col.names = c("Sample"))


# GBS: reading data ------------------------------------------------------------

topw <- vcf2pooldata(vcf.file = "GBS/POPOOLATION2/gbs_chr1.filtered.vcf",
                     poolsizes = rep(30,42))  


topw.fst <- computeFST(topw, sliding.window.size = 50)
topw.fst$FST

topw.fst$sliding.windows.fst$CumulatedPosition/1e6 %>% head(3)
topw.fst$sliding.windows.fst$MultiLocusFst %>% head(3)
plot(topw.fst$sliding.windows.fst$CumulatedPosition/1e6,
     topw.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
     col='black',pch=16)
abline(h=topw.fst$FST,lty=2)

topw.pairwisefst <- compute.pairwiseFST(topw, verbose = F)
pdf("chr1_gbs_heatmap_fst.pdf")
heatmap(topw.pairwisefst)
dev.off()


# GBS: PCA ----------------------------------------------------------------
frq <- topw@refallele.readcount/topw@readcoverage 

colnames(frq) <- topw@poolnames

frq %>% head(3)

res.gbs.pca <- prcomp(frq %>% t(), scale = TRUE)


p1 <- fviz_pca_ind(res.gbs.pca,
             col.ind = c(
               rep("EB", 3), rep("LT", 3), rep("M1-50", 3), rep("M1-75", 3),  rep("M2-50", 3),  rep("M2-75", 3),  rep("M3-50", 3),  rep("M3-75", 3), 
               # "Alg", "Ara", "Ari", "Art", "Rep", 
               rep("Alg", 3),rep("Ara", 3),rep("Arc", 3),rep("Ari", 3),rep("Art", 3),rep("Rep", 3)
               
             ) , # color by groups
             #  palette = groups$col,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = TRUE
             
) + ggtitle("GBS", subtitle = sprintf("PCA with %d SNPs", topw@nsnp))# + ggp.cust

ggsave(p1, filename = "GBS_chr1_PCA.pdf", dpi = 'print')
# AMPSEQ: reading data ----------------------------------------------------


as.topw <- vcf2pooldata(vcf.file = "AMPLICON_SEQUENCING/POPOOLATION2/ampseq_30_all.filtered.DP60.vcf",
                     poolsizes = rep(30,48))  

as.topw@poolnames <- samples.df[c(1:48),] 
as.topw.subset<-pooldata.subset(as.topw,pool.index=c(1:18, 20:24, 31:48), min.cov.per.pool = 60, verbose=FALSE )
as.topw.pure<-pooldata.subset(as.topw,pool.index=c(25,26,28:48), min.cov.per.pool = 60, verbose=FALSE )


as.topw.pure.pairwisefst <- compute.pairwiseFST(as.topw.pure, verbose = F)
heatmap(as.topw.pure.pairwisefst)





as.topw.fst <- computeFST(as.topw.subset, sliding.window.size = 50)
as.topw.fst$FST

plot(as.topw.fst$sliding.windows.fst$CumulatedPosition/1e6,
     as.topw.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst",
     col='black',pch=16)
abline(h=as.topw.fst$FST,lty=2)



# AMPSEQ: PCA --------------------------------------------------------------

as.frq <- as.topw.subset@refallele.readcount/as.topw.subset@readcoverage 
colnames(as.frq) <- as.topw.subset@poolnames

as.frq %>% head(3)

res.pca <- prcomp(as.frq %>% t(), scale = TRUE)

as.topw.subset@nsnp
p2 <- fviz_pca_ind(res.pca,
             col.ind = c(
               rep("EB", 3), rep("LT", 3), rep("M1-50", 3), rep("M1-75", 3),  rep("M2-50", 3),  rep("M2-75", 3),  rep("M3-50", 2),  rep("M3-75", 3), 
               # "Alg", "Ara", "Ari", "Art", "Rep", 
               rep("Alg", 3),rep("Ara", 3),rep("Arc", 3),rep("Ari", 3),rep("Art", 3),rep("Rep", 3)
               
             ) , # color by groups
          #  palette = groups$col,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Groups",
             repel = TRUE
             
) + ggtitle("Amplicon sequencing", subtitle = sprintf("PCA with %d SNPs", as.topw.subset@nsnp))# + ggp.cust

ggsave(p2, filename = "AMPSEQ_PCA.pdf", dpi = 'print')
