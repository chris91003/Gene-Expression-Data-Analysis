setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week9")

#1
library(fibroEset)
data(fibroEset)
head(fibroEset$species)
#Levels: b g h

#2
dat <- exprs(fibroEset)
ann <-fibroEset$species

df_sample <- dat[sample(nrow(dat), 50), ]
df_sample <- as.matrix(df_sample)
colnames(df_sample) <- ann

#3
man_dist <- dist(t(df_sample), method = "manhattan")
hc <- hclust(man_dist, method = "median")



plot(
  hc, labels = colnames((df_sample)), sub="", xlab="", 
  main = "Hierarchical Clustering  fibroEset data"
) 

#4
hm.rg <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")
heatmap(as.matrix(df_sample),col=hm.rg, main = "Heatmap Genes Vs Samples")


heatmap(df_sample, distfun = function(x) dist(x, method ="manhattan"), hclustfun=function(d) 
  hclust(d, method="median"), main = "Heatmap of 50 random genes from fibroEset data", xlab = "sample", ylab = "genes")

#5


# Kernel PCA
library(kernlab)
# calculate PCA for first 2 PCs

dat.pca <- prcomp(t(df_sample), scale= T)
dat.2 <- dat.pca$x[,1:2]
k_pca <- kmeans(dat.2, centers = 3, iter.max = 20)

dat.3 <- dat.pca$rotation[,1:2]

#6
plot(dat.2, col= k_pca$cluster, xlab= "PC1",ylab="PC2",pch=16,cex=1,main="Kernel PCA of FibroEset Data")
points(k_pca$centers, col = 1:2, pch = 8, cex = 2)

