setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week8")

#1
dat = read.table(file = "Sotiriou.txt", row.names = 1, header = TRUE)

ann = read.table(file= "Sotiriou_annotations.txt", row.names = 1, header = TRUE)

#2

dat.pca <- prcomp(t(dat),cor=F)

# plot all 3 components in 3 plots
dat.loadings <- dat.pca$x[,1:2]
plot(range(dat.loadings[,1]),range(dat.loadings[,2]),type="n",xlab='p1',ylab='p2',main='PCA plot of Breast Cancer Data KIU vs OXF')
points(dat.loadings[,1][ann$site== "KIU"], dat.loadings[,2][ann$site == "KIU"],col=1,bg='red',pch=21,cex=1.5)
points(dat.loadings[,1][ann$site == "OXF"], dat.loadings[,2][ann$site == "OXF"],col=1,bg='blue',pch=21,cex=1.5)
legend(-32,-17,c("KIU","OXF"),col=c("red","blue"),pch=15,cex=.9,horiz=F)

#3

dat.pca.var <- round(dat.pca$sdev^2 / sum(dat.pca$sdev^2)*100,2)
plot(c(1:length(dat.pca.var)),dat.pca.var,type="b",xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1.5)
title("Scree plot showing % variability explained by each eigenvalue\n Cancer Dataset")

#4
library(MASS)
library(multtest)

dat.dist <- dist(t(dat))


# classical metric MDS on samples (no stress value provided)
dat.loc <- cmdscale(dat.dist)
plot(dat.loc, type = "n")
points(dat.loc[,1][ann$site== "KIU"], dat.loc[,2][ann$site== "KIU"],col= "red",pch=16,cex=1.5)
points(dat.loc[,1][ann$site == "OXF"], dat.loc[,2][ann$site == "OXF"],col="blue",pch=16,cex=1.5)
title(main="Classic MDS plot of Breast Cancer data ")
legend(-32,15,c("KIU","OXF"),col=c("red","blue"),pch=15,cex=.7,horiz=F)



dat.mds <- isoMDS(dat.dist)
plot(dat.mds$points, type = "n")
points(dat.mds$points[,1][ann$site== "KIU"], dat.mds$points[,2][ann$site== "KIU"],col= "red",pch=16,cex=1.5)
points(dat.mds$points[,1][ann$site == "OXF"], dat.mds$points[,2][ann$site == "OXF"],col="blue",pch=16,cex=1.5)
title(main="Nonmetric MDS plot of Breast Cancer data-stress=20% ")
legend(30,15,c("KIU","OXF"),col=c("red","blue"),pch=15,cex=.7,horiz=F)


#5

temp <- t(dat)
temp <- scale(temp,center=T,scale=T) 


# The weighted graph Laplacian
k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    n <- attr(dis, "Size")
    full <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}
phi <- k.speClust2(t(temp),qnt=0.005)
plot(range(phi[,1]),range(phi[,2]),xlab="phi1",ylab="phi2",main="Weighted Graph Laplacian plot of Breast Cancer Data\nepsilon=0.005")
points(phi[,1][ann$site== "KIU"],phi[,2][ann$site== "KIU"],col="red",pch=16,cex=1.5)
points(phi[,1][ann$site == "OXF"],phi[,2][ann$site == "OXF"],col="blue",pch=16,cex=1.5)
legend("bottomleft",c("KIU","OXF"),col=c("red","blue"),pch=15,cex=.7,horiz=F)




