library(ggplot2)
library(gplots)
library(ggcorrplot)
library(impute)

#Hw 1

#1
setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Hw")

renal <- read.table("renal_cell_carcinoma.txt", header= T,na.strings = "NA", 
                    blank.lines.skip = F,row.names = 1)

dim(renal)

dat <- read.table("renal_cell_carcinoma.txt", header= T,na.strings = "NA", 
                    blank.lines.skip = F,row.names = 1)


#2
renal_labels <- read.table(text = gsub(",", "\t", readLines("renal_carcinoma_annotation.txt")))

tumor_id <- renal_labels$V9

#Get Original colnames
renal_names <- colnames(renal)

colnames(renal) <- paste(colnames(renal), tumor_id, sep="_")

head(renal)

#4
cor_matrix <- cor(renal)

ggcorrplot::ggcorrplot(cor(renal))

library(gplots)
cor_renal <- cor(renal)

#Correlation Plot
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(cor_renal,na.rm=T),max(cor_renal,na.rm=T),length=10)
par(mar = rep(2, 4))

image(cor_renal,main="Correlation plot Normal/Tumor data",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(cor_renal)),label=dimnames(cor_renal)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(cor_renal)),label=dimnames(cor_renal)[[2]],cex.axis=0.9,las=2)
image(as.matrix(leg), col = cx, axes = T) 
tmp <- round(leg, 2) 
axis(1, at = seq(0, 1, length = length(leg)), labels = tmp, cex.axis = 1) 



#Cluster Dendogram
dat <- t(renal)
dat.dist <- dist(dat, method= "euclidean")
dat.clust <- hclust(dat.dist, method= "single")
plot(dat.clust, labels=names(dat), cex= 0.75, xlab = "Samples")

#Cv vs Mean PLot
renal_mean <- apply(log2(renal),2,mean)		# calculate mean for each sample
renal.sd <- sqrt(apply(log2(renal),2,var))		# calculate st.deviation for each sample
renal.cv <- renal.sd/renal_mean			#calculate cv

plot(renal_mean,renal.cv,main="Renal dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.1,type="n")
points(renal_mean,renal.cv,bg="lightblue",col=1,pch=21)
text(renal_mean,renal.cv,label=dimnames(renal)[[2]],pos=1,cex=0.5)


# average correlation plot
renal.avg <- apply(cor_matrix,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(renal.avg)),range(renal.avg),type="n",xlab="Samples",ylab="Avg r",main="Avg correlation of Tumor/Normal samples",axes=F)
points(renal.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(renal.avg)),labels=dimnames(renal)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


#6
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("impute")


#7 Remove Outliers
#GSM146798 GSM146799

new_renal <- subset(renal, select= -c(GSM146798_Normal, GSM146799_Tumor))

#8

#Load in Data
KNG1 <- read.csv("KNG1.csv")
AQP2 <- read.csv("AQP2.csv")

#Combine to one dataset
KNG1 <- new_renal["206054_at", ] 
KNG2 <- new_renal["217512_at", ] 
AQP  <- new_renal["206672_at", ] 
combined_probes <- rbind(KNG1, KNG2, AQP) 

#Get profile plot for each probeset

at_206054 <- subset(KNG1, select= -c(X217512_at))

at_217512 <- subset(KNG1, select= -c(X206054_at))

genes_combined <- cbind(KNG1, AQP2$X206672_at)

plot(c(1,ncol(new_renal)),range(combined_probes),type='n',main="Profile plot of Samples",xlab="Samples",ylab="Expression",axes=F)
axis(side=1,at=c(1:ncol(combined_probes)),labels=dimnames(new_renal)[[2]],cex.axis=0.4,las=2)
axis(side=2)
for(i in 1:nrow(combined_probes)) {
  dat.y <- as.numeric(combined_probes[i,])
  lines(c(1:ncol(new_renal)),dat.y,col=i,lwd=2)
}

#Tumor samples have much lower expression than the normal samples


#9
new_renal <- as.matrix(new_renal)
renal2 <- new_renal["206054_at", ]
actual_value <- new_renal["206054_at", "GSM146784_Normal"]

new_renal["206054_at", "GSM146784_Normal"] <- NA 


#10
impute_renal <-impute.knn(new_renal, k= 6) 
k <- impute_renal$data["206054_at",]
head(impute_renal)

#11
impute_value <- impute_renal$data["206054_at", "GSM146784_Normal"]

error <- abs(actual_value-impute_value)/(actual_value)

error

#12

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("pcaMethods")

library(pcaMethods)

SVD_renal <- pca(new_renal, method= "svdImpute", nPcs= 9)
SVD_output <- completeObs(SVD_renal) 

summary(SVD_renal)

SVD_output["206054_at", "GSM146784_Normal"]

SVD <- SVD_output["206054_at", ] 
head(SVD)

#13

combined_values <- rbind(renal2, k, SVD) 

plot(c(1,ncol(new_renal)),range(combined_values, na.rm = TRUE),type='n',main="Profile plot of Renal Samples",xlab="Samples",ylab="Expression",axes=F)
axis(side=1,at=c(1:ncol(new_renal)),labels=dimnames(combined_values)[[2]],cex.axis=0.4,las=2)
axis(side=2)
for(i in 1:nrow(combined_values)) {
  dat.y <- as.numeric(combined_values[i,])
  lines(c(1:ncol(combined_values)),dat.y,col=i,lwd=2)
}



results <- list(original= as.numeric(renal2), 
                knn= as.numeric(impute_renal$data["206054_at",]), svd= as.numeric(SVD_output["206054_at",]))

plot(c(1,length(renal2)),range(renal2, na.rm = TRUE),type='n',main="Profile plot of 8 random genes",xlab="Samples",ylab="Expression",axes=F)
axis(side=1,at=c(1:length(renal2)),labels=dimnames(renal2)[[2]],cex.axis=0.4,las=2)
axis(side=2)

for(i in 1:length(results)) {
  dat.y <- results[[i]]
  lines(c(1:length(renal2)),dat.y,col=i,lwd=2)
}

