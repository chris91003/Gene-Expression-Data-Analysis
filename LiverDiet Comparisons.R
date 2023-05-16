library(GEOquery)

#"Analysis of liver from C57Bl/6J males fed isocaloric diets containing sucrose or palatinose. Sucrose is a high glycemic index (GI) sugar; 
#palatinose is a low GI sugar. Results provide insight into molecular mechanisms underlying the role of sucrose in the development of 
#non-alcoholic fatty liver."



#Column 1-7 Sucrose Containing Diet, Column 8-14 Palatinose Containing Diet

dat= getGEO("GDS5435")	
data= Table(dat)

gen_names <- data$IDENTIFIER
dat_mat <- data[, -c(1:2)]


#Mean vs CV Plot

dat.mean <- apply(log2(dat_mat),2,mean)		# calculate mean for each sample
dat.sd <- sqrt(apply(log2(dat_mat),2,var))		# calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean			#calculate cv

plot(dat.mean,dat.cv,main="Diabetes dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat_mat)[[2]],pos=1,cex=0.5)

#Average Correlation PLot
dat.cor <- cor(dat_mat)
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of sucrose/palatinose samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat_mat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

#PCA Plot

dat.pca <- prcomp(t(dat_mat))
dat.loads <- dat.pca$x[,1:2]
plot(dat.loads[,1],dat.loads[,2],main="Sample PCA plot",xlab="p1",ylab="p2",col='red',cex=1,pch=15)
text(dat.loads,label=dimnames(dat_mat)[[2]],pos=1,cex=0.5)


#Outlier Removal Based on Plots
dat <- dat_mat[ , -which(names(dat_mat) %in% c("GSM1322815", "GSM1322818"))]


#Filter genes using a criterion (Here we filtered out expression < 25%)
quantile(log2(rowMeans(dat)))

filtered_dat <- subset(dat, log2(rowMeans(dat)) > 3.866249)


#Ttest
filtered_dat <- log2(filtered_dat)
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

pv <- apply(filtered_dat,1,t.test.all.genes,s1=c(1:6),s2=c(7:12))


par(mfrow=c(1,2))
hist(pv,col="lightblue",xlab="p-values",main="P-value distance between Sucrose and Palatinose Groups",cex.main=0.9)
abline(v=.05,col=2,lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)", main="-log10(pv) dist'n between\nSucrose and Palatinose Groups",cex.main=0.9)
abline(v= -log10(.05),col=2,lwd=2)



#Get P-values < threshold and subset
new_pv <- pv[pv <.01]
length(new_pv)
hist(new_pv, xlab= "P-Values", main= "Histogram of Log P-Values after Threshold")



dat <- filtered_dat[pv<.01,]
hist(as.numeric(unlist(dat)), xlab= "P-Values", main= "Histogram of P-Values after Threshold")
hist(as.numeric(unlist(logdat)), xlab= "P-Values", main= "Histogram of Log P-Values after Threshold")


#Run PCA
dat.pca <- prcomp(t(dat),cor=F)

# plot all 3 components in 3 plots
dat.loadings <- dat.pca$x[,1:3]
colors= c("red", "blue")
plot(dat.loadings, col = colors) 

#Perform Kmeans

# kmeans clustering with random data
library(mva)
cl <- kmeans(dat, centers=4, iter.max=20)
plot(dat, col = cl$cluster,cex=1)
points(cl$centers, col = 1:4, pch = "*",cex=2.5)


#Get Minimum and Maximum Values
head(sort(new_pv, decreasing=FALSE), 5)
