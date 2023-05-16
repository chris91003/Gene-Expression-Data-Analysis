setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week4")


#2
eisen <- read.table("eisen.txt", header= T,na.strings = "NA", blank.lines.skip = F,row.names = 1)
head(eisen)


#3
eisen_class <- read.table("eisenClasses.txt", header= T)
head(eisen_class)

#4

colchar <- as.character(eisen_class[,2])

eisen <- eisen[, colchar]

#Subset
eisen1 <- colchar[1:19]
eisen2 <- colchar[20:39]

#5 Separate X and Y variables and plot


x <- as.numeric(eisen["100", eisen1])

x <- x[!is.na(x)]
x

y <- as.numeric(eisen["100", eisen2])
y <- y[!is.na(y)]
y

boxplot(x, y, col=c("red", "blue"), xlab= "Class 1 vs Class 2", ylab= "Expression Values" )

#Histogram
par(mfrow=c(2,1)) 
hist(x, xlab= "Class 1", ylab= "Intensity", col = ("red"))
hist(y, xlab= "Class 2", ylab= "Intensity", col= ("blue"))


#6

# size of each group
nx <- length(x)
ny <- length(y)

# pooled variance
pool.var <- (((nx-1)*var(x)) + ((ny-1)*var(y)))/(nx+ny-2)
pool.var

# sample size calculation based on a 1.5-fold difference between means
dif.15fold <- log2(1.5)/sqrt(pool.var)
pl.ss3 <- power.t.test(d=dif.15fold,sig.level=.01,power=0.8,type="two.sample")
pl.ss3



#7 delta = 2.088773

gene_size <- power.t.test(n= NULL, delta = 2.088773, sig.level=.01,power=0.8,type="two.sample")
gene_size



#8

library(ssize)
library(gdata) 

sd <- apply(eisen,2, sd, na.rm= T)

new_sd<- transform(eisen, apply(eisen,2, sd, na.rm= T) )
hist(sd,n= 15, xlab= "Standard Deviation for data on log2 scale", ylab= "Counts", main= "Histogram of SD Data")

#9

fold.change= 3
all.size <- ssize(sd=sd, delta=log2(3), sig.level=.05, power=0.8) 
ssize.plot(all.size, lwd=2, col="magenta", xlim=c(1,20)) 
xmax <- par("usr")[2]-1; 
ymin <- par("usr")[3] + 0.05 
title("Sample Size to Detect 3-Fold Change")  
