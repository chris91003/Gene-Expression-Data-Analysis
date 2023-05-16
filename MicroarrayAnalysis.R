setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week5")

library("marray")
data(swirl)

#2
par(mfrow=c(1,2))

maPlot(swirl[,3], lines.func = NULL, legend.func = NULL, main= "Regular Plot") 


#3

mnorm <- maNorm(swirl[,3], norm="median")
mnorm

#4

maPlot(mnorm,lines.func = NULL, legend.func = NULL, main= "Normalized Plot")

#5 The difference is that the middle line (0 line) changes for the normalized plot as there are 
#some points that go to 3 on y axis 


#6
mnorm_loess <- maNorm(swirl[,3], norm="loess")
maPlot(mnorm_loess,lines.func = NULL, legend.func = NULL)

#7 
#It appears that the loess normalization is better which is more centered at 0

#8
dir.path <- "C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week5"
a.cdna <- read.GenePix(path=dir.path,name.Gf = "F532 Median",name.Gb ="B532 Median", 
                       name.Rf = "F635 Median", name.Rb = "B635 Median",name.W ="Flags")

#9
a_1 <- maNorm(a.cdna[,1], norm="none")
a_12 <- maNorm(a.cdna[,1], norm="printTipLoess")
a_13 <- maNorm(a.cdna[,1], norm="scalePrintTipMAD")

par(mfrow=c(3,1))
maPlot(a_1, lines.func = NULL, legend.func = NULL) 
maPlot(a_12, lines.func = NULL, legend.func = NULL) 
maPlot(a_13, lines.func = NULL, legend.func = NULL) 


#9-2
a_2 <- maNorm(a.cdna[,2], norm="none")
a_22 <- maNorm(a.cdna[,2], norm="printTipLoess")
a_23 <- maNorm(a.cdna[,2], norm="scalePrintTipMAD")

par(mfrow=c(3,1))
maPlot(a_2, lines.func = NULL, legend.func = NULL) 
maPlot(a_22, lines.func = NULL, legend.func = NULL) 
maPlot(a_23, lines.func = NULL, legend.func = NULL) 

#10
df <- data.frame(maGnames(a_1)@maLabels, a_1@maM, a_2@maM)
df
#11
library(affy)
library(limma)
library(affydata)
library(affyPLM)
library(fpc)

#12
dir.path <- "C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week5"
fns <- sort(list.celfiles(path=dir.path,full.names=TRUE))
data.affy <- ReadAffy(filenames=fns,phenoData=NULL)


#13

rma <- expresso(data.affy, bgcorrect.method = "rma", normalize.method = "quantiles", 
                summary.method = "medianpolish", pmcorrect.method = "pmonly")

MAS <- expresso(data.affy, bgcorrect.method = "rma", normalize.method = "quantiles", 
                summary.method = "medianpolish", pmcorrect.method = "mas")


cor_rma <- cor(exprs(rma))
cor_mas <- cor(exprs(MAS))

#cor rma has a higher overall correlation structure