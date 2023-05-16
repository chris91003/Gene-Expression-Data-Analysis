setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week7")


#2
affy = read.table(file = "agingStudy11FCortexAffy.txt", row.names = 1, header = TRUE)

info = read.table(file= "agingStudy1FCortexAffyAnn.txt", row.names = 1, header = TRUE)


#3
sorted_age = info[order(info$Age),]
new_info <- paste(dimnames(sorted_age)[[1]], sorted_age[, 2], sorted_age[, 1], sep = ".") 	

data = affy[, new_info]

# gender comparison gene vector
g.g <- c(1394,  1474,  1917,  2099,  2367,  2428, 2625,  3168,  3181,  3641,  3832,  4526,
         4731,  4863,  6062,  6356,  6684,  6787,  6900,  7223,  7244,  7299,  8086,  8652,
         8959,  9073,  9145,  9389, 10219, 11238, 11669, 11674, 11793)

# age comparison gene vector
g.a <- c(25, 302,  1847,  2324,  246,  2757, 3222, 3675,  4429,  4430,  4912,  5640, 5835, 5856,  6803,  7229,  7833,  8133, 8579,  8822,  8994, 10101, 11433, 12039, 12353,
         12404, 12442, 67, 88, 100)


#4


t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative= c("two.sided"), var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}


# s1 and s2 are dimensions of the two samples
# run function on each gene in the data frame
gender <- apply(affy,1,t.test.all.genes,s1=c(1:18),s2=c(19:30))
age <- apply(data,1,t.test.all.genes,s1=c(1:12),s2=c(13:30))

gender <- apply(affy[g.g, ],1,t.test.all.genes,s1=info$Gender == "M",s2=info$Gender == "F")
age <- apply(data[g.a,],1,t.test.all.genes,s1= sorted_age$Age < 50, s2= sorted_age$Age> 50)



#Apply holms method
library(base)
p.gender <- p.adjust(gender,method= "holm")
p.age <- p.adjust(age,method="holm")


#5

sort.gen = sort(p.gender)
sort.age = sort(p.age)
gender1= sort(gender)
age1= sort(age)

gender_sort = cbind(gender1, sort.gen)
matplot(gender_sort, type= "l", pch = 1, col= 1:2, main= "Gender Plot Holm", ylab= "P- Values")

age_sort= cbind(age1, sort.age)
matplot(age_sort, type= "l", pch = 1, col= 1:2, main= "Age Plot Holm", ylab= "P- Values" )




#6

#Apply Bonfernni method


library(base)

b.gender <- p.adjust(gender,method= "bonferroni")
b.age <- p.adjust(age,method="bonferroni")

sort.gen = sort(p.gender)
sort.age = sort(p.age)
gender2= sort(gender)
age2= sort(age)

gender_sortb = cbind(gender2, b.gender)
matplot(gender_sortb, type= "l", pch = 1, col= 1:2, main= "Gender Plot Bonferonni", ylab= "P- Values")

age_sortb= cbind(age2, b.age)
matplot(age_sortb, type= "l", pch = 1, col= 1:2, main= "Age Plot Bonferroni", ylab= "P- Values")


#7
tcga = read.table(file = "tcga_brca_fpkm.txt", row.names = 1, header = TRUE)

sam = read.table(file= "tcga_brca_fpkm_sam.txt", row.names = 1, header = TRUE, fill = TRUE)

#8
gata = grep("GATA3", tcga)

mat1= tcga[grep(pattern="GATA3", rownames(tcga)), ]
mat1= as.numeric(mat1)

#get percentile
quantile(mat1)

#9
##vector where patients are upper 25% 
group = ifelse(mat1 > 12.94155, 1, 0 )

#10
sam$group = group
head(sam)

#11

#Set status variable 
sam <- ifelse(sam$vital_status == "Living", 0, 1)

sam$status[sam$vital_status == "LIVING"] <- 0
sam$status[sam$vital_status == "DECEASED"] <- 1

surv_gata <- survdiff(formula= Surv(months_to_event, status) ~ group, data= sam)

#Extract p-value
library(broom)
diff_pvalue <- broom::glance(surv_gata)$p.value

#12
fit <- coxph(Surv(months_to_event, status)~group,data =sam)
summary(fit)
plot( survfit( fit),xlab="Months",ylab="S(t)",main="Cox PH Curve")

#Extract p-value and hazard ratio
cox_pvalue <- summary(fit)[['coefficients']][,'Pr(>|z|)']
cox_pvalue

cox_pvalue <-  0.08162176
cox_ratio <- broom::tidy(fit)

cox_ratio <- exp(fit$coefficients)
cox_ratio <- as.numeric(cox_ratio)

#13
survival2 <- survfit(Surv(months_to_event, status)~group,data =sam)
summary(survival2)
plot(survival2, lty = 1:2,col = 1:2, xlab="Weeks",ylab="S(t)")
legend(85, .9, c("Low Expression", "High Expression"), col = 1:2, lty= 1:2)
legend(1, .2, legend=rev(mylegend), pch = 15,cex=.6)
title("Kaplan-Meier Curves\nfor GATA3 Study")

group0 <- 86
group1 <- 30

#Proive ratios and p-values

mylegend <- c(paste0("Surv Diff PValue", diff_pvalue),
              paste0("Cox PValue ", cox_pvalue),
              paste0("Cox Ratio ", cox_ratio), 
              paste0(" Group 0 Size ", group0),
              paste0("Group 1 size ", group1))

#14
#This Result does match the study because high GATA3 resulted in higher surival

