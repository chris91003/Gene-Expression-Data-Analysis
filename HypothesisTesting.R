setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week6")

#1 and 2
rat = read.table("rat_KD.txt", header= T, row.names = 1)
rat

#3
new_rat <- log2(rat)
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

names(rat)
pv <- apply(new_rat,1,t.test.all.genes,s1=c(1:6),s2=c(7:11))


#4

pv_df <- as.data.frame(pv)
hist(pv)

pv1= sum(pv_df$pv < 0.05)
print(pv1)

pv2= sum(pv_df$pv < .01)
print(pv2)

threshold = .05/length(pv)

#Get p values <  3.140112e-06
pv3= sum(pv < threshold)


#5
#Get mean of group
control = apply(new_rat[,1:6], 1, mean)

#get mean of test group
keto = apply(new_rat[, 7:11], 1, mean) 

#6
foldchange <- control - keto 
hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

linear= 2^(foldchange)
min_foldchange= min(linear)
max_foldchange= max(linear)

pv[pv < threshold & foldchange >2]



#7

#1367553_x_at   1370239_at 1370240_x_at 1371102_x_at 1388608_x_at 
#hemoglobin subunit. They all seem to be subunits of hemoglobin


#8
# transpose p-values
t.test.run <- apply(new_rat,1,t.test.all.genes,s1=c(1:6),s2=c(7:11))
fold <- apply(log(rat[,c(1:6)]),1,mean) - apply(log(rat[,c(7:11)]),1,mean)
p.trans <- -1 * log(t.test.run)


# volcano plot #1
plot(range(p.trans),range(foldchange),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot Keto and Control differences')
points(p.trans,foldchange,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&foldchange>log2(2))],fold[(p.trans> -log10(.05)&foldchange>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&foldchange< -log2(2))],foldchange[(p.trans> -log10(.05)&foldchange< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))


