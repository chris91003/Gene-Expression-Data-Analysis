setwd("C:/Users/cdima/OneDrive/Desktop/Gene Expression Data Analysis Visualization/Week10")

#1
d = read.table(file = "lung_cancer_data.txt", row.names = 1, header = TRUE)

#2
library(MASS)

new_name= colnames(d)[1:length(colnames(d))]
names <- gsub("\\d+", "", new_name)

data  <- data.frame(names, t(d))
dim(data)



#3

data1 <- data[1:6,]
data2 <- data[11:16,]
data3 <- data[20:22,]

#Create Train dataset
train <- rbind(data1, data2, data3)

test  <- rbind(data[7:10,], data[17:19, ], data[23:24,])

#Get Names
test_names <- test$names

#Remove first column

test = subset(test, select = -c(names) )

#4

#Get train names
clas <- train$names
#clas  <- factor(gsub('[[:digit:]]+', '', clas))





#Get first two training samples
train_2 <- train[, c(2,3)]

#Run Lda
dat.lda <- lda(clas~ .,train_2)
dat.pred <- predict(dat.lda,test[,c(1,2)])


#Confusion Matrix
table(dat.pred$class,test_names)

#5 labels were misclassified

#5

plot(range(dat.pred$x[, 1]), range(dat.pred$x[, 2]), type = "n", xlab = "LD1",ylab = "LD2",
  main = "LDA 2 Genes",)
points(dat.pred$x, col = c(rep("Green", 4), rep("Blue", 3), rep("Red", 2)),)
legend("bottomright", c("Adeno", "SCLC", "Normal"), 
       col = c("Green", "Blue", "Red"),
       pch = c(1:3))


plot(dat.pred$x,bg=as.numeric(factor(clas)),pch=21,col=1,ylab="Discriminant function",axes=F,xlab="Score",main="Discriminant function for Lung dataset 2 Samples")
axis(1,at=c(1:24),names(d),las=3,cex.axis=0.5)
axis(2)


plot(dat.pred$x,col= dat.pred$class, cex= 1,axes=T,xlab="LD1", ylab= "LD2",main="Discriminant function for Lung dataset")
legend("bottomright", c("Adeno", "SCLC", "Normal"), col = c("Green", "Blue", "Red"), pch = c(1:3))


#6 and 7
train = subset(train, select = -c(names) )

dat.lda2 <- lda(clas~ ., train)
dat.pred2 <- predict(dat.lda2,test)

#Confusion Matrix
table(dat.pred2$class,test_names)


#No misclassified samples


plot(range(dat.pred2$x[, 1]), range(dat.pred2$x[, 2]), type = "n", xlab = "LD1",ylab = "LD2",
     main = "LDA 2 Genes",)
points(dat.pred2$x, col = c(rep("Green", 4), rep("Blue", 3), rep("Red", 2)),)
legend("bottomright", c("Adeno", "SCLC", "Normal"), 
       col = c("Green", "Blue", "Red"),
       pch = c(16:18))


plot(dat.pred2$x,col= dat.pred2$class, cex= 1,axes=T,xlab="LD1", ylab= "LD2",main="Discriminant function for Lung dataset All Data")
legend("bottomright", c("Adeno", "SCLC", "Normal"), col = c("Green", "Blue", "Red"), pch = c(1:3))















