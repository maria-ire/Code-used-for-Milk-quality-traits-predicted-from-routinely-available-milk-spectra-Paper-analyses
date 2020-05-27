rm(list=ls())
setwd("C:/Users/direction to the file")
list.files()

# Read in data

dataset = read.csv("file name.csv")

# Clean the dataset #load the package of interest

library("pls")
library(glmnet)  # for ridge regression
library(dplyr)   # for data cleaning
library(psych)   # for function tr() to compute trace of a matrix
library(penalized)

trait = "name of the trait to be predicted"
dataset <- dataset[complete.cases(dataset[ , trait]),]

#division of the data in 4 groups for 4-fold cross-validation
library(groupdata2)

train.fold <- fold(
  dataset,
  k = 4,
  num_col = trait,
  method = "n_dist")


fold1<- which(train.fold$.folds==1)
fold2<- which(train.fold$.folds==2)
fold3<- which(train.fold$.folds==3)
fold4<- which(train.fold$.folds==4)

# Definition of brinary trait
dataset <- mutate(dataset, binary = (ifelse(Heat_stability <= 4, 1, 0))) #insert the threshold for the division
dataset$binary <- as.factor(dataset$binary)


## creation of a table where report the results
R <- 4
out <- matrix(NA, R, 60)
colnames(out) <- c("pls.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc",
                   "dt.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc",
                   "rf.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc",
                   "bo.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc",
                   "svm.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc",
                   "lr.n","sen", "spec", "ari", "acc","pls.n","sen", "spec", "ari", "acc")


###############################################################################
#############################dataset1#########################################
# Division in train and test

training<-data.frame(dataset[-fold1,])
testing<-data.frame(dataset[fold1,])


# Generic code set up
y.train = training$binary
x.train = as.matrix(training[,56:589])
train = cbind(y.train, x.train)
y.test = testing$binary
x.test = as.data.frame(testing[,56:589])
test = cbind(y.test, x.test)
train <- as.data.frame(train)
test <- as.data.frame(test)
names(train)[1] <- "binary"
names(test)[1] <- "binary"
y.train <- as.factor(y.train)
y.test <- as.factor(y.test)
train$binary <- as.factor(train$binary)
test$binary <- as.factor(test$binary)

library(mclust)
library(caret)
library(pROC)

train <- mutate(train, binary = (ifelse(binary == 1, 0, 1)))
#######################################################################################
########################## Apply PLS DA#########################################

plsda.model <- plsda(x.train, y.train, ncomp = 2, probMethod = "bayes")
yhat_plsda <- predict (plsda.model , newdata= x.train)   #predict yhat
yhat_plsda_test <- predict (plsda.model , newdata= x.test)   #predict yhat


##measures of accuracy calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_plsda #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[1,1] <- length(x)
out[1,2] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,3] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,4] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,5] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_plsda_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,6] <- length(x)
out[1,7] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,8] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,9] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,10] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


#######################################################################################
########################## Apply Decision Tree#########################################

library(C50)

##train the model
ynew1<-as.factor(y.train)
tree_mod <- C5.0(x = x.train, y = ynew1)     #estimate decision tree model
str(tree_mod)

##predict y
yhat_dt <- predict (tree_mod , newdata= x.train)   #predict yhat
yhat_dt <- as.matrix(yhat_dt)
yhat_dt <- as.numeric(yhat_dt)
yhat_test_dt <- predict (tree_mod , newdata = x.test) #predict yhat from the test dataset
yhat_test_dt <- as.matrix(yhat_test_dt)
yhat_test_dt <- as.numeric(yhat_test_dt)


##measures of accuracy calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_dt #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[1,11] <- length(x)
out[1,12] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,13] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,14] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,15] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_dt #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,16] <- length(x)
out[1,17] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,18] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,19] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,20] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

#inser in mtry the number of variables to consider at each split
# suggestion: mtry = root square of the total number of variables

##train the model
y.train <- as.factor(y.train)
RF_mod <- randomForest (x= x.train, y= y.train, ntree = 500, mtry =23)   #model code

## predict y
yhat_RF <- predict (RF_mod , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.vector(yhat_RF)
yhat_RF_test <- predict (RF_mod , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.vector(yhat_RF_test)

##measures of accuracy in calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[1,21] <- length(x)
out[1,22] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,23] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,24] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,25] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy in validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,26] <- length(x)
out[1,27] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,28] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,29] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,30] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

############################################################################################
############################## Boosting ####################################################

library(gbm)
##train the model
boost2_mod <- gbm(binary ~.,data = train,n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

# the algorithm recognize the trait as binomial

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict.gbm(boost2_mod,Xnew1,n.trees = 500, type="response")             #predict y hat
yhat1_test_boosting<-predict.gbm(boost2_mod,Xnew1_test,n.trees = 500, type="response")   #predict y hat from test dataset

#define the threshold
roc_objbo <- roc(y.train, yhat1_boosting)
plot(roc_objbo)
coorbo <- coords(roc_objbo, "best", "threshold", transpose = FALSE)
taubo <- coorbo$threshold

# application of the threshold
yhat1_boosting <- ifelse(yhat1_boosting < taubo, 0, 1)
yhat1_test_boosting <- ifelse(yhat1_test_boosting < taubo, 0, 1)

##measures of accuracy in calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,31] <- length(x)
out[1,32] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,33] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,34] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,35] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy in validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,36] <- length(x)
out[1,37] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,38] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,39] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,40] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

###########################################################################################
############################Support Vector Machine SVM#####################################
library(e1071)

classifier = svm(formula = binary ~.,  data = train, 
                 type = 'C-classification', 
                 kernel = 'linear')
y_pred = predict(classifier, newdata = x.train)
y_pred1 = predict(classifier, newdata = x.test)


##measures of accuracy in calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_pred #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,41] <- length(x)
out[1,42] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,43] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,44] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,45] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy in validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_pred1 #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,46] <- length(x)
out[1,47] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,48] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,49] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,50] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


###########################################################################################
############################Logistic Regression############################################
train <- as.data.frame(train)
logfit <- glm(binary~., data=train, family=binomial("logit"), maxit=100)
x.train <- as.data.frame(x.train)
x.test <- as.data.frame(x.test)
yhat_lr <- predict.glm(logfit, newdata = x.train, type = c("response"))
yhat_lr_test <- predict(logfit, newdata = x.test, type = c("response"))


#define the threshold
roc_objlr <- roc(y.train, yhat_lr)
plot(roc_objlr)
coorlr <- coords(roc_objlr, "best", "threshold", transpose = FALSE)
taulr <- coorlr$threshold

#apply the threshold
yhat_lr <- ifelse(yhat_lr < taulr, 0, 1)
yhat_lr_test <- ifelse(yhat_lr_test < taulr, 0, 1)

##measures of accuracy in calibration##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_lr #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,51] <- length(x)
out[1,52] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,53] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,54] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,55] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_lr_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[1,56] <- length(x)
out[1,57] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[1,58] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[1,59] <- ARI.pls <- adjustedRandIndex(x, y)
out[1,60] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


#######################################################################################
#############################dataset2#########################################
# All the steps presented for the first dataset are repeated also for the next datasets

# Division in train and test

training<-data.frame(dataset[-fold2,])
testing<-data.frame(dataset[fold2,])


# Generic code set up
y.train = training$binary
x.train = as.matrix(training[,56:589])
train = cbind(y.train, x.train)
y.test = testing$binary
x.test = as.data.frame(testing[,56:589])
test = cbind(y.test, x.test)
train <- as.data.frame(train)
test <- as.data.frame(test)
names(train)[1] <- "binary"
names(test)[1] <- "binary"
y.train <- as.factor(y.train)
y.test <- as.factor(y.test)
train$binary <- as.factor(train$binary)
test$binary <- as.factor(test$binary)

library(mclust)
library(caret)
library(pROC)

train <- mutate(train, binary = (ifelse(binary == 1, 0, 1)))
#######################################################################################
########################## Apply PLS DA#########################################

plsda.model <- plsda(x.train, y.train, ncomp = 2, probMethod = "bayes")
yhat_plsda <- predict (plsda.model , newdata= x.train)   #predict yhat
yhat_plsda_test <- predict (plsda.model , newdata= x.test)   #predict yhat


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_plsda #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[2,1] <- length(x)
out[2,2] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,3] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,4] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,5] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


pls.cal <- data.frame(ARI.pls, accuracy.pls)

##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_plsda_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,6] <- length(x)
out[2,7] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,8] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,9] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,10] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

pls.test <- data.frame(ARI.pls.test, accuracy.pls.test)



#######################################################################################
########################## Apply Decision Tree#########################################

library(C50)

##train the model
ynew1<-as.factor(y.train)
tree_mod <- C5.0(x = x.train, y = ynew1)     #estimate decision tree model
str(tree_mod)

##predict y
yhat_dt <- predict (tree_mod , newdata= x.train)   #predict yhat
yhat_dt <- as.matrix(yhat_dt)
yhat_dt <- as.numeric(yhat_dt)
yhat_test_dt <- predict (tree_mod , newdata = x.test) #predict yhat from the test dataset
yhat_test_dt <- as.matrix(yhat_test_dt)
yhat_test_dt <- as.numeric(yhat_test_dt)


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_dt #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[2,11] <- length(x)
out[2,12] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,13] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,14] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,15] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_dt #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,16] <- length(x)
out[2,17] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,18] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,19] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,20] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.factor(y.train)
RF_mod <- randomForest (x= x.train, y= y.train, ntree = 500, mtry =23)   #model code

## predict y
yhat_RF <- predict (RF_mod , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.vector(yhat_RF)
yhat_RF_test <- predict (RF_mod , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.vector(yhat_RF_test)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[2,21] <- length(x)
out[2,22] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,23] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,24] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,25] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,26] <- length(x)
out[2,27] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,28] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,29] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,30] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

############################################################################################
############################## Boosting ####################################################

##################

#Algorithm 2
library(gbm)
##train the model
boost2_mod <- gbm(binary ~.,data = train,n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict.gbm(boost2_mod,Xnew1,n.trees = 500, type="response")             #predict y hat
yhat1_test_boosting<-predict.gbm(boost2_mod,Xnew1_test,n.trees = 500, type="response")   #predict y hat from test dataset

#define the threshold
roc_objbo <- roc(y.train, yhat1_boosting)
plot(roc_objbo)
coorbo <- coords(roc_objbo, "best", "threshold", transpose = FALSE)
taubo <- coorbo$threshold


yhat1_boosting <- ifelse(yhat1_boosting < taubo, 0, 1)
yhat1_test_boosting <- ifelse(yhat1_test_boosting < taubo, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[2,31] <- length(x)
out[2,32] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,33] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,34] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,35] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,36] <- length(x)
out[2,37] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,38] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,39] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,40] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

###########################################################################################
############################Support Vector Machine SVM#####################################
library(e1071)

classifier = svm(formula = binary ~.,  data = train, 
                 type = 'C-classification', 
                 kernel = 'linear')
y_pred = predict(classifier, newdata = x.train)
y_pred1 = predict(classifier, newdata = x.test)


##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_pred #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[2,41] <- length(x)
out[2,42] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,43] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,44] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,45] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_pred1 #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,46] <- length(x)
out[2,47] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,48] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,49] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,50] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


###########################################################################################
############################Logistic Regression############################################
train <- as.data.frame(train)
logfit <- glm(binary~., data=train, family=binomial("logit"), maxit=100)
x.train <- as.data.frame(x.train)
x.test <- as.data.frame(x.test)
yhat_lr <- predict.glm(logfit, newdata = x.train, type = c("response"))
yhat_lr_test <- predict(logfit, newdata = x.test, type = c("response"))


#define the threshold
roc_objlr <- roc(y.train, yhat_lr)
plot(roc_objlr)
coorlr <- coords(roc_objlr, "best", "threshold", transpose = FALSE)
taulr <- coorlr$threshold

yhat_lr <- ifelse(yhat_lr < taulr, 0, 1)
yhat_lr_test <- ifelse(yhat_lr_test < taulr, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_lr #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,51] <- length(x)
out[2,52] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,53] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,54] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,55] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_lr_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[2,56] <- length(x)
out[2,57] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[2,58] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[2,59] <- ARI.pls <- adjustedRandIndex(x, y)
out[2,60] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#######################################################################################
#############################dataset3#########################################
# Division in train and test

training<-data.frame(dataset[-fold3,])
testing<-data.frame(dataset[fold3,])


# Generic code set up
y.train = training$binary
x.train = as.matrix(training[,56:589])
train = cbind(y.train, x.train)
y.test = testing$binary
x.test = as.data.frame(testing[,56:589])
test = cbind(y.test, x.test)
train <- as.data.frame(train)
test <- as.data.frame(test)
names(train)[1] <- "binary"
names(test)[1] <- "binary"
y.train <- as.factor(y.train)
y.test <- as.factor(y.test)
train$binary <- as.factor(train$binary)
test$binary <- as.factor(test$binary)

library(mclust)
library(caret)
library(pROC)

train <- mutate(train, binary = (ifelse(binary == 1, 0, 1)))
#######################################################################################
########################## Apply PLS DA#########################################

plsda.model <- plsda(x.train, y.train, ncomp = 2, probMethod = "bayes")
yhat_plsda <- predict (plsda.model , newdata= x.train)   #predict yhat
yhat_plsda_test <- predict (plsda.model , newdata= x.test)   #predict yhat


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_plsda #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[3,1] <- length(x)
out[3,2] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,3] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,4] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,5] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


pls.cal <- data.frame(ARI.pls, accuracy.pls)

##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_plsda_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,6] <- length(x)
out[3,7] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,8] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,9] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,10] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

pls.test <- data.frame(ARI.pls.test, accuracy.pls.test)



#######################################################################################
########################## Apply Decision Tree#########################################

library(C50)

##train the model
ynew1<-as.factor(y.train)
tree_mod <- C5.0(x = x.train, y = ynew1)     #estimate decision tree model
str(tree_mod)

##predict y
yhat_dt <- predict (tree_mod , newdata= x.train)   #predict yhat
yhat_dt <- as.matrix(yhat_dt)
yhat_dt <- as.numeric(yhat_dt)
yhat_test_dt <- predict (tree_mod , newdata = x.test) #predict yhat from the test dataset
yhat_test_dt <- as.matrix(yhat_test_dt)
yhat_test_dt <- as.numeric(yhat_test_dt)


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_dt #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[3,11] <- length(x)
out[3,12] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,13] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,14] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,15] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_dt #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,16] <- length(x)
out[3,17] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,18] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,19] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,20] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.factor(y.train)
RF_mod <- randomForest (x= x.train, y= y.train, ntree = 500, mtry =23)   #model code

## predict y
yhat_RF <- predict (RF_mod , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.vector(yhat_RF)
yhat_RF_test <- predict (RF_mod , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.vector(yhat_RF_test)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[3,21] <- length(x)
out[3,22] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,23] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,24] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,25] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,26] <- length(x)
out[3,27] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,28] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,29] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,30] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

############################################################################################
############################## Boosting ####################################################

##################

#Algorithm 2
library(gbm)
##train the model
boost2_mod <- gbm(binary ~.,data = train,n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict.gbm(boost2_mod,Xnew1,n.trees = 500, type="response")             #predict y hat
yhat1_test_boosting<-predict.gbm(boost2_mod,Xnew1_test,n.trees = 500, type="response")   #predict y hat from test dataset

#define the threshold
roc_objbo <- roc(y.train, yhat1_boosting)
plot(roc_objbo)
coorbo <- coords(roc_objbo, "best", "threshold", transpose = FALSE)
taubo <- coorbo$threshold


yhat1_boosting <- ifelse(yhat1_boosting < taubo, 0, 1)
yhat1_test_boosting <- ifelse(yhat1_test_boosting < taubo, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[3,31] <- length(x)
out[3,32] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,33] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,34] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,35] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,36] <- length(x)
out[3,37] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,38] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,39] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,40] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

###########################################################################################
############################Support Vector Machine SVM#####################################
library(e1071)

classifier = svm(formula = binary ~.,  data = train, 
                 type = 'C-classification', 
                 kernel = 'linear')
y_pred = predict(classifier, newdata = x.train)
y_pred1 = predict(classifier, newdata = x.test)


##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_pred #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[3,41] <- length(x)
out[3,42] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,43] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,44] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,45] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_pred1 #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,46] <- length(x)
out[3,47] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,48] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,49] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,50] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


###########################################################################################
############################Logistic Regression############################################
train <- as.data.frame(train)
logfit <- glm(binary~., data=train, family=binomial("logit"), maxit=100)
x.train <- as.data.frame(x.train)
x.test <- as.data.frame(x.test)
yhat_lr <- predict.glm(logfit, newdata = x.train, type = c("response"))
yhat_lr_test <- predict(logfit, newdata = x.test, type = c("response"))


#define the threshold
roc_objlr <- roc(y.train, yhat_lr)
plot(roc_objlr)
coorlr <- coords(roc_objlr, "best", "threshold", transpose = FALSE)
taulr <- coorlr$threshold

yhat_lr <- ifelse(yhat_lr < taulr, 0, 1)
yhat_lr_test <- ifelse(yhat_lr_test < taulr, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_lr #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,51] <- length(x)
out[3,52] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,53] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,54] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,55] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_lr_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[3,56] <- length(x)
out[3,57] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[3,58] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[3,59] <- ARI.pls <- adjustedRandIndex(x, y)
out[3,60] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


#######################################################################################
#############################dataset4#########################################
# Division in train and test

training<-data.frame(dataset[-fold4,])
testing<-data.frame(dataset[fold4,])


# Generic code set up
y.train = training$binary
x.train = as.matrix(training[,56:589])
train = cbind(y.train, x.train)
y.test = testing$binary
x.test = as.data.frame(testing[,56:589])
test = cbind(y.test, x.test)
train <- as.data.frame(train)
test <- as.data.frame(test)
names(train)[1] <- "binary"
names(test)[1] <- "binary"
y.train <- as.factor(y.train)
y.test <- as.factor(y.test)
train$binary <- as.factor(train$binary)
test$binary <- as.factor(test$binary)

library(mclust)
library(caret)
library(pROC)

train <- mutate(train, binary = (ifelse(binary == 1, 0, 1)))
#######################################################################################
########################## Apply PLS DA#########################################

plsda.model <- plsda(x.train, y.train, ncomp = 2, probMethod = "bayes")
yhat_plsda <- predict (plsda.model , newdata= x.train)   #predict yhat
yhat_plsda_test <- predict (plsda.model , newdata= x.test)   #predict yhat


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_plsda #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[4,1] <- length(x)
out[4,2] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,3] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,4] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,5] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


pls.cal <- data.frame(ARI.pls, accuracy.pls)

##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_plsda_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,6] <- length(x)
out[4,7] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,8] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,9] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,10] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

pls.test <- data.frame(ARI.pls.test, accuracy.pls.test)



#######################################################################################
########################## Apply Decision Tree#########################################

library(C50)

##train the model
ynew1<-as.factor(y.train)
tree_mod <- C5.0(x = x.train, y = ynew1)     #estimate decision tree model
str(tree_mod)

##predict y
yhat_dt <- predict (tree_mod , newdata= x.train)   #predict yhat
yhat_dt <- as.matrix(yhat_dt)
yhat_dt <- as.numeric(yhat_dt)
yhat_test_dt <- predict (tree_mod , newdata = x.test) #predict yhat from the test dataset
yhat_test_dt <- as.matrix(yhat_test_dt)
yhat_test_dt <- as.numeric(yhat_test_dt)


##measures of accuracy cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_dt #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[4,11] <- length(x)
out[4,12] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,13] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,14] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,15] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_dt #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,16] <- length(x)
out[4,17] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,18] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,19] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,20] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.factor(y.train)
RF_mod <- randomForest (x= x.train, y= y.train, ntree = 500, mtry =23)   #model code

## predict y
yhat_RF <- predict (RF_mod , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.vector(yhat_RF)
yhat_RF_test <- predict (RF_mod , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.vector(yhat_RF_test)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[4,21] <- length(x)
out[4,22] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,23] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,24] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,25] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,26] <- length(x)
out[4,27] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,28] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,29] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,30] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

############################################################################################
############################## Boosting ####################################################

##################

#Algorithm 2
library(gbm)
##train the model
boost2_mod <- gbm(binary ~.,data = train,n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict.gbm(boost2_mod,Xnew1,n.trees = 500, type="response")             #predict y hat
yhat1_test_boosting<-predict.gbm(boost2_mod,Xnew1_test,n.trees = 500, type="response")   #predict y hat from test dataset

#define the threshold
roc_objbo <- roc(y.train, yhat1_boosting)
plot(roc_objbo)
coorbo <- coords(roc_objbo, "best", "threshold", transpose = FALSE)
taubo <- coorbo$threshold


yhat1_boosting <- ifelse(yhat1_boosting < taubo, 0, 1)
yhat1_test_boosting <- ifelse(yhat1_test_boosting < taubo, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##


table.pls <- as.data.frame(table(x,y))

out[4,31] <- length(x)
out[4,32] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,33] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,34] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,35] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,36] <- length(x)
out[4,37] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,38] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,39] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,40] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])

###########################################################################################
############################Support Vector Machine SVM#####################################
library(e1071)

classifier = svm(formula = binary ~.,  data = train, 
                 type = 'C-classification', 
                 kernel = 'linear')
y_pred = predict(classifier, newdata = x.train)
y_pred1 = predict(classifier, newdata = x.test)


##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_pred #values of the predicted eg fat_predicted#

##

table.pls <- as.data.frame(table(x,y))

out[4,41] <- length(x)
out[4,42] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,43] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,44] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,45] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_pred1 #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,46] <- length(x)
out[4,47] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,48] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,49] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,50] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


###########################################################################################
############################Logistic Regression############################################
train <- as.data.frame(train)
logfit <- glm(binary~., data=train, family=binomial("logit"), maxit=100)
x.train <- as.data.frame(x.train)
x.test <- as.data.frame(x.test)
yhat_lr <- predict.glm(logfit, newdata = x.train, type = c("response"))
yhat_lr_test <- predict(logfit, newdata = x.test, type = c("response"))


#define the threshold
roc_objlr <- roc(y.train, yhat_lr)
plot(roc_objlr)
coorlr <- coords(roc_objlr, "best", "threshold", transpose = FALSE)
taulr <- coorlr$threshold

yhat_lr <- ifelse(yhat_lr < taulr, 0, 1)
yhat_lr_test <- ifelse(yhat_lr_test < taulr, 0, 1)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_lr #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,51] <- length(x)
out[4,52] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,53] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,54] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,55] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])


##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_lr_test #values of the predicted eg fat_predicted#

##
table.pls <- as.data.frame(table(x,y))

out[4,56] <- length(x)
out[4,57] <- (table.pls[1, 3]/(table.pls[1, 3]+table.pls[2, 3]))
out[4,58] <- (table.pls[4, 3]/(table.pls[4, 3]+table.pls[3, 3]))
out[4,59] <- ARI.pls <- adjustedRandIndex(x, y)
out[4,60] <-accuracy.pls <- (table.pls[1, 3]+table.pls[4, 3])/(table.pls[1, 3]+table.pls[2, 3]+table.pls[3, 3]+table.pls[4, 3])



#########################################################################################
######################results#############################################################
## Look at the results

apply(out, 2, summary)  #mean, min, max
boxplot(out)   #boxplot of the results
apply(out, 2, sd)   #standard deviation of the results

#save the output

out <- as.data.frame(out)

library("writexl")
write_xlsx(out,"C:/Users/direction where to save the file/file name.xlsx")
