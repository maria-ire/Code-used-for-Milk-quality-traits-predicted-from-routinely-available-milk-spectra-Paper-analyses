# Set the direction where to find the dataset to be imported

rm(list=ls())
setwd("C:/Users/direction to the file")
list.files()

# Read in data

dataset = read.csv("file name.csv")

#load the useful packages

library("pls")
library(glmnet)  # for ridge regression
library(dplyr)   # for data cleaning
library(psych)   # for function tr() to compute trace of a matrix
library(penalized)

##define the trait of interest (name present in the dataset) to analyse and remove the Na

trait = "name of the trait"
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

# Generic code set up
y.train = dataset$beta_casein             #vector of the reference data
x.train = as.matrix(dataset[,56:589])     #matrix of the spectra
train = cbind(x.train, y.train)
train <- as.data.frame(train)
names(train)[1] <- "hs"

## definition of a table where results are going to be reported

R <- 4                               #number of n-fold cross validation
out <- matrix(NA, R, 77)             #number of output columns
colnames(out) <- c("plsr.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "ridge.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "lasso.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "en.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "average.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "pcr.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "ppr.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "ssr.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "rf.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "boosting.n", "RMSE", "r", "n", "RMSECV", "Slope", "r",
                   "brnn.n", "RMSE", "r", "n", "RMSECV", "Slope", "r")


##Following the algorithms are repeated 4 times, and every time a different 
# fold as test dataset was considered

##############################################################################
############################################################################
###############define x.train y.train x.test y.test based on the fold############
x.train = as.matrix(train[-fold1,-1])    #matrix of spectra of the train dataset
y.train = train[-fold1,1]                #vector of the trait of interest train dataset
x.test = as.matrix(train[fold1,-1])      #matrix od spectra of the test dataset
y.test= train [fold1,1]                  #vector of the trait of interest test dataset

############################################################################################
###################################### Apply PLS############################################
library("pls")
library(magrittr)
library("parallel")

# Generic code set up: put training and test x and y data into data frames
train.pls = data.frame(y = y.train)
train.pls$x = x.train
test.pls = data.frame(y = y.test)
test.pls$x = x.test

# Fit a PLS regression to the training data
# Use multiple cores for the CV if available
pls.options(parallel = makeCluster(detectCores(), type = "PSOCK"))
res.pls = plsr(y ~ x, ncomp=20, data = train.pls, validation = "LOO", scale = TRUE)
stopCluster(pls.options()$parallel)


# Examine results through plots
summary(res.pls)   
plot(RMSEP(res.pls), legendpos = "topright",main="Fat Variable PLS")
axis(side = 1, at=1:10)
plot(res.pls, ncomp = 4, asp = 1, line = TRUE, main="4  Component Fit")
plot(res.pls, plottype = "scores", comps = 1:4, main="Scores")
plot(res.pls, "loadings", comps = 1:4, legendpos = "bottomright", xlab = "spectrum", main="Fat Dataset 1: PLS Components")
abline(h = 0)

# Predict y hat in training
y_hat_pls = predict(res.pls, ncomp = 10, newdata = train.pls$x)

# Predict y hat in testing
y_hat_pls_test = predict(res.pls, ncomp = 10, newdata = test.pls$x)

##measures of accuracy
##define X and Y

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_pls #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,1] <- n
n1
bias
out[1,2] <- MPE
RPE
out[1,3] <- R
B1

##measures of accuracy in validation
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_pls_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,4] <- n
n1
bias
out[1,5] <- MPE
RPE
out[1,7] <- R
out[1,6] <- B1
############################################################################################
###############################Apply Ridge Regression ######################################
library(glmnet)
#set lambda
set.seed(123)
lambdas_ridge <- 10^seq(-2, 2, length.out = 1000) 
#TRAIN model
ridge_cv <- cv.glmnet(x.train, y.train, alpha = 0, lambda = lambdas_ridge,
                      standardize = TRUE, nfolds = 10)
lambda_ridge_cv <- ridge_cv$lambda.min
model_ridge_cv <- glmnet(x.train, y.train, alpha = 0, lambda = lambda_ridge_cv, standardize = TRUE)

#Y hat in training 
y_hat_ridge_cv <- predict(model_ridge_cv, x.train)

# number of wavelength selected
beta_ridge<-model_ridge_cv$beta
Bot_ridge <- beta_ridge[which(model_ridge_cv$beta !=0),]
Bot_ridge<-as.matrix(Bot_ridge)
dim(Bot_ridge)
nrow(Bot_ridge)

#testing
y_hat_ridge_cv_test_A<-predict(model_ridge_cv, x.test)


##measures of accuracy in calibration
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,8] <- n
n1
bias
out[1,9] <- MPE
RPE
out[1,10] <- R
B1


##measures of accuracy in validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv_test_A #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,11] <- n
n1
bias
out[1,12] <- MPE
RPE
out[1,14] <- R
out[1,13] <- B1

##################################################################################
###################################LASSO########################################## 
#cross-validation to select lambda
lambdas_lasso <- 10^seq(-3, 3, length.out = 1000)
#TRAIN the lasso model
lasso_cv <- cv.glmnet(x.train, y.train, alpha = 1, lambda = lambdas_lasso,
                      standardize = TRUE, nfolds = 10)
lambda_lasso_cv <- lasso_cv$lambda.min
model_lasso_cv <- glmnet(x.train, y.train, alpha = 1, lambda = lambda_lasso_cv, standardize = TRUE)
# ppredict y hat in calibration
y_hat_lasso_cv <- predict(model_lasso_cv, x.train)
# number of wavelengths selected
beta_lasso<-model_lasso_cv$beta
Bot_lasso <- beta_lasso[which(model_lasso_cv$beta !=0),] 
Bot_lasso1<-as.matrix(Bot_lasso)
dim(Bot_lasso)
nrow(Bot_lasso)

#test 
y_hat_lasso_cv_test <- predict(model_lasso_cv, x.test)


##measures of accuracy in calibration##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,15] <- n
n1
bias
out[1,16] <- MPE
RPE
out[1,17] <- R
B1


##measures of accuracy in validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,18] <- n
n1
bias
out[1,19] <- MPE
RPE
out[1,21] <- R
out[1,20] <- B1


####################################################################################
################################Apply Elastic Net###################################

library(glmnet)
# train the en model
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,intercept=FALSE)
fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")
fit.elnetLambda<-fit.elnet.cv$lambda.min
fit.elnet2 <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda =fit.elnetLambda, 
                     intercept=FALSE)
beta<-fit.elnet2$beta
# number of wavelengths selected
Bot_en <- beta[which(fit.elnet2$beta !=0),] 
Bot_en1<-as.matrix(Bot_en)
dim(Bot_en)
nrow(Bot_en)
# y hat in calibration
y_hat_en<-x.train%*%beta

# y hat in  validation
y_hat_en_test<-x.test%*%beta


##measures of accuracy in cross validation
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_en #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,22] <- n
n1
bias
out[1,23] <- MPE
RPE
out[1,24] <- R
B1


##measures of accuracy in validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_en_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,25] <- n
n1
bias
out[1,26] <- MPE
RPE
out[1,28] <- R
out[1,27] <- B1

##################################################################################
############################### averaging model###################################

library(broom)
library(tidyverse)

#evaluate y hat as average of regression models in calibration
yhat_cv <- cbind(y_hat_lasso_cv, y_hat_en, y_hat_pls, y_hat_ridge_cv)
yhat_cv <- as.matrix(yhat_cv)
yhat_cv <- as.data.frame(yhat_cv)
names(yhat_cv)[1] <- "lasso"
names(yhat_cv)[2] <- "en"
names(yhat_cv)[3] <- "pls"
names(yhat_cv)[4] <- "ridge"
yhat_cv <- transform(yhat_cv, lasso = as.numeric(lasso))
yhat_cv <- transform(yhat_cv, en = as.numeric(en))
yhat_cv <- transform(yhat_cv, pls = as.numeric(pls))
yhat_cv <- transform(yhat_cv, ridge = as.numeric(ridge))
yhat_cv1 <- apply(yhat_cv,1,mean)

#evaluate y hat as average of regression models in validation
yhat_cv_test <- cbind(y_hat_lasso_cv_test, y_hat_en_test, y_hat_pls_test, y_hat_ridge_cv_test_A)
yhat_cv_test <- as.matrix(yhat_cv_test)
yhat_cv_test <- as.data.frame(yhat_cv_test)
names(yhat_cv_test)[1] <- "lasso"
names(yhat_cv_test)[2] <- "en"
names(yhat_cv_test)[3] <- "pls"
names(yhat_cv_test)[4] <- "ridge"
yhat_cv_test <- transform(yhat_cv_test, lasso = as.numeric(lasso))
yhat_cv_test <- transform(yhat_cv_test, en = as.numeric(en))
yhat_cv_test <- transform(yhat_cv_test, pls = as.numeric(pls))
yhat_cv_test <- transform(yhat_cv_test, ridge = as.numeric(ridge))
yhat_cv_test1 <- apply(yhat_cv_test,1,mean)


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_cv1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,29] <- n
n1
bias
out[1,30] <- MPE
RPE
out[1,31] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_cv_test1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,32] <- n
n1
bias
out[1,33] <- MPE
RPE
out[1,35] <- R
out[1,34] <- B1

#########################################################################################
############## Principal Component Analysis #############################################
library (pls)

##train the model
set.seed(123)
pcr_mod <- pcr(hs ~., data=train, ncomp=20, scale=TRUE)
summary(pcr_mod)
##predict y hat in calibration
yhat_pcr <- pcr_mod %>% predict(x.train, ncomp=20)               #predict y hat
##predict y hat in validation
yhat_test_pcr <- pcr_mod %>% predict(x.test, ncomp=20)     # predict y hat from test dataset


##measures of accuracy cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,36] <- n
n1
bias
out[1,37] <- MPE
RPE
out[1,38] <- R
B1

##measures of accuracy validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,39] <- n
n1
bias
out[1,40] <- MPE
RPE
out[1,42] <- R
out[1,41] <- B1

#########################################################################################
############## Projection Pursuit Regression ############################################
library (stats)

##trai the model
set.seed(123)
ppr_mod <- ppr(x.train, y.train, nterms=1, max.terms=2)
ppr_plot <- plot(ppr_mod)
## predict y calibration
yhat_ppr <- ppr_mod %>% predict(x.train)               #predict y hat
## predict y validation
yhat_test_ppr <- ppr_mod %>% predict(x.test)     # predict y hat from test dataset



##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,43] <- n
n1
bias
out[1,44] <- MPE
RPE
out[1,45] <- R
B1

##measures of accuracy in validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,46] <- n
n1
bias
out[1,47] <- MPE
RPE
out[1,49] <- R
out[1,48] <- B1

#########################################################################################
############## Spike and Slab Regression ################################################
library (spikeslab)
library (plyr)

###train the model
spikeslab_mod <- spikeslab (hs ~., data = train, max.var = 1000)
## predict y in calibration
yhat_spikeslab <- spikeslab_mod %>% predict(x.train)               #predict y hat
yhat_spikeslab <- yhat_spikeslab$yhat.bma
## predict y in validation
yhat_test_spikeslab <- spikeslab_mod %>% predict(x.test)     # predict y hat from test dataset
yhat_test_spikeslab <- yhat_test_spikeslab$yhat.bma

# wavelengths selected by spike and slab
attributes(spikeslab_mod)
spikeslab_mod$gnet
length(spikeslab_mod$gnet)
plot(spikeslab_mod$gnet, type="h")
sum(spikeslab_mod$gnet == 0)
ssr_bot1 <- as.data.frame(which(spikeslab_mod$gnet != 0))
dim(ssr_bot)


##define X and Y in cross validation##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,50] <- n
n1
bias
out[1,51] <- MPE
RPE
out[1,52] <- R
B1

##define X and Y in external validation##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,53] <- n
n1
bias
out[1,54] <- MPE
RPE
out[1,56] <- R
out[1,55] <- B1

#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.vector(train[-fold1,1])
RF_mod1 <- randomForest (x= train[-fold1,-1], y= train[-fold1,1], ntree = 500, mtry = 40)   #model code

## predict y
yhat_RF <- predict (RF_mod1 , newdata= train[-fold1,-1])           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.numeric(yhat_RF)
yhat_RF_test <- predict (RF_mod1 , newdata = train[fold1,-1]) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.numeric(yhat_RF_test)

##measures of accuracy in cross validation##

##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,57] <- as.numeric(n)
n1
bias
out[1,58] <-MPE
RPE
out[1,59] <-R
B1

rf.cv <- data.frame(n, n1, bias, MPE, RPE, R, B1)

##measures of accuracy in external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,60] <-n
n1
bias
out[1,61] <-MPE
RPE
out[1,63] <-R
out[1,62] <-B1

############################################################################################
############################## Boosting ####################################################

library(gbm)
##train the model
boost2_mod <- gbm(hs ~.,data = train,distribution = "gaussian",n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict(boost2_mod,Xnew1,n.trees = 500)             #predict y hat
yhat1_test_boosting<-predict(boost2_mod,Xnew1_test,n.trees = 500)   #predict y hat from test dataset


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,64] <- n
n1
bias
out[1,65] <- MPE
RPE
out[1,66] <- R
B1


##measures of accuracy in external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,67] <- n
n1
bias
out[1,68] <- MPE
RPE
out[1,70] <- R
out[1,69] <- B1


library(caret)

#########################################################################################
##############Neural Network #####################################
library(brnn)
##train model
set.seed(123)
brnn_mod <- brnn(x.train, y.train)

###prediction
yhat_brnn <- brnn_mod %>% predict(x.train)               #predict y hat
yhat_test_brnn <- brnn_mod %>% predict(x.test)     # predict y hat from test dataset


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,71] <- n
n1
bias
out[1,72] <- MPE
RPE
out[1,73] <- R
B1

##measures of accuracy in cross validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[1,74] <- n
n1
bias
out[1,75] <- MPE
RPE
out[1,77] <- R
out[1,76] <- B1

##############################################################################
############################################################################
###############define x.train y.train x.test y.test#######################
x.train = as.matrix(train[-fold2,-1])
y.train = train[-fold2,1]
x.test = as.matrix(train[fold2,-1])
y.test= train [fold2,1]

############################################################################################
###################################### Apply PLS############################################
library("pls")
library(magrittr)
library("parallel")

# Generic code set up: put training and test x and y data into data frames
train.pls = data.frame(y = y.train)
train.pls$x = x.train
test.pls = data.frame(y = y.test)
test.pls$x = x.test

# Fit a PLS regression to the training data
# Use multiple cores for the CV if available
pls.options(parallel = makeCluster(detectCores(), type = "PSOCK"))
res.pls = plsr(y ~ x, ncomp=20, data = train.pls, validation = "LOO", scale = TRUE)
stopCluster(pls.options()$parallel)


# Examine results through plots
summary(res.pls)   
plot(RMSEP(res.pls), legendpos = "topright",main="Fat Variable PLS")
axis(side = 1, at=1:10)
plot(res.pls, ncomp = 4, asp = 1, line = TRUE, main="4  Component Fit")
plot(res.pls, plottype = "scores", comps = 1:4, main="Scores")
plot(res.pls, "loadings", comps = 1:4, legendpos = "bottomright", xlab = "spectrum", main="Fat Dataset 1: PLS Components")
abline(h = 0)

# Predict y hat in cross validation
y_hat_pls = predict(res.pls, ncomp = 10, newdata = train.pls$x)

# Predict y hat in external validation
y_hat_pls_test = predict(res.pls, ncomp = 10, newdata = test.pls$x)

##measures of accuracy in cross validation
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_pls #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,1] <- n
n1
bias
out[2,2] <- MPE
RPE
out[2,3] <- R
B1
fit$coefficients


##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_pls_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,4] <- n
n1
bias
out[2,5] <- MPE
RPE
out[2,7] <- R
out[2,6] <- B1
fit$coefficients

############################################################################################
###############################Apply Ridge Regression ######################################
library(glmnet)
#set lambda
set.seed(123)
lambdas_ridge <- 10^seq(-2, 2, length.out = 1000) 
#TRAIN model
ridge_cv <- cv.glmnet(x.train, y.train, alpha = 0, lambda = lambdas_ridge,
                      standardize = TRUE, nfolds = 10)
lambda_ridge_cv <- ridge_cv$lambda.min
model_ridge_cv <- glmnet(x.train, y.train, alpha = 0, lambda = lambda_ridge_cv, standardize = TRUE)

#Y hat in cross validation
y_hat_ridge_cv <- predict(model_ridge_cv, x.train)

# number of wavelength selected
beta_ridge<-model_ridge_cv$beta
Bot_ridge <- beta_ridge[which(model_ridge_cv$beta !=0),]
Bot_ridge<-as.matrix(Bot_ridge)
dim(Bot_ridge)
nrow(Bot_ridge)

#test
y_hat_ridge_cv_test_A<-predict(model_ridge_cv, x.test)


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,8] <- n
n1
bias
out[2,9] <- MPE
RPE
out[2,10] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv_test_A #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,11] <- n
n1
bias
out[2,12] <- MPE
RPE
out[2,14] <- R
out[2,13] <- B1

##################################################################################
###################################LASSO########################################## 
#cross-validation to select lambda
lambdas_lasso <- 10^seq(-3, 3, length.out = 1000)
#TRAIN the lasso model
lasso_cv <- cv.glmnet(x.train, y.train, alpha = 1, lambda = lambdas_lasso,
                      standardize = TRUE, nfolds = 10)
lambda_lasso_cv <- lasso_cv$lambda.min
model_lasso_cv <- glmnet(x.train, y.train, alpha = 1, lambda = lambda_lasso_cv, standardize = TRUE)
# ppredict y hat in cross validation
y_hat_lasso_cv <- predict(model_lasso_cv, x.train)
# number of wavelength selected
beta_lasso<-model_lasso_cv$beta
Bot_lasso <- beta_lasso[which(model_lasso_cv$beta !=0),] 
Bot_lasso2<-as.matrix(Bot_lasso)
dim(Bot_lasso)
nrow(Bot_lasso)

#test 
y_hat_lasso_cv_test <- predict(model_lasso_cv, x.test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,15] <- n
n1
bias
out[2,16] <- MPE
RPE
out[2,17] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,18] <- n
n1
bias
out[2,19] <- MPE
RPE
out[2,21] <- R
out[2,20] <- B1

####################################################################################
################################Apply Elastic Net###################################

library(glmnet)
# train the en model
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,intercept=FALSE)
fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")
fit.elnetLambda<-fit.elnet.cv$lambda.min
fit.elnet2 <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda =fit.elnetLambda, 
                     intercept=FALSE)
beta<-fit.elnet2$beta
# number of wavelength selected
Bot_en <- beta[which(fit.elnet2$beta !=0),] 
Bot_en2<-as.matrix(Bot_en)
dim(Bot_en)
nrow(Bot_en)
# y hat in cross validation
y_hat_en<-x.train%*%beta

# y hat in external validation
y_hat_en_test<-x.test%*%beta


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_en #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,22] <- n
n1
bias
out[2,23] <- MPE
RPE
out[2,24] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_en_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,25] <- n
n1
bias
out[2,26] <- MPE
RPE
out[2,28] <- R
out[2,27] <- B1

##################################################################################
############################### averaging model###################################

library(broom)
library(tidyverse)

#evaluate y hat as average of regression models in cross validation
yhat_cv <- cbind(y_hat_lasso_cv, y_hat_en, y_hat_pls, y_hat_ridge_cv)
yhat_cv <- as.matrix(yhat_cv)
yhat_cv <- as.data.frame(yhat_cv)
names(yhat_cv)[1] <- "lasso"
names(yhat_cv)[2] <- "en"
names(yhat_cv)[3] <- "pls"
names(yhat_cv)[4] <- "ridge"
yhat_cv <- transform(yhat_cv, lasso = as.numeric(lasso))
yhat_cv <- transform(yhat_cv, en = as.numeric(en))
yhat_cv <- transform(yhat_cv, pls = as.numeric(pls))
yhat_cv <- transform(yhat_cv, ridge = as.numeric(ridge))
yhat_cv1 <- apply(yhat_cv,1,mean)

#evaluate y hat as average of regression models in external validation
yhat_cv_test <- cbind(y_hat_lasso_cv_test, y_hat_en_test, y_hat_pls_test, y_hat_ridge_cv_test_A)
yhat_cv_test <- as.matrix(yhat_cv_test)
yhat_cv_test <- as.data.frame(yhat_cv_test)
names(yhat_cv_test)[1] <- "lasso"
names(yhat_cv_test)[2] <- "en"
names(yhat_cv_test)[3] <- "pls"
names(yhat_cv_test)[4] <- "ridge"
yhat_cv_test <- transform(yhat_cv_test, lasso = as.numeric(lasso))
yhat_cv_test <- transform(yhat_cv_test, en = as.numeric(en))
yhat_cv_test <- transform(yhat_cv_test, pls = as.numeric(pls))
yhat_cv_test <- transform(yhat_cv_test, ridge = as.numeric(ridge))
yhat_cv_test1 <- apply(yhat_cv_test,1,mean)


##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_cv1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,29] <- n
n1
bias
out[2,30] <- MPE
RPE
out[2,31] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_cv_test1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,32] <- n
n1
bias
out[2,33] <- MPE
RPE
out[2,35] <- R
out[2,34] <- B1

#########################################################################################
############## Principal Component Analysis #############################################
library (pls)

##train the model
set.seed(123)
pcr_mod <- pcr(hs ~., data=train, ncomp=20, scale=TRUE)
summary(pcr_mod)
##predict y hat in cross validation
yhat_pcr <- pcr_mod %>% predict(x.train, ncomp=20)               #predict y hat
##predict y hat in validation
yhat_test_pcr <- pcr_mod %>% predict(x.test, ncomp=20)     # predict y hat from test dataset

##measures of accuracy cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,36] <- n
n1
bias
out[2,37] <- MPE
RPE
out[2,38] <- R
B1

##measures of accuracy external validation##

##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,39] <- n
n1
bias
out[2,40] <- MPE
RPE
out[2,42] <- R
out[2,41] <- B1

#########################################################################################
############## Projection Pursuit Regression ############################################
library (stats)

##trai the model
set.seed(123)
ppr_mod <- ppr(x.train, y.train, nterms=1, max.terms=2)
ppr_plot <- plot(ppr_mod)
## predict y cross validation
yhat_ppr <- ppr_mod %>% predict(x.train)               #predict y hat
## predict y external validation
yhat_test_ppr <- ppr_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,43] <- n
n1
bias
out[2,44] <- MPE
RPE
out[2,45] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,46] <- n
n1
bias
out[2,47] <- MPE
RPE
out[2,49] <- R
out[2,48] <- B1

#########################################################################################
############## Spike and Slab Regression ################################################
library (spikeslab)
library (plyr)

###train the model
spikeslab_mod <- spikeslab (hs ~., data = train, max.var = 1000)
## predict y in corss validation
yhat_spikeslab <- spikeslab_mod %>% predict(x.train)               #predict y hat
yhat_spikeslab <- yhat_spikeslab$yhat.bma
## predict y in external validation
yhat_test_spikeslab <- spikeslab_mod %>% predict(x.test)     # predict y hat from test dataset
yhat_test_spikeslab <- yhat_test_spikeslab$yhat.bma

# wavelength selected by spike and slab
attributes(spikeslab_mod)
spikeslab_mod$gnet
length(spikeslab_mod$gnet)
plot(spikeslab_mod$gnet, type="h")
sum(spikeslab_mod$gnet == 0)
ssr_bot2 <- as.data.frame(which(spikeslab_mod$gnet != 0))
dim(ssr_bot)

##define X and Y in cross validation##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,50] <- n
n1
bias
out[2,51] <- MPE
RPE
out[2,52] <- R
B1

##define X and Y in external validation##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,53] <- n
n1
bias
out[2,54] <- MPE
RPE
out[2,56] <- R
out[2,55] <- B1

#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.vector(y.train)
RF_mod2 <- randomForest (x= x.train, y= y.train, ntree = 500, mtry = 40)   #model code

## predict y
yhat_RF <- predict (RF_mod2 , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.numeric(yhat_RF)
yhat_RF_test <- predict (RF_mod2 , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.numeric(yhat_RF_test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,57] <- as.numeric(n)
n1
bias
out[2,58] <-MPE
RPE
out[2,59] <-R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,60] <-n
n1
bias
out[2,61] <-MPE
RPE
out[2,63] <-R
out[2,62] <-B1

############################################################################################
############################## Boosting ####################################################

library(gbm)
##train the model
boost2_mod <- gbm(hs ~.,data = train,distribution = "gaussian",n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict(boost2_mod,Xnew1,n.trees = 500)             #predict y hat
yhat1_test_boosting<-predict(boost2_mod,Xnew1_test,n.trees = 500)   #predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,64] <- n
n1
bias
out[2,65] <- MPE
RPE
out[2,66] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,67] <- n
n1
bias
out[2,68] <- MPE
RPE
out[2,70] <- R
out[2,69] <- B1

#########################################################################################
############## Neural Network #####################################
library(brnn)
##train model
set.seed(123)
brnn_mod <- brnn(x.train, y.train)

###prediction
yhat_brnn <- brnn_mod %>% predict(x.train)               #predict y hat
yhat_test_brnn <- brnn_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,71] <- n
n1
bias
out[2,72] <- MPE
RPE
out[2,73] <- R
B1

##measures of accuracy in cross validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[2,74] <- n
n1
bias
out[2,75] <- MPE
RPE
out[2,77] <- R
out[2,76] <- B1

##############################################################################
############################################################################
###############define x.train y.train x.test y.test#######################
x.train = as.matrix(train[-fold3,-1])
y.train = train[-fold3,1]
x.test = as.matrix(train[fold3,-1])
y.test= train [fold3,1]

############################################################################################
###################################### Apply PLS############################################
library("pls")
library(magrittr)
library("parallel")

# Generic code set up: put training and test x and y data into data frames
train.pls = data.frame(y = y.train)
train.pls$x = x.train
test.pls = data.frame(y = y.test)
test.pls$x = x.test

# Fit a PLS regression to the training data
# Use multiple cores for the CV if available
pls.options(parallel = makeCluster(detectCores(), type = "PSOCK"))
res.pls = plsr(y ~ x, ncomp=20, data = train.pls, validation = "LOO", scale = TRUE)
stopCluster(pls.options()$parallel)

# Examine results through plots
summary(res.pls)   
plot(RMSEP(res.pls), legendpos = "topright",main="Fat Variable PLS")
axis(side = 1, at=1:10)
plot(res.pls, ncomp = 4, asp = 1, line = TRUE, main="4  Component Fit")
plot(res.pls, plottype = "scores", comps = 1:4, main="Scores")
plot(res.pls, "loadings", comps = 1:4, legendpos = "bottomright", xlab = "spectrum", main="Fat Dataset 1: PLS Components")
abline(h = 0)

# Predict y hat in cross validation
y_hat_pls = predict(res.pls, ncomp = 10, newdata = train.pls$x)

# Predict y hat in external validation
y_hat_pls_test = predict(res.pls, ncomp = 10, newdata = test.pls$x)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_pls #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,1] <- n
n1
bias
out[3,2] <- MPE
RPE
out[3,3] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_pls_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,4] <- n
n1
bias
out[3,5] <- MPE
RPE
out[3,7] <- R
out[3,6] <- B1


############################################################################################
###############################Apply Ridge Regression ######################################
library(glmnet)
#set lambda
set.seed(123)
lambdas_ridge <- 10^seq(-2, 2, length.out = 1000) 
#TRAIN model
ridge_cv <- cv.glmnet(x.train, y.train, alpha = 0, lambda = lambdas_ridge,
                      standardize = TRUE, nfolds = 10)
lambda_ridge_cv <- ridge_cv$lambda.min
model_ridge_cv <- glmnet(x.train, y.train, alpha = 0, lambda = lambda_ridge_cv, standardize = TRUE)

#Y hat in cross validation
y_hat_ridge_cv <- predict(model_ridge_cv, x.train)

# number of wavelength selected
beta_ridge<-model_ridge_cv$beta
Bot_ridge <- beta_ridge[which(model_ridge_cv$beta !=0),]
Bot_ridge<-as.matrix(Bot_ridge)
dim(Bot_ridge)
nrow(Bot_ridge)

#test
y_hat_ridge_cv_test_A<-predict(model_ridge_cv, x.test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,8] <- n
n1
bias
out[3,9] <- MPE
RPE
out[3,10] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv_test_A #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,11] <- n
n1
bias
out[3,12] <- MPE
RPE
out[3,14] <- R
out[3,13] <- B1

##################################################################################
###################################LASSO########################################## 
#cross-validation to select lambda
lambdas_lasso <- 10^seq(-3, 3, length.out = 1000)
#TRAIN the lasso model
lasso_cv <- cv.glmnet(x.train, y.train, alpha = 1, lambda = lambdas_lasso,
                      standardize = TRUE, nfolds = 10)
lambda_lasso_cv <- lasso_cv$lambda.min
model_lasso_cv <- glmnet(x.train, y.train, alpha = 1, lambda = lambda_lasso_cv, standardize = TRUE)
# ppredict y hat in cross validation
y_hat_lasso_cv <- predict(model_lasso_cv, x.train)
# number of wavelength selected
beta_lasso<-model_lasso_cv$beta
Bot_lasso <- beta_lasso[which(model_lasso_cv$beta !=0),] 
Bot_lasso3<-as.matrix(Bot_lasso)
dim(Bot_lasso)
nrow(Bot_lasso)

#test 
y_hat_lasso_cv_test <- predict(model_lasso_cv, x.test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,15] <- n
n1
bias
out[3,16] <- MPE
RPE
out[3,17] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,18] <- n
n1
bias
out[3,19] <- MPE
RPE
out[3,21] <- R
out[3,20] <- B1

####################################################################################
################################Apply Elastic Net###################################

library(glmnet)
# train the en model
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,intercept=FALSE)
fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")
fit.elnetLambda<-fit.elnet.cv$lambda.min
fit.elnet2 <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda =fit.elnetLambda, 
                     intercept=FALSE)
beta<-fit.elnet2$beta
# number of wavelength selected
Bot_en <- beta[which(fit.elnet2$beta !=0),] 
Bot_en3<-as.matrix(Bot_en)
dim(Bot_en)
nrow(Bot_en)
# y hat in cross validation
y_hat_en<-x.train%*%beta

# y hat in external validation
y_hat_en_test<-x.test%*%beta

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_en #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,22] <- n
n1
bias
out[3,23] <- MPE
RPE
out[3,24] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_en_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,25] <- n
n1
bias
out[3,26] <- MPE
RPE
out[3,28] <- R
out[3,27] <- B1

##################################################################################
############################### averaging model###################################

library(broom)
library(tidyverse)

#evaluate y hat as average of regression models in cross validation
yhat_cv <- cbind(y_hat_lasso_cv, y_hat_en, y_hat_pls, y_hat_ridge_cv)
yhat_cv <- as.matrix(yhat_cv)
yhat_cv <- as.data.frame(yhat_cv)
names(yhat_cv)[1] <- "lasso"
names(yhat_cv)[2] <- "en"
names(yhat_cv)[3] <- "pls"
names(yhat_cv)[4] <- "ridge"
yhat_cv <- transform(yhat_cv, lasso = as.numeric(lasso))
yhat_cv <- transform(yhat_cv, en = as.numeric(en))
yhat_cv <- transform(yhat_cv, pls = as.numeric(pls))
yhat_cv <- transform(yhat_cv, ridge = as.numeric(ridge))
yhat_cv1 <- apply(yhat_cv,1,mean)

#evaluate y hat as average of regression models in external validation
yhat_cv_test <- cbind(y_hat_lasso_cv_test, y_hat_en_test, y_hat_pls_test, y_hat_ridge_cv_test_A)
yhat_cv_test <- as.matrix(yhat_cv_test)
yhat_cv_test <- as.data.frame(yhat_cv_test)
names(yhat_cv_test)[1] <- "lasso"
names(yhat_cv_test)[2] <- "en"
names(yhat_cv_test)[3] <- "pls"
names(yhat_cv_test)[4] <- "ridge"
yhat_cv_test <- transform(yhat_cv_test, lasso = as.numeric(lasso))
yhat_cv_test <- transform(yhat_cv_test, en = as.numeric(en))
yhat_cv_test <- transform(yhat_cv_test, pls = as.numeric(pls))
yhat_cv_test <- transform(yhat_cv_test, ridge = as.numeric(ridge))
yhat_cv_test1 <- apply(yhat_cv_test,1,mean)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_cv1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,29] <- n
n1
bias
out[3,30] <- MPE
RPE
out[3,31] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_cv_test1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,32] <- n
n1
bias
out[3,33] <- MPE
RPE
out[3,35] <- R
out[3,34] <- B1

#########################################################################################
############## Principal Component Analysis #############################################
library (pls)

##train the model
set.seed(123)
pcr_mod <- pcr(hs ~., data=train, ncomp=20, scale=TRUE)
summary(pcr_mod)
##predict y hat in cross validation
yhat_pcr <- pcr_mod %>% predict(x.train, ncomp=20)               #predict y hat
##predict y hat in external validation
yhat_test_pcr <- pcr_mod %>% predict(x.test, ncomp=20)     # predict y hat from test dataset

##measures of accuracy cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,36] <- n
n1
bias
out[3,37] <- MPE
RPE
out[3,38] <- R
B1

##measures of accuracy external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,39] <- n
n1
bias
out[3,40] <- MPE
RPE
out[3,42] <- R
out[3,41] <- B1

#########################################################################################
############## Projection Pursuit Regression ############################################
library (stats)

##trai the model
set.seed(123)
ppr_mod <- ppr(x.train, y.train, nterms=1, max.terms=2)
ppr_plot <- plot(ppr_mod)
## predict y cross validation
yhat_ppr <- ppr_mod %>% predict(x.train)               #predict y hat
## predict y external validation
yhat_test_ppr <- ppr_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,43] <- n
n1
bias
out[3,44] <- MPE
RPE
out[3,45] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,46] <- n
n1
bias
out[3,47] <- MPE
RPE
out[3,49] <- R
out[3,48] <- B1

#########################################################################################
############## Spike and Slab Regression ################################################
library (spikeslab)
library (plyr)

###train the model
spikeslab_mod <- spikeslab (hs ~., data = train, max.var = 1000)
## predict y in corss validation
yhat_spikeslab <- spikeslab_mod %>% predict(x.train)               #predict y hat
yhat_spikeslab <- yhat_spikeslab$yhat.bma
## predict y in external validation
yhat_test_spikeslab <- spikeslab_mod %>% predict(x.test)     # predict y hat from test dataset
yhat_test_spikeslab <- yhat_test_spikeslab$yhat.bma

# wavelength selected by spike and slab
attributes(spikeslab_mod)
spikeslab_mod$gnet
length(spikeslab_mod$gnet)
plot(spikeslab_mod$gnet, type="h")
sum(spikeslab_mod$gnet == 0)
ssr_bot3 <- as.data.frame(which(spikeslab_mod$gnet != 0))
dim(ssr_bot)

##define X and Y in cross validation##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,50] <- n
n1
bias
out[3,51] <- MPE
RPE
out[3,52] <- R
B1

##define X and Y in external validation##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,53] <- n
n1
bias
out[3,54] <- MPE
RPE
out[3,56] <- R
out[3,55] <- B1

#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.vector(y.train)
RF_mod3 <- randomForest (x= x.train, y= y.train, ntree = 500, mtry = 40)   #model code

## predict y
yhat_RF <- predict (RF_mod3 , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.numeric(yhat_RF)
yhat_RF_test <- predict (RF_mod3 , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.numeric(yhat_RF_test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,57] <- as.numeric(n)
n1
bias
out[3,58] <-MPE
RPE
out[3,59] <-R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,60] <-n
n1
bias
out[3,61] <-MPE
RPE
out[3,63] <-R
out[3,62] <-B1

############################################################################################
############################## Boosting ####################################################

library(gbm)
##train the model
boost2_mod <- gbm(hs ~.,data = train,distribution = "gaussian",n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict(boost2_mod,Xnew1,n.trees = 500)             #predict y hat
yhat1_test_boosting<-predict(boost2_mod,Xnew1_test,n.trees = 500)   #predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,64] <- n
n1
bias
out[3,65] <- MPE
RPE
out[3,66] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,67] <- n
n1
bias
out[3,68] <- MPE
RPE
out[3,70] <- R
out[3,69] <- B1

#########################################################################################
############## Neural Network #####################################
library(brnn)
##train model
set.seed(123)
brnn_mod <- brnn(x.train, y.train)

###prediction
yhat_brnn <- brnn_mod %>% predict(x.train)               #predict y hat
yhat_test_brnn <- brnn_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,71] <- n
n1
bias
out[3,72] <- MPE
RPE
out[3,73] <- R
B1


##measures of accuracy in cross validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[3,74] <- n
n1
bias
out[3,75] <- MPE
RPE
out[3,77] <- R
out[3,76] <- B1

##############################################################################
############################################################################
###############define x.train y.train x.test y.test#######################
x.train = as.matrix(train[-fold4,-1])
y.train = train[-fold4,1]
x.test = as.matrix(train[fold4,-1])
y.test= train [fold4,1]

############################################################################################
###################################### Apply PLS############################################
library("pls")
library(magrittr)
library("parallel")

# Generic code set up: put training and test x and y data into data frames
train.pls = data.frame(y = y.train)
train.pls$x = x.train
test.pls = data.frame(y = y.test)
test.pls$x = x.test

# Fit a PLS regression to the training data
# Use multiple cores for the CV if available
pls.options(parallel = makeCluster(detectCores(), type = "PSOCK"))
res.pls = plsr(y ~ x, ncomp=20, data = train.pls, validation = "LOO", scale = TRUE)
stopCluster(pls.options()$parallel)

# Examine results through plots
summary(res.pls)   
plot(RMSEP(res.pls), legendpos = "topright",main="Fat Variable PLS")
axis(side = 1, at=1:10)
plot(res.pls, ncomp = 4, asp = 1, line = TRUE, main="4  Component Fit")
plot(res.pls, plottype = "scores", comps = 1:4, main="Scores")
plot(res.pls, "loadings", comps = 1:4, legendpos = "bottomright", xlab = "spectrum", main="Fat Dataset 1: PLS Components")
abline(h = 0)

# Predict y hat in cross validation
y_hat_pls = predict(res.pls, ncomp = 10, newdata = train.pls$x)

# Predict y hat in external validation
y_hat_pls_test = predict(res.pls, ncomp = 10, newdata = test.pls$x)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_pls #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,1] <- n
n1
bias
out[4,2] <- MPE
RPE
out[4,3] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_pls_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,4] <- n
n1
bias
out[4,5] <- MPE
RPE
out[4,7] <- R
out[4,6] <- B1

############################################################################################
###############################Apply Ridge Regression ######################################
library(glmnet)
#set lambda
set.seed(123)
lambdas_ridge <- 10^seq(-2, 2, length.out = 1000) 
#TRAIN model
ridge_cv <- cv.glmnet(x.train, y.train, alpha = 0, lambda = lambdas_ridge,
                      standardize = TRUE, nfolds = 10)
lambda_ridge_cv <- ridge_cv$lambda.min
model_ridge_cv <- glmnet(x.train, y.train, alpha = 0, lambda = lambda_ridge_cv, standardize = TRUE)

#Y hat in cross validation
y_hat_ridge_cv <- predict(model_ridge_cv, x.train)

# number of wavelength selected
beta_ridge<-model_ridge_cv$beta
Bot_ridge <- beta_ridge[which(model_ridge_cv$beta !=0),]
Bot_ridge<-as.matrix(Bot_ridge)
dim(Bot_ridge)
nrow(Bot_ridge)

#test
y_hat_ridge_cv_test_A<-predict(model_ridge_cv, x.test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,8] <- n
n1
bias
out[4,9] <- MPE
RPE
out[4,10] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_ridge_cv_test_A #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,11] <- n
n1
bias
out[4,12] <- MPE
RPE
out[4,14] <- R
out[4,13] <- B1

##################################################################################
###################################LASSO########################################## 
#cross-validation to select lambda
lambdas_lasso <- 10^seq(-3, 3, length.out = 1000)
#TRAIN the lasso model
lasso_cv <- cv.glmnet(x.train, y.train, alpha = 1, lambda = lambdas_lasso,
                      standardize = TRUE, nfolds = 10)
lambda_lasso_cv <- lasso_cv$lambda.min
model_lasso_cv <- glmnet(x.train, y.train, alpha = 1, lambda = lambda_lasso_cv, standardize = TRUE)
# ppredict y hat in cross validation
y_hat_lasso_cv <- predict(model_lasso_cv, x.train)
# number of wavelength selected
beta_lasso<-model_lasso_cv$beta
Bot_lasso <- beta_lasso[which(model_lasso_cv$beta !=0),] 
Bot_lasso4<-as.matrix(Bot_lasso)
dim(Bot_lasso)
nrow(Bot_lasso)

#test 
y_hat_lasso_cv_test <- predict(model_lasso_cv, x.test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,15] <- n
n1
bias
out[4,16] <- MPE
RPE
out[4,17] <- R
B1


##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_lasso_cv_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,18] <- n
n1
bias
out[4,19] <- MPE
RPE
out[4,21] <- R
out[4,20] <- B1

####################################################################################
################################Apply Elastic Net###################################

library(glmnet)
# train the en model
fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5,intercept=FALSE)
fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
                          family="gaussian")
fit.elnetLambda<-fit.elnet.cv$lambda.min
fit.elnet2 <- glmnet(x.train, y.train, family="gaussian", alpha=.5, lambda =fit.elnetLambda, 
                     intercept=FALSE)
beta<-fit.elnet2$beta
# number of wavelength selected
Bot_en <- beta[which(fit.elnet2$beta !=0),] 
Bot_en4<-as.matrix(Bot_en)
dim(Bot_en)
nrow(Bot_en)
# y hat in cross validation
y_hat_en<-x.train%*%beta

# y hat in external validation
y_hat_en_test<-x.test%*%beta

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= y_hat_en #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,22] <- n
n1
bias
out[4,23] <- MPE
RPE
out[4,24] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= y_hat_en_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,25] <- n
n1
bias
out[4,26] <- MPE
RPE
out[4,28] <- R
out[4,27] <- B1

##################################################################################
############################### averaging model###################################

library(broom)
library(tidyverse)

#evaluate y hat as average of regression models in cross validation
yhat_cv <- cbind(y_hat_lasso_cv, y_hat_en, y_hat_pls, y_hat_ridge_cv)
yhat_cv <- as.matrix(yhat_cv)
yhat_cv <- as.data.frame(yhat_cv)
names(yhat_cv)[1] <- "lasso"
names(yhat_cv)[2] <- "en"
names(yhat_cv)[3] <- "pls"
names(yhat_cv)[4] <- "ridge"
yhat_cv <- transform(yhat_cv, lasso = as.numeric(lasso))
yhat_cv <- transform(yhat_cv, en = as.numeric(en))
yhat_cv <- transform(yhat_cv, pls = as.numeric(pls))
yhat_cv <- transform(yhat_cv, ridge = as.numeric(ridge))
yhat_cv1 <- apply(yhat_cv,1,mean)

#evaluate y hat as average of regression models in external validation
yhat_cv_test <- cbind(y_hat_lasso_cv_test, y_hat_en_test, y_hat_pls_test, y_hat_ridge_cv_test_A)
yhat_cv_test <- as.matrix(yhat_cv_test)
yhat_cv_test <- as.data.frame(yhat_cv_test)
names(yhat_cv_test)[1] <- "lasso"
names(yhat_cv_test)[2] <- "en"
names(yhat_cv_test)[3] <- "pls"
names(yhat_cv_test)[4] <- "ridge"
yhat_cv_test <- transform(yhat_cv_test, lasso = as.numeric(lasso))
yhat_cv_test <- transform(yhat_cv_test, en = as.numeric(en))
yhat_cv_test <- transform(yhat_cv_test, pls = as.numeric(pls))
yhat_cv_test <- transform(yhat_cv_test, ridge = as.numeric(ridge))
yhat_cv_test1 <- apply(yhat_cv_test,1,mean)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_cv1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,29] <- n
n1
bias
out[4,30] <- MPE
RPE
out[4,31] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_cv_test1 #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,32] <- n
n1
bias
out[4,33] <- MPE
RPE
out[4,35] <- R
out[4,34] <- B1

#########################################################################################
############## Principal Component Analysis #############################################
library (pls)

##train the model
set.seed(123)
pcr_mod <- pcr(hs ~., data=train, ncomp=20, scale=TRUE)
summary(pcr_mod)
##predict y hat in cross validation
yhat_pcr <- pcr_mod %>% predict(x.train, ncomp=20)               #predict y hat
##predict y hat in external validation
yhat_test_pcr <- pcr_mod %>% predict(x.test, ncomp=20)     # predict y hat from test dataset

##measures of accuracy cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)

xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,36] <- n
n1
bias
out[4,37] <- MPE
RPE
out[4,38] <- R
B1

##measures of accuracy external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_pcr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,39] <- n
n1
bias
out[4,40] <- MPE
RPE
out[4,42] <- R
out[4,41] <- B1

#########################################################################################
############## Projection Pursuit Regression ############################################
library (stats)

##trai the model
set.seed(123)
ppr_mod <- ppr(x.train, y.train, nterms=1, max.terms=2)
ppr_plot <- plot(ppr_mod)
## predict y cross validation
yhat_ppr <- ppr_mod %>% predict(x.train)               #predict y hat
## predict y external validation
yhat_test_ppr <- ppr_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,43] <- n
n1
bias
out[4,44] <- MPE
RPE
out[4,45] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_ppr #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,46] <- n
n1
bias
out[4,47] <- MPE
RPE
out[4,49] <- R
out[4,48] <- B1

#########################################################################################
############## Spike and Slab Regression ################################################
library (spikeslab)
library (plyr)

###train the model
spikeslab_mod <- spikeslab (hs ~., data = train, max.var = 1000)
## predict y in corss validation
yhat_spikeslab <- spikeslab_mod %>% predict(x.train)               #predict y hat
yhat_spikeslab <- yhat_spikeslab$yhat.bma
## predict y in external validation
yhat_test_spikeslab <- spikeslab_mod %>% predict(x.test)     # predict y hat from test dataset
yhat_test_spikeslab <- yhat_test_spikeslab$yhat.bma

# wavelength selected by spike and slab
attributes(spikeslab_mod)
spikeslab_mod$gnet
length(spikeslab_mod$gnet)
plot(spikeslab_mod$gnet, type="h")
sum(spikeslab_mod$gnet == 0)
ssr_bot4 <- as.data.frame(which(spikeslab_mod$gnet != 0))
dim(ssr_bot)

##define X and Y in cross validation##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,50] <- n
n1
bias
out[4,51] <- MPE
RPE
out[4,52] <- R
B1

##define X and Y in external validation##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_spikeslab #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,53] <- n
n1
bias
out[4,54] <- MPE
RPE
out[4,56] <- R
out[4,55] <- B1

#############################################################################
########################## Apply Random Forest###############################
library(randomForest)

##train the model
y.train <- as.vector(y.train)
RF_mod4 <- randomForest (x= x.train, y= y.train, ntree = 500, mtry = 40)   #model code

## predict y
yhat_RF <- predict (RF_mod4 , newdata= x.train)           #predict yhat CV
yhat_RF <- as.matrix(yhat_RF)
yhat_RF <- as.numeric(yhat_RF)
yhat_RF_test <- predict (RF_mod4 , newdata = x.test) #predict yhat EV
yhat_RF_test <- as.matrix(yhat_RF_test)
yhat_RF_test <- as.numeric(yhat_RF_test)

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_RF #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,57] <- as.numeric(n)
n1
bias
out[4,58] <-MPE
RPE
out[4,59] <-R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_RF_test #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,60] <-n
n1
bias
out[4,61] <-MPE
RPE
out[4,63] <-R
out[4,62] <-B1

############################################################################################
############################## Boosting ####################################################

library(gbm)
##train the model
boost2_mod <- gbm(hs ~.,data = train,distribution = "gaussian",n.trees = 500,
                  shrinkage = 0.01, interaction.depth = 4)              #model code

## predic y
Xnew1 <- as.data.frame(x.train)
Xnew1_test <- as.data.frame(x.test)
yhat1_boosting<-predict(boost2_mod,Xnew1,n.trees = 500)             #predict y hat
yhat1_test_boosting<-predict(boost2_mod,Xnew1_test,n.trees = 500)   #predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat1_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,64] <- n
n1
bias
out[4,65] <- MPE
RPE
out[4,66] <- R
B1

##measures of accuracy in external validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat1_test_boosting #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,67] <- n
n1
bias
out[4,68] <- MPE
RPE
out[4,70] <- R
out[4,69] <- B1

#########################################################################################
############## Neural Network #####################################
library(brnn)
##train model
set.seed(123)
brnn_mod <- brnn(x.train, y.train)

###prediction
yhat_brnn <- brnn_mod %>% predict(x.train)               #predict y hat
yhat_test_brnn <- brnn_mod %>% predict(x.test)     # predict y hat from test dataset

##measures of accuracy in cross validation##
##define X and Y##

x= y.train  #values of the character analysed eg fat_content#
y= yhat_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,71] <- n
n1
bias
out[4,72] <- MPE
RPE
out[4,73] <- R
B1

##measures of accuracy in cross validation##
##define X and Y##

x= y.test  #values of the character analysed eg fat_content#
y= yhat_test_brnn #values of the predicted eg fat_predicted#

##
n <- length(x)
n1 = length(y)
xbar = mean(x)
covariance = (sum((x-mean(x))*(y-mean(y))))/(n-1)  #coavariance x,y
varx = (sum((x-mean(x))**2))/(n-1)                 #variance x
vary = (sum((y-mean(y))**2))/(n-1)                 #variance y
bias = ((sum(x))/n)-((sum(y))/n1)
resid = sum(x-y)
SSE = sum((x-y)**2)
MPSE = SSE / (n-1)                                 #Mean Square Prediction Error; Fuentes-Pila 1995
MPE = sqrt (MPSE)                                  #Mean Prediction Error; Fuentes-Pila 1995
R = covariance / (sqrt(varx*vary))                 
B1 = covariance / vary                             #Slope
B0 = (mean(y)) - (B1* mean(y))                       #Intercept

out[4,74] <- n
n1
bias
out[4,75] <- MPE
RPE
out[4,77] <- R
out[4,76] <- B1

#########################################################################################
########################################################################################
########################results##########################################################

apply(out, 2, summary)   #mean, max, min, and median of each predictors for each model
boxplot(out)             #boxplot of the results
apply(out, 2, sd)        #standard deviation of the results

out <- as.data.frame(out)   #save the results table as data frame for export it

library("writexl")
write_xlsx(out,"C:/Users/Utente/Desktop/lavoro teagasc/out beta casein.xlsx")  #export the results in excel file


#other interesting results

#number of wavelenghts selected by the models
#Lasso
dim(Bot_lasso1)
dim(Bot_lasso2)
dim(Bot_lasso3)
dim(Bot_lasso4)

#EN
dim(Bot_en1)
dim(Bot_en2)
dim(Bot_en3)
dim(Bot_en4)

#Spike and Slab
dim(ssr_bot1)
dim(ssr_bot2)
dim(ssr_bot3)
dim(ssr_bot4)

#Variable importance in Random Forest
importance(RF_mod1)
importance(RF_mod2)
importance(RF_mod3)
importance(RF_mod4)