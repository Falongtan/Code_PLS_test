
###library packages
library(glmnet)
library(psych)
library(LassoSIR)
library(VariableScreening)
library(pROC)
####Predict whether the new molecule is a musk or non-musk
real_data <- read.csv("E:/桌面/local's code/code of simulation and real data/clean1.csv",header=T)
#response variable
Y <- real_data$Y
Y <- as.numeric(Y)
#X
X_0 <- real_data[,1:(ncol(real_data)-1)]
X <- scale(X_0)
n <- nrow(X)
p <- ncol(X)


#####test
set.seed(222)
#data spliying
index_x<-sample(1:n,floor(n/2),replace = FALSE)  
#data1
x1 <- as.matrix(X[-index_x,]) 
x2 <- as.matrix(X[index_x,])
n1 <- nrow(x1)
#data2
y1 <- as.matrix(Y[-index_x])
y2 <- as.matrix(Y[index_x]) 
n2 <- nrow(x2)

#residuals and projections based on x1 y1 and dimension reduction
#beta_hat projection and residuals based on lasso and data x1 y1
lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
len_beta1_non0 <- length(index_beta1_non0)
if(len_beta1_non0 == 0){
  pred1 <- mean(y1)
  U1 <- y1 - pred1                                                         # residuals based on x1 y1
}else{
  x1_sec <- x1[,index_beta1_non0]                                          # second estimation
  y1_x1_sec <- as.data.frame(cbind(y1,x1_sec))
  colnames(y1_x1_sec) <- c("y1",colnames(x1_sec))
  sec_model1 <- glm(y1~.,data=y1_x1_sec, family = binomial(link = "logit"))
  pred1 = predict(sec_model1, newx = x2_sec, type="response")
  pred1 = matrix(unname(pred1), ncol = 1) 
  U1 <- y1 - pred1                                                           # residuals based on x1 y1
}

# projections based on screening and lassosir for x1 U1
screen_num1 <- floor(n1/log(n1))
sdr_Upro1 <- tryCatch({
  if(p <= screen_num1){
    # sir projection without screening
    sir_U1 <- LassoSIR(x1, U1, H=10, choosing.d="automatic",
                       solution.path=FALSE, categorical=FALSE, nfolds=5,
                       screening=FALSE)
    sir_Ubeta1<- sir_U1$beta
    sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
  }else{
    #screening
    rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
    index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
    x1_Uscreen <- x1[,index_U1]
    sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=10, choosing.d="automatic",
                            solution.path=FALSE, categorical=FALSE, nfolds=5,
                            screening=FALSE)
    sir_Ubeta1<- sir_Ulasso1$beta
    sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
    sir_Upro1[index_U1,] <- sir_Ubeta1
    sir_Upro1/sqrt(colSums(sir_Upro1^2))
  }
}, error =function(e){ 
  NA })

# projections based on screening and lassosir for x1 y1
sdr_ypro1 <-  tryCatch({
  if(p <= screen_num1){
    # sir projection without screening
    sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic",
                            solution.path=FALSE, categorical=TRUE, nfolds=5,
                            screening=TRUE)
    sir_ybeta1<- sir_ylasso1$beta
    sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
  }else{
    #screening
    rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                     # screening method based on distance correlation
    index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
    x1_yscreen <- x1[,index_y1]
    sir_ylasso1 <- LassoSIR(x1_yscreen, y1, H=2, choosing.d="automatic",
                            solution.path=FALSE, categorical=TRUE, nfolds=5,
                            screening=TRUE)
    sir_ybeta1<- sir_ylasso1$beta
    sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
    sir_ypro1[index_y1,] <- sir_ybeta1
    sir_ypro1/sqrt(colSums(sir_ypro1^2))
  }
}, error =function(e){ 
  NA })

# residuals and projections based on x2 y2 and dimension reduction
# beta_hat projection and residuals based on lasso and data x2 y2
lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso)!=0]
len_beta2_non0 <- length(index_beta2_non0)
if(len_beta2_non0 == 0){
  pred2 <- mean(y2)
  U2 <- y2-pred2                                                          # residuals based on x2 y2
}else{
  x2_sec <- x2[,index_beta2_non0]
  y2_x2_sec <- as.data.frame(cbind(y2,x2_sec))
  colnames(y2_x2_sec) <- c("y2",colnames(x2_sec))
  sec_model2 <- glm(y2~.,data=y2_x2_sec, family = binomial(link = "logit"))
  # sec_model2 <- glm(y2~x2_sec,family = binomial(link = "logit"))
  pred2 = predict(sec_model2, newx = x2_sec,type="response")
  pred2 = matrix(unname(pred2), ncol = 1)
  U2 <- y2-pred2
}
#sir projections based on screening and x2 U2
screen_num2 <- floor(n2/log(n2))
sdr_Upro2 <- tryCatch({
  if(p <= screen_num2){
    # sir projection without screen
    sir_Ulasso2 <- LassoSIR(x2, U2, H=10, choosing.d="automatic",
                            solution.path=FALSE, categorical=FALSE, nfolds=5,
                            screening=FALSE)
    sir_Ubeta2<- sir_Ulasso2$beta
    sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
  }else{
    #screening
    rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
    index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
    x2_Uscreen <- x2[,index_U2]
    sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=10, choosing.d="automatic",
                            solution.path=FALSE, categorical=FALSE, nfolds=5,
                            screening=FALSE)
    sir_Ubeta2 <- sir_Ulasso2$beta
    sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
    sir_Upro2[index_U2,] <- sir_Ubeta2
    sir_Upro2/sqrt(colSums(sir_Upro2^2))
  }
}, error =function(e){ 
  NA })


#sir projections based on screening and x2 y2
sdr_ypro2 <-tryCatch({
  if(p <= screen_num2){
    # sir projection without screen
    sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic",
                            solution.path=FALSE, categorical=TRUE, nfolds=5,
                            screening=TRUE)
    sir_ybeta2<- sir_ylasso2$beta
    sir_ybeta2/sqrt(colSums(sir_ybeta2^2))
  }else{
    #screening
    rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                     # screening method based on distance correlation
    index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
    x2_yscreen <- x2[,index_y2]
    sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=2, choosing.d="automatic",
                            solution.path=FALSE, categorical=TRUE, nfolds=5,
                            screening=TRUE)
    sir_ybeta2 <- sir_ylasso2$beta
    sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
    sir_ypro2[index_y2,] <- sir_ybeta2
    sir_ypro2/sqrt(colSums(sir_ypro2^2))
  }
}, error =function(e){ 
  NA })
#construct PLS test statistic based on x1 y1 and projections based on x2 y2
pro2 <- cbind(sdr_Upro2, sdr_ypro2)                                           # projections based on x2 y2
###To ensure that the projections are linearly independent
pro2 <- tryCatch({
  qr_obj2 <- qr(pro2)
  if( qr_obj2$rank == ncol(pro2) ){
    pro2
  }else{
    pivot2 <- qr_obj2$pivot[1:qr_obj2$rank]
    pro2[, pivot2, drop = FALSE]
  }
}, error =function(e){
  NA })
pro2 <- as.matrix(pro2)
pro2_num <- ncol(pro2)                                                        # the number of projections
h1 <- 1*(n1^(-2/9))                                                         # bandwidth

pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) # p-value matrix
errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
for(q in 1:pro2_num){
  if(any(is.na(pro2[,q])) ){
    pval_matrix1[,q] <- NA
  }else{
    x_pro1 <- x1%*%pro2[,q]
    x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
    indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
    kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel 
    #PLS test statistics
    Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
    pval1 <- 1-pnorm(Tn1)
    pval_matrix1[,q] <- pval1
  }
  
}

#construct PLS test statistics based on x2 y2 and projections based on x1 y1
pro1 <- cbind(sdr_Upro1, sdr_ypro1)                                           # projections based on x1 y1
###To ensure that the projections are linearly independent
pro1 <- tryCatch({                                                         
  qr_obj1 <- qr(pro1)
  if( qr_obj1$rank == ncol(pro1) ){
    pro1
  }else{
    pivot1 <- qr_obj1$pivot[1:qr_obj1$rank]
    pro1[, pivot1, drop = FALSE]
  }
}, error =function(e){
  NA })
pro1 <- as.matrix(pro1)
pro1_num <- ncol(pro1)                                                        # the number of projections
h2 <- 1*(n2^(-2/9))                                                         # bandwidth 

pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) # p-value matrix
errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
for(l in 1:pro1_num){
  if(any(is.na(pro1[,l])) ){
    pval_matrix2[,l] <- NA
  }else{
    x_pro2 <- x2%*%pro1[,l]
    x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2   # kernel function matrix
    indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
    kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                 # Epanechnikov kernel 
    # PLS test statistics
    Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
    pval2 <- 1-pnorm(Tn2)
    pval_matrix2[,l] <- pval2
  }
  
}
#p-value matrix
pval_matrix <- cbind(pval_matrix1, pval_matrix2)
pval_num <- pro1_num + pro2_num
#fihser
pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
#min
pval_min1 <-   min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
pval_min2 <-  min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
#result
logistic_result <- c(pval_fisher,pval_min)
print(logistic_result)





#####prediction
sss <- 50
result_prediction <- matrix(nrow=sss,ncol=5)

suppressMessages({
  for(i in 1:sss){
    set.seed(i)
    data_y_x <- as.data.frame(cbind(Y,X))
    # stratified 75/25 split using base R (replaces caTools::sample.split)
    split <- logical(nrow(data_y_x))
    for(cls in unique(data_y_x$Y)){
      idx <- which(data_y_x$Y == cls)
      split[sample(idx, size = floor(0.75*length(idx)))] <- TRUE
    }
    data_train <- subset(data_y_x, split == TRUE)
    data_test <- subset(data_y_x, split == FALSE)

    p=ncol(X)
    train_x <- data_train[,-1]
    colnames(train_x) <- colnames(X)

    train_y <- data_train[,1]
    test_x <- data_test[,-1]
    colnames(test_x) <- colnames(X)
    test_y <- data_test[,1]

  
    ##fit model
    train_lasso_model <- cv.glmnet(as.matrix(train_x), train_y, family="binomial", intercept = T)
    train_lasso_beta <- coef(train_lasso_model, s = "lambda.min")[-1]
    train_index_beta_non0 <- seq(1:p)[as.numeric(train_lasso_beta) != 0]
  
    train_len_beta_non0 <- length(train_index_beta_non0)
    if(train_len_beta_non0 == 0){
      result_prediction[i,1] <- NA
      result_prediction[i,2] <- NA
      result_prediction[i,3] <- NA
      result_prediction[i,4] <- NA
      result_prediction[i,5] <- NA
    }else{
      train_x_sec <- train_x[,train_index_beta_non0]                                          # second estimation
      train_y_x_sec <- as.data.frame(cbind(train_y,train_x_sec))
      train_sec_model <- glm(train_y~.,data=train_y_x_sec, family = binomial(link = "logit"))
      train_sec_beta = unname(train_sec_model$coefficients)[-1]
      #  probability of prediction
      train_probs <- predict(train_sec_model,newdata=as.data.frame(train_x),type = "response")
      #roc of train
      roc_train <- roc(as.vector(train_y), train_probs)
      #Youden’s Rule
      bestp <- roc_train$thresholds[
        which.max(roc_train$sensitivities+roc_train$specificities-1)
      ]

      predictions <- predict(train_sec_model, newdata = as.data.frame(test_x[, train_index_beta_non0]), type="response")
      Prob=predictions
      # Prob=exp(predictions)/(1+exp(predictions))
      # Prob[which(is.nan(Prob))]=1
      predicted_classes <- ifelse(Prob>bestp,1,0)
      confusion_matrix <- table(as.vector(test_y),predicted_classes)
      if(ncol(confusion_matrix)==1){
        if(colnames(confusion_matrix)=="1"){
          result_prediction[i,4] <- confusion_matrix[2,1]/sum(confusion_matrix)
          result_prediction[i,3] <- 1
          result_prediction[i,1] <- 0
          result_prediction[i,2] <- 0
        }else{
          result_prediction[i,4] <- confusion_matrix[1,1]/sum(confusion_matrix)
          result_prediction[i,1] <- 1
          result_prediction[i,3] <- 0
          result_prediction[i,2] <- 0
        }
      }else{
        #True Negative Rate
        result_prediction[i,1] <- (confusion_matrix[1,1])/(confusion_matrix[1,1]+confusion_matrix[1,2])
        #F1
        result_prediction[i,2] <- 2*((confusion_matrix[2,2])/(confusion_matrix[2,1]+confusion_matrix[2,2]))*((confusion_matrix[1,1])/(confusion_matrix[1,1]+confusion_matrix[1,2]))/
          ((confusion_matrix[2,2])/(confusion_matrix[2,1]+confusion_matrix[2,2])+(confusion_matrix[1,1])/(confusion_matrix[1,1]+confusion_matrix[1,2]))
        #True Positive Rate
        result_prediction[i,3] <- (confusion_matrix[2,2])/(confusion_matrix[2,1]+confusion_matrix[2,2])
        #Accuracy
        result_prediction[i,4] <- (confusion_matrix[1,1]+confusion_matrix[2,2])/sum(confusion_matrix)
      }
    
      # ROC
      roc_obj <- roc(as.vector(test_y), as.vector(Prob))
      #plot(roc_obj, main = "ROC Curve", col = "blue")
      auc_value <- auc(roc_obj)
      result_prediction[i,5] <- auc_value
    }
  
  
  
  }
})


###col1:True Negative Rate   col2:F1    col3:True Positive Rate
###col4: Accuracy   col5:ROC
result_prediction1 = apply(result_prediction,2,mean,na.rm=TRUE)
print(result_prediction1)
