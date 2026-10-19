

###generation squared terms and interaction terms
xpand_features <- function(X) {
  # to ensure X is a matrix
  if (!is.data.frame(X)) {
    X <- as.data.frame(X)
  }
  
  # colnames
  col_names <- names(X)
  n_cols <- ncol(X)
  
  # Create squared terms
  squares <- lapply(col_names, function(col) {
    X[[col]]^2
  })
  names(squares) <- paste0(col_names, "_sq")
  squares_df <- as.data.frame(squares)
  
  # create interaction terms
  cross_products <- list()
  cp_names <- character()
  
  for (i in 1:(n_cols-1)) {
    for (j in (i+1):n_cols) {
      col_i <- col_names[i]
      col_j <- col_names[j]
      cp_name <- paste0(col_i, "_x_", col_j)
      
      cross_products[[cp_name]] <- X[[col_i]] * X[[col_j]]
      cp_names <- c(cp_names, cp_name)
    }
  }
  
  cross_products_df <- as.data.frame(cross_products)
  names(cross_products_df) <- cp_names
  
  # cbind X, squared terms and interaction terms
  result <- cbind(X, squares_df, cross_products_df)
  return(result)
}

##library packages
library(readxl)
library(glmnet)
library(psych)
library(LassoSIR)
library(VariableScreening)


###read data
linear_data <- read_excel("E:/桌面/local's code/code of simulation and real data/CommViolPredUnnormalizedData.xlsx", sheet = 1, col_names = TRUE)
###drop missing values
linear_data_nona <- linear_data[which(rowSums(is.na(linear_data))==0),]
linear_data_final <- matrix(nrow=0,ncol=120)
for (i in 1:2197){
  if(any(linear_data_nona[i,]=="?")){
    linear_data_final <- linear_data_final
  }else{
    linear_data_final <- rbind(linear_data_final,linear_data_nona[i,])
  }
}
#response value
Y_matrix <- linear_data_final[,103:120]
Y_matrix <- as.data.frame(Y_matrix)
colnames(Y_matrix) <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9","Y10","Y11","Y12","Y13","Y14","Y15","Y16","Y17","Y18")
for(i in 1:8){
  o <- c(3,4,9,10,15,16,17,18)[i]
  Y_matrix[,o] <- as.numeric(Y_matrix[,o])
}
###X,Y
Y <- Y_matrix$Y17
Y <- as.numeric(Y)
#scale
Y <- scale(Y)
X_0 <- linear_data_final[,-c(103:120)]
X <- scale(X_0)
##sample size and dimension
n <- nrow(X)
p <- ncol(X)
##seed
set.seed(222)
###data splitting
index_x<-sample(1:n, floor(n/2),replace = FALSE)  
index_x<-sort(index_x)
p <- ncol(X)
#data1
x1 <- X[-index_x,] 
y1 <- as.matrix(Y[-index_x]) 
n1 <- nrow(x1)
#data2
x2 <- X[index_x,]
y2 <- as.matrix(Y[index_x])
n2 <- nrow(x2)
# residuals and projections based on x1 y1 and dimension reduction
# beta_hat projection and residuals based on lasso and data x1 y1
#Since the data has been standardized, the intercept is set to F.
lasso_model1 <- cv.glmnet(as.matrix(x1), y1, family="gaussian", intercept = F)
lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
len_beta1_non0 <- length(index_beta1_non0)
if(len_beta1_non0 == 0){
  U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
}else{
  x1_sec <- x1[,index_beta1_non0]                                            # second estimation
  sec_model1 <- glm(y1~x1_sec-1, family = gaussian)
  pred1 = predict(sec_model1, newx = x1_sec, type="response")
  U1 <- y1-pred1   # residual based on x1 y1
}
#screening number
screen_num1 <- floor(n1/log(n1))
# projections based on screening and lassosir or seas for x1 U1
sdr_Upro1 <- tryCatch({
  if(p <= screen_num1){
    # projection based on lassosir without screening
    sir_U1 <- LassoSIR(x1, U1, H=10, choosing.d="automatic", solution.path=FALSE, 
                       categorical=FALSE, nfolds=5, screening=FALSE)
    sir_Ubeta1<- sir_U1$beta
    sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
  }else{
    #screening
    rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
    index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
    x1_Uscreen <- x1[,index_U1]
    
    # projection based on lassosir after screening
    sir_U1 <- LassoSIR(x1_Uscreen, U1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta1<- sir_U1$beta
    sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
    sir_Upro1[index_U1,] <- sir_Ubeta1
    sir_Upro1/sqrt(colSums(sir_Upro1^2))
  }
},error=function(e){
  NA
})
sdr_ypro1 <- tryCatch({
  if(p <= screen_num1){
    # laosssir projection without screening
    sir_y1 <- LassoSIR(x1, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
  }else{
    #screening
    rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                     # screening method 
    index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
    x1_yscreen <- x1[,index_y1]
    
    # projection based on lassosir after screening
    sir_y1 <- LassoSIR(x1_yscreen, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
    sir_ypro1[index_y1,] <- sir_ybeta1
    sir_ypro1/sqrt(colSums(sir_ypro1^2))
  }
},error = function(e){
  NA
})



# residuals and projections based on x2 y2 and dimension reduction
# residuals based on lasso and data x2 y2
#Since the data has been standardized, the intercept is set to F.
lasso_model2 <- cv.glmnet(as.matrix(x2), y2, family="gaussian", intercept = F)
lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                   # estimated beta0 without intercept 
index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  # index of beta_n with nonzero components
len_beta2_non0 <- length(index_beta2_non0)
if(len_beta2_non0 == 0){
  U2 <- y2-mean(y2)                                                       # residual based on x2 y2
}else{
  x2_sec <- x2[,index_beta2_non0]
  sec_model2 = glm(y2~x2_sec-1,family = gaussian)                           # second estimation
  pred2 = predict(sec_model2, newx = x2_sec, type="response")
  U2 <- y2-pred2   # residual based on x1 y1                                                  # residual based on x2 y2
}

#sir projections based on screening and x2 U2
screen_num2 <- floor(n2/log(n2))
sdr_Upro2 <- tryCatch({
  if(p <= screen_num2){
    # lassosir projection without screen
    sir_U2 <- LassoSIR(x2, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta2<- sir_U2$beta
    sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))                     # sdr projection based on lassosir and x2 U2
  }else{
    #screening
    rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
    index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
    x2_Uscreen <- x2[,index_U2]
    
    # projection based on lassosir after screening
    sir_U2 <- LassoSIR(x2_Uscreen, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta2 <- sir_U2$beta
    sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
    sir_Upro2[index_U2,] <- sir_Ubeta2
    sir_Upro2/sqrt(colSums(sir_Upro2^2))                       # sdr projection based on lassosir and x2 U2
  }
},error=function(e){
  NA
})

#sir projections based on screening and x2 y2
sdr_ypro2 <- tryCatch({
  if(p <= screen_num2){
    # Lassosir projection without screen
    sir_y2 <- LassoSIR(x2, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta2<- sir_y2$beta
    sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                     # sdr projection based on lassosir and x2 y2
  }else{ 
    #screening
    rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                     # screening method based on distance correlation
    index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
    x2_yscreen <- x2[,index_y2]
    
    # projection based on lassosir after screening
    sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                            nfolds=5, screening=FALSE)
    sir_ybeta2 <- sir_ylasso2$beta
    sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
    sir_ypro2[index_y2,] <- sir_ybeta2
    sir_ypro2/sqrt(colSums(sir_ypro2^2))                       # sdr projection based on lassosir and x2 y2
  }
},error=function(e){
  NA
})

#construct test statistic based on x1 y1 and projections based on x2 y2
pro2 <- cbind(sdr_Upro2, sdr_ypro2)  
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

pval_matrix <- cbind(pval_matrix1, pval_matrix2)
pval_num <- pro1_num + pro2_num


pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  

pval_min1 <-   min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
pval_min2 <-  min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
linear_result_lasso <- c(pval_fisher,pval_min)
print(linear_result_lasso)


######scatterplot_1
p <- ncol(X)
#Since the data has been standardized, the intercept is set to F.
lasso_model <- cv.glmnet(as.matrix(X), Y, family="gaussian", intercept = F)
beta_fit <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric( beta_fit) != 0]
x_sec <- X[,index_beta_non0]                                          
sec_model = glm(Y~x_sec-1,family = gaussian)                           

pred = predict(sec_model, newx = x_sec, type="response")
U <- Y-pred  
xbeta <- x_sec%*%sec_model$coefficients   #[-1]+sec_model$coefficients[1]
#####plot_1
pdf("linear-fit.pdf", width = 14, height = 14)   # auto-save figure (14 x 14 in = 1008 pt)
plot_range <- range(c(Y, xbeta))
plot(xbeta, Y, 
     xlab = "",  
     ylab = "Y",
     pch = 16,
     ylim = plot_range,
     cex.lab = 2.5,
     cex.axis = 2,
     mgp = c(2, 0.8, 0)  
)
mtext(expression(hat(beta)^T*X), side = 1,  
      line = 3,  
      cex = 2.5)
dev.off()




###############Polynomial Regression
p <- ncol(X)
#select variable
set.seed(222)
#Since the data has been standardized, the intercept is set to F.
lasso_model <- cv.glmnet(as.matrix(X), Y, family="gaussian", intercept = F)
beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
x_sec <- X[,index_beta_non0]   
###generate quadratic term and interaction term
X_quad_inter <- xpand_features(x_sec)
p <- ncol(X_quad_inter)
#colnames
new_names <- paste("col", 1:p, sep = "_")
colnames(X_quad_inter) <- new_names
#scale
X_quad_inter[,(ncol(x_sec)+1):ncol(X_quad_inter)] <- scale(X_quad_inter[,(ncol(x_sec)+1):ncol(X_quad_inter)])
####data splitting
####using the same index
x1 <- X_quad_inter[-index_x,] 
x2 <- X_quad_inter[index_x,]
x1 <- as.matrix(x1)
x2 <- as.matrix(x2)
n1 <- nrow(x1)
n2 <- nrow(x2)
# residuals and projections based on x1 y1 and dimension reduction
# beta_hat projection and residuals based on lasso and data x1 y1
#Since the data has been standardized, the intercept is set to F.
lasso_model1 <- cv.glmnet(as.matrix(x1), y1, family="gaussian", intercept = F)
lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
len_beta1_non0 <- length(index_beta1_non0)
if(len_beta1_non0 == 0){
  U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
}else{
  x1_sec <- as.matrix(x1[,index_beta1_non0])                                            # second estimation
  sec_model1 <- glm(y1~x1_sec-1, family = gaussian)
  pred1 = predict(sec_model1, newx = x1_sec, type="response")
  U1 <- y1-pred1
}

# projections based on screening and lassosir or seas for x1 U1
screen_num1 <- floor(n1/log(n1))
sdr_Upro1 <- tryCatch({
  if(p <= screen_num1){
    # projection based on lassosir without screening
    sir_U1 <- LassoSIR(x1, U1, H=10, choosing.d="automatic", solution.path=FALSE, 
                       categorical=FALSE, nfolds=5, screening=FALSE)
    sir_Ubeta1<- sir_U1$beta
    sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
  }else{
    #screening
    rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
    index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
    x1_Uscreen <- x1[,index_U1]
    
    # projection based on lassosir after screening
    sir_U1 <- LassoSIR(x1_Uscreen, U1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta1<- sir_U1$beta
    sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
    sir_Upro1[index_U1,] <- sir_Ubeta1
    sir_Upro1/sqrt(colSums(sir_Upro1^2))
  }
},error=function(e){
  NA
})

# projections based on screening and lassosir or seas for x1 y1
sdr_ypro1 <- tryCatch({
  if(p <= screen_num1){
    # laosssir projection without screening
    sir_y1 <- LassoSIR(x1, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
  }else{
    #screening
    rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                     # screening method 
    index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
    x1_yscreen <- x1[,index_y1]
    
    # projection based on lassosir after screening
    sir_y1 <- LassoSIR(x1_yscreen, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
    sir_ypro1[index_y1,] <- sir_ybeta1
    sir_ypro1/sqrt(colSums(sir_ypro1^2))
  }
},error = function(e){
  NA
})

# residuals and projections based on x2 y2 and dimension reduction
# residuals based on lasso and data x2 y2
#Since the data has been standardized, the intercept is set to F.
lasso_model2 <- cv.glmnet(as.matrix(x2), y2, family="gaussian", intercept = F)
lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                   # estimated beta0 without intercept 
index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  # index of beta_n with nonzero components
len_beta2_non0 <- length(index_beta2_non0)
if(len_beta2_non0 == 0){
  U2 <- y2-mean(y2)                                                       # residual based on x2 y2
}else{
  x2_sec <- as.matrix(x2[,index_beta2_non0])
  sec_model2 = glm(y2~x2_sec-1,family = gaussian)                           # second estimation
  pred2 = predict(sec_model2, newx = x2_sec, type="response")
  U2 <- y2-pred2
}

#sir projections based on screening and x2 U2
screen_num2 <- floor(n2/log(n2))
sdr_Upro2 <- tryCatch({
  if(p <= screen_num2){
    # lassosir projection without screen
    sir_U2 <- LassoSIR(x2, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta2<- sir_U2$beta
    sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))                     # sdr projection based on lassosir and x2 U2
  }else{
    #screening
    rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
    index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
    x2_Uscreen <- x2[,index_U2]
    
    # projection based on lassosir after screening
    sir_U2 <- LassoSIR(x2_Uscreen, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_Ubeta2 <- sir_U2$beta
    sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
    sir_Upro2[index_U2,] <- sir_Ubeta2
    sir_Upro2/sqrt(colSums(sir_Upro2^2))                       # sdr projection based on lassosir and x2 U2
  }
},error=function(e){
  NA
})

#sir projections based on screening and x2 y2
sdr_ypro2 <- tryCatch({
  if(p <= screen_num2){
    # Lassosir projection without screen
    sir_y2 <- LassoSIR(x2, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                       nfolds=5, screening=FALSE)
    sir_ybeta2<- sir_y2$beta
    sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                     # sdr projection based on lassosir and x2 y2
  }else{ 
    #screening
    rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                     # screening method based on distance correlation
    index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
    x2_yscreen <- x2[,index_y2]
    
    # projection based on lassosir after screening
    sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                            nfolds=5, screening=FALSE)
    sir_ybeta2 <- sir_ylasso2$beta
    sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
    sir_ypro2[index_y2,] <- sir_ybeta2
    sir_ypro2/sqrt(colSums(sir_ypro2^2))                       # sdr projection based on lassosir and x2 y2
  }
},error=function(e){
  NA
})

#construct test statistic based on x1 y1 and projections based on x2 y2
pro2 <- cbind(sdr_Upro2, sdr_ypro2)              # projections based on x1 y1
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

#construct test statistics based on x2 y2 and projections based on x1 y1
pro1 <- cbind(sdr_Upro1, sdr_ypro1)                           # projections based on x1 y1
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

pval_matrix <- cbind(pval_matrix1, pval_matrix2)
pval_num <- pro1_num + pro2_num

###fisher 
pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
###min
pval_min1 <-   min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
pval_min2 <-  min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
#result
nonlinear_result_lasso <- c(pval_fisher,pval_min)
print(nonlinear_result_lasso)



######scatterplot_2
p <- ncol(X_quad_inter)
#Since the data has been standardized, the intercept is set to F.
lasso_model <- cv.glmnet(as.matrix(X_quad_inter), Y, family="gaussian", intercept = F)
beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
x_sec <- as.matrix(X_quad_inter[,index_beta_non0] )                                         
sec_model = glm(Y~x_sec-1,family = gaussian)                           # second estimation
pred = predict(sec_model, newx = x_sec, type="response")
U <- Y-pred
xbeta <- x_sec%*%sec_model$coefficients   
plot_range_x <- range(xbeta)  
plot_range_y <- range(U)      
y_ticks <- seq(from = -2, to = 4, by = 1)
pdf("nonlinear-residual.pdf", width = 14, height = 14)   # auto-save figure (14 x 14 in = 1008 pt)
plot(xbeta, U,
     xlab = "",
     ylab = "Residuals",
     pch = 16,
     axes = FALSE,       
     xlim = plot_range_x,
     ylim = plot_range_y,
     cex.lab = 2.5,
     cex.axis = 2,
     mgp = c(2, 0.8, 0)
)
axis(1, cex.axis = 2, mgp = c(2, 0.8, 0))
axis(2, at = y_ticks, cex.axis = 2, mgp = c(2, 0.8, 0))
box()
mtext(expression(hat(Y)), side = 1,
      line = 3,
      cex = 2.5)
dev.off()








