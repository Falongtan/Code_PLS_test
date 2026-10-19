# PLS_logistic
###Check if the logistic model is properly formulated
###use min and fisher combination method to combine the p-values on two data, respectively
###use Cauchy combination method to combine the p-values from two data
###compare with grp test
#library packages
library(MASS)
library(glmnet)
library(ncvreg)
library(LassoSIR) 
library(VariableScreening)
library(psych)
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(doRNG)   
library(GRPtests)



####define function
####parameter explanation
#### p : dimension   
#### n : sample size 
#### a : type of hypothesis, 0: the null hypothesis, 1:the alternative hypothesis
#### example: type of linear model, 1:H11  2:H12  3:H13  4:H14
#### c_h : bandwidth parameteer
#### pho : correlation coefficient between X
#### rpgrp_switch : the switch to run RP and GRP.  TRUE:run RP and GRP   FALSE: don't run rp and grp
#### estimate_type: type of estimation , estimate_type = 1: Lasso ; 2: SCAD ; others:MCP
#### estimate_num: number of estimation,   estimate_num = 1: once estimation; others: twice estimation
PLS_logit_parallel_1 <- function(p,n,a,c_h,pho,example,rpgrp_switch,estimate_type,estimate_num,seed){
  
  
  s <- 1000
  mu <- rep(0, p)
  if(pho == 0){
    sigma <- diag(rep(1,p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }   
  beta0 <- c(rep(1,5),rep(0,p-5))
  
  #parallel
  # tic()
  if (detectCores(logical=F) >= 16){
    cores <- 14
  }else if(detectCores(logical=F)<=6){
    cores <- detectCores(logical=F)
  }else{
    cores <- detectCores(logical=F)-6
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  set.seed(seed)                     # seed the master process (reviewer comment #1)
  registerDoRNG(seed)                # reproducible parallel RNG streams per iteration
  PLS_power <- foreach(k = 1:s, .combine='rbind', .errorhandling = "remove",
                       .packages = c('MASS','glmnet','ncvreg','LassoSIR','VariableScreening','psych',
                                     'foreach','parallel','iterators','doParallel','GRPtests')) %dopar%  
    {
      #generate x and y 
      x <- mvrnorm(n, mu, sigma)
      if(example==1){
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                              # H21
      }else{
        z <- x %*% beta0 + a *(x[,1]*x[,2] + x[,2]*x[,3] + x[,3]*x[,4] + x[,4]*x[,5])    # H22
      }
      
      pr <- 1/(1 + exp(-z)) 
      y <- matrix(rbinom(n, 1, pr),ncol = 1)
      
      ###grp
      if(rpgrp_switch && c_h==1){
        pval_grp <- GRPtest(x, y, fam = "binomial", nsplits = 1,RP_function = NULL,output_all=FALSE,penalize=TRUE)
      }else{
        pval_grp <- NA
      }
      
      
      #data splitting  
      index_x<-sample(1:n,floor(n/2),replace = FALSE)  
      index_x<-sort(index_x)
      #x
      x1 <- x[-index_x,] 
      x2 <- x[index_x,]
      #y
      y1 <- as.matrix(y[-index_x,])
      y2 <- as.matrix(y[index_x,]) 
      #sample size
      n1 <- nrow(x1)
      n2 <- nrow(x2)
      
      #residuals and projections based on x1 y1 and dimension reduction
      if(estimate_type == 1){                                                        #Lasso
        lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
        beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
        }else{
          if(estimate_num == 1){                                                     # once estimation
            pred1 =predict(lasso_model1, newx = x1, s = "lambda.min", type = "response")
          }else{
            x1_sec <- x1[,index_beta1_non0]                                          # twice estimation
            sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
            pred1 = predict(sec_model1, newx = x1_sec, type="response")
            pred1 = matrix(unname(pred1), ncol = 1) 
          }
        }
      }else if(estimate_type == 2){                                                                        #scad
        scad_model1 <- cv.ncvreg(x1, y1, penalty = "SCAD", family = "binomial")
        scad_beta1 <- coef(scad_model1, s = "lambda.min")[-1]  
        index_beta1_non0 <- seq(1:p)[which(as.numeric(scad_beta1) != 0)]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
        } else {
          if(estimate_num == 1){                                                     # once estimation
            pred1 <- predict(object = scad_model1,X = x1, type = "response",lambda = scad_model1$lambda.min )
          }else{
            x1_sec <- x1[,index_beta1_non0]                                          # twice estimation
            sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
            pred1 = predict(sec_model1, newx = x1_sec, type="response")
            pred1 = matrix(unname(pred1), ncol = 1) 
          }
        }
      }else{
        mcp_model1 <- cv.ncvreg(x1, y1, penalty = "MCP", family = "binomial")
        mcp_beta1 <- coef(mcp_model1, s = "lambda.min")[-1]  
        index_beta1_non0 <- seq(1:p)[which(as.numeric(mcp_beta1) != 0)]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
        } else {
          if(estimate_num == 1){                                                     # once estimation
            pred1 <- predict(object = mcp_model1,X = x1, type = "response",lambda = mcp_model1$lambda.min )
          }else{
            x1_sec <- x1[,index_beta1_non0]                                          # twice estimation
            sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
            pred1 = predict(sec_model1, newx = x1_sec, type="response")
            pred1 = matrix(unname(pred1), ncol = 1) 
          }
        }
      }
      U1 <- y1 - pred1                                                         # residuals based on x1 y1
      
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
      if(estimate_type == 1){                                                        #Lasso
        lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
        beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso) != 0]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
        }else{
          if(estimate_num == 1){                                                     # once estimation
            pred2 =predict(lasso_model2, newx = x2, s = "lambda.min", type = "response")
          }else{
            x2_sec <- x2[,index_beta2_non0]                                          # twice estimation
            sec_model2 <- glm(y2~x2_sec, family = binomial(link = "logit"))
            pred2 = predict(sec_model2, newx = x2_sec, type="response")
            pred2 = matrix(unname(pred2), ncol = 1) 
          }
        }
      }else if(estimate_type==2){                                                                        #scad
        scad_model2 <- cv.ncvreg(x2, y2, penalty = "SCAD", family = "binomial")
        scad_beta2 <- coef(scad_model2, s = "lambda.min")[-1]  
        index_beta2_non0 <- seq(1:p)[which(as.numeric(scad_beta2) != 0)]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
        } else {
          if(estimate_num == 1){                                                     # once estimation
            pred2 <- predict(object = scad_model2,X = x2, type = "response",lambda = scad_model2$lambda.min )
          }else{
            x2_sec <- x2[,index_beta2_non0]                                          # twice estimation
            sec_model2 <- glm(y2~x2_sec, family = binomial(link = "logit"))
            pred2 = predict(sec_model2, newx = x2_sec, type="response")
            pred2 = matrix(unname(pred2), ncol = 1) 
          }
        }
      }else{
        mcp_model2 <- cv.ncvreg(x2, y2, penalty = "MCP", family = "binomial")
        mcp_beta2 <- coef(mcp_model2, s = "lambda.min")[-1]  
        index_beta2_non0 <- seq(1:p)[which(as.numeric(mcp_beta2) != 0)]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
        } else {
          if(estimate_num == 1){                                                     # once estimation
            pred2 <- predict(object = mcp_model2,X = x2, type = "response",lambda = mcp_model2$lambda.min )
          }else{
            x2_sec <- x2[,index_beta2_non0]                                          # twice estimation
            sec_model2 <- glm(y2~x2_sec, family = binomial(link = "logit"))
            pred2 = predict(sec_model2, newx = x2_sec, type="response")
            pred2 = matrix(unname(pred2), ncol = 1) 
          }
        }
      }
      U2 <- y2 - pred2                                                         # residuals based on x2 y2
      
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
      h1 <- c_h*(n1^(-2/9))                                                         # bandwidth
      
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
      h2 <- c_h*(n2^(-2/9))                                                         # bandwidth 
      
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
      
      #result
      result <- c(pval_fisher1,pval_fisher2,pval_fisher,
                  pval_min1,pval_min2,pval_min,pval_grp)
      
    }
  #end parallel
  stopCluster(cl)
  
  
  # toc()
  #power
  PLS_pval_fisher1 <- PLS_power[,1]
  PLS_pval_fisher2 <- PLS_power[,2]
  PLS_pval_fisher  <- PLS_power[,3]
  PLS_pval_min1  <- PLS_power[,4]
  PLS_pval_min2  <- PLS_power[,5]
  PLS_pval_min  <- PLS_power[,6]
  #grp
  PLS_pval_grp <- PLS_power[,7]
  
  PLS_fisher1_power <- mean(PLS_pval_fisher1<=0.05,na.rm=TRUE)
  PLS_fisher2_power <- mean(PLS_pval_fisher2<=0.05,na.rm=TRUE)
  PLS_fisher_power <- mean(PLS_pval_fisher<=0.05,na.rm=TRUE)
  PLS_min1_power <- mean(PLS_pval_min1,na.rm=TRUE)
  PLS_min2_power <- mean(PLS_pval_min2,na.rm=TRUE)
  PLS_min_power <- mean(PLS_pval_min<=0.05,na.rm=TRUE)
  PLS_grp_power <- mean(PLS_pval_grp<=0.05,na.rm=TRUE)
  
  #return result
  return(list(PLS_fisher1_power=PLS_fisher1_power,
              PLS_fisher2_power=PLS_fisher2_power,
              PLS_fisher_power=PLS_fisher_power,
              PLS_min1_power=PLS_min1_power,
              PLS_min2_power=PLS_min2_power,
              PLS_min_power=PLS_min_power,
              PLS_grp_power=PLS_grp_power))
}





#####analysis of Empirical sizes and powers of our tests with different bandwidths in logistic models
example <- c(1,2) 
a <- c(0,1) #null and alternative
n <- 600                      #sample size
p<-  c(50,100,300,600, 900,1200)   #dimension
c_h <-  c(0.6,0.8, 1,1.2, 1.4)    #bandwidths
pho <- c(0,0.4,0.8)           #correlation
estimate_type <- c(1,2,3)
rpgrp_run <- FALSE           ###don't run  grp
estimate_num=2


r <- 7
compar_logit_result <- matrix(0,nrow = r*length(estimate_type)*length(example)*length(pho)*length(a)*length(c_h),ncol = length(p)+6)
for(nn in 1:length(estimate_type)){
  for(ll in 1:length(example)){
    for (kk in 1:length(pho)){
      for (ii in 1:length(a)){
        for (jj in 1:length(p)){
          for (cc in 1:length(c_h)){
            
            t3 <- Sys.time()
            seed <- (nn-1)*length(example)*length(pho)*length(a)*length(p)*length(c_h) +
                    (ll-1)*length(pho)*length(a)*length(p)*length(c_h) +
                    (kk-1)*length(a)*length(p)*length(c_h) +
                    (ii-1)*length(p)*length(c_h) +
                    (jj-1)*length(c_h) + cc    # unique seed per parameter combination
            result_pls <- PLS_logit_parallel_1(p[jj], n, a[ii],c_h[cc],pho[kk],example[ll],rpgrp_switch=rpgrp_run,estimate_type[nn],estimate_num,seed=seed)
            #index
            base_index <- r * ((nn-1)*length(example)*length(pho)*length(a)*length(c_h) + (ll-1)*length(pho)*length(a)*length(c_h) 
                               + (kk-1)*length(a)*length(c_h) + (ii-1)*length(c_h) + (cc-1) ) 
            #result
            compar_logit_result[base_index + 1, jj]  <- result_pls$PLS_fisher1_power
            compar_logit_result[base_index + 2, jj]  <- result_pls$PLS_fisher2_power
            compar_logit_result[base_index + 3, jj]  <- result_pls$PLS_fisher_power
            compar_logit_result[base_index + 4, jj]  <- result_pls$PLS_min1_power
            compar_logit_result[base_index + 5, jj]  <- result_pls$PLS_min2_power
            compar_logit_result[base_index + 6, jj]  <- result_pls$PLS_min_power
            compar_logit_result[base_index + 7, jj]  <- result_pls$PLS_grp_power
            
            t4 <- Sys.time()
            cat("example = ",example[ll],"p = ", p[jj],", n = ", n, ", a = ", a[ii] ,", c_h = ", c_h[cc], ", pho = ", pho[kk],
                ",estimate_type=",estimate_type[nn],
                ", PLS_fisher1_power = ", result_pls$PLS_fisher1_power,
                ", PLS_fisher2_power= ", result_pls$PLS_fisher2_power,
                ", PLS_fisher_power = ", result_pls$PLS_fisher_power,
                ", PLS_min1_power = ", result_pls$PLS_min1_power, 
                ", PLS_min2_power = ", result_pls$PLS_min2_power,
                ", PLS_min_power= ", result_pls$PLS_min_power, 
                ", grp_power = ",  result_pls$PLS_grp_power, 
                ", time:",t4-t3, "\n")
            
            ##the label with respect to result
            compar_logit_result[base_index + 1:r, (length(p)+1):(length(p)+6)] <- 
              matrix(c(example[ll],n, a[ii], pho[kk], c_h[cc],estimate_type[nn]), nrow = r, ncol = 6, byrow = TRUE)
          }
        }
      }
    }
  }
}

outcome.name <- paste("PLS_logistic_H2_result", ".csv", sep="")
write.csv(compar_logit_result, file = outcome.name, row.names = FALSE)



#compare with grp
example <- c(1,2) 
a <- c(0,1) #null and alternative
n <- 600                      #sample size
p<-  c(50,100,300,600, 900,1200)   #dimension
c_h <-  c(1)                          #bandwidths
pho <- c(0,0.4,0.8)           #correlation
estimate_type <- c(1,2,3)
rpgrp_run <- TRUE           ### run grp
estimate_num=2


r <- 7
compar_logit_result <- matrix(0,nrow = r*length(estimate_type)*length(example)*length(pho)*length(a)*length(c_h),ncol = length(p)+6)
for(nn in 1:length(estimate_type)){
  for(ll in 1:length(example)){
    for (kk in 1:length(pho)){
      for (ii in 1:length(a)){
        for (jj in 1:length(p)){
          for (cc in 1:length(c_h)){
            
            t3 <- Sys.time()
            seed <- (nn-1)*length(example)*length(pho)*length(a)*length(p)*length(c_h) +
                    (ll-1)*length(pho)*length(a)*length(p)*length(c_h) +
                    (kk-1)*length(a)*length(p)*length(c_h) +
                    (ii-1)*length(p)*length(c_h) +
                    (jj-1)*length(c_h) + cc    # unique seed per parameter combination
            result_pls <- PLS_logit_parallel_1(p[jj], n, a[ii],c_h[cc],pho[kk],example[ll],rpgrp_switch=rpgrp_run,estimate_type[nn],estimate_num,seed=seed)
            #index
            base_index <- r * ((nn-1)*length(example)*length(pho)*length(a)*length(c_h) + (ll-1)*length(pho)*length(a)*length(c_h) 
                               + (kk-1)*length(a)*length(c_h) + (ii-1)*length(c_h) + (cc-1) ) 
            #result
            compar_logit_result[base_index + 1, jj]  <- result_pls$PLS_fisher1_power
            compar_logit_result[base_index + 2, jj]  <- result_pls$PLS_fisher2_power
            compar_logit_result[base_index + 3, jj]  <- result_pls$PLS_fisher_power
            compar_logit_result[base_index + 4, jj]  <- result_pls$PLS_min1_power
            compar_logit_result[base_index + 5, jj]  <- result_pls$PLS_min2_power
            compar_logit_result[base_index + 6, jj]  <- result_pls$PLS_min_power
            compar_logit_result[base_index + 7, jj]  <- result_pls$PLS_grp_power
            
            t4 <- Sys.time()
            cat("example = ",example[ll],"p = ", p[jj],", n = ", n, ", a = ", a[ii] ,", c_h = ", c_h[cc], ", pho = ", pho[kk],
                ",estimate_type=",estimate_type[nn],
                ", PLS_fisher1_power = ", result_pls$PLS_fisher1_power,
                ", PLS_fisher2_power= ", result_pls$PLS_fisher2_power,
                ", PLS_fisher_power = ", result_pls$PLS_fisher_power,
                ", PLS_min1_power = ", result_pls$PLS_min1_power, 
                ", PLS_min2_power = ", result_pls$PLS_min2_power,
                ", PLS_min_power= ", result_pls$PLS_min_power, 
                ", grp_power = ",  result_pls$PLS_grp_power, 
                ", time:",t4-t3, "\n")
            
            ##the label with respect to result
            compar_logit_result[base_index + 1:r, (length(p)+1):(length(p)+6)] <- 
              matrix(c(example[ll],n, a[ii], pho[kk], c_h[cc],estimate_type[nn]), nrow = r, ncol = 6, byrow = TRUE)
          }
        }
      }
    }
  }
}

outcome.name <- paste("PLS_logistic_H2_grp_result", ".csv", sep="")
write.csv(compar_logit_result, file = outcome.name, row.names = FALSE)



