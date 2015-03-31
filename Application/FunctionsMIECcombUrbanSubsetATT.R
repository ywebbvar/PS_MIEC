IPTWest <- function(txmodel, data, outcomemodel){
    data<-data
    model <- glm(formula=txmodel, data=data, family = "binomial")
    data$ps<- predict(model, data=data, type="response")
    data$iptw<-ifelse(data$tertscore==1, 1, data$ps/(1-data$ps))
    outmodel<-lm(formula=outcomemodel, weights=data$iptw, data=data)

    return(c(summary(outmodel)$coef[2], sqrt(diag(vcovHC(outmodel))[2]), summary(outmodel)$coef[2] - (1.96*sqrt(diag(vcovHC(outmodel))[2])), summary(outmodel)$coef[2]+ (1.96*sqrt(diag(vcovHC(outmodel))[2])) ))

}

multiple_imputation_EC   <- function(data_main, data_calib, outcomevector){
  MIEC_data <- MIEC(data_main,data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
  imputed_cols <- 14:ncol(MIEC_data)
                         
  gamma_hat <- var_gamma_hat <- c()
                         
  Y <- outcomevector
  A <- MIEC_data[,"tertscore"]
  Z1 <- MIEC_data[,"SEXF"]
  Z2 <- MIEC_data[,"age_cent"]
  Z3 <- MIEC_data[,"moth"]
  Z4 <- MIEC_data[,"fath"]
  Z5 <- MIEC_data[,"black"]
  Z6 <- MIEC_data[,"latino"]
  Z7 <- MIEC_data[,"other"]
  Z8 <- MIEC_data[,"midwest"]
  Z9 <- MIEC_data[,"south"]
  Z10 <- MIEC_data[,"west"]
  Z11 <- MIEC_data[,"cinc"]
                      
  for(k in imputed_cols){
    X <- MIEC_data[,k]
                         #  X <- MIEC_data[,14]
    txmodel <- glm(A ~ poly(X,2) + Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 + Z11 , family = "binomial")
                        
    ps<- predict(txmodel, type="response")
    iptw<-ifelse(A==1, 1, ps/(1-ps))
    outmodel<-lm(Y ~ A, weights=iptw)
    gamma_hat[k-13]<-summary(outmodel)$coef[2]
    var_gamma_hat[k-13] <- diag(vcovHC(outmodel))[2]
                            }
                         
combine_foo <- function(coef, vars){
  gamma_hat_MI <- mean(coef)
  gamma_matrix <- matrix(coef, ncol=3, byrow=T)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
                             
  W <- sum(sapply(1:m, function(x){
    sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)}))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
  U <- mean(vars)
  
  # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
  
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
        ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2)))
  if(T_MI < 0) df <- m-1
                             
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
                             
    return(c(coef=gamma_hat_MI,se=sqrt(T_MI), CI_low=CI_low, CI_upp=CI_upp))
                        }
                         
    return(t(combine_foo(gamma_hat, var_gamma_hat)))
                         }

multiple_imputation_EC_congenial   <- function(data_main, data_calib, outcomevector){
  MIEC_data <- MIEC(data_main,data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
  imputed_cols <- 14:ncol(MIEC_data)
                         
  gamma_hat <- var_gamma_hat <- c()
                         
  Y <- outcomevector
  A <- MIEC_data[,"tertscore"]
  Z1 <- MIEC_data[,"SEXF"]
  Z2 <- MIEC_data[,"age_cent"]
  Z3 <- MIEC_data[,"moth"]
  Z4 <- MIEC_data[,"fath"]
  Z5 <- MIEC_data[,"black"]
  Z6 <- MIEC_data[,"latino"]
  Z7 <- MIEC_data[,"other"]
  Z8 <- MIEC_data[,"midwest"]
  Z9 <- MIEC_data[,"south"]
  Z10 <- MIEC_data[,"west"]
  Z11 <- MIEC_data[,"cinc"]
                      
  for(k in imputed_cols){
    X <- MIEC_data[,k]
                         #  X <- MIEC_data[,14]
    txmodel <- glm(A ~ X + Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 + Z11 , family = "binomial")
                        
    ps<- predict(txmodel, type="response")
    iptw<-ifelse(A==1, 1, ps/(1-ps))
    outmodel<-lm(Y ~ A, weights=iptw)
    gamma_hat[k-13]<-summary(outmodel)$coef[2]
    var_gamma_hat[k-13] <- diag(vcovHC(outmodel))[2]
                            }
                         
combine_foo <- function(coef, vars){
  gamma_hat_MI <- mean(coef)
  gamma_matrix <- matrix(coef, ncol=3, byrow=T)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
                             
  W <- sum(sapply(1:m, function(x){
    sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)}))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
  U <- mean(vars)
  
  # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
  
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
        ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2)))
  if(T_MI < 0) df <- m-1
                             
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
                             
    return(c(coef=gamma_hat_MI,se=sqrt(T_MI), CI_low=CI_low, CI_upp=CI_upp))
                        }
                         
    return(t(combine_foo(gamma_hat, var_gamma_hat)))
                         }


multiple_imputation_EC_noY   <- function(data_main, data_calib, outcomevector){
  MIEC_data <- MIEC(data_main,data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
  imputed_cols <- 13:ncol(MIEC_data)
                         
  gamma_hat <- var_gamma_hat <- c()
                         
  Y <- outcomevector
  A <- MIEC_data[,"tertscore"]
  Z1 <- MIEC_data[,"SEXF"]
  Z2 <- MIEC_data[,"age_cent"]
  Z3 <- MIEC_data[,"moth"]
  Z4 <- MIEC_data[,"fath"]
  Z5 <- MIEC_data[,"black"]
  Z6 <- MIEC_data[,"latino"]
  Z7 <- MIEC_data[,"other"]
  Z8 <- MIEC_data[,"midwest"]
  Z9 <- MIEC_data[,"south"]
  Z10 <- MIEC_data[,"west"]
  Z11 <- MIEC_data[,"cinc"]          

  for(k in imputed_cols){
    X <- MIEC_data[,k]
                      
    txmodel <- glm(A ~ X + Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 + Z11 , family = "binomial")
    ps<- predict(txmodel, type="response")
    iptw<-ifelse(A==1, 1, ps/(1-ps))
    outmodel<-lm(Y ~ A,  weights=iptw)
    gamma_hat[k-12]<-summary(outmodel)$coef[2]
    var_gamma_hat[k-12] <- diag(vcovHC(outmodel))[2]
                            }
                         
combine_foo <- function(coef, vars){
# coefficient of interest (gamma_x_hat)
  gamma_hat_MI <- mean(coef)
                             
  gamma_matrix <- matrix(coef, ncol=3, byrow=T)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
                             
  W <- sum(sapply(1:m, function(x){
       sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)
          }))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
                             
  U <- mean(vars)
                             
 # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
                             
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
        ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2)))
  if(T_MI < 0) df <- m-1
                             
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
                             
    return(c(coef=gamma_hat_MI,se=sqrt(T_MI), CI_low=CI_low, CI_upp=CI_upp))
                        }
                         
    return(t(combine_foo(gamma_hat, var_gamma_hat)))
                         }

multiple_imputation_EC_noTY   <- function(data_main, data_calib, treatmentvector, outcomevector){
  MIEC_data <- MIEC(data_main,data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
  imputed_cols <- 12:ncol(MIEC_data)
                         
  gamma_hat <- var_gamma_hat <- c()
                         
  Y <- outcomevector
  A <- treatmentvector
  Z1 <- MIEC_data[,"SEXF"]
  Z2 <- MIEC_data[,"age_cent"]
  Z3 <- MIEC_data[,"moth"]
  Z4 <- MIEC_data[,"fath"]
  Z5 <- MIEC_data[,"black"]
  Z6 <- MIEC_data[,"latino"]
  Z7 <- MIEC_data[,"other"]
  Z8 <- MIEC_data[,"midwest"]
  Z9 <- MIEC_data[,"south"]
  Z10 <- MIEC_data[,"west"]
  Z11 <- MIEC_data[,"cinc"]          

  for(k in imputed_cols){
    X <- MIEC_data[,k]
                      
    txmodel <- glm(A ~ X + Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 + Z9 + Z10 + Z11 , family = "binomial")
    ps<- predict(txmodel, type="response")
    iptw<-ifelse(A==1, 1, ps/(1-ps))
    outmodel<-lm(Y ~ A,  weights=iptw)
    gamma_hat[k-11]<-summary(outmodel)$coef[2]
    var_gamma_hat[k-11] <- diag(vcovHC(outmodel))[2]
                            }
                         
combine_foo <- function(coef, vars){
# coefficient of interest (gamma_x_hat)
  gamma_hat_MI <- mean(coef)
                             
  gamma_matrix <- matrix(coef, ncol=3, byrow=T)
  mean_n_gamma_hat <- apply(gamma_matrix, 1, mean)
                             
  W <- sum(sapply(1:m, function(x){
       sum((gamma_matrix[x,] - mean_n_gamma_hat[x])^2)
          }))/(m*(n-1))
  B <- sum((mean_n_gamma_hat - gamma_hat_MI)^2)/(m-1)
                             
  U <- mean(vars)
                             
 # Calculating variance of our coefficient of interest
  T_MI <- U - W + (1+1/m)*B - W/n
  if(T_MI < 0) T_MI <- (1+1/m)*B
                             
  # Confidence intervals are done with t-distribution with these
  # degrees of freedom
  df <- 1/(((((1+1/m)*B)^2)/((m-1)*T_MI^2)) + 
        ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2)))
  if(T_MI < 0) df <- m-1
                             
  CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
  CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
                             
    return(c(coef=gamma_hat_MI,se=sqrt(T_MI), CI_low=CI_low, CI_upp=CI_upp))
                        }
                         
    return(t(combine_foo(gamma_hat, var_gamma_hat)))
                         }
