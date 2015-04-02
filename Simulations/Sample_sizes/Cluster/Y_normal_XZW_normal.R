################################ Simulation of propensity score with covariates
##measured with error
##
##  This simulation includes the following changes: 
##  - Changes size of calibration sample
##  method. 
##
##  It also does:
##   - Implements simple imputation method
##   - Removes reliabilities < 0.5
##   - Uses m=12, n=3  (as all the previous simulations)
##   - Runs reliabilities 0.9 and 0.999
##   - Modified the 'sampling' function to also run an 'extra small' variance (with reliability = 0.99)
##   - Eliminates CC, RP. Only calculates naive. 
##   - Uses 4 types of MIEC
##   - Sets reliabilities instead of variances
##   - Uses both a correct model (uses X and Z in the 
##      estimation model) as well as IPTW

##  We use this code to evaluate a simple imputation procedure

cluster=TRUE

library(mvtnorm)
library(mitools)
library(MCMCpack)
library(xtable)
library(survey)
library(ecodist)

if(cluster){source("/home/bst/student/ywebbvar/Github/PS_MIEC/MIEC_original/MI-EC_algorithm.r")
}else{source("~/GitHub/PS_MIEC/MIEC_original/MI-EC_algorithm.r")}


## ----expitlogit----------------------------------------------------------
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(p) log(p/(1-p))


## ----, eval=TRUE---------------------------------------------------------
if(cluster){
  temp <- commandArgs(TRUE)
  i    <-as.numeric(temp[1])
}else{
  i = 1
}


## ----def_params----------------------------------------------------------
# Sample size and number of simulations
N_main  <- c(2500,10000)
N_calib <- c(100,500,1000)
N_sim   <- ifelse(cluster, 10, 2)

# Distribution of T | X,Z
gamma_0  <- 0 #This will determine the proportion receiving treatment
gamma_Z  <- 0.4
gamma_Xs <- 0.4
gamma_Xl <- 1.2

# Distribution of (X,Z)
rho_l    <- 0.3
rho_m    <- 0.6
rho_h    <- 0.9


# Distribution of W | Y,X,Z
beta_0   <- 0
beta_1   <- 1

#The variance of X is obtain by doing

rels      <- c(0.5,0.7,0.9,0.999)
sigmas    <- ((1/rels) - beta_1^(2))
sigma2_xs <- sigmas[4]
sigma2_s  <- sigmas[3]
sigma2_m  <- sigmas[2]
sigma2_l  <- sigmas[1]

# Distribution of Y(T) | T, X, Z
Delta    <- 2
delta_X  <- 0.5
delta_Z  <- 0.1
tau2     <- 1


## ----sampling_fun--------------------------------------------------------
sampling <- function(cor_level, X_effect, m_error, n_main, n_calib){
  
  # Sampling ($X,Z$) from a multivariate distribution, according to a set level of correlation.
  if(cor_level == "low"){
    sigma_X_Z <- matrix(c(1,rho_l,rho_l,1), ncol=2)
  }else if(cor_level == "med") {
    sigma_X_Z = matrix(c(1,rho_m,rho_m,1), ncol=2)
  }else if(cor_level == "high") {
    sigma_X_Z = matrix(c(1,rho_h,rho_h,1), ncol=2)
  }else {stop("Correlation level must be 'high', 'med', or 'low'")}
  
  X_Z <- rmvnorm(n_main+n_calib, mean = c(0,0), sigma = sigma_X_Z)
  colnames(X_Z) <- c("X", "Z")
  
  
  # Sample Y(0), Y(1) \mid X,Z
  
  Y0 <- rnorm(n_main+n_calib, mean=(delta_X*X_Z[,"X"] + delta_Z*X_Z[,"Z"]), 
              sqrt(tau2))
  Y1 <- rnorm(n_main+n_calib, mean=(Delta + delta_X*X_Z[,"X"] + delta_Z*X_Z[,"Z"]), sqrt(tau2))
  
  #Y1 <- Y0 + Delta     # rank preserving
  
  
  # Sample from the distribution of $T \mid X,Z,\psi$
  
  if(X_effect == "small"){
    logit_T <- gamma_0 + gamma_Xs*X_Z[,"X"] + gamma_Z*X_Z[,"Z"]
  } else if(X_effect == "large") {
    logit_T <- gamma_0 + gamma_Xl*X_Z[,"X"] + gamma_Z*X_Z[,"Z"]
  } else {stop("Effect of X must be 'small' or 'large'")}
  
  T <- rbinom(n_main+n_calib, 1, p= expit(logit_T))
  
  
  # Creating Y_obs given Y(0), Y(1), T
  
  Y_obs <- Y1*T + Y0*(1-T)
  
  
  #Sample from the distribution of $W \mid T,X,Z$
  
  if(m_error == "small"){
    W <- rnorm(n_main+n_calib, mean= (beta_0+beta_1*X_Z[,"X"]), sd=sqrt(sigma2_s))
  } else if(m_error == "moderate") {
    W <- rnorm(n_main+n_calib, mean= (beta_0+beta_1*X_Z[,"X"]), sd=sqrt(sigma2_m))
  } else if(m_error == "large") {
    W <- rnorm(n_main+n_calib, mean= (beta_0+beta_1*X_Z[,"X"]), sd=sqrt(sigma2_l))
  } else if(m_error == "extra small") {
    W <- rnorm(n_main+n_calib, mean= (beta_0+beta_1*X_Z[,"X"]), sd=sqrt(sigma2_xs))
  } else {stop("Measurement error must be 'extra small', 'small', 'moderate', or 'large'")}
  
  
  #Setting datasets for calibration and for main inference
  
  i.calib <- 1:n_calib
  i.main  <- (n_calib+1):(n_calib+n_main)
  
  data.main  <- data.frame(W = W[i.main], T = T[i.main], Z=X_Z[i.main,"Z"], 
                           Y_obs = Y_obs[i.main])
  data.cause <- data.frame(p= expit(logit_T[i.main]), X= X_Z[i.main,"X"], Y0 = Y0[i.main], Y1 = Y1[i.main])
  data.calib <- data.frame(X= X_Z[i.calib,"X"], W = W[i.calib])
  
  return(list(main=data.main, calib=data.calib, cause=data.cause))
}




## ----Correction_functions------------------------------------------------
true_regression <- function(data_main, Xtrue){
  data_main$Xtrue <- Xtrue
  
  model <- glm(T~Xtrue+Z, data=data_main, family=binomial)
  p_hat <- predict(model,type="response")
  
  data_main$wt <- ifelse(data_main$T==1, 1/p_hat, 1/(1-p_hat))
  
  design.ate   <- svydesign(ids=~1, weights=~wt, data=data_main)
  survey.model <- svyglm(Y_obs~T, design=design.ate)
  coef         <- summary(survey.model)$coeff["T",c("Estimate", "Std. Error")]
  
  CI_low <- coef["Estimate"] + qnorm(0.025)*coef["Std. Error"]
  CI_upp <- coef["Estimate"] + qnorm(0.975)*coef["Std. Error"]
  
  coef <- cbind(coef["Estimate"],coef["Std. Error"],CI_low, CI_upp)
  colnames(coef) <- c("results", "se", "(lower", "upper)")
  
  return(coef)
}


naive_regression <- function(data_main){
  model <- glm(T ~ W + Z, data=data_main, family=binomial)
  p_hat <- predict(model,type="response")
  
  data_main$wt <- ifelse(data_main$T==1, 1/p_hat, 1/(1-p_hat))
  
  design.ate   <- svydesign(ids=~1, weights=~wt, data=data_main)
  survey.model <- svyglm(Y_obs~T, design=design.ate)
  coef         <- summary(survey.model)$coeff["T",c("Estimate", "Std. Error")]
  
  CI_low <- coef["Estimate"] + qnorm(0.025)*coef["Std. Error"]
  CI_upp <- coef["Estimate"] + qnorm(0.975)*coef["Std. Error"]
  
  coef <- cbind(coef["Estimate"],coef["Std. Error"],CI_low, CI_upp)
  colnames(coef) <- c("results", "se", "(lower", "upper)")
  
  return(coef)
}


simple_imputation <- function(data_main, data_calib, N_data=12*3){
  #1) estimate the variance of the error sigma^2 in the calibration sample
  
  sd_XW <- sd(lm(W ~ X-1, data=data_calib)$res)
  
  #2) generate an error variable in a multiple imputation setting following 
  #   N(0,sigma^2) to be subtracted (or added) to the observed variable X, over 
  #   N datasets (following their approach to combine across datasets).
  
  errors <- sapply(1:N_data, function(n) rnorm(nrow(data_main), mean=0, sd=sd_XW))
  X_imp  <- data_main$W - errors
  
  delta_MI <- matrix(NA,ncol=ncol(X_imp),nrow=2)
  
  MI_data  <- data.frame(
    T = data_main[,"T"],
    Z = data_main[,"Z"],
    Y_obs = data_main[,"Y_obs"])
  
  for(k in 1:ncol(X_imp)){
    MI_data$X <- X_imp[,k]
    model <- glm(T ~ X + Z, data=MI_data, family=binomial)
    p_hat <- predict(model,type="response")
    
    MI_data$wt <- ifelse(data_main$T==1, 1/p_hat, 1/(1-p_hat))
    
    design.ate   <- svydesign(ids=~1, weights=~wt, data=MI_data)
    survey.model <- svyglm(Y_obs~T, design=design.ate)
    delta_MI[,k] <- summary(survey.model)$coeff["T",c("Estimate", "Std. Error")]
  }
  
  
  estimate.mitools <- summary(MIcombine(results = as.list(delta_MI[1,]), 
                                        variances = as.list(delta_MI[2,]^2)))[,-5]
  return(estimate.mitools)   
}

## ----MIEC_funs-----------------------------------------------------------
multiple_imputation_EC   <- function(data_main, data_calib, option="Yout",n_main, n_calib){
  #Option: select 'Ycov', 'noT', 'noY', 'noT'
  # Yout (includes Y) uses $T$ as outcome, $Z$ and $Y_{obs}$ as helpful covariates
  # noY (no Y) uses $T$ as outcome, $Z$ as helpful covariate
  # noT (no T) uses $Y_{obs}$ as outcome, $Z$ as helpful covariate
  # noTY (no T, nor Y) uses no outcome, $Z$ as helpful covariate
  
  # Other parameters for MIEC/ two-stage imputation procedure
  m <- 12 #Number of draws from parameter distribution
  n <- 3  #Number of samples (and imputations) for each m
  
  # Generating Multiple Imputations:
  
  if(option=="Yout"){ # Y as a helpful covariate
    q <- 2  #Dimension of T. T = (T,Y_obs)
    r <- 1  #Dimension of Z. Z = (Z)
    MIEC_data <- MIEC(data_main[,c("W","T","Y_obs","Z")],data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
    
  }else if(option=="noY"){
    q <- 1  #Dimension of T. T = T
    r <- 1  #Dimension of Z. Z = Z
    MIEC_data <- MIEC(data_main[,c("W","T","Z")],data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
    
  }else if(option=="noT"){
    q <- 1  #Dimension of T. T = Y
    r <- 1  #Dimension of Z. Z = Z
    MIEC_data <- MIEC(data_main[,c("W","Y_obs","Z")],data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
    
  }else if(option=="noTY"){
    q <- 0  #Dimension of T. T = NULL
    r <- 1  #Dimension of Z. Z = Z
    MIEC_data <- MIEC(data_main[,c("W","Z")],data_calib,n_calib,n_main,M=m,N=n,K=q,S=r)
    
  }else {stop("Only options 'Yout', 'noT', 'noY', 'noT' are accepted")}
  
  
  imputed_cols <- (q+r+1):ncol(MIEC_data)
  
  delta_MI <- matrix(NA,ncol=length(imputed_cols),nrow=2)
  
  MI_data  <- data.frame(
    T = data_main[,"T"],
    Z = data_main[,"Z"],
    Y_obs = data_main[,"Y_obs"])
  
  for(k in imputed_cols){
    MI_data$X <- MIEC_data[,k]
    model <- glm(T ~ X + Z, data=MI_data, family=binomial)
    p_hat <- predict(model,type="response")
    
    MI_data$wt <- ifelse(data_main$T==1, 1/p_hat, 1/(1-p_hat))
    
    design.ate   <- svydesign(ids=~1, weights=~wt, data=MI_data)
    survey.model <- svyglm(Y_obs~T, design=design.ate)
    delta_MI[,k-(q+r)] <- summary(survey.model)$coeff["T",c("Estimate", "Std. Error")]
  }
  
  
  estimate.mitools <- summary(MIcombine(results = as.list(delta_MI[1,]), 
                                        variances = as.list(delta_MI[2,]^2)))[,-5]
  
  
  combine_foo <- function(coef, vars){
    
    # coefficient of interest (gamma_x_hat)
    gamma_hat_MI <- mean(coef)
    
    # The following code gives non-sensical results
    #Calculating W,B,U
    gamma_matrix <- matrix(coef, ncol=n, byrow=TRUE)
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
               ((((1+1/n)*W)^2)/(m*(n-1)*T_MI^2))) #Note that here was the mistake
    if(T_MI < 0) df <- m-1
    
    CI_low <- gamma_hat_MI + qt(0.025, df=df)*sqrt(T_MI)
    CI_upp <- gamma_hat_MI + qt(0.975, df=df)*sqrt(T_MI)
    
    return(c(coef=gamma_hat_MI,se=sqrt(T_MI), CI_low=CI_low, CI_upp=CI_upp))
  }
  
  
  estimate.reiter  <- combine_foo(delta_MI[1,], delta_MI[2,]^2)
  
  return(rbind(estimate.mitools,estimate.reiter))
}




## ----simulation_fun------------------------------------------------------
full_simulation <- function(n_main=2500, n_calib=500,cor_level = "low", X_effect = "large", m_error = "large"){
  data <- sampling(cor_level = cor_level, X_effect = X_effect, m_error = m_error,n_main=n_main, n_calib=n_calib)
  
  Xtrue     <- true_regression(data$main, data$cause$X)
  naive     <- naive_regression(data$main)
  MIEC_Yout <- multiple_imputation_EC(data$main, data$calib, "Yout",n_main=n_main, n_calib=n_calib)
  MIEC_noY  <- multiple_imputation_EC(data$main, data$calib, "noY",n_main=n_main, n_calib=n_calib)
  MIEC_noT  <- multiple_imputation_EC(data$main, data$calib, "noT",n_main=n_main, n_calib=n_calib)
  MIEC_noTY <- multiple_imputation_EC(data$main, data$calib, "noTY",n_main=n_main, n_calib=n_calib)
  MI_simple <- simple_imputation(data$main, data$calib)
  
  result_table      <- rbind(Xtrue, naive,MI_simple,MIEC_Yout, MIEC_noY, MIEC_noT, MIEC_noTY)
  mean_insample     <- mean(data$cause$Y1 - data$cause$Y0) 
  result_string     <- c(mean_insample,t(result_table))
  
  names(result_string) <- c("insample mean", 
                            paste0("Xtrue.", colnames(Xtrue)),
                            paste0("naive.", colnames(Xtrue)),
                            paste0("MI_simple.", colnames(Xtrue)),
                            paste0("MIEC_Yout_mi.", colnames(Xtrue)),
                            paste0("MIEC_Yout_re.", colnames(Xtrue)),
                            paste0("MIEC_noY_mi.", colnames(Xtrue)),
                            paste0("MIEC_noY_re.", colnames(Xtrue)),
                            paste0("MIEC_noT_mi.", colnames(Xtrue)),
                            paste0("MIEC_noT_re.", colnames(Xtrue)),
                            paste0("MIEC_noTY_mi.", colnames(Xtrue)),
                            paste0("MIEC_noTY_re.", colnames(Xtrue))
  )
  
  return(result_string)
}


## ----simulations,cache=TRUE, eval=TRUE-----------------------------------
time1 <- Sys.time()
set.seed(954054+i) 

results <- list()

results[["n_main_m"]][["sigma2_xs"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "extra small"))

results[["n_main_m"]][["sigma2_xs"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "extra small"))

results[["n_main_m"]][["sigma2_xs"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "extra small"))


results[["n_main_l"]][["sigma2_xs"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "extra small"))

results[["n_main_l"]][["sigma2_xs"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "extra small"))

results[["n_main_l"]][["sigma2_xs"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "extra small"))




results[["n_main_m"]][["sigma2_s"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "small"))

results[["n_main_m"]][["sigma2_s"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "small"))

results[["n_main_m"]][["sigma2_s"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "small"))


results[["n_main_l"]][["sigma2_s"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "small"))

results[["n_main_l"]][["sigma2_s"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "small"))

results[["n_main_l"]][["sigma2_s"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "small"))




results[["n_main_m"]][["sigma2_m"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "moderate"))

results[["n_main_m"]][["sigma2_m"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "moderate"))

results[["n_main_m"]][["sigma2_m"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "moderate"))


results[["n_main_l"]][["sigma2_m"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "moderate"))

results[["n_main_l"]][["sigma2_m"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "moderate"))

results[["n_main_l"]][["sigma2_m"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "moderate"))




results[["n_main_m"]][["sigma2_l"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "large"))

results[["n_main_m"]][["sigma2_l"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "large"))

results[["n_main_m"]][["sigma2_l"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[1], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "large"))


results[["n_main_l"]][["sigma2_l"]][["n_calib_s"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[1], cor_level = "med", X_effect = "large", m_error = "large"))

results[["n_main_l"]][["sigma2_l"]][["n_calib_m"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[2], cor_level = "med", X_effect = "large", m_error = "large"))

results[["n_main_l"]][["sigma2_l"]][["n_calib_l"]] <- sapply(1:N_sim, function(x) full_simulation(n_main=N_main[2], n_calib=N_calib[3], cor_level = "med", X_effect = "large", m_error = "large"))
time2 <- Sys.time()


## ----save_results, eval=TRUE---------------------------------------------
time2-time1

if(cluster) save(results,file=paste("Simulation_sample_size_", i,"_", format(Sys.time(), "%Y%m%d-%H%M"),".Rdata", sep=""))
