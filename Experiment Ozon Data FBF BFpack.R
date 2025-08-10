
library(BayesVarSel)
library(mice)
library(BFpack)

#X <- rmvnorm(178,mean=rep(0,7))
#y <- rnorm(178)
X <- as.matrix(Ozone35[,2:8])
colnames(X) <- NULL
y <- Ozone35$y

# get inclusion probs using FBFs in BFpack based on the full data set
# step 1. fit a basic lm as input of BF() function
lm1 <- lm(y ~ 1 + X)
# names of coefficients that can be constrained in the models
get_estimates(lm1)

# Define variable names
vars <- names(get_estimates(lm1)$estimate[-1])

# Initialize vector to collect combinations
all_combinations <- c()

# Create character strings for all 2^K constrained hypotheses for variable selection
for (k in 1:length(vars)) {
  cmb <- combn(vars, k, simplify = FALSE)
  # For each combination, create a string like "X1=0 & X2=0"
  formatted <- sapply(cmb, function(x) paste0(x, "=0", collapse = " & "))
  all_combinations <- c(all_combinations, formatted)
}
# Collapse into a single string separated by "; "
all_hypotheses <- paste(all_combinations, collapse = "; ")
# get PMPs based on FBF
BF1 <- BF(lm1,hypothesis=all_hypotheses)

# Obtain inclusion probabilities for each predictor by summing the posterior model
# probabilities which includes that predictor
posterior_probs <- BF1$PHP_confirmatory  # PMPs
# Initialize a named vector to store inclusion probabilities
inclusion_probs <- setNames(numeric(length(vars)), vars)
# Loop through each variable
for (var in vars) {
  # Find models that do not constrain this variable to zero
  not_constrained <- !grepl(paste0(var, "=0"), names(posterior_probs))
  # Sum the posterior probabilities of these models
  inclusion_probs[var] <- sum(posterior_probs[not_constrained])
}
# View the result of full data set.
print(inclusion_probs)

# data
X <- as.matrix(Ozone35[,2:8])
colnames(X) <- NULL
y <- Ozone35$y

# create missing data
N <- 1000
data10 <- list()
data20 <- list()
data30 <- list()

for (rem.prop in c(0.1,0.20,0.30)){
  cat("rem.prop:", rem.prop, "\n")
  n.missing <- round(n * rem.prop)
  for (J in 1:N){
    set.seed(rem.prop * 100 + J)
    X.miss <- X
    for(j in 3:7){
      col.miss=rep(1,7)
      col.miss[j]=0
      aux=ampute(X, prop=rem.prop, patterns=col.miss, mech = "MAR")
      X.miss[,j]=aux$amp[,j]
    }
    if (rem.prop==0.1) data10[[J]]<- data.frame(cbind(y,X.miss)) #JM
    if (rem.prop==0.2) data20[[J]]<- data.frame(cbind(y,X.miss)) #JM
    if (rem.prop==0.3) data30[[J]]<- data.frame(cbind(y,X.miss)) #JM
  }
}

# Obtain expected values of the quantities of interest based on 500 imputations.
M <- 250
inclusion_data10 <- inclusion_data10_listwise <- matrix(NA,ncol=ncol(X),nrow=length(data10))
for(r in 1:length(data10)){
  mice_miss <- mice :: mice(data = data10[[r]], m = M,
                            meth = rep("norm",8),
                            diagnostics = F, printFlag = F)
  relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
    y <- mice::complete(mice_miss, m)[,1] # same as 'y'
    X <- as.matrix(mice::complete(mice_miss, m)[,-1])
    colnames(X) <- NULL
    lm1_miss <- lm(y ~ 1 + X)
    # get PMPs based on FBF
    BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
    #extract quantities of interest
    c(BF1_miss$BFtable_confirmatory[, 1:4])
  },mc.cores=5)),ncol = M)
  while(is.character(relmeas_all[1,1])){ # in the (rare) case the imputed data gives errors Then redo imputation.
    mice_miss <- mice :: mice(data = data10[[r]], m = M,
                              meth = rep("norm",8),
                              diagnostics = F, printFlag = F)
    relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
      y <- mice::complete(mice_miss, m)[,1] # same as 'y'
      X <- as.matrix(mice::complete(mice_miss, m)[,-1])
      colnames(X) <- NULL
      lm1_miss <- lm(y ~ 1 + X)
      # get PMPs based on FBF
      BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
      #extract quantities of interest
      c(BF1_miss$BFtable_confirmatory[, 1:4])
    },mc.cores=5)),ncol = M)
  }
  # get expected values of quantities of interest
  relmeas <- matrix(apply(relmeas_all, 1, mean),ncol = 4)
  #get FBFs
  FBFs <- relmeas[,3] / relmeas[,1] * relmeas[,4] / relmeas[,2]
  #get PMPs
  PMPs <- FBFs / sum(FBFs)
  names(PMPs) <- names(posterior_probs)
  #get inclusion probabilities (same as above)
  inclusion_probs_miss <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained <- !grepl(paste0(var, "=0"), names(PMPs))
    # Sum the posterior probabilities of these models
    inclusion_probs_miss[var] <- sum(PMPs[not_constrained])
  }
  inclusion_data10[r,] <- inclusion_probs_miss
  
  # inclusion probs based on FBF after list-wise deletion
  data10_listwise_r <- data10[[r]][!is.na(apply(data10[[r]],1,sum)),]
  y <- data10_listwise_r[,1]
  X <- as.matrix(data10_listwise_r[,-1])
  colnames(X) <- NULL
  lm1_listwise_r <- lm(y ~ 1 + X)
  BF1_listwise_r <- BF(lm1_listwise_r,hypothesis=all_hypotheses)
  # Obtain inclusion probabilities for each predictor by summing the posterior model
  # probabilities which includes that predictor
  posterior_probs_r <- BF1_listwise_r$PHP_confirmatory  # PMPs
  # Initialize a named vector to store inclusion probabilities
  inclusion_probs_r <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained_r <- !grepl(paste0(var, "=0"), names(posterior_probs_r))
    # Sum the posterior probabilities of these models
    inclusion_probs_r[var] <- sum(posterior_probs_r[not_constrained_r])
  }
  inclusion_data10_listwise[r,] <- inclusion_probs_r
  
  print(c(10,r/length(data10)))
}

inclusion_data20 <- inclusion_data20_listwise <- matrix(NA,ncol=ncol(X),nrow=length(data10))
for(r in 1:length(data20)){
  mice_miss <- mice :: mice(data = data20[[r]], m = M,
                            meth = rep("norm",8),
                            diagnostics = F, printFlag = F)
  relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
    y <- mice::complete(mice_miss, m)[,1] # same as 'y'
    X <- as.matrix(mice::complete(mice_miss, m)[,-1])
    colnames(X) <- NULL
    lm1_miss <- lm(y ~ 1 + X)
    # get PMPs based on FBF
    BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
    #extract quantities of interest
    c(BF1_miss$BFtable_confirmatory[, 1:4])
  },mc.cores=5)),ncol = M)
  while(is.character(relmeas_all[1,1])){ # in the (rare) case the imputed data gives errors Then redo imputation.
    mice_miss <- mice :: mice(data = data20[[r]], m = M,
                              meth = rep("norm",8),
                              diagnostics = F, printFlag = F)
    relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
      y <- mice::complete(mice_miss, m)[,1] # same as 'y'
      X <- as.matrix(mice::complete(mice_miss, m)[,-1])
      colnames(X) <- NULL
      lm1_miss <- lm(y ~ 1 + X)
      # get PMPs based on FBF
      BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
      #extract quantities of interest
      c(BF1_miss$BFtable_confirmatory[, 1:4])
    },mc.cores=5)),ncol = M)
  }
  # get expected values of quantities of interest
  relmeas <- matrix(apply(relmeas_all, 1, mean),ncol = 4)
  #get FBFs
  FBFs <- relmeas[,3]/relmeas[,1] * relmeas[,4]/relmeas[,2]
  #get PMPs
  PMPs <- FBFs / sum(FBFs)
  names(PMPs) <- names(posterior_probs)
  #get inclusion probabilities (same as above)
  inclusion_probs_miss <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained <- !grepl(paste0(var, "=0"), names(PMPs))
    # Sum the posterior probabilities of these models
    inclusion_probs_miss[var] <- sum(PMPs[not_constrained])
  }
  inclusion_data20[r,] <- inclusion_probs_miss
  
  # inclusion probs based on FBF after list-wise deletion
  data20_listwise_r <- data20[[r]][!is.na(apply(data20[[r]],1,sum)),]
  y <- data20_listwise_r[,1]
  X <- as.matrix(data20_listwise_r[,-1])
  colnames(X) <- NULL
  lm1_listwise_r <- lm(y ~ 1 + X)
  BF1_listwise_r <- BF(lm1_listwise_r,hypothesis=all_hypotheses)
  # Obtain inclusion probabilities for each predictor by summing the posterior model
  # probabilities which includes that predictor
  posterior_probs_r <- BF1_listwise_r$PHP_confirmatory  # PMPs
  # Initialize a named vector to store inclusion probabilities
  inclusion_probs_r <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained_r <- !grepl(paste0(var, "=0"), names(posterior_probs_r))
    # Sum the posterior probabilities of these models
    inclusion_probs_r[var] <- sum(posterior_probs_r[not_constrained_r])
  }
  inclusion_data20_listwise[r,] <- inclusion_probs_r
  
  print(c(20,r/length(data10)))
}

inclusion_data30 <- inclusion_data30_listwise <- matrix(NA,ncol=ncol(X),nrow=length(data10))
for(r in 1:length(data30)){
  mice_miss <- mice :: mice(data = data30[[r]], m = M,
                            meth = rep("norm",8),
                            diagnostics = F, printFlag = F)
  relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
    y <- mice::complete(mice_miss, m)[,1] # same as 'y'
    X <- as.matrix(mice::complete(mice_miss, m)[,-1])
    colnames(X) <- NULL
    lm1_miss <- lm(y ~ 1 + X)
    # get PMPs based on FBF
    BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
    #extract quantities of interest
    c(BF1_miss$BFtable_confirmatory[, 1:4])
  },mc.cores=5)),ncol = M)
  while(is.character(relmeas_all[1,1])){ # in the (rare) case the imputed data gives errors Then redo imputation.
    mice_miss <- mice :: mice(data = data30[[r]], m = M,
                              meth = rep("norm",8),
                              diagnostics = F, printFlag = F)
    relmeas_all <- matrix(unlist(mclapply(1:M,function(m){
      y <- mice::complete(mice_miss, m)[,1] # same as 'y'
      X <- as.matrix(mice::complete(mice_miss, m)[,-1])
      colnames(X) <- NULL
      lm1_miss <- lm(y ~ 1 + X)
      # get PMPs based on FBF
      BF1_miss <- BF(lm1_miss,hypothesis=all_hypotheses)
      #extract quantities of interest
      c(BF1_miss$BFtable_confirmatory[, 1:4])
    },mc.cores=5)),ncol = M)
  }
  # get expected values of quantities of interest
  relmeas <- matrix(apply(relmeas_all, 1, mean),ncol = 4)
  #get FBFs
  FBFs <- relmeas[,3]/relmeas[,1] * relmeas[,4]/relmeas[,2]
  #get PMPs
  PMPs <- FBFs / sum(FBFs)
  names(PMPs) <- names(posterior_probs)
  #get inclusion probabilities (same as above)
  inclusion_probs_miss <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained <- !grepl(paste0(var, "=0"), names(PMPs))
    # Sum the posterior probabilities of these models
    inclusion_probs_miss[var] <- sum(PMPs[not_constrained])
  }
  inclusion_data30[r,] <- inclusion_probs_miss
  
  # inclusion probs based on FBF after list-wise deletion
  data30_listwise_r <- data30[[r]][!is.na(apply(data30[[r]],1,sum)),]
  y <- data30_listwise_r[,1]
  X <- as.matrix(data30_listwise_r[,-1])
  colnames(X) <- NULL
  lm1_listwise_r <- lm(y ~ 1 + X)
  BF1_listwise_r <- BF(lm1_listwise_r,hypothesis=all_hypotheses)
  # Obtain inclusion probabilities for each predictor by summing the posterior model
  # probabilities which includes that predictor
  posterior_probs_r <- BF1_listwise_r$PHP_confirmatory  # PMPs
  # Initialize a named vector to store inclusion probabilities
  inclusion_probs_r <- setNames(numeric(length(vars)), vars)
  # Loop through each variable
  for (var in vars) {
    # Find models that do not constrain this variable to zero
    not_constrained_r <- !grepl(paste0(var, "=0"), names(posterior_probs_r))
    # Sum the posterior probabilities of these models
    inclusion_probs_r[var] <- sum(posterior_probs_r[not_constrained_r])
  }
  inclusion_data30_listwise[r,] <- inclusion_probs_r
  
  print(c(30,r/length(data10)))
}
#combine results based on MI and list-wise delation
inclusion_data10_combined <- matrix(NA,nrow=nrow(inclusion_data10),ncol=3*ncol(inclusion_data10))
inclusion_data10_combined[,1+3*(0:6)] <- inclusion_data10_listwise
inclusion_data10_combined[,2+3*(0:6)] <- 1
inclusion_data10_combined[,3+3*(0:6)] <- inclusion_data10
inclusion_data20_combined <- matrix(NA,nrow=nrow(inclusion_data20),ncol=3*ncol(inclusion_data20))
inclusion_data20_combined[,1+3*(0:6)] <- inclusion_data20_listwise
inclusion_data20_combined[,2+3*(0:6)] <- 1
inclusion_data20_combined[,3+3*(0:6)] <- inclusion_data20
inclusion_data30_combined <- matrix(NA,nrow=nrow(inclusion_data30),ncol=3*ncol(inclusion_data30))
inclusion_data30_combined[,1+3*(0:6)] <- inclusion_data30_listwise
inclusion_data30_combined[,2+3*(0:6)] <- 1
inclusion_data30_combined[,3+3*(0:6)] <- inclusion_data30

colnames(inclusion_data10_combined) <- colnames(inclusion_data20_combined) <-
  colnames(inclusion_data30_combined) <- rep(paste0("x",4:10),each=3)
boxplot(inclusion_data10_combined,yaxt = "n",border=rep(c(3,0,4),length=21),col="white",
        boxlwd=2,staplelwd=2,whisklwd=2,whisklty=1)
axis(2, yaxp = c(0, 1, 4))
for(yy in 1:4){
  abline(h=yy/4,col="grey",lty=1)
}
for(yy in 1:8){
  abline(h=yy/8,col="grey",lty=3)
}
points(2+3*(0:6),inclusion_probs,col="red",pch=15)
boxplot(inclusion_data20_combined,yaxt = "n",border=rep(c(3,0,4),length=21),col="white",
        boxlwd=2,staplelwd=2,whisklwd=2,whisklty=1)
axis(2, yaxp = c(0, 1, 4))
for(yy in 1:4){
  abline(h=yy/4,col="grey",lty=1)
}
for(yy in 1:8){
  abline(h=yy/8,col="grey",lty=3)
}
points(2+3*(0:6),inclusion_probs,col="red",pch=15)
boxplot(inclusion_data30_combined,yaxt = "n",border=rep(c(3,0,4),length=21),col="white",
        boxlwd=2,staplelwd=2,whisklwd=2,whisklty=1)
axis(2, yaxp = c(0, 1, 4))
for(yy in 1:4){
  abline(h=yy/4,col="grey",lty=1)
}
for(yy in 1:8){
  abline(h=yy/8,col="grey",lty=3)
}
points(2+3*(0:6),inclusion_probs,col="red",pch=15)
