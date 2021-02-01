#### Naive Bayes Classifier ####

library(tidyr)
UpdateNB <- function(df, target = NULL, alphas, beta1, beta2) {
  if (!is.null(target)) {
    df <- cbind(target,df)
  }
  
  names(df)[1] <- "target"
  
  # Extract summary statistics
  Nc <- df %>%
    group_by(target) %>%
    summarise(Nc = n())
  targets <- Nc$target
  num_targets <- length(unique(targets))
  alphas <- alphas + Nc$Nc
  Nc <- matrix(data = Nc$Nc, nrow=num_targets, ncol=ncol(df)-1, byrow=F)
  
  Njc <- df %>%
    group_by(target) %>%
    summarise_all(sum)
  
  beta1s <- beta1 + Njc[,-1]
  beta2s <- beta2 + Nc - Njc[,-1]
  
  post_mean_pi <- alphas / sum(alphas)
  
  post_mean_thetas <- beta1s / (Nc + beta1 + beta2)
  
  # Change names
  row.names(post_mean_pi) <- targets
  row.names(post_mean_thetas) <- targets
  row.names(beta1s) <- targets
  row.names(beta2s) <- targets
  row.names(alphas) <- targets
  
  out <- list(post_pi = post_mean_pi, 
              post_theta = post_mean_thetas,
              beta1 = beta1s,
              beta2 = beta2s,
              alphas = alphas)
  return(out)
}

predictNB <- function(X, post_thetas, post_pi) {
  X_comp <- abs(X-1)
  post_ys <- matrix(data=NA, nrow=nrow(X), ncol=nrow(post_pi))
  for (i in 1:nrow(X)) { # Make prediction for each row of X
    post_y <- NULL
    for (p in 1:nrow(post_pi)) { # Calculate Posterior probs of each category.
      post_y[p] <- prod(post_thetas[p,] ** X[i,]) * prod((1-post_thetas[p,])**X_comp[i,]) * post_pi[p]
    }
    post_ys[i,] <- post_y/sum(post_y) # Append Predicted values to post_ys
  }
  colnames(post_ys) <- row.names(post_thetas)
  return(post_ys)
}


post_theta_calc <- function(mat) { #input is beta[[j]]
  mat <- as.matrix(mat)
  out <- matrix(data = NA, nrow=nrow(mat),ncol=ncol(mat)) # same dimensions
  for (i in 1:nrow(mat)) {
    out[i, ] <- mat[i,] / sum(mat[i,]) # this is being converted to a list!!!!!!! 
  }
  row.names(out) <- row.names(mat)
  colnames(out) <- colnames(mat)
  return(out)
}


UpdateMNNB <- function(feats, target = NULL, alphas = NULL, betas = NULL) {
  
  if(!is.null(target)) { # adjust if features passed in seperate from data
    feats <- cbind(target,feats)
  }
  
  names(feats)[1] <- "target"  # renames column 1 = target
  num_targ <- length(unique(feats$target))
  d <- ncol(feats) - 1
  
  if(is.null(alphas)) { # initialize alphas if not passed in
    alphas <- matrix(data= 1, nrow=num_targ)
  }
  
  if(is.null(betas)) { # initialize betas if none passed in
    betas <- list()
    kjs <- NULL
    for (j in 1:d) {
      kjs[j] <- length(unique(feats[,j+1]))
      betas[[j]] <- matrix(data=1, nrow=num_targ, ncol=kjs[j])
    }
  }
  
  # compute summary statistics
  Nc <- feats %>%
    group_by(target) %>%
    summarise(Nc = n())
  targets <- Nc$target
  alphas <- alphas + Nc$Nc
  
  # update betas
  for (j in 1:d) {
    Ncjk <- 
      feats %>% select(1,feat = j+1) %>%
      group_by(target, feat) %>%
      summarise(Ncjk = n()) %>%
      pivot_wider(names_from = feat, values_from = Ncjk)
    
    betas[[j]] <- Ncjk[-1] + betas[[j]]
    row.names(betas[[j]]) <- targets
    colnames(betas[[j]]) <- colnames(Ncjk)[-1]
    names(betas)[j] <- names(feats)[j+1]
  }
  
  post_pi <- alphas / sum(alphas)
  post_thetas <- lapply(betas, FUN = post_theta_calc)
  
  # Change names
  row.names(post_pi) <- targets
  row.names(alphas) <- targets
  
  out <- list(post_pi = post_pi, 
              post_theta = post_thetas,
              betas = betas,
              alphas = alphas)
  return(out)
}

predictMNNB <- function(Xnew, post_pi, post_theta) {
  post_ps <- matrix(nrow=nrow(Xnew), ncol=nrow(post_pi))
  Pc <- NULL
  for (i in 1:nrow(Xnew)) {
    X <- as.vector(Xnew[i,])
    for (c in 1:nrow(post_pi)) {
      Lc <- log(post_pi[c]) 
      for (j in 1:length(X)) {
        Lc <- Lc + log(post_theta[[j]][c, colnames(post_theta[[j]]) == as.character(X[j])])
      }
      Pc[c] <- exp(Lc)
    }
    post_ps[i,] <- Pc / sum(Pc)
  }
  return(post_ps)
}





