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


library(tidyr)
archer1 <- 
  scrape_statcast_savant(
    start_date = "2016-03-25", 
    end_date = "2016-10-01", 
    playerid = 502042, # Chris Archer, primarily a 2-3 pitch guy. 
    player_type='pitcher')

archer_clean <- archer1 %>%
  mutate(pitch_type = ifelse(pitch_type == "FF", pitch_type, "OS")) %>%
  map_df(rev) %>%
  mutate(previous_pitch = ifelse(pitch_number == 1, "null", lag(pitch_type)),
         previous_des =   ifelse(pitch_number == 1, "null", lag(description)))

archer_clean$count <- paste0(archer_clean$balls, "-", archer_clean$strikes)
archer_clean$ahead <- ifelse(archer_clean$count %in% c("0-0", "0-1", "0-2",
                                                       "1-2", "2-2", "3-2"), 1, 0)
archer_clean$strikes2 <- ifelse(archer_clean$strikes == 2, 1, 0)

archer_clean$scoring_pos <- ifelse(is.na(archer_clean$on_2b) & is.na(archer_clean$on_3b), 1, 0)
archer_clean$left <- ifelse(archer_clean$stand == "L", 1, 0) 


# subset
df_bin <- archer_clean[archer_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up")]
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)

# Normalize before so we don't have to worry about it in the model.
#df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
#names(df_bin)[1] <- "pitch_type"
df_test <- df_bin[3001:nrow(df_bin), ]
df_bin <- df_bin[1:3000, ]
N <- nrow(df_bin)
J <- length(pitches_map)
d <- ncol(df_bin)


feats <- ncol(df_bin) - 1
feats

alphas <- matrix(data = 1, nrow=2)
alphas

beta1 <- matrix(data=1, nrow=2, ncol=4)
beta2 <- matrix(data=1, nrow=2, ncol=4)


m0 <- UpdateNB(df_bin[,], alphas = alphas, beta1 = beta1, beta2=beta2)
m1 <- UpdateMNNB(feats = df_bin[,-6])
predsnb <- predictNB(df_test[,c(-1,-6)], m0$post_theta, m0$post_pi)
nrow(predsnb)
nrow(df_test[,1])
predsnb <- cbind(predsnb, apply(predsnb, 1, which.max), ifelse(df_test[,1] == "FF",1,2))
mean(predsnb[,3] == predsnb[,4]) 
table(predsnb[,3],df_test$pitch_type) 

Xnew <- df_test[,c(-1,-6)]
post_pi <- m1$post_pi
post_theta <- m1$post_theta

post_ps <- matrix(nrow=nrow(Xnew), ncol=nrow(post_pi))
Pc <- NULL
for (i in 1:nrow(Xnew)) {
  X <- as.vector(Xnew[i,])
  for (c in 1:nrow(post_pi)) {
    Lc <- log(post_pi[c]) 
    for (j in 1:length(X)) {
      Lc <- Lc + log(post_theta[[j]][c, as.integer(X[j])+1])
    }
    Pc[c] <- exp(Lc)
  }
  post_ps[i,] <- Pc / sum(Pc)
}
post_ps

preds <- cbind(post_ps, apply(post_ps, 1, which.max), ifelse(df_test[,1] == "FF",1,2))
head(preds)
mean(preds[,3] == preds[,4])


Lc <- log(pi[c]) 
for (j in 1:length(X)) {
  Lc <- Lc + log(theta[[j]][c, as.integer(X[j])+1])
}
Pc <- exp(Lc)
post_ps

archer_clean$previous_des[archer_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
archer_clean$previous_des[archer_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip")] <- 1
archer_clean$previous_des[archer_clean$previous_des %in% c("called_strike")] <- 2
archer_clean$previous_des[archer_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
archer_clean$previous_des <- as.integer(archer_clean$previous_des)

df_bin <- archer_clean[archer_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up",
                                                          "previous_des", "balls", "strikes", "pitch_number")]

df_test <- df_bin[3001:nrow(df_bin), ]
df_bin <- df_bin[1:3000, ]

m2 <- UpdateMNNB(feats = df_bin)
Xnew <- df_test[,c(-1)]
post_pi <- m2$post_pi
post_theta <- m2$post_theta

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
post_ps

preds <- cbind(post_ps, apply(post_ps, 1, which.max), ifelse(df_test[,1] == "FF",1,2))
head(preds)
mean(preds[,3] == preds[,4])


max_clean <- max1 %>%
  filter(pitch_type != "null") %>%
  mutate(pitch_type = ifelse(pitch_type %in% c("FF", "FC", "FT"), "FF", 
                             ifelse(pitch_type %in% c("SL", "CH"), "SL-CH", "CU"))) %>%
  map_df(rev) %>%
  mutate(previous_pitch = ifelse(pitch_number == 1, "null", lag(pitch_type)),
         previous_des =   ifelse(pitch_number == 1, "null", lag(description)))

max_clean$count <- paste0(max_clean$balls, "-", max_clean$strikes)
max_clean$ahead <- ifelse(max_clean$count %in% c("0-0", "0-1", "0-2",
                                                 "1-2", "2-2", "3-2"), 1, 0)
max_clean$strikes2 <- ifelse(max_clean$strikes == 2, 1, 0)

max_clean$scoring_pos <- ifelse(is.na(max_clean$on_2b) & is.na(max_clean$on_3b), 1, 0)
max_clean$left <- ifelse(max_clean$stand == "L", 1, 0) 

# subset

# df_bin <- max_clean[max_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up",
#                                                           "previous_des", "balls", "strikes", "at_bat_number", "pitch_number")]
# head(df_bin)
# df_bin$pitch_type <- as.factor(df_bin$pitch_type)
# pitches_map <- levels(df_bin$pitch_type)
# # Normalize before so we don't have to worry about it in the model.
# df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
# names(df_bin)[1] <- "pitch_type"
# df_test <- df_bin[3001:nrow(df_bin), ]
# df_bin <- df_bin[1:3000, ]
# N <- nrow(df_bin)
# J <- length(pitches_map)
# d <- ncol(df_bin)


max_clean$previous_des[max_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
max_clean$previous_des[max_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip", "missed_bunt")] <- 1
max_clean$previous_des[max_clean$previous_des %in% c("called_strike")] <- 2
max_clean$previous_des[max_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
max_clean$previous_des <- as.integer(max_clean$previous_des)

table(max_clean$pitch_type)
pitches_map <- levels(max_clean$pitch_type)
#max_clean$pitch_type[max_clean$pitch_type %in% c("FF", "FT", "FC")] <- "FF"
max_clean$pitch_number[max_clean$pitch_number > 5] <- 5
max_clean$pitch_type <- as.factor(max_clean$pitch_type)
as.integer(max_clean$pitch_type)

df_bin <- max_clean[max_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "outs_when_up",
                                                          "previous_des", "pitch_number", "previous_pitch")]

max_df <- df_bin
train_ind <- sample(1:nrow(max_df), size = nrow(max_df)*.75)
max_test <- max_df[-train_ind, ]
max_train <- max_df[train_ind, ]

m2 <- UpdateMNNB(feats = max_train)
Xnew <- max_test[,c(-1)]
post_pi <- m2$post_pi
post_theta <- m2$post_theta

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
post_ps

preds <- cbind(post_ps, apply(post_ps, 1, which.max), as.integer(df_test$pitch_type))
head(preds)
mean(preds[,5] == preds[,6])
table(preds[,5], preds[,6])

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

phats <- predictMNNB(Xnew, m2$post_pi, m2$post_theta)
preds <- as.data.frame(cbind(phats, apply(phats, 1, which.max), as.integer(max_test$pitch_type)))
names(preds) <- c(pitches_map,"pred_type", "pitch_type")
preds$pred_type <- pitches_map[preds$pred_type]
preds$pitch_type <- pitches_map[preds$pitch_type]

head(preds)
write.csv(preds, "MaxPredsNB.csv", row.names=F)

set.seed(420)
archer_df <- df_bin
archer_df$pitch_number[archer_df$pitch_number > 5] <- 6
train_ind <- sample(1:nrow(archer_df), size = nrow(archer_df)*.75)
archer_train <- archer_df[train_ind,]
archer_test <- archer_df[-train_ind, ]

archMNNB <- UpdateMNNB(feats=archer_train)
arch_phats <- predictMNNB(archer_test[,-1], archMNNB$post_pi, archMNNB$post_theta)
arch_preds <- as.data.frame(arch_phats)
names(arch_preds) <- pitches_map
arch_preds$pred_class <- pitches_map[apply(arch_phats,1,which.max)]
arch_preds$actual <- archer_test$pitch_type
head(arch_preds)

mean(arch_preds$pred_class == arch_preds$actual)
table(arch_preds$pred_class, arch_preds$actual)

write.csv(arch_preds, "ArchPredsNB.csv", row.names=F)
