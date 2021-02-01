setwd("~/URI/bayes-545/Project")

# read in Data
library(baseballr)
library(ggplot2)
library(dplyr)
library(purrr)
library(mcmcplots)
library(rjags)
library(runjags)

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
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
df_bin <- df_bin[,-5]
train_ind <- sample(1:nrow(df_bin), size=nrow(df_bin)*.75)
df_test <- df_bin[-train_ind, ]
df_train <- df_bin[train_ind, ]
N <- nrow(df_train)
J <- length(pitches_map)
d <- ncol(df_train)

# fit models 

#### Model 1 left + ahead + scoring_pos + strikes2 + outs_when_up #### 
model_string1 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * strikes2[i] + 
        b[6,j] * outs_when_up[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:6){
    b[k,1] <- 0          ## Constrain to 0 for base comparison
    for(j in 2:J){
      b[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}
"

model_string3 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * outs_when_up[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:5){
    b[k,1] <- 0          ## Constrain to 0 for base comparison
    for(j in 2:J){
      b[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}
"


inits=function(){ list(b = matrix(data = c(rep(NA,d), rep(0,d * (J-1))),
                                  nrow=d,ncol=J))}

data <- list(left = df_train$left,
             ahead = df_train$ahead,
             scoring_pos = df_train$scoring_pos,
             outs_when_up = df_train$outs_when_up,
             pitch_type = df_train$pitch_type,
             N = N,
             J = J)

archer.fit3=run.jags(model=model_string3,
                    monitor= c("b"),
                    burnin=1000, sample=10000,
                    data=data, n.chains=2, method="rjags", inits=inits)

write.jagsfile(archer.fit3, file="archer3.txt")

plot(archer.fit3, plot.type = c("trace", "ecdf", "histogram",
                               "autocorr"), vars = NA,
     layout = runjags.getOption("plot.layout"),
     new.windows = runjags.getOption("new.windows"), file = "",
     mutate = NULL, col = NA, trace.iters = NA, separate.chains = NA,
     trace.options = NA, density.options = NA, histogram.options = NA,
     ecdfplot.options = NA, acplot.options = NA)

mod3 <- as.mcmc(archer.fit3)
write.csv(mod3, "Archer3mcmc.csv", row.names=F)
autocorr.diag(mod1)
autocorr.plot(mod1)
effectiveSize(mod1)


par(mfrow=c(2,3))
densplot(mod1[,7:12], xlim=c(-0.5,0.5))

dic3 <- dic.samples(as.jags(archer.fit3), n.iter=1e3, "pD")
dic3
# Mean deviance:  3968 
# penalty 12.5 
# Penalized deviance: 3981 

# Mean deviance:  3496 
# penalty 6.992 
# Penalized deviance: 3503
names(mod3)[6:10] <- c("intercept", "Left", "ahead", "scoring", "outs")
par(mfrow=c(2,3))
for (i in 6:10) {
  if (prod(quantile(mod3[,i], c(0.025, 0.975))) > 0) { # same sign Significant!!
    densplot(mod3[,i], main=names(mod3)[i],xlim=c(-0.40,0.40),col="red")
    abline(v=0,col='darkblue')
  } else {
    densplot(mod3[,i], main=names(mod3)[i],xlim=c(-0.40,0.40))
    abline(v=0,col='darkblue')
  }
}


pm_coef <- colMeans(mod3)
pm_coef
# going to need to change this to get multiclass predictions. 
pm_Xb = pm_coef[6] + as.matrix(df_test[,-1]) %*% pm_coef[7:10]
phat <- 1.0 / (1 + exp(-pm_Xb))
head(phat)

preds <- cbind(phat, (phat > 0.5) + 1, df_test$pitch_type)
head(preds)
mean(preds[,2] == preds[,3])

# Make prediction
pm_coef <- colMeans(mod1)
pm_coef
# going to need to change this to get multiclass predictions. 
pm_Xb = pm_coef[7] + as.matrix(df_test[,-1]) %*% pm_coef[8:12]
phat <- 1.0 / (1 + exp(-pm_Xb))
head(phat)

preds <- cbind(phat, (phat > 0.5) + 1, df_test$pitch_type)
head(preds)
mean(preds[,2] == preds[,3])

# .55 accuracy LAME

#### Model 2 ####

# subset
# fix up description column
archer_clean$previous_des[archer_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
archer_clean$previous_des[archer_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip")] <- 1
archer_clean$previous_des[archer_clean$previous_des %in% c("called_strike")] <- 2
archer_clean$previous_des[archer_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
archer_clean$previous_des <- as.integer(archer_clean$previous_des)

df_bin <- archer_clean[archer_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up",
                                                          "previous_des", "balls", "strikes", "at_bat_number", "pitch_number")]
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)

# Normalize before so we don't have to worry about it in the model.
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
df_test <- df_bin[3001:nrow(df_bin), ]
df_bin <- df_bin[1:3000, ]
N <- nrow(df_bin)
J <- length(pitches_map)
d <- ncol(df_bin)

# fit models 

#### Model 1 left + ahead + scoring_pos + strikes2 + outs_when_up #### 
model_string2 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * strikes2[i] + 
        b[6,j] * outs_when_up[i] + 
        b[7,j] * previous_des[i] + 
        b[8,j] * balls[i] + 
        b[9,j] * strikes[i] + 
        b[10,j] * at_bat_number[i] + 
        b[11,j] * pitch_number[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:11){
    b[k,1] <- 0          ## Constrain to 0 for base comparison
    for(j in 2:J){
      b[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}
"



inits=function(){ list(b = matrix(data = c(rep(NA,d), rep(0,d * (J-1))),
                                  nrow=d,ncol=J))}

data <- list(left = df_bin$left,
             ahead = df_bin$ahead,
             scoring_pos = df_bin$scoring_pos,
             strikes2 = df_bin$strikes2,
             outs_when_up = df_bin$outs_when_up,
             previous_des = df_bin$previous_des,
             balls = df_bin$balls,
             strikes = df_bin$strikes,
             at_bat_number = df_bin$at_bat_number,
             pitch_number = df_bin$pitch_number,
             pitch_type = df_bin$pitch_type,
             N = N,
             J = J)

archer.fit2=run.jags(model=model_string2,
                     monitor= c("b"),burnin=1000, sample=10000,
                     data=data, n.chains=2, method="rjags", inits=inits)
                     

write.jagsfile(archer.fit2, file="archer2.txt")

plot(archer.fit2, plot.type = c("trace", "ecdf", "histogram",
                                "autocorr"), vars = NA,
     layout = runjags.getOption("plot.layout"),
     new.windows = runjags.getOption("new.windows"), file = "",
     mutate = NULL, col = NA, trace.iters = NA, separate.chains = NA,
     trace.options = NA, density.options = NA, histogram.options = NA,
     ecdfplot.options = NA, acplot.options = NA)

mod2 <- as.mcmc(archer.fit2)

autocorr.diag(mod2)
autocorr.plot(mod2)
effectiveSize(mod2)


par(mfrow=c(2,3))
densplot(mod2[,12:22], xlim=c(-0.5,0.5))

dic2 <- dic.samples(as.jags(archer.fit2), n.iter=1e3, "pD")
dic2
#
#
#

# Make prediction
pm_coef <- colMeans(mod2)
pm_coef
# going to need to change this to get multiclass predictions. 
pm_Xb = pm_coef[12] + as.matrix(df_test[,-1]) %*% pm_coef[13:22]
phat <- 1.0 / (1 + exp(-pm_Xb))
head(phat)

preds <- cbind(phat, (phat > 0.5) + 1, df_test$pitch_type)
head(preds)
mean(preds[,2] == preds[,3])

# Accuracy



#### Jake Arrieta 453562 ####
# try stepwise regression to identify best variables, then use jags



#### Max Scherzer 453286 ####
max1 <- 
  scrape_statcast_savant(
    start_date = "2016-03-25", 
    end_date = "2016-10-01", 
    playerid = 453286, # Chris Archer, primarily a 2-3 pitch guy. 
    player_type='pitcher')

max_clean <- max1 %>%
  filter(pitch_type != "null") %>%
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

df_bin <- max_clean[, c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up")]
head(df_bin)
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)
# Normalize before so we don't have to worry about it in the model.
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
df_test <- df_bin[3001:nrow(df_bin), ]
df_bin <- df_bin[1:3000, ]
N <- nrow(df_bin)
J <- length(pitches_map)
d <- ncol(df_bin)


model_string1 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * strikes2[i] + 
        b[6,j] * outs_when_up[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:6){
    b[k,1] <- 0          ## Constrain to 0 for base comparison
    for(j in 2:J){
      b[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}
"


inits=function(){ list(b = matrix(data = c(rep(NA,d), rep(0,d * (J-1))),
                                  nrow=d,ncol=J))}


data <- list(left = df_bin$left,
             ahead = df_bin$ahead,
             scoring_pos = df_bin$scoring_pos,
             strikes2 = df_bin$strikes2,
             outs_when_up = df_bin$outs_when_up,
             pitch_type = df_bin$pitch_type,
             N = N,
             J = J)

max.fit1=run.jags(model=model_string1,
                     monitor= c("b"),
                     burnin=1000, sample=10000,
                     data=data, n.chains=1, method="rjags", inits=inits)
# write.jagsfile(max.fit1, file="max1.txt")

max_fit <- read.jagsfile("max1.txt")

max2 <- read.csv("max2_mcmc.csv")

qm_coef <- colMeans(max2)
qm_coef
Xnew <- df_test[,-1]
qhats <- matrix(nrow=nrow(Xnew), ncol=6)
phats <- matrix(nrow=nrow(Xnew), ncol=6)
# going to need to change this to get multiclass predictions. 
for ( i in 1:6 ) {
  inds <- i + (i-1)*10
  qhats[,i] <- qm_coef[inds] + as.matrix(Xnew) %*% qm_coef[(inds+1):(inds+10)]
}

for (i in 1:nrow(qhats)) {
  phats[i,] <- exp(qhats[i,]) / sum(exp(qhats[i,]))
}

head(phats)
sum(phats[1,])

###### Things to redo:
# Redo new max_clean columns
# Copy code from model 2 to get right code in place for max.fit
# Fix df_test etc to get predictions based on the better model.

#### Model 2 ####

# subset
# fix up description column
max_clean$previous_des[max_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
max_clean$previous_des[max_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip", "missed_bunt")] <- 1
max_clean$previous_des[max_clean$previous_des %in% c("called_strike")] <- 2
max_clean$previous_des[max_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
max_clean$previous_des <- as.integer(max_clean$previous_des)

df_bin <- max_clean[max_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "strikes2", "outs_when_up",
                                                          "previous_des", "balls", "strikes", "at_bat_number", "pitch_number")]
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)

# Normalize before so we don't have to worry about it in the model.
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
df_test <- df_bin[3001:nrow(df_bin), ]
df_bin <- df_bin[1:3000, ]
N <- nrow(df_bin)
J <- length(pitches_map)
d <- ncol(df_bin)

# fit models 

#### Model 1 left + ahead + scoring_pos + strikes2 + outs_when_up #### 
model_string2 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * strikes2[i] + 
        b[6,j] * outs_when_up[i] + 
        b[7,j] * previous_des[i] + 
        b[8,j] * balls[i] + 
        b[9,j] * strikes[i] + 
        b[10,j] * at_bat_number[i] + 
        b[11,j] * pitch_number[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:11){
    b[k,1] <- 0          ## Constrain to 0 for base comparison
    for(j in 2:J){
      b[k,j] ~ dnorm(0, 0.1)
    }  # close J loop
  }  # close K loop
}
"



inits=function(){ list(b = matrix(data = c(rep(NA,d), rep(0,d * (J-1))),
                                  nrow=d,ncol=J))}

data <- list(left = df_bin$left,
             ahead = df_bin$ahead,
             scoring_pos = df_bin$scoring_pos,
             strikes2 = df_bin$strikes2,
             outs_when_up = df_bin$outs_when_up,
             previous_des = df_bin$previous_des,
             balls = df_bin$balls,
             strikes = df_bin$strikes,
             at_bat_number = df_bin$at_bat_number,
             pitch_number = df_bin$pitch_number,
             pitch_type = df_bin$pitch_type,
             N = N,
             J = J)



