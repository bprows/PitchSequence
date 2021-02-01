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

dic3 <- dic.samples(as.jags(archer.fit3), n.iter=1e3, "pD")
dic3
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
# Accuracy = 57%


#  subset
# fix up description column
archer_clean$previous_des[archer_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
archer_clean$previous_des[archer_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip")] <- 1
archer_clean$previous_des[archer_clean$previous_des %in% c("called_strike")] <- 2
archer_clean$previous_des[archer_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
archer_clean$previous_des <- as.integer(archer_clean$previous_des)
archer_clean$previous_pitch <- ifelse(archer_clean$previous_pitch == "FF", 1, 2)

df_bin <- archer_clean[archer_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "outs_when_up",
                                                          "previous_des", "pitch_number", "previous_pitch")]
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)

# Normalize before so we don't have to worry about it in the model.
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
train_ind <- sample(1:nrow(df_bin), size=nrow(df_bin)*.75)
df_test <- df_bin[-train_ind, ]
df_train <- df_bin[train_ind, ]
N <- nrow(df_train)
J <- length(pitches_map)
d <- ncol(df_train)

model_string4 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * outs_when_up[i] + 
        b[6,j] * previous_des[i] +
        b[7,j] * pitch_number[i] +
        b[8,j] * previous_pitch[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:8){
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
             previous_des = df_train$previous_des,
             pitch_number = df_train$pitch_number,
             previous_pitch = df_train$previous_pitch,
             pitch_type = df_train$pitch_type,
             N = N,
             J = J)

archer.fit4=run.jags(model=model_string4,
                     monitor= c("b"),
                     burnin=1000, sample=10000,
                     data=data, n.chains=2, method="rjags", inits=inits)

write.jagsfile(archer.fit4, file="archer4.txt")

plot(archer.fit4, plot.type = c("trace", "ecdf", "histogram",
                                "autocorr"), vars = NA,
     layout = runjags.getOption("plot.layout"),
     new.windows = runjags.getOption("new.windows"), file = "",
     mutate = NULL, col = NA, trace.iters = NA, separate.chains = NA,
     trace.options = NA, density.options = NA, histogram.options = NA,
     ecdfplot.options = NA, acplot.options = NA)

mod4 <- as.mcmc(archer.fit4)
write.csv(mod4, "Archer4mcmc.csv", row.names=F)

dic4 <- dic.samples(as.jags(archer.fit4), n.iter=1e3, "pD")
dic4
# Mean deviance:  3369 
# penalty 2.687 
# Penalized deviance: 3372 

names(mod4)[9:16] <- c("intercept", "Left", "ahead", "scoring", "outs",
                       "previous_des", "pitch_number", "prev_pitch_type")
par(mfrow=c(2,4))
for (i in 9:16) {
  if (prod(quantile(mod4[,i], c(0.025, 0.975))) > 0) { # same sign Significant!!
    densplot(mod4[,i],main=names(mod4)[i],xlab="",xlim=c(-0.60,0.60),col="red")
    abline(v=0,col='darkblue')
  } else {
    densplot(mod4[,i],main=names(mod4)[i],xlab="",xlim=c(-0.60,0.60))
    abline(v=0,col='darkblue')
  }
}

mod4_table <- NULL
mod4_table$Feature <- names(mod4)[9:16]
mod4_table$CILB <- apply(mod4[,9:16], 2, function(x) quantile(x, 0.025))
mod4_table$CIUB <- apply(mod4[,9:16], 2, function(x) quantile(x, 0.975))
mod4_table$Significance <- ifelse(mod4_table$CILB * mod4_table$CIUB > 0,"*","") 

as.data.frame(mod4_table, row.names=mod4_table$Feature)[,-1]

pm_coef <- colMeans(mod4)
pm_coef
# going to need to change this to get multiclass predictions. 
pm_Xb = pm_coef[9] + as.matrix(df_test[,-1]) %*% pm_coef[10:16]
phat <- 1.0 / (1 + exp(-pm_Xb))
head(phat)

preds <- cbind(phat, (phat > 0.5) + 1, df_test$pitch_type)
head(preds)
mean(preds[,2] == preds[,3])
preds <- as.data.frame(preds)
names(preds) <- c("ProbOS", "Predicted", "Actual")
preds$Predicted <- pitches_map[preds$Predicted]
preds$Actual <- pitches_map[preds$Actual]

write.csv(preds, "Archer_PredsLR.csv", row.names=F)


#### Scherzer ####
max1 <- 
  scrape_statcast_savant(
    start_date = "2016-03-25", 
    end_date = "2016-10-01", 
    playerid = 453286, # Chris Archer, primarily a 2-3 pitch guy. 
    player_type='pitcher')

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

max_clean$previous_des[max_clean$previous_des %in% c("null", "ball", "blocked_ball", "pitchout")] <- 0
max_clean$previous_des[max_clean$previous_des %in% c("foul", "foul_bunt", "foul_tip", "missed_bunt")] <- 1
max_clean$previous_des[max_clean$previous_des %in% c("called_strike")] <- 2
max_clean$previous_des[max_clean$previous_des %in% c("swinging_strike", "swinging_strike_blocked")] <- 3
max_clean$previous_des <- as.integer(max_clean$previous_des)
max_clean$previous_pitch <- as.integer(as.factor(max_clean$previous_pitch))

df_bin <- max_clean[max_clean$pitch_type != "PO", c("pitch_type", "left", "ahead", "scoring_pos", "outs_when_up",
                                                    "previous_des", "pitch_number", "previous_pitch")]
df_bin$pitch_type <- as.factor(df_bin$pitch_type)
pitches_map <- levels(df_bin$pitch_type)

# Normalize before so we don't have to worry about it in the model.
df_bin <- as.data.frame(cbind(df_bin$pitch_type, scale(df_bin[,-1])))
names(df_bin)[1] <- "pitch_type"
train_ind <- sample(1:nrow(df_bin), size=nrow(df_bin)*.75)
df_test <- df_bin[-train_ind, ]
df_train <- df_bin[train_ind, ]
N <- nrow(df_train)
J <- length(pitches_map)
d <- ncol(df_train)

# fit models 

#### Model 1 left + ahead + scoring_pos + strikes2 + outs_when_up #### 
model_string5 <- "model
{
for(i in 1:N){
    pitch_type[i] ~ dcat(p[i, 1:J])
    
    for (j in 1:J){
      log(q[i,j]) <-  b[1,j] + 
        b[2,j] * left[i] + 
        b[3,j] * ahead[i] + 
        b[4,j] * scoring_pos[i]  + 
        b[5,j] * outs_when_up[i] + 
        b[6,j] * previous_des[i] +
        b[7,j] * pitch_number[i] +
        b[8,j] * previous_pitch[i]
      p[i,j] <- q[i,j]/(sum(q[i,1:J]))  ## should be familiar from MLE notes: q is exp(Xb)
    }   # close J loop
  }  # close N loop
  
  for(k in 1:8){
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
             previous_des = df_train$previous_des,
             pitch_number = df_train$pitch_number,
             previous_pitch = df_train$previous_pitch,
             pitch_type = df_train$pitch_type,
             N = N,
             J = J)

max.fit5=run.jags(model=model_string5,
                     monitor= c("b"),
                     burnin=1000, sample=10000,
                     data=data, n.chains=2, method="rjags", inits=inits)

write.jagsfile(max.fit5, file="max.fit5.txt")
mod5 <- as.mcmc(max.fit5)
write.csv(mod5, "Max5mcmc.csv", row.names=F)

dic5 <- dic.samples(as.jags(max.fit5), n.iter=1e3, "pD")
dic5
# Mean deviance:  4294 
# penalty 19.49 
# Penalized deviance: 4314 

#Densplots


# Predictions
qm_coef <- colMeans(mod5)
qm_coef
Xnew <- df_test[,-1]
qhats <- matrix(nrow=nrow(Xnew), ncol=3)
phats <- matrix(nrow=nrow(Xnew), ncol=3)
# going to need to change this to get multiclass predictions.
for ( i in 1:3 ) {
  inds <- i + (i-1)*7
  qhats[,i] <- qm_coef[inds] + as.matrix(Xnew) %*% qm_coef[(inds+1):(inds+7)]
}

for (i in 1:nrow(qhats)) {
  phats[i,] <- exp(qhats[i,]) / sum(exp(qhats[i,]))
}

head(phats)
preds <- data.frame(phats,
               pred_type = apply(phats,1,which.max),
               pitch_type = df_test$pitch_type)
preds$pred_type <- pitches_map[preds$pred_type]
preds$pitch_type <- pitches_map[preds$pitch_type]

mean(preds$pred_type==preds$pitch_type)
table(preds$pred_type, preds$pitch_type)

write.csv(preds, "Scherzer_Preds.csv", row.names = F)
