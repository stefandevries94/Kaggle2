setwd("~/Documents/Behavioral Data Science/BDA/Kaggle 2")

library(tidyverse)

dir("physical-activity-recognition/RawData/Train/", pattern = "^acc")

act_labels <- read_delim("physical-activity-recognition/activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% dplyr::select("X1","X2")

labels <- read_delim("physical-activity-recognition/RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')

labels <- labels %>% mutate(activity = act_labels$X2[activity])

train_labels = labels

# Entropy function --------------------------------------------------------

entropy  <- function(x, nbreaks = nclass.Sturges(x)) {
  r = range(x)
  x_binned = findInterval(x, seq(r[1], r[2], len= nbreaks))
  h = tabulate(x_binned, nbins = nbreaks) # fast histogram
  p = h/sum(h)
  -sum(p[p>0] * log(p[p>0]))
}

# Spectrum function -------------------------------------------------------

spect <- function(x) {
  s <- spectrum(x, log='n', plot=F)
  s$freq[which.max(s$spec)]
}



extractTimeDomainFeatures <- 
  function(filename, labels = train_labels) {
    # extract user and experimental run ID's from file name
    username = gsub(".+user(\\d+).+", "\\1", filename)
    expname =  gsub(".+exp(\\d+).+", "\\1", filename)
    
    # import the data from the file
    user01 <- read_delim(filename, " ", col_names = F, progress = TRUE, col_types = "ddd")
    
    # select this users activity labels from the `labels` frame for this experimental run
    user_labels <- 
      labels %>% 
      dplyr::filter(userid==as.numeric(username) & trial == as.numeric(expname)) %>% # select rows pertaining to current signals
      mutate(segment = row_number()) %>%                                             # segment identifies different WALKING_UPSTAIRS etc
      gather(start_end, vec, -trial, -userid, -activity, -segment) %>%               # stack start and end on top of each other
      arrange(vec) %>%                                                               # arrange rows in natural order
      mutate(activity = ifelse(start_end == 'end', NA, activity), activity_id = row_number()) # remove activity label `end` time rows
    
    # add activity labels to each sample
    user <- 
      user01 %>% 
      mutate(sample = row_number()-1) %>%
      mutate(activity_id = findInterval(sample, user_labels$vec)) %>%
      left_join(user_labels, by = "activity_id") 
    
    
    # split in epochs of 128 samples and compute features per epoch
    usertimedom <- 
      user %>%
      mutate(epoch = sample %/% 128) %>% # epoch = 2.56 sec
      group_by(epoch) %>%
      summarise(
        user_id = username, # added to identify user in data frame rows
        exp_id = expname,   # added to identify experimental run in data frame rows
        activity = names(which.max(table(c("-", activity)))),
        sample = sample[1],
        p1 = mean(X1)^2,
        p2 = mean(X2)^2,
        p3 = mean(X3)^2,
        p12 = mean(sqrt(X1^2 + X2^2))^2,
        p13 = mean(sqrt(X1^2 + X3^2))^2,
        p23 = mean(sqrt(X2^2 + X3^2))^2,
        p123 = mean(sqrt(X1^2 + X2^2 + X3^2))^2,
        sd1 = sd(X1),
        sd2 = sd(X2),
        sd3 = sd(X3),
        min1 = min(X1),
        min2 = min(X2),
        min3 = min(X3),
        max1 = max(X1),
        max2 = max(X2),
        max3 = max(X3),
        q1_25 = quantile(X1, .25),
        q1_50 = quantile(X1, .50),
        q1_75 = quantile(X1, .75),
        q2_25 = quantile(X2, .25),
        q2_50 = quantile(X2, .50),
        q2_75 = quantile(X2, .75),
        q3_25 = quantile(X3, .25),
        q3_50 = quantile(X3, .50),
        q3_75 = quantile(X3, .75),
        skew1 = e1071::skewness(X1),
        skew2 = e1071::skewness(X2),
        skew3 = e1071::skewness(X3),
        AR1.1 = cor(X1, lag(X1), use = "pairwise"),
        AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
        AR2.1 = cor(X2, lag(X2), use = "pairwise"),
        AR2.2 = cor(X2, lag(X2, n = 2), use = "pairwise"),
        AR3.1 = cor(X3, lag(X3), use = "pairwise"),
        AR3.2 = cor(X3, lag(X3, n = 2), use = "pairwise"),
        AR12 = cor(X1, X2, use = "pairwise"),
        AR13 = cor(X1, X3, use = "pairwise"),
        AR23 = cor(X2, X3, use = "pairwise"),
        AR12.1 = cor(X1, lag(X2), use = "pairwise"),
        AR13.1 = cor(X1, lag(X3), use = "pairwise"),
        AR23.1 = cor(X2, lag(X3), use = "pairwise"),
        AR21.1 = cor(X2, lag(X1), use = "pairwise"),
        AR31.1 = cor(X3, lag(X1), use = "pairwise"),
        AR32.1 = cor(X3, lag(X2), use = "pairwise"),
        entr.1 = entropy(X1),
        entr.2 = entropy(X2),
        entr.3 = entropy(X3),
        spec.1 = spect(X1),
        spec.2 = spect(X2),
        spec.3 = spect(X3),
        n=n()
      ) 
    
    usertimedom 
  }

filenamesACC <- dir("physical-activity-recognition/RawData/Train/", "^acc", full.names = TRUE)
filenamesGYRO <- dir("physical-activity-recognition/RawData/Train/", "^gyro", full.names = TRUE)

Acc_data =  filenamesACC %>%
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Acc"))

Gyro_data = filenamesGYRO %>% 
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Gyro"))

Full_Train <- inner_join(Acc_data, Gyro_data) %>% 
  filter(activity != "-")

filenamesACC <- dir("physical-activity-recognition/RawData/Test/", "^acc", full.names = TRUE)
filenamesGYRO <- dir("physical-activity-recognition/RawData/Test/", "^gyro", full.names = TRUE)

Acc_data =  filenamesACC %>%
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Acc"))

Gyro_data = filenamesGYRO %>% 
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Gyro"))

Full_Test <- inner_join(Acc_data, Gyro_data)

Full_Train %>%
  ggplot(aes(AR12_Acc)) + 
  geom_histogram(bins=40, fill=1, alpha=0.5) + 
  geom_histogram(aes(AR13_Acc), bins=40, fill = 2, alpha=0.5) + 
  geom_histogram(aes(AR23_Acc), bins=40, fill = 4, alpha=0.5) +
  facet_wrap(~activity, scales = "free_y")

#################################################################################
# Models
library(MASS)
lda_test <- lda(activity ~ . + p1_Acc * p1_Gyro + p2_Acc * p2_Gyro + p3_Acc * p3_Gyro, data = Full_Train[, -c(1:3,5,55)], CV = T)

# Assess the accuracy of the prediction
# percent correct for each category of G
ct <- table(Full_Train$activity, lda_test$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))


lda_test <- lda(activity ~ . + p1_Acc * p1_Gyro + p2_Acc * p2_Gyro + p3_Acc * p3_Gyro, data = Full_Train[, -c(1:3,5,55)])
predicted <- predict(lda_test, Full_Test)
Predicted <- as.vector(predicted$class)

Full_Test$activity <-  as.vector(predicted$class)
Full_Test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions.csv")

file.show("test_set_predictions.csv")



