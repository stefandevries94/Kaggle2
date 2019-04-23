library(tidyverse)

dir("physical-activity-recognition/RawData/Train/", pattern = "^acc")

act_labels <- read_delim("physical-activity-recognition/activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% dplyr::select("X1","X2")

labels <- read_delim("physical-activity-recognition/RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')

labels <- labels %>% mutate(activity = act_labels$X2[activity])

train_labels = labels

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
        m1 = mean(X1), 
        m2 = mean(X2),
        m3 = mean(X3),
        sd1 = sd(X1),
        sd2 = sd(X2),
        sd3 = sd(X3),
        #q1_25 = quantile(X1, .25),
        skew1 = e1071::skewness(X1),
        skew2 = e1071::skewness(X2),
        skew3 = e1071::skewness(X3),
        AR1.1 = cor(X1, lag(X1), use = "pairwise"),
        AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
        AR2.1 = cor(X2, lag(X2), use = "pairwise"),
        AR2.2 = cor(X2, lag(X2, n = 2), use = "pairwise"),
        AR3.1 = cor(X3, lag(X3), use = "pairwise"),
        AR3.2 = cor(X3, lag(X3, n = 2), use = "pairwise"),
        n=n()
      ) 
    
    usertimedom 
  }

filename = "physical-activity-recognition/RawData/Train/acc_exp01_user01.txt"
extractTimeDomainFeatures(filename)

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
  ggplot(aes(m1_Gyro)) + 
  geom_histogram(bins=40, fill=1, alpha=0.5) + 
  geom_histogram(aes(m2_Gyro), bins=40, fill = 2, alpha=0.5) + 
  geom_histogram(aes(m3_Gyro), bins=40, fill = 4, alpha=0.5) +
  facet_wrap(~activity, scales = "free_y")



#################################################################################
# Models
library(MASS)
lda_test <- lda(activity ~ ., data = Full_Train[, -c(1:3,5,21)])
predicted <- predict(lda_test, Full_Test)
summary(lda_test)
Predicted <- as.vector(predicted$class)


Full_Test$activity <- as.vector(predicted$class)
Full_Test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions.csv")

file.show("test_set_predictions.csv")


##### QDA
qda_test <- qda(activity ~ ., data = Full_Train[, -c(1:3,5,21)])
predicted <- predict(qda_test, Full_Test)
summary(qda_test)
Predicted <- as.vector(predicted$class)


Full_Test$activity <- as.vector(predicted$class)
Full_Test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions_qda.csv")

file.show("test_set_predictions.csv")

