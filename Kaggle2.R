library(tidyverse)
dir("physical-activity-recognition/RawData/Train/", pattern = "^acc")


# Data import -------------------------------------------------------------

# import the activity labels
act_labels <- read_delim("physical-activity-recognition/activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% dplyr::select("X1","X2")

# import the data labels
labels <- read_delim("physical-activity-recognition/RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')

# name activity labels
labels <- labels %>% mutate(activity = act_labels$X2[activity])


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

# Function for extracting time domain features ----------------------------

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
        p1 = mean(X1^2),
        p2 = mean(X2^2),
        p3 = mean(X3^2),
        sd1 = sd(X1),
        sd2 = sd(X2),
        sd3 = sd(X3),
        min1 = min(X1),
        min2 = min(X2),
        min3 = min(X3),
        max1 = max(X1),
        max2 = max(X2),
        max3 = max(X3),
        skew1 = e1071::skewness(X1),
        skew2 = e1071::skewness(X2),
        skew3 = e1071::skewness(X3),
        AR1.1 = cor(X1, lag(X1), use = "pairwise"),
        AR1.2 = cor(X1, lag(X1, n = 2), use = "pairwise"),
        AR2.1 = cor(X2, lag(X2), use = "pairwise"),
        AR2.2 = cor(X2, lag(X2, n = 2), use = "pairwise"),
        AR3.1 = cor(X3, lag(X3), use = "pairwise"),
        AR3.2 = cor(X3, lag(X3, n = 2), use = "pairwise"),
        AR12.1 = cor(X1, lag(X2), use = "pairwise"),
        AR12.2 = cor(X1, lag(X2, n = 2), use = "pairwise"),
        AR13.1 = cor(X1, lag(X3), use = "pairwise"),
        AR13.2 = cor(X1, lag(X3, n = 2), use = "pairwise"),
        AR23.1 = cor(X2, lag(X3), use = "pairwise"),
        AR23.2 = cor(X2, lag(X3, n = 2), use = "pairwise"),
        n=n(),
        entr.1 = entropy(X1),
        entr.2 = entropy(X2),
        entr.3 = entropy(X3),
        spec.1 = spect(X1),
        spec.2 = spect(X2),
        spec.3 = spect(X3)
      ) 
    
    usertimedom 
  }


# Importing data with the function ----------------------------------------

# get all filenames for acceleration and gyroscope train data
filenames_train_acc <- dir("physical-activity-recognition/RawData/Train/", "^acc", full.names = TRUE)
filenames_train_gyro <- dir("physical-activity-recognition/RawData/Train/", "^gyro", full.names = TRUE)

# run data function over all files and store in tibble
acc_train =  filenames_train_acc %>%
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Acc"))

gyro_train = filenames_train_gyro %>% 
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Gyro"))

# join tibbles
data_train <- inner_join(acc_train, gyro_train) %>% 
  filter(activity != "-")


## repeat for test data
filenames_test_acc <- dir("physical-activity-recognition/RawData/Test/", "^acc", full.names = TRUE)
filenames_test_gyro <- dir("physical-activity-recognition/RawData/Test/", "^gyro", full.names = TRUE)

# extract data
acc_test =  filenames_test_acc %>%
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Acc"))

gyro_test = filenames_test_gyro %>% 
  map_dfr(extractTimeDomainFeatures) %>% 
  rename_at(vars(-c(epoch, user_id,exp_id, activity, sample, n)), ~ paste0(., "_Gyro"))

# join tibbles
data_test <- inner_join(acc_test, gyro_test)


# Visualisations ----------------------------------------------------------

data_train %>%
  ggplot(aes(m1_Gyro)) + 
  geom_histogram(bins=40, fill=1, alpha=0.5) + 
  geom_histogram(aes(m2_Gyro), bins=40, fill = 2, alpha=0.5) + 
  geom_histogram(aes(m3_Gyro), bins=40, fill = 4, alpha=0.5) +
  facet_wrap(~activity, scales = "free_y")


# Modelling ---------------------------------------------------------------
library(MASS)
library(caret)

# recode user and experiment id to factor
data_train <- data_train %>% 
  mutate_if(is.character, as.factor)

data_test <- data_test %>% 
  mutate_if(is.character, as.factor)

data_train[,-1]

### caret
qda_fit <- train(activity ~ ., data = data_train[,-1], method = "stepQDA")












##### QDA
qda_test <- qda(activity ~ ., data = data_train[, -1])
predicted <- predict(qda_test, data_test)
summary(qda_test)
Predicted <- as.vector(predicted$class)


Full_Test$activity <- as.vector(predicted$class)
Full_Test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions_qda.csv")

file.show("test_set_predictions.csv")

