#-----------------------------------------------------------------------------#
# Project: Microbial clusters
# Script:  ml_predictor.R
# Purpose: This script is intended to generate a set of ML predictors of immunotherapy outcome
# Author:  Miguel Parra-Martinez
# Date:    2026-02-04
#-----------------------------------------------------------------------------#


library(caret)
library(pROC)


# Definition of function for selecting hyperparameters during crossvalidation:
batch_mean_selector <- function(x, metric, maximize) {
  # x is a data frame where each row is a hyperparameter combination
  # columns contain performance metrics from each fold

  # Calculate mean across folds for each parameter combination
  means <- rowMeans(x[, grep(metric, names(x))], na.rm = TRUE)

  if (maximize) {
    best_idx <- which.max(means)
  } else {
    best_idx <- which.min(means)
  }
  return(best_idx)
}


# Based on bacteria ------------------------------------------------------------

loso_sgb = read_rds("output/data/metagenomics/microbeICB/meta-analysis-nocovariates/sgb.mmuphin.loso.rds")


## Elastic net PA ------------------------------------------------------------
metadata <- coltibble(tse_int.sgb.fil.pretreatment) %>% 
  mutate(Response_fct = as.factor(ifelse(Response == "Responder", "R", "NR")))

# Loop through all the Studies.  to performe meta-analysis. 
Studies <- unique(metadata$Study)

# Storage for results
outer_results <- list()

# Outer CV loop
for (Study_loop in Studies) {
  # Split data into outer train and test sets
  tse_int.sgb.fil.pretreatment_x <- filter_se(
    tse_int.sgb.fil.pretreatment,
    Study != Study_loop
  ) %>%
    transformAssay(method = "relabundance") %>% 
    transformAssay(method = "pa")

  tse_int.sgb.fil.pretreatment_y <- filter_se(
    tse_int.sgb.fil.pretreatment,
    Study == Study_loop
  ) %>%
    transformAssay(method = "relabundance") %>% 
    transformAssay(method = "pa")

  metadata_x = metadata %>%
    filter(Study != Study_loop)

  metadata_y = metadata %>%
    filter(Study == Study_loop)


  # Load LOSO meta-analysis results
  top_sgb = loso_sgb[[Study_loop]]$meta_fits %>%
    arrange(pval) %>%
    slice_head(n = 200) %>%
    pull(feature)


  data_train <- tse_int.sgb.fil.pretreatment_x[top_sgb, ] %>% 
    assay("pa") %>% 
    t()
  
  data_test <-tse_int.sgb.fil.pretreatment_y[top_sgb, ] %>% 
    assay("pa") %>% 
    t()

  train_studies <- unique(coltibble(tse_int.sgb.fil.pretreatment_x)$Study)

  # Custom indices for inner CV
  inner_index <- list()
  for (j in 1:length(train_studies)) {
    # Training indices = all batches except j-th batch
    train_idx <- which(
      coltibble(tse_int.sgb.fil.pretreatment_x)$Study != train_studies[j]
    )
    inner_index[[j]] <- train_idx
  }

  inner_control <- trainControl(
    method = "cv",
    number = length(train_studies),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    selectionFunction = batch_mean_selector,
    returnResamp = "all" # Important: keep all resampling results
  )

  #elastic net
  model <- train(
    x = data_train,
    y = factor(metadata_x$Response_fct),
    method = "glmnet", # Regularized logistic regression (elastic net)
    trControl = inner_control,
    # preProcess = c("center", "scale"),
    tuneGrid = expand.grid(
      alpha = c(0, 0.1, 0.5, 0.7, 0.9, 1), # Favor EN (1) for feature selection
      lambda = 10^seq(-3, 1, length.out = 30) # Wide range of regularization
    ),
    metric = "ROC",
    family = "binomial", 
    verboseIter = TRUE 
  )

  # Get predicted probabilities on outer test set
  pred_probs_train <- predict(model, data_train, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_train <- roc(
    factor(metadata_x$Response_fct),
    pred_probs_train[, "R"]
  )
  auc_value_train <- auc(roc_obj_train)
  # print(paste0("auc train: ", auc_value_train))

  # Get predicted probabilities on outer test set
  pred_probs_test <- predict(model, data_test, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)
  # print(paste0(Study_loop, " auc test: ", auc_value_test))
  pred_class_test <- predict(model, data_test, type = "raw")

  cm <- confusionMatrix(
    pred_class_test,
    factor(metadata_y$Response_fct),
    positive = "R" # Specify which class is "positive"
  )

  # 2. Extract specific metrics
  accuracy <- cm$overall["Accuracy"]
  balanced_acc <- cm$byClass["Balanced Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  f1 <- cm$byClass["F1"]
  ppv <- cm$byClass["Pos Pred Value"]
  npv <- cm$byClass["Neg Pred Value"]

  # 3. AUC (you're already doing this)
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)

  # 4. Find optimal threshold using Youden's Index
  coords_obj <- coords(roc_obj_test, "best", best.method = "youden")
  optimal_threshold <- coords_obj$threshold

  outer_results[[Study_loop]] <- list(
    fold = Study_loop,

    # AUC
    auc_test = as.numeric(auc_value_test),
    auc_train = as.numeric(auc_value_train),

    # Classification metrics at default threshold (0.5)
    accuracy = as.numeric(accuracy),
    balanced_accuracy = as.numeric(balanced_acc),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    f1 = as.numeric(f1),
    ppv = as.numeric(ppv),
    npv = as.numeric(npv),

    # Additional info
    confusion_matrix = cm$table,
    optimal_threshold = coords_obj$threshold,
    n_test = nrow(data_test),
    class_distribution = table(metadata_y$Response_fct),

    # Models and predictions
    model = model,
    roc_obj = roc_obj_test,
    predictions = data.frame(
      true_class = metadata_y$Response_fct,
      pred_prob = pred_probs_test[, "R"],
      pred_class = pred_class_test
    )
  )

  cat(sprintf(
    "\nOuter fold %s:\n  AUC = %.3f\n  F1 = %.3f\n  Bal Acc = %.3f\n  Sens = %.3f\n  Spec = %.3f\n",
    Study_loop,
    as.numeric(auc_value_test),
    as.numeric(f1),
    as.numeric(balanced_acc),
    as.numeric(sensitivity),
    as.numeric(specificity)
  ))
}

# Save results

outer_results %>% 
  write_rds("output/data/metagenomics/microbeICB/ml/EN_pa_200f.rds")

# outer_results <- read_rds("output/data/metagenomics/microbeICB/ml/EN_pa_200f.rds")

sapply( names(outer_results), function(batch){
  outer_results[[batch]]$auc_test
})
sapply( names(outer_results), function(batch){
  outer_results[[batch]]$auc_train
})

# After your loop completes, create a summary
results_EN_pa_df <- do.call(rbind, lapply(outer_results, function(x) {
  data.frame(
    Study = x$fold,
    AUC = x$auc_test,
    F1 = x$f1,
    Balanced_Acc = x$balanced_accuracy,
    Sensitivity = x$sensitivity,
    Specificity = x$specificity,
    N_test = x$n_test
  )
})) %>% 
  mutate_if(is.numeric,~round(., 3))

# Add mean and SD
results_EN_pa_df <- rbind(
  results_EN_pa_df,
  data.frame(
    Study = "Mean ± SD",
    AUC = sprintf("%.3f ± %.3f", mean(results_EN_pa_df$AUC), sd(results_EN_pa_df$AUC)),
    F1 = sprintf("%.3f ± %.3f", mean(results_EN_pa_df$F1, na.rm=TRUE), sd(results_EN_pa_df$F1, na.rm=TRUE)),
    Balanced_Acc = sprintf("%.3f ± %.3f", mean(results_EN_pa_df$Balanced_Acc), sd(results_EN_pa_df$Balanced_Acc)),
    Sensitivity = sprintf("%.3f ± %.3f", mean(results_EN_pa_df$Sensitivity), sd(results_EN_pa_df$Sensitivity)),
    Specificity = sprintf("%.3f ± %.3f", mean(results_EN_pa_df$Specificity), sd(results_EN_pa_df$Specificity)),
    N_test = sum(as.numeric(results_EN_pa_df$N_test[-nrow(results_EN_pa_df)]))
  )
)

print(results_EN_pa_df)



## RF PA ------------------------------------------------------------

outer_results_rf <- list()

# Outer CV loop
for(Study_loop in Studies) {
  
  tse_int.sgb.fil.pretreatment_x <- filter_se(
    tse_int.sgb.fil.pretreatment,
    Study != Study_loop
  ) %>%
    transformAssay(method = "relabundance") %>% 
    transformAssay(method = "pa")

  tse_int.sgb.fil.pretreatment_y <- filter_se(
    tse_int.sgb.fil.pretreatment,
    Study == Study_loop
  ) %>%
    transformAssay(method = "relabundance") %>% 
    transformAssay(method = "pa")

  metadata_x = metadata %>%
    filter(Study != Study_loop)

  metadata_y = metadata %>%
    filter(Study == Study_loop)

  # Load LOSO meta-analysis results
  top_sgb = loso_sgb[[Study_loop]]$meta_fits %>%
    arrange(pval) %>%
    slice_head(n = 200) %>%
    pull(feature)

  # # Use SGB based on clusters.
  # top_sgb = metadata_icbm_hm %>%
  #   filter(cluster_small %in% c(1, 2)) %>%
  #   pull(SGB)


  data_train <- tse_int.sgb.fil.pretreatment_x[top_sgb, ] %>% 
    assay("pa") %>% 
    t()
  
  data_test <-tse_int.sgb.fil.pretreatment_y[top_sgb, ] %>% 
    assay("pa") %>% 
    t()

  train_studies <- unique(coltibble(tse_int.sgb.fil.pretreatment_x)$Study)
  train_studies <- unique(coltibble(tse_int.sgb.fil.pretreatment_x)$Study)
  
  # Custom indices for inner CV
  inner_index <- list()
  for(j in 1:length(train_studies)) {
    # Training indices = all batches except j-th batch
    train_idx <- which(coltibble(tse_int.sgb.fil.pretreatment_x)$Study != train_studies[j])
    inner_index[[j]] <- train_idx
  }

  inner_control <- trainControl(
    method = "cv",
    number = length(train_studies),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    selectionFunction = batch_mean_selector,
    returnResamp = "all"  # Important: keep all resampling results
  )

  n_features <- ncol(data_train)

  #RF
  model  <- train(
  x = data_train,
  y = factor(metadata_x$Response_fct),
  method = "ranger",  # Fast random forest
  trControl = inner_control,
  tuneGrid = expand.grid(
    mtry = unique(round(c(2, 4, sqrt(n_features), n_features/10, n_features/5))),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(5, 10, 20)
  ),
  metric = "ROC"
)

  # Get predicted probabilities on outer test set
  pred_probs_train <- predict(model, data_train, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_train <- roc(
    factor(metadata_x$Response_fct),
    pred_probs_train[, "R"]
  )
  auc_value_train <- auc(roc_obj_train)
  # print(paste0("auc train: ", auc_value_train))

  # Get predicted probabilities on outer test set
  pred_probs_test <- predict(model, data_test, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)
  # print(paste0(Study_loop, " auc test: ", auc_value_test))
  pred_class_test <- predict(model, data_test, type = "raw")

  cm <- confusionMatrix(
    pred_class_test,
    factor(metadata_y$Response_fct),
    positive = "R" # Specify which class is "positive"
  )

  # 2. Extract specific metrics
  accuracy <- cm$overall["Accuracy"]
  balanced_acc <- cm$byClass["Balanced Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  f1 <- cm$byClass["F1"]
  ppv <- cm$byClass["Pos Pred Value"]
  npv <- cm$byClass["Neg Pred Value"]

  # 3. AUC (you're already doing this)
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)

  # 4. Find optimal threshold using Youden's Index
  coords_obj <- coords(roc_obj_test, "best", best.method = "youden")
  optimal_threshold <- coords_obj$threshold

  outer_results_rf[[Study_loop]] <- list(
    fold = Study_loop,

    # AUC
    auc_test = as.numeric(auc_value_test),
    auc_train = as.numeric(auc_value_train),

    # Classification metrics at default threshold (0.5)
    accuracy = as.numeric(accuracy),
    balanced_accuracy = as.numeric(balanced_acc),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    f1 = as.numeric(f1),
    ppv = as.numeric(ppv),
    npv = as.numeric(npv),

    # Additional info
    confusion_matrix = cm$table,
    optimal_threshold = coords_obj$threshold,
    n_test = nrow(data_test),
    class_distribution = table(metadata_y$Response_fct),

    # Models and predictions
    model = model,
    roc_obj = roc_obj_test,
    predictions = data.frame(
      true_class = metadata_y$Response_fct,
      pred_prob = pred_probs_test[, "R"],
      pred_class = pred_class_test
    )
  )
  
  cat(sprintf("\nOuter fold %s:\n  AUC = %.3f\n  F1 = %.3f\n  Bal Acc = %.3f\n  Sens = %.3f\n  Spec = %.3f\n", 
              Study_loop, 
              as.numeric(auc_value_test),
              as.numeric(f1),
              as.numeric(balanced_acc),
              as.numeric(sensitivity),
              as.numeric(specificity)))
              
}


# Save results
outer_results_rf %>% 
  write_rds("output/data/metagenomics/microbeICB/ml/rf_pa_200f.rds")

# After your loop completes, create a summary
results_rf_df <- do.call(rbind, lapply(outer_results_rf, function(x) {
  data.frame(
    Study = x$fold,
    AUC = x$auc_test,
    F1 = x$f1,
    Balanced_Acc = x$balanced_accuracy,
    Sensitivity = x$sensitivity,
    Specificity = x$specificity,
    N_test = x$n_test
  )
})) %>% 
  mutate_if(is.numeric,~round(., 3))

# Add mean and SD
results_rf_df <- rbind(
  results_rf_df,
  data.frame(
    Study = "Mean ± SD",
    AUC = sprintf("%.3f ± %.3f", mean(results_rf_df$AUC), sd(results_rf_df$AUC)),
    F1 = sprintf("%.3f ± %.3f", mean(results_rf_df$F1, na.rm=TRUE), sd(results_rf_df$F1, na.rm=TRUE)),
    Balanced_Acc = sprintf("%.3f ± %.3f", mean(results_rf_df$Balanced_Acc), sd(results_rf_df$Balanced_Acc)),
    Sensitivity = sprintf("%.3f ± %.3f", mean(results_rf_df$Sensitivity), sd(results_rf_df$Sensitivity)),
    Specificity = sprintf("%.3f ± %.3f", mean(results_rf_df$Specificity), sd(results_rf_df$Specificity)),
    N_test = sum(as.numeric(results_rf_df$N_test[-nrow(results_rf_df)]))
  )
)

print(results_rf_df)



# Based on metacyc pathways -----------------------------------------------------

loso_metacyc = read_rds("output/data/metagenomics/microbeICB/meta-analysis-nocovariates/metacyc.mmuphin.loso.rds")

# Load data:

tse_int_metacyc <- read_rds("data/metagenomics/int.metacyc.rda") %>% 
  convertFromPhyloseq() %>% 
  transformAssay(method = "relabundance") %>% 
  transformAssay(method = "clr", pseudocount = T) %>% 
  transformAssay(method = "pa")

metacyc.names = data.table::fread("data/resources/map_metacyc-pwy_name.txt.gz", header = F) %>% 
  rename(ID = V1, name = V2) %>% 
  as_tibble()

# Join names pathways
rowData(tse_int_metacyc) <- rowtibble(tse_int_metacyc) %>%
  rename(ID = Species) %>%
  left_join(metacyc.names) %>%
  DataFrame()

# Filter to only samples that are also in Metaphlan pre-treatment 
tse_int_metacyc_fil <- tse_int_metacyc %>% 
  filter_se(Sample %in% colData(tse_int.sgb.fil.pretreatment)$Sample)

# Add community score and merge Kim batches. 
colData(tse_int_metacyc_fil) <- coldf(tse_int_metacyc_fil) %>%
  inner_join(
    coldf(tse_int.sgb.fil.pretreatment) %>% select(Sample, ratio_abun)
  ) %>%
  mutate(
    Study = ifelse(str_detect(Study, "Kim"), "Kim", Study)
  ) %>%
  mutate(
    dysbiosis_level = case_when(
      ratio_abun < quantile(ratio_abun, 0.4, na.rm = T) ~ "Eubiotic",
      ratio_abun > quantile(ratio_abun, 0.6, na.rm = T) ~ "Dysbiotic",
      .default = "other"
    ),
    .by = Study
  ) %>%
  DataFrame()
colnames(tse_int_metacyc_fil) <- colData(tse_int_metacyc_fil)$Sample


##  Elastic Net PA ------------------------------------------------------------

metadata <- coltibble(tse_int_metacyc_fil) %>% 
  mutate(Response_fct = as.factor(ifelse(Response == "Responder", "R", "NR")))

# Loop through all the Studies.  to performe meta-analysis. 
Studies <- unique(metadata$Study)

# Storage for results
outer_results <- list()

# Create tse object used in the training loop
tse = tse_int_metacyc_fil

# Outer CV loop
for (Study_loop in Studies) {
  # Split data into outer train and test sets
  tse_x <- filter_se(
    tse,
    Study != Study_loop
  )

  tse_y <- filter_se(
    tse,
    Study == Study_loop
  )

  metadata_x = metadata %>%
    filter(Study != Study_loop)

  metadata_y = metadata %>%
    filter(Study == Study_loop)

  # Filter features based on meta-analysis

  # Load LOSO meta-analysis results
  top_metacyc = loso_metacyc[[Study_loop]]$meta_fits %>%
    arrange(pval) %>%
    slice_head(n = 200) %>%
    pull(feature)

  data_train <- tse_x[top_metacyc, ] %>%
    assay("pa") %>%
    t()

  data_test <- tse_y[top_metacyc, ] %>%
    assay("pa") %>%
    t()

  train_studies <- unique(coltibble(tse_x)$Study)

  # Custom indices for inner CV
  inner_index <- list()
  for (j in 1:length(train_studies)) {
    # Training indices = all batches except j-th batch
    train_idx <- which(
      coltibble(tse_x)$Study != train_studies[j]
    )
    inner_index[[j]] <- train_idx
  }

  inner_control <- trainControl(
    method = "cv",
    number = length(train_studies),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    selectionFunction = batch_mean_selector,
    returnResamp = "all" # Important: keep all resampling results
  )

  #elastic net
  model <- train(
    x = data_train,
    y = factor(metadata_x$Response_fct),
    method = "glmnet", # Regularized logistic regression (elastic net)
    trControl = inner_control,
    # preProcess = c("center", "scale"),
    tuneGrid = expand.grid(
      alpha = c(0, 0.1, 0.5, 0.7, 0.9, 1), # Favor EN (1) for feature selection
      lambda = 10^seq(-3, 1, length.out = 30) # Wide range of regularization
    ),
    metric = "ROC",
    family = "binomial",
    verboseIter = TRUE
  )

  # Get predicted probabilities on outer test set
  pred_probs_train <- predict(model, data_train, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_train <- roc(
    factor(metadata_x$Response_fct),
    pred_probs_train[, "R"]
  )
  auc_value_train <- auc(roc_obj_train)
  # print(paste0("auc train: ", auc_value_train))

  # Get predicted probabilities on outer test set
  pred_probs_test <- predict(model, data_test, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)
  # print(paste0(Study_loop, " auc test: ", auc_value_test))
  pred_class_test <- predict(model, data_test, type = "raw")

  cm <- confusionMatrix(
    pred_class_test,
    factor(metadata_y$Response_fct),
    positive = "R" # Specify which class is "positive"
  )

  # 2. Extract specific metrics
  accuracy <- cm$overall["Accuracy"]
  balanced_acc <- cm$byClass["Balanced Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  f1 <- cm$byClass["F1"]
  ppv <- cm$byClass["Pos Pred Value"]
  npv <- cm$byClass["Neg Pred Value"]

  # 3. AUC (you're already doing this)
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)

  # 4. Find optimal threshold using Youden's Index
  coords_obj <- coords(roc_obj_test, "best", best.method = "youden")
  optimal_threshold <- coords_obj$threshold

  outer_results[[Study_loop]] <- list(
    fold = Study_loop,

    # AUC
    auc_test = as.numeric(auc_value_test),
    auc_train = as.numeric(auc_value_train),

    # Classification metrics at default threshold (0.5)
    accuracy = as.numeric(accuracy),
    balanced_accuracy = as.numeric(balanced_acc),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    f1 = as.numeric(f1),
    ppv = as.numeric(ppv),
    npv = as.numeric(npv),

    # Additional info
    confusion_matrix = cm$table,
    optimal_threshold = coords_obj$threshold,
    n_test = nrow(data_test),
    class_distribution = table(metadata_y$Response_fct),

    # Models and predictions
    model = model,
    roc_obj = roc_obj_test,
    predictions = data.frame(
      true_class = metadata_y$Response_fct,
      pred_prob = pred_probs_test[, "R"],
      pred_class = pred_class_test
    )
  )

  cat(sprintf(
    "\nOuter fold %s:\n  AUC = %.3f\n  F1 = %.3f\n  Bal Acc = %.3f\n  Sens = %.3f\n  Spec = %.3f\n",
    Study_loop,
    as.numeric(auc_value_test),
    as.numeric(f1),
    as.numeric(balanced_acc),
    as.numeric(sensitivity),
    as.numeric(specificity)
  ))
}

# Save results

outer_results %>% 
  write_rds("output/data/metagenomics/microbeICB/ml/metacyc_EN_pa_200f.rds")

sapply( names(outer_results), function(batch){
  outer_results[[batch]]$auc_test
})
sapply( names(outer_results), function(batch){
  outer_results[[batch]]$auc_train
})

# After your loop completes, create a summary
results_df <- do.call(rbind, lapply(outer_results, function(x) {
  data.frame(
    Study = x$fold,
    AUC = x$auc_test,
    F1 = x$f1,
    Balanced_Acc = x$balanced_accuracy,
    Sensitivity = x$sensitivity,
    Specificity = x$specificity,
    N_test = x$n_test
  )
})) %>% 
  mutate_if(is.numeric,~round(., 3))

# Add mean and SD
results_df <- rbind(
  results_df,
  data.frame(
    Study = "Mean ± SD",
    AUC = sprintf("%.3f ± %.3f", mean(results_df$AUC), sd(results_df$AUC)),
    F1 = sprintf("%.3f ± %.3f", mean(results_df$F1, na.rm=TRUE), sd(results_df$F1, na.rm=TRUE)),
    Balanced_Acc = sprintf("%.3f ± %.3f", mean(results_df$Balanced_Acc), sd(results_df$Balanced_Acc)),
    Sensitivity = sprintf("%.3f ± %.3f", mean(results_df$Sensitivity), sd(results_df$Sensitivity)),
    Specificity = sprintf("%.3f ± %.3f", mean(results_df$Specificity), sd(results_df$Specificity)),
    N_test = sum(as.numeric(results_df$N_test[-nrow(results_df)]))
  )
)

print(results_df)

## RF PA ------------------------------------------------------------

metadata <- coltibble(tse_int_metacyc_fil) %>% 
  mutate(Response_fct = as.factor(ifelse(Response == "Responder", "R", "NR")))

outer_results_rf <- list()

# Create tse object used in the training loop
tse = tse_int_metacyc_fil

Studies = unique(tse_int_metacyc_fil$Study)

# Outer CV loop
for(Study_loop in Studies) {
  
  tse_x <- filter_se(
    tse,
    Study != Study_loop
  )

  tse_y <- filter_se(
    tse,
    Study == Study_loop
  )

  metadata_x = metadata %>%
    filter(Study != Study_loop)

  metadata_y = metadata %>%
    filter(Study == Study_loop)

  # Load LOSO meta-analysis results
  top_metacyc = loso_metacyc[[Study_loop]]$meta_fits %>%
    arrange(pval) %>%
    slice_head(n = 200) %>%
    pull(feature)

  data_train <- tse_x[top_metacyc, ] %>%
    assay("pa") %>%
    t()

  data_test <- tse_y[top_metacyc, ] %>%
    assay("pa") %>%
    t()

  train_studies <- unique(coltibble(tse_x)$Study)
  train_studies <- unique(coltibble(tse_x)$Study)
  
  # Custom indices for inner CV
  inner_index <- list()
  for(j in 1:length(train_studies)) {
    # Training indices = all batches except j-th batch
    train_idx <- which(coltibble(tse_x)$Study != train_studies[j])
    inner_index[[j]] <- train_idx
  }

  inner_control <- trainControl(
    method = "cv",
    number = length(train_studies),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    selectionFunction = batch_mean_selector,
    returnResamp = "all"  # Important: keep all resampling results
  )

  n_features <- ncol(data_train)

  #RF
  model <- train(
    x = data_train,
    y = factor(metadata_x$Response_fct),
    method = "ranger", # Fast random forest
    trControl = inner_control,
    tuneGrid = expand.grid(
      mtry = unique(round(c(
        2,
        4,
        sqrt(n_features),
        n_features / 10,
        n_features / 5
      ))),
      splitrule = c("gini", "extratrees"),
      min.node.size = c(5, 10, 20)
    ),
    metric = "ROC"
  )

  # Get predicted probabilities on outer test set
  pred_probs_train <- predict(model, data_train, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_train <- roc(
    factor(metadata_x$Response_fct),
    pred_probs_train[, "R"]
  )
  auc_value_train <- auc(roc_obj_train)
  # print(paste0("auc train: ", auc_value_train))

  # Get predicted probabilities on outer test set
  pred_probs_test <- predict(model, data_test, type = "prob")
  # Calculate AUC on outer test set
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)
  # print(paste0(Study_loop, " auc test: ", auc_value_test))
  pred_class_test <- predict(model, data_test, type = "raw")

  cm <- confusionMatrix(
    pred_class_test,
    factor(metadata_y$Response_fct),
    positive = "R" # Specify which class is "positive"
  )

  # 2. Extract specific metrics
  accuracy <- cm$overall["Accuracy"]
  balanced_acc <- cm$byClass["Balanced Accuracy"]
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  f1 <- cm$byClass["F1"]
  ppv <- cm$byClass["Pos Pred Value"]
  npv <- cm$byClass["Neg Pred Value"]

  # 3. AUC (you're already doing this)
  roc_obj_test <- roc(factor(metadata_y$Response_fct), pred_probs_test[, "R"])
  auc_value_test <- auc(roc_obj_test)

  # 4. Find optimal threshold using Youden's Index
  coords_obj <- coords(roc_obj_test, "best", best.method = "youden")
  optimal_threshold <- coords_obj$threshold

  outer_results_rf[[Study_loop]] <- list(
    fold = Study_loop,

    # AUC
    auc_test = as.numeric(auc_value_test),
    auc_train = as.numeric(auc_value_train),

    # Classification metrics at default threshold (0.5)
    accuracy = as.numeric(accuracy),
    balanced_accuracy = as.numeric(balanced_acc),
    sensitivity = as.numeric(sensitivity),
    specificity = as.numeric(specificity),
    f1 = as.numeric(f1),
    ppv = as.numeric(ppv),
    npv = as.numeric(npv),

    # Additional info
    confusion_matrix = cm$table,
    optimal_threshold = coords_obj$threshold,
    n_test = nrow(data_test),
    class_distribution = table(metadata_y$Response_fct),

    # Models and predictions
    model = model,
    roc_obj = roc_obj_test,
    predictions = data.frame(
      true_class = metadata_y$Response_fct,
      pred_prob = pred_probs_test[, "R"],
      pred_class = pred_class_test
    )
  )
  
  cat(sprintf("\nOuter fold %s:\n  AUC = %.3f\n  F1 = %.3f\n  Bal Acc = %.3f\n  Sens = %.3f\n  Spec = %.3f\n", 
              Study_loop, 
              as.numeric(auc_value_test),
              as.numeric(f1),
              as.numeric(balanced_acc),
              as.numeric(sensitivity),
              as.numeric(specificity)))
              
}

# Save results
outer_results_rf %>% 
  write_rds("output/data/metagenomics/microbeICB/ml/metacyc_rf_pa_200f.rds")

# After your loop completes, create a summary
results_rf_df <- do.call(rbind, lapply(outer_results_rf, function(x) {
  data.frame(
    Study = x$fold,
    AUC = x$auc_test,
    F1 = x$f1,
    Balanced_Acc = x$balanced_accuracy,
    Sensitivity = x$sensitivity,
    Specificity = x$specificity,
    N_test = x$n_test
  )
})) %>% 
  mutate_if(is.numeric,~round(., 3))

# Add mean and SD
results_rf_df <- rbind(
  results_rf_df,
  data.frame(
    Study = "Mean ± SD",
    AUC = sprintf("%.3f ± %.3f", mean(results_rf_df$AUC), sd(results_rf_df$AUC)),
    F1 = sprintf("%.3f ± %.3f", mean(results_rf_df$F1, na.rm=TRUE), sd(results_rf_df$F1, na.rm=TRUE)),
    Balanced_Acc = sprintf("%.3f ± %.3f", mean(results_rf_df$Balanced_Acc), sd(results_rf_df$Balanced_Acc)),
    Sensitivity = sprintf("%.3f ± %.3f", mean(results_rf_df$Sensitivity), sd(results_rf_df$Sensitivity)),
    Specificity = sprintf("%.3f ± %.3f", mean(results_rf_df$Specificity), sd(results_rf_df$Specificity)),
    N_test = sum(as.numeric(results_rf_df$N_test[-nrow(results_rf_df)]))
  )
)

print(results_rf_df)
