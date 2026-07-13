library(tidyverse) # Collection of packages for data wrangling + plotting etc. 
library(tidymodels) # Collection of packages for machine learning
library(themis) # Balancing data
library(doParallel) # Parallel processing
library(ranger) # Quantile regression forests
library(furrr) # Parallel processing (with map)
library(tictoc)


## Set working directory

setwd("...")


## Load regression models

load("./new_template_models.RData")


## Split data into 60% training and 40% test set
## Re-run this so the materials and species are present in the dataset
set.seed(123)

data_split_reg <- initial_split(data_analysis, strata=strata,
                                prop=0.6)


data_train_reg <- training(data_split_reg)
data_test_reg <- testing(data_split_reg)




## Model performance (performance metrics)

# Function to extract cross validation results
cross_val_performance <- function(x, y) {
  collect_metrics(x) %>% inner_join(y)
}

# Function to extract cross validation performance metrics
predictions <- function(x, y) {
  collect_predictions(x) %>% inner_join(y) %>% select(id, .pred, .row, endpoint_value)
}

# Combine all individual models into object
models_reg <- tibble(model=list(bart_final, cubist_final, knn_final, linear_final, mars_final, nnet_final, rf_final, svm_final, xgboost_final), 
                     tune=list(bart_tune, cubist_tune, knn_tune, linear_tune, mars_tune, nnet_tune, rf_tune, svm_tune, xgboost_tune), 
                     parms=list(bart_best_parms, cubist_best_parms, knn_best_parms, linear_best_parms, mars_best_parms, nnet_best_parms,
                                rf_best_parms, svm_best_parms, xgboost_best_parms),
                     recipe=list(bart_final_rec, cubist_final_rec, knn_final_rec, linear_final_rec, mars_final_rec, nnet_final_rec, 
                                 rf_final_rec, svm_final_rec, xgboost_final_rec),
                     algorithm=c("bart", "cubist", "knn", "linear", "mars", "nnet", "rf", "svm", "xgboost"))

# Extract performance metrics for each model
models_reg <- models_reg %>% mutate(ext_validation=map(model, collect_metrics),
                                    int_validation=map2(.x=tune, .y=parms, cross_val_performance),
                                    ext_pred=map(model, collect_predictions),
                                    int_pred=map2(.x=tune, .y=parms, predictions))



################################################################################################################################


## Remove all objects in environments except the tibble with all models

data_reg <- data_analysis

rm(list=setdiff(ls(), c("data_split_reg", "data_train_reg", "data_test_reg", "models_reg", "data_reg")))






################################################################################################################################

## Load classification models

load("./new_template_classification_models.RData")



## Split data into 60% training and 40% test set
## Re-run this so the materials and species are present in the dataset
set.seed(123)

data_split_class <- initial_split(data_analysis, strata=strata,
                                  prop=0.6)


data_train_class <- training(data_split_class)
data_test_class <- testing(data_split_class)



## Model performance (performance metrics)

# Function to extract cross validation results
cross_val_performance <- function(x, y) {
  collect_metrics(x) %>% inner_join(y)
}

# Function to extract cross validation performance metrics
predictions <- function(x, y) {
  collect_predictions(x) %>% inner_join(y) %>% select(id, .row, .pred_class, `.pred_very toxic`, .pred_toxic,
                                                      `.pred_not harmful`, endpoint_category)
}

# Combine all individual models into object
models_class <- tibble(model=list(knn_final, multinom_final, nnet_final, rf_final, svm_final, xgboost_final), 
                       tune=list(knn_tune, multinom_tune, nnet_tune, rf_tune, svm_tune, xgboost_tune), 
                       parms=list(knn_best_parms, multinom_best_parms, nnet_best_parms, 
                                  rf_best_parms, svm_best_parms, xgboost_best_parms),
                       recipe=list(knn_final_rec, multinom_final_rec, nnet_final_rec, rf_final_rec, 
                                   svm_final_rec, xgboost_final_rec),
                       algorithm=c("knn", "multinom", "nnet", "rf", "svm", "xgboost"))

# Extract performance metrics for each model
models_class <- models_class %>% mutate(ext_validation=map(model, collect_metrics),
                                        int_validation=map2(.x=tune, .y=parms, cross_val_performance),
                                        ext_pred=map(model, collect_predictions),
                                        int_pred=map2(.x=tune, .y=parms, predictions))





################################################################################################################################


## Remove all objects in environments except the tibble with all models

data_class <- data_analysis

rm(list=setdiff(ls(), c("data_split_reg", "data_train_reg", "data_test_reg", "models_reg", "data_reg",
                        "data_split_class", "data_train_class", "data_test_class", "models_class", "data_class")))






################################################################################################################################


## Stacked classification models ##



## Collect model predictions to stack (the cross-validation results) and create training set for stacking

stack_train_class <- models_class %>% 
  select(algorithm, int_pred) %>% 
  unnest(int_pred) %>% 
  select(-c(.pred_class)) %>% 
  pivot_wider(names_from=algorithm, values_from=c(`.pred_very toxic`,
                                                  `.pred_toxic`,
                                                  `.pred_not harmful`),
              names_glue = "{algorithm}_{.value}") %>%
  arrange(.row) %>%
  select(-c(id, .row))

# Rename column names
# Order columns alphabetically
colnames(stack_train_class) <- colnames(stack_train_class) %>% str_remove("\\.pred_") %>% str_replace("\\s", "_")
stack_train_class <- stack_train_class %>% select(order(colnames(stack_train_class)))




## Collect model predictions to stack and create testing set for stacking

stack_test_class <- models_class %>% 
  select(algorithm, ext_pred) %>% 
  unnest(ext_pred) %>% 
  select(-c(.pred_class)) %>% 
  pivot_wider(names_from=algorithm, values_from=c(`.pred_very toxic`,
                                                  `.pred_toxic`,
                                                  `.pred_not harmful`),
              names_glue = "{algorithm}_{.value}") %>% 
  select(-c(id, .row, .config))


# Rename column names
# Order columns alphabetically
colnames(stack_test_class) <- colnames(stack_test_class) %>% str_remove("\\.pred_") %>% str_replace("\\s", "_")
stack_test_class <- stack_test_class %>% select(order(colnames(stack_test_class)))



## Build stacked models


# Create split object from training and testing data
stack_split_class <- make_splits(stack_train_class, assessment=stack_test_class)

# Create ensemble data

stack_train_class <- training(stack_split_class)
stack_test_class <- testing(stack_split_class)


# 10-fold cross validation of training set
set.seed(345)
stack_folds_class <- vfold_cv(stack_train_class, v=10, repeats=1)


stack_rec_class <- recipe(endpoint_category ~ ., data=stack_train_class)

# Specify model and engine
stack_spec_class <- multinom_reg(penalty=tune(), mixture=tune()) %>% 
  set_mode("classification") %>%
  set_engine("glmnet")

# Create workflow with pre-processing recipe and model specification
stack_workflow_class <- workflow() %>% 
  add_recipe(stack_rec_class) %>% 
  add_model(stack_spec_class) 


# Set hyperparameter ranges
param <- extract_parameter_set_dials(stack_workflow_class) %>% 
  update(penalty = penalty(c(0, 2), NULL),
         mixture = mixture(c(0, 1), NULL))


# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Set metrics for evaluation

metrics <- metric_set(precision, accuracy, roc_auc, sens, spec, f_meas, bal_accuracy,
                      mcc, mn_log_loss, gain_capture)

# Initial tuning
set.seed(456)
stack_tune_class <- tune_grid(stack_workflow_class, 
                              resamples=stack_folds_class, 
                              grid=100,
                              control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                              metrics=metrics,
                              param_info=param)

stopCluster(cl) # Shut down cluster


stack_tune_class %>% show_best("roc_auc", n=10)
stack_tune_class %>% show_best("accuracy", n=10)
stack_tune_class %>% show_best("bal_accuracy", n=10)
stack_tune_class %>% show_best("mcc", n=10)


stack_best_parms_class <- stack_tune_class %>% select_best("mcc") # Select best parameters based on MCC


## Final model


# Update hyperparameters in recipe and model 
stack_final_rec_class <- finalize_recipe(stack_rec_class, stack_best_parms_class)
stack_final_model_class <- finalize_model(stack_spec_class, stack_best_parms_class)

# Create new workflow for final model
stack_final_wf_class <- workflow() %>% 
  add_recipe(stack_final_rec_class) %>%
  add_model(stack_final_model_class)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
stack_final_class <- stack_final_wf_class %>%
  last_fit(stack_split_class,
           metrics=metrics)


collect_metrics(stack_final_class) # Metrics of stacked model (test set)
models_class %>% filter(algorithm=="xgboost") %>% select(ext_validation) %>% unnest(ext_validation) # Metrics of XGBoost model






################################################################################################################################


## Stacked regression model ##


# At the moment tidymodels doesn't fully support quantile regression forests (assesses models as if they're a normal
# randomforest instead of quantile regression forest)




# Collect model predictions to stack (the cross-validation results) and create training set for stacking

stack_train_reg <- models_reg %>% 
  select(algorithm, int_pred) %>% 
  unnest(int_pred) %>%  
  pivot_wider(names_from=algorithm, values_from=.pred) %>% 
  arrange(.row) %>%
  select(-c(id, .row))



# Collect model predictions to stack and create testing set for stacking

stack_test_reg <- models_reg %>% 
  select(algorithm, ext_pred) %>% 
  unnest(ext_pred) %>%  
  pivot_wider(names_from=algorithm, values_from=.pred) %>% 
  select(-c(id, .config, .row))


# Predict class probabilities using XGBoost model (best classification model) for training and test set
# Rename column names
class_train_pred <- models_class %>% filter(algorithm=="xgboost") %>% pull(model) %>% pluck(1) %>% extract_workflow() %>% 
  predict(data_train_reg, type="prob")

colnames(class_train_pred) <- colnames(class_train_pred) %>% str_replace("\\.pred", "prob") %>% str_replace("\\s", "_")


class_test_pred <- models_class %>% filter(algorithm=="xgboost") %>% pull(model) %>% pluck(1) %>% extract_workflow() %>% 
  predict(data_test_reg, type="prob")

colnames(class_test_pred) <- colnames(class_test_pred) %>% str_replace("\\.pred", "prob") %>% str_replace("\\s", "_")

# Combine classification probabilities with regression stack dataset
stack_train_reg <- stack_train_reg %>% bind_cols(class_train_pred)
stack_test_reg <- stack_test_reg %>% bind_cols(class_test_pred)

# Remove linear model results (to improve results)
stack_train_reg <- stack_train_reg %>% select(-c(linear))
stack_test_reg <- stack_test_reg %>% select(-c(linear))



## Build stacked models



# Create split object from training and testing data
stack_split_reg <- make_splits(stack_train_reg, assessment=stack_test_reg)

# Create ensemble data

stack_train_reg <- training(stack_split_reg)
stack_test_reg <- testing(stack_split_reg)


# 10-fold cross validation of training set
set.seed(345)
stack_folds_reg <- vfold_cv(stack_train_reg, v=10, repeats=1)



################################################################################################################################


## Quantile regression forest tuning ##


## Function to calculate standard error
std.error <- function(x) sd(x)/sqrt(length(x))


## Function for cross validation and extracting performance metrics + predictions
cross_val <- function(i, workflow, data_train) {
  
  # Get row ids from split object
  train_ids <- stack_folds_reg %>% get_rsplit(i) %>% 
    pluck(2)
  
  # Split into train and test set based on row ids
  train <- data_train %>% mutate(id=1:nrow(data_train)) %>% filter(id %in% train_ids)
  test <- data_train %>% mutate(id=1:nrow(data_train)) %>% filter(!id %in% train_ids)
  
  set.seed(456) # Set seed for reproducibility
  
  fit <- workflow %>% fit(data=train) # Fit model
  
  preds <- fit %>% extract_fit_engine() %>% predict(test, type="quantiles") # Predict using quantile regression forest
  
  # Put predictions into tibble
  predictions <- preds$predictions %>% as_tibble() %>% 
    rename(.pred=2) %>% select(.pred) %>% bind_cols(test %>% select(id)) %>% 
    rename(row_id=id)
  
  # Calculate performance metrics and put into tibble
  metrics <- preds$predictions %>% as_tibble() %>% 
    rename(.pred=2) %>% select(.pred) %>%
    bind_cols(test %>% select(endpoint_value)) %>%
    metrics(truth=endpoint_value,
            estimate=.pred)
  
  # Put all results into tibble
  res <- tibble(fold=i,
                trees=fit$fit$fit$fit$num.trees,
                mtry=fit$fit$fit$fit$mtry,
                min_n=fit$fit$fit$fit$min.node.size,
                predictions=list(predictions),
                metrics=list(metrics))
  
  
  
  return(res)
  
}


## Function for grid search
tune_grid <- function(parameters) {
  
  mtry <- parameters$mtry
  trees <- parameters$trees
  min_n <- parameters$min_n
  
  spec <- rand_forest(mtry=mtry,
                      trees=trees,
                      min_n=min_n) %>% 
    set_mode("regression") %>%
    set_engine("ranger", quantreg=TRUE, keep.inbag=TRUE)
  
  # Normalize (center and scale) the data
  recipe <- recipe(endpoint_value ~ ., data=stack_train_reg) %>% 
    step_normalize(all_predictors())
  
  # Create workflow with pre-processing recipe and model specification
  workflow <- workflow() %>% 
    add_recipe(recipe) %>% 
    add_model(spec) 
  
  # Apply recipe (normalization) to training data
  data_train <- recipe %>% prep() %>% juice()
  
  # Fit models to all folds
  cross_val_results <- map(.x=c(1:10), .f=\(x) cross_val(x, workflow, data_train))
  
  
  return(cross_val_results)
  
}




# Set performance metrics to assess
metrics <- metric_set(rmse, mae, rsq, ccc, huber_loss, rpiq, rpd, iic)

# Set range for parameters
# Initial tuning with grid_random and then grid_regular is used (grid_random doesn't allow to choose just one value within the range)
parm <- grid_random(mtry(range=c(3, 4), NULL), 
                    trees(rang=c(1000, 1500), NULL),
                    min_n(range=c(1, 15), NULL),
                    size=50) 

parm <- grid_regular(mtry(range=c(4, 4), NULL), 
                     trees(rang=c(1100, 1300), NULL),
                     min_n(range=c(6, 6), NULL),
                     levels=50)

parm


# Set backend for parallel processing
plan(multisession, workers = 4)


# Loop for cross-validation on entire parameter range
tic()
quant_tune <- future_map(.x=c(1:nrow(parm)), .f=\(x) tune_grid(parm[x, ]),
                         .options=furrr_options(packages=c("tidymodels"), 
                                                seed=TRUE))
toc()

# RMSE for all parameters
bind_rows(quant_tune) %>% select(-predictions) %>% unnest(metrics) %>%
  filter(.metric=="rmse") %>%
  group_by(trees, mtry, min_n) %>% 
  summarise(mean=mean(.estimate),
            std_error=std.error(.estimate)) 


bind_rows(quant_tune) %>% select(-predictions) %>% unnest(metrics) %>%
  filter(.metric=="rmse") %>%
  group_by(trees, mtry, min_n) %>% 
  summarise(mean=mean(.estimate),
            std_error=std.error(.estimate)) %>% 
  pivot_longer(trees:min_n, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")



# Select best parameters to build final model (based on RMSE)
stack_best_parms_reg <- bind_rows(quant_tune) %>% select(-predictions) %>% unnest(metrics) %>%
  filter(.metric=="rmse") %>%
  group_by(trees, mtry, min_n) %>% 
  summarise(mean=mean(.estimate),
            std_error=std.error(.estimate)) %>% 
  ungroup() %>% arrange(mean) %>% slice(1) %>% select(mtry, trees, min_n)


bind_rows(quant_tune) %>% select(-metrics) %>% inner_join(stack_best_parms_reg) %>% unnest(predictions) %>% 
  arrange(row_id) %>% bind_cols(stack_train_reg %>% select(endpoint_value)) %>% 
  ggplot(aes(endpoint_value, .pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + 
  geom_abline(intercept=0.5, linetype=2) + 
  geom_abline(intercept=-0.5, linetype=2)

################################################################################################################################

## Fit same model with tidymodels using best parameters


# Normalize (center and scale) the data
stack_rec_reg <- recipe(endpoint_value ~ ., data=stack_train_reg) %>% 
  step_normalize(all_predictors())


# Specify model and engine
# Add additional settings for quantile regression forests and for calculation of confidence intervals
stack_spec_reg <- rand_forest(mtry=tune(),
                              trees=tune(),
                              min_n=tune()) %>% 
  set_mode("regression") %>%
  set_engine("ranger", quantreg=TRUE, keep.inbag=TRUE)

# Create workflow with pre-processing recipe and model specification
stack_workflow_reg <- workflow() %>% 
  add_recipe(stack_rec_reg) %>% 
  add_model(stack_spec_reg) 


# Set hyperparameter range for mtry (otherwise there are errors)
param <- extract_parameter_set_dials(stack_workflow_reg) %>% 
  update(trees=trees(c(stack_best_parms_reg$trees, stack_best_parms_reg$trees)),
         mtry=mtry(c(stack_best_parms_reg$mtry, stack_best_parms_reg$mtry), NULL),
         min_n=min_n(c(stack_best_parms_reg$min_n, stack_best_parms_reg$min_n), NULL))


# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Set metrics for evaluation

metrics <- metric_set(rmse, mae, rsq, ccc, huber_loss, rpiq, rpd, iic)

# Initial tuning
set.seed(456)
stack_tune_reg <- tune::tune_grid(stack_workflow_reg, 
                                  resamples=stack_folds_reg, 
                                  grid=1,
                                  control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                                  metrics=metrics,
                                  param_info=param)

stopCluster(cl) # Shut down cluster

stack_best_parms_reg <- stack_tune_reg %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
stack_final_rec_reg <- finalize_recipe(stack_rec_reg, stack_best_parms_reg)
stack_final_model_reg <- finalize_model(stack_spec_reg, stack_best_parms_reg)

# Create new workflow for final model
stack_final_wf_reg <- workflow() %>% 
  add_recipe(stack_final_rec_reg) %>%
  add_model(stack_final_model_reg)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
stack_final_reg <- stack_final_wf_reg %>%
  last_fit(stack_split_reg,
           metrics=metrics)


# Modify tune results to add the predictions and performance metrics from the quantile regression forest model
# Reconstruct object
stack_tune_reg <- stack_tune_reg %>% as_tibble() %>% 
  bind_cols(bind_rows(quant_tune) %>% 
              inner_join(stack_best_parms_reg) %>% 
              select(predictions, metrics)) %>% 
  rename(.predictions_quantreg=predictions,
         .metrics_quantreg=metrics) %>% 
  rset_reconstruct(., stack_tune_reg)




# Get predictions on test set for quantile regression forest model
quant_reg_preds <- stack_final_reg %>% 
  extract_fit_engine() %>% predict(stack_final_reg %>% extract_recipe() %>% bake(stack_test_reg), 
                                   quantiles=c(0.05, 0.5, 0.95), type="quantiles")


collect_metrics(stack_final_reg) # Metrics of model (test set)
quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% 
  bind_cols(data_test_reg) %>% metrics(., truth=endpoint_value, estimate=.pred) # Metrics of quantile regression forest
models_reg %>% filter(algorithm=="xgboost") %>% select(ext_validation) %>% unnest(ext_validation) # Metrics of XGBoost
models_reg %>% filter(algorithm=="rf") %>% select(ext_validation) %>% unnest(ext_validation) # Metrics of Randomforest


collect_metrics(stack_tune_reg) %>% inner_join(stack_best_parms_reg) %>% select(.metric, mean, n) # Metrics of model (training set)
as_tibble(stack_tune_reg) %>% select(id, .metrics_quantreg) %>% unnest(cols=.metrics_quantreg) %>% 
  group_by(.metric) %>% summarise(.estimate=mean(.estimate),
                                  n=n())
models_reg %>% filter(algorithm=="xgboost") %>% 
  select(int_validation) %>% 
  unnest(int_validation) %>% 
  select(.metric, mean, n) # Metrics of XGBoost
models_reg %>% filter(algorithm=="rf") %>% 
  select(int_validation) %>% 
  unnest(int_validation) %>% 
  select(.metric, mean, n) # Metrics of Randomforest

# Actual vs predicted plot (randomforest)
collect_predictions(stack_final_reg) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)


collect_predictions(stack_final_reg) %>% mutate(offset=endpoint_value-.pred,
                                          freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(stack_test_reg))

# Actual vs predicted (quantile regression forest with prediction intervals)

quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% bind_cols(data_test_reg) %>% 
  ggplot(aes(endpoint_value, .pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)


quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% bind_cols(data_test_reg) %>% 
  ggplot(aes(endpoint_value, .pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  geom_errorbar(aes(ymin=.pred_lower,
                    ymax=.pred_upper)) + 
  facet_wrap(.~material)

quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% bind_cols(data_test_reg) %>% 
  ggplot(aes(endpoint_value, .pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  geom_errorbar(aes(ymin=.pred_lower,
                    ymax=.pred_upper)) + 
  facet_wrap(.~species_group)

quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% 
  bind_cols(data_test_reg) %>% 
  mutate(offset=endpoint_value-.pred,
         freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test_reg))


## Percent coverage
quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% 
  bind_cols(data_test_reg %>% select(endpoint_value)) %>% 
  mutate(range=case_when(endpoint_value <= .pred_upper & endpoint_value >= .pred_lower ~ "inside",
                        .default="outside")) %>% 
  group_by(range) %>% summarise(count=n()) %>% 
  mutate(perc=(count/(nrow(quant_reg_preds$predictions)))*100)


# Actual vs predicted plot (XGBoost)
models_reg %>% filter(algorithm=="xgboost") %>% select(model) %>% pull(model) %>% pluck(1) %>% 
  collect_predictions() %>% 
  ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)


models_reg %>% filter(algorithm=="xgboost") %>% select(model) %>% pull(model) %>% pluck(1) %>% 
  collect_predictions() %>% mutate(offset=endpoint_value-.pred,
                                                freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test_reg))















































