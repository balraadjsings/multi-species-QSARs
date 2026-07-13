library(tidyverse) # Collection of packages for data wrangling + plotting etc. 
library(tidymodels) # Collection of packages for machine learning
library(themis) # Balancing data
library(vip) # Variable importance
library(doParallel) # Parallel processing
library(doRNG) # Random seed generator during parallel processing
library(tictoc) # Timing script
library(furrr) # Parallel processing for purrr functions



# Lime for local explanation with package lime

# Pdp for partial dependence profiles with packages iml
## Problem with pdp is that it cannot give pdp's per group like DALEX. When done it gives same results for all groups
## Package pdp seems broken and is slower compared to iml. iml and pdp give similar results, go with iml

# Variable importance with package vip or iml. iml seems faster?

# DALEX seems easiest to use but output vs plots make no sense. Also offers no parallel processing

# For model performance show:
# actual vs predicted (cross-validation results + on test set)
# actual vs residuals (cross-validation results + on test set) --> Also density plot of residuals



## Set working directory

setwd("...")


## Load data

load("./new_stack.RData")



## Add stacked classification model to tibble with individual models

# Function to extract cross validation results
cross_val_performance <- function(x, y) {
  collect_metrics(x) %>% inner_join(y)
}

# Function to extract cross validation performance metrics
predictions <- function(x, y) {
  collect_predictions(x) %>% inner_join(y) %>% select(id, .row, .pred_class, `.pred_very toxic`, .pred_toxic,
                                                      `.pred_not harmful`, endpoint_category)
}

# Put stacked model and results into object
stack_class <- tibble(model=list(stack_final_class), 
                      tune=list(stack_tune_class), 
                      parms=list(stack_best_parms_class),
                      recipe=list(stack_final_rec_class),
                      algorithm=c("stack"))

stack_class <- stack_class %>% mutate(ext_validation=map(model, collect_metrics),
                                      int_validation=map2(.x=tune, .y=parms, cross_val_performance),
                                      ext_pred=map(model, collect_predictions),
                                      int_pred=map2(.x=tune, .y=parms, predictions))

# Add stacked model to object with individual models
models_class <- models_class %>% bind_rows(stack_class)



## Add stacked regression model to tibble with individual models

# Function to extract cross validation results
cross_val_performance <- function(x, y) {
  collect_metrics(x) %>% inner_join(y)
}

# Function to extract cross validation performance metrics
predictions <- function(x, y) {
  collect_predictions(x) %>% inner_join(y) %>% select(id, .pred, .row, endpoint_value)
}

# Put stacked model and results into object
stack_reg <- tibble(model=list(stack_final_reg), 
                    tune=list(stack_tune_reg),
                    parms=list(stack_best_parms_reg),
                    recipe=list(stack_final_rec_reg),
                    algorithm=c("stack"))


stack_reg <- stack_reg %>% 
  mutate(ext_validation=list(quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% 
                               bind_cols(data_test_reg) %>% metrics(., truth=endpoint_value, estimate=.pred)),
         int_validation=list(stack_tune_reg %>% as_tibble() %>% select(.metrics_quantreg) %>% 
                               unnest(cols=.metrics_quantreg) %>% group_by(.metric) %>% 
                               summarise(mean=mean(.estimate),
                                         std_err=std.error(.estimate))),
         ext_pred=list(quant_reg_preds$predictions %>% as_tibble() %>% rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% 
                         bind_cols(data_test_reg %>% select(endpoint_value)) %>% select(-c(.pred_lower, .pred_upper))),
         int_pred=list(stack_tune_reg %>% as_tibble() %>% select(.predictions_quantreg) %>% 
                         unnest(cols=.predictions_quantreg) %>% rename(.row=row_id)))


# Add stacked model to object with individual models
models_reg <- models_reg %>% bind_rows(stack_reg)


################################################################################################################################


## Remove all objects in environments except the tibble with all models



rm(list=setdiff(ls(), c("data_split_reg", "data_train_reg", "data_test_reg", "models_reg", "data_reg",
                        "data_split_class", "data_train_class", "data_test_class", "models_class", "data_class",
                        "stack_train_class", "stack_test_class", "stack_train_reg", "stack_test_reg")))



################################################################################################################################


## Variable importance (permutation) (regression models)

## knn works (2 sim: 90-120sec)
## rf works (2 sim: 20sec)
## bart works (2 sim: 140-160sec)
## cubist works (2 sim: 20sec)
## linear works (2 sim: 1-2sec)
## mars works (2 sim: 2-3sec)
## nnet works (2 sim: 1-2sec)
## svm works (2 sim: 5-8sec)
## xgboost works (2 sim: 2-3sec)

var_imp_reg <- function(type) {
# Prepare training set
dat <- models_reg %>% filter(algorithm==type) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% prep() %>% juice()

# Create function for pred_wrapper
# For metric="rmse" to work $.pred is needed
pred <- function(object, newdata) {
  stats::predict(object, newdata)$.pred
}


# Set up parallel processing
# Because parallel processing is used, seed must be set through doRNG because foreach will ignore set.seed()
cl <- makeCluster(4)
registerDoParallel(cl)
registerDoRNG(seed = 123456789)

# Calculate permutation variable importance
# .packages must be passed on to foreach, otherwise it won't work because of issues with foreach and loading packages
var_imp <- models_reg %>% filter(algorithm==type) %>% select(model) %>% pull(model) %>% pluck(1) %>% 
  extract_fit_parsnip() %>% 
  vi_permute(target="endpoint_value",
             metric="rmse",
             train=dat,
             pred_wrapper=pred,
             nsim=1000,
             sample_size=NULL,
             parallel=TRUE,
             .packages="parsnip")

return(var_imp)

stopCluster(cl)

}

var_imp_reg_results <- list()

set.seed(123)

tic()
#var_imp_reg_results[["linear"]] <- var_imp_reg("linear")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["nnet"]] <- var_imp_reg("nnet")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["xgboost"]] <- var_imp_reg("xgboost")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["mars"]] <- var_imp_reg("mars")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["svm"]] <- var_imp_reg("svm")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["rf"]] <- var_imp_reg("rf")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["cubist"]] <- var_imp_reg("cubist")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["knn"]] <- var_imp_reg("knn")
toc()

set.seed(123)

tic()
#var_imp_reg_results[["bart"]] <- var_imp_reg("bart")
toc()

var_imp_reg_results

showConnections()



# Plot variable importance results
ggplot(var_imp_reg_results[["xgboost"]], aes(Importance, fct_reorder(Variable, Importance))) + 
  geom_col(fill="#F8766D") +
  geom_errorbar(aes(xmin=Importance-StDev, xmax=Importance+StDev), col="black", width=0.5, linewidth=0.75, alpha=0.5) +
  labs(x="Importance (permutation)", y="Variable") +
  theme(axis.title = element_text(face="bold"))









## Variable importance (permutation) (classification models)


## discrim works (2 sim: 2-12sec)
## knn works (2 sim: 150-160sec)
## multinom works (2 sim: 2-3sec)
## nnet works (2 sim: 2sec)
## rf works (2 sim: 11-13sec)
## svm works (2 sim: 17-21sec)
## xgboost works (2 sim: 3-5sec)

var_imp_class <- function(type) {

# Prepare training set
dat <- models_class %>% filter(algorithm==type) %>% select(recipe) %>% 
  pull(recipe) %>% pluck(1) %>% prep() %>% juice()

# Create function for pred_wrapper
# For metric="mcc" to work $.pred_class is needed
pred <- function(object, newdata) {
  stats::predict(object, newdata)$.pred_class
}


# Set up parallel processing
# Because parallel processing is used, seed must be set through doRNG because foreach will ignore set.seed()
cl <- makeCluster(4)
registerDoParallel(cl)
registerDoRNG(seed = 123456789)

# Calculate permutation variable importance
# .packages must be passed on to foreach, otherwise it won't work because of issues with foreach and loading packages

var_imp <- models_class %>% filter(algorithm==type) %>% select(model) %>% pull(model) %>% pluck(1) %>%  
  extract_fit_parsnip() %>% 
  vi_permute(target="endpoint_category",
             metric=mcc_vec,
             smaller_is_better = FALSE,
             train=dat,
             pred_wrapper=pred,
             nsim=1000,
             sample_size=NULL,
             parallel=TRUE,
             .packages="parsnip")

return(var_imp)

stopCluster(cl)


}

var_imp_class_results <- list()

set.seed(123)

tic()
#var_imp_class_results[["multinom"]] <- var_imp_class("multinom")
toc()

set.seed(123)

tic()
#var_imp_class_results[["nnet"]] <- var_imp_class("nnet")
toc()

set.seed(123)

tic()
#var_imp_class_results[["xgboost"]] <- var_imp_class("xgboost")
toc()

set.seed(123)

tic()
#var_imp_class_results[["rf"]] <- var_imp_class("rf")
toc()

set.seed(123)

tic()
#var_imp_class_results[["svm"]] <- var_imp_class("svm")
toc()

set.seed(123)

tic()
#var_imp_class_results[["knn"]] <- var_imp_class("knn")
toc()


var_imp_class_results


# Plot variable importance results
ggplot(var_imp_class_results[["multinom"]], aes(Importance, fct_reorder(Variable, Importance))) + 
  geom_col(fill="#F8766D") +
  geom_errorbar(aes(xmin=Importance-StDev, xmax=Importance+StDev), col="black", width=0.5, linewidth=0.75, alpha=0.5) +
  labs(x="Importance (permutation)", y="Variable") +
  theme(axis.title = element_text(face="bold"))




#######################################################################################################

## Remove all objects in environments (because results VIP results and models will be loaded)

rm(list=ls())

#######################################################################################################



#######################################################################################################


## Variable importance XGBoost (regression)

var_imp_reg_results$xgboost %>% select(Variable) %>% print(n=1000)

var_imp_reg_results$xgboost %>% 
  mutate(Group=case_when(Variable %in% c("crystallinity_anatase", "crystallinity_rutile", "surface_area", 
                                         "primary_length", "primary_diameter", "shape_spheroid", "shape_irregular",
                                         "shape_rod", "shape_other", "coating_coated", "coating_uncoated", 
                                         "solubility_group_soluble", "solubility_group_slow", "material_group_metal",
                                         "material_group_metal.oxide") ~ "P-chem properties",
                         Variable %in% c("test_duration", "illumination_light", "water_hardness", "ph", 
                                         "temperature", "nat_org_matter_no", "nat_org_matter_yes", 
                                         "shaking_during_experiment_no", "shaking_during_experiment_yes", 
                                         "dispersion_sonicated_no", "dispersion_sonicated_yes", 
                                         "dispersion_stirred_no", "dispersion_stirred_yes",
                                         "dispersion_shaked_no", "dispersion_shaked_other", 
                                         "dispersion_vortexed_no", "dispersion_vortexed_yes", 
                                         "dispersion_homogenized_no", "dispersion_homogenized_other") ~ "Exposure conditions", 
                         Variable %in% c("MW", "Sv", "Se", "Sp", "Mv", "Me", "GD", "NssO", 
                                         "MRcons", "Vx", "VvdwZAZ", "SAscore") ~ "Molecular descriptors",
                         Variable %in% c("endpoint_mortality", "endpoint_growth.inhibition") ~ "Other",
                         .default="Species traits")) %>% 
  ggplot(aes(Importance, fct_reorder(Variable, Importance), fill=Group)) + 
  geom_col() +
  geom_errorbar(aes(xmin=Importance-StDev, xmax=Importance+StDev), col="black", width=0.5, linewidth=0.75, alpha=0.5)



## Variable importance XGBoost (classification)


var_imp_class_results$xgboost %>% 
  mutate(Group=case_when(Variable %in% c("crystallinity_anatase", "crystallinity_rutile", "surface_area", 
                                         "primary_length", "primary_diameter", "shape_spheroid", "shape_irregular",
                                         "shape_rod", "shape_other", "coating_coated", "coating_uncoated", 
                                         "solubility_group_soluble", "solubility_group_slow", "material_group_metal",
                                         "material_group_metal.oxide") ~ "P-chem properties",
                         Variable %in% c("test_duration", "illumination_light", "water_hardness", "ph", 
                                         "temperature", "nat_org_matter_no", "nat_org_matter_yes", 
                                         "shaking_during_experiment_no", "shaking_during_experiment_yes", 
                                         "dispersion_sonicated_no", "dispersion_sonicated_yes", 
                                         "dispersion_stirred_no", "dispersion_stirred_yes",
                                         "dispersion_shaked_no", "dispersion_shaked_other", 
                                         "dispersion_vortexed_no", "dispersion_vortexed_yes", 
                                         "dispersion_homogenized_no", "dispersion_homogenized_other") ~ "Exposure conditions", 
                         Variable %in% c("MW", "Sv", "Se", "Sp", "Mv", "Me", "GD", "NssO", 
                                         "MRcons", "Vx", "VvdwZAZ", "SAscore") ~ "Molecular descriptors",
                         Variable %in% c("endpoint_mortality", "endpoint_growth.inhibition") ~ "Other",
                         .default="Species traits")) %>% 
  ggplot(aes(Importance, fct_reorder(Variable, Importance), fill=Group)) + 
  geom_col() +
  geom_errorbar(aes(xmin=Importance-StDev, xmax=Importance+StDev), col="black", width=0.5, linewidth=0.75, alpha=0.5)


# Calculate ranking (regression models)


var_imp_ranking <- list()

for (i in 1:length(var_imp_reg_results)) {
  
  ## Arrange importance score from high to low and give ranking based on that
  ## Rename column names
  ## Add variable group
  
  var_imp_ranking[[i]] <- var_imp_reg_results[[i]] %>% arrange(desc(Importance)) %>% mutate(Position=seq(1:93)) %>% 
    rename(!! paste(names(var_imp_reg_results)[i], "_rank", sep="") := Position) %>% 
    select(-c(Importance, StDev)) %>% arrange(Variable) %>% 
    mutate(Group=case_when(Variable %in% c("crystallinity_anatase", "crystallinity_rutile", "surface_area", 
                                           "primary_length", "primary_diameter", "shape_spheroid", "shape_irregular",
                                           "shape_rod", "shape_other", "coating_coated", "coating_uncoated", 
                                           "solubility_group_soluble", "solubility_group_slow", "material_group_metal",
                                           "material_group_metal.oxide") ~ "P-chem properties",
                           Variable %in% c("test_duration", "illumination_light", "water_hardness", "ph", 
                                           "temperature", "nat_org_matter_no", "nat_org_matter_yes", 
                                           "shaking_during_experiment_no", "shaking_during_experiment_yes", 
                                           "dispersion_sonicated_no", "dispersion_sonicated_yes", 
                                           "dispersion_stirred_no", "dispersion_stirred_yes",
                                           "dispersion_shaked_no", "dispersion_shaked_other", 
                                           "dispersion_vortexed_no", "dispersion_vortexed_yes", 
                                           "dispersion_homogenized_no", "dispersion_homogenized_other") ~ "Exposure conditions", 
                           Variable %in% c("MW", "Sv", "Se", "Sp", "Mv", "Me", "GD", "NssO", 
                                           "MRcons", "Vx", "VvdwZAZ", "SAscore") ~ "Molecular descriptors",
                           Variable %in% c("endpoint_mortality", "endpoint_growth.inhibition") ~ "Other",
                           .default="Species traits"))
  
  
}

## Bind everything together

var_imp_ranking <- bind_cols(var_imp_ranking) %>% rename(variable=`Variable...1`,
                                                         group=`Group...3`) %>% select(!starts_with("Variable...")) %>% 
  select(!starts_with("Group..."))


## Relocate columns from best to worst model

var_imp_ranking <- var_imp_ranking %>% 
  select(variable, group, xgboost_rank, rf_rank, cubist_rank, knn_rank, bart_rank, svm_rank, 
         nnet_rank, mars_rank, linear_rank)


## Categorize the ranking into groups
## 4 categories are created based by evenly dividing the 53 variables into 4 parts

var_imp_category <- var_imp_ranking %>% select(variable, group, xgboost_rank, rf_rank, cubist_rank, knn_rank, bart_rank, svm_rank, 
                                               nnet_rank, mars_rank, linear_rank) %>% 
  mutate(xgboost_rank=case_when(xgboost_rank <= 14 ~ "very important", 
                                xgboost_rank > 14 & xgboost_rank <= 28 ~ "important", 
                                xgboost_rank > 28 & xgboost_rank <= 42 ~ "less important", 
                                xgboost_rank > 42 ~ "least important"), 
         rf_rank=case_when(rf_rank <= 14 ~ "very important", 
                           rf_rank > 14 & rf_rank <= 28 ~ "important", 
                           rf_rank > 28 & rf_rank <= 42 ~ "less important", 
                           rf_rank > 42 ~ "least important"),
         cubist_rank=case_when(cubist_rank <= 14 ~ "very important", 
                               cubist_rank > 14 & cubist_rank <= 28 ~ "important", 
                               cubist_rank > 28 & cubist_rank <= 42 ~ "less important", 
                               cubist_rank > 42 ~ "least important"),
         knn_rank=case_when(knn_rank <= 14 ~ "very important", 
                            knn_rank > 14 & knn_rank <= 28 ~ "important", 
                            knn_rank > 28 & knn_rank <= 42 ~ "less important", 
                            knn_rank > 42 ~ "least important"),
         bart_rank=case_when(bart_rank <= 14 ~ "very important", 
                             bart_rank > 14 & bart_rank <= 28 ~ "important", 
                             bart_rank > 28 & bart_rank <= 42 ~ "less important", 
                             bart_rank > 42 ~ "least important"),
         svm_rank=case_when(svm_rank <= 14 ~ "very important", 
                            svm_rank > 14 & svm_rank <= 28 ~ "important", 
                            svm_rank > 28 & svm_rank <= 42 ~ "less important", 
                            svm_rank > 42 ~ "least important"),
         nnet_rank=case_when(nnet_rank <= 14 ~ "very important", 
                             nnet_rank > 14 & nnet_rank <= 28 ~ "important", 
                             nnet_rank > 28 & nnet_rank <= 42 ~ "less important", 
                             nnet_rank > 42 ~ "least important"),
         mars_rank=case_when(mars_rank <= 14 ~ "very important", 
                             mars_rank > 14 & mars_rank <= 28 ~ "important", 
                             mars_rank > 28 & mars_rank <= 42 ~ "less important", 
                             mars_rank > 42 ~ "least important"),
         linear_rank=case_when(linear_rank <= 14 ~ "very important", 
                               linear_rank > 14 & linear_rank <= 28 ~ "important", 
                               linear_rank > 28 & linear_rank <= 42 ~ "less important", 
                               linear_rank > 42 ~ "least important"))

## Consensus of the importance groups

var_imp_category_consensus <- var_imp_category %>% mutate(very_important = rowSums(var_imp_category=="very important"), 
                                                          important = rowSums(var_imp_category=="important"), 
                                                          less_important = rowSums(var_imp_category=="less important"),
                                                          least_important = rowSums(var_imp_category=="least important")) %>% 
  select(variable, group, very_important:least_important)


var_imp_category_consensus %>% mutate(consensus = case_when(very_important > 4 ~ "very important", 
                                                            important > 4 ~ "important",
                                                            less_important > 4 ~ "less important",
                                                            least_important > 4 ~ "least important",
                                                            .default="mixed"))











## Calculate ranking (classification models)


var_imp_ranking_class <- list()

for (i in 1:length(var_imp_class_results)) {
  
  ## Arrange importance score from high to low and give ranking based on that
  ## Rename column names
  ## Add variable group
  
  var_imp_ranking_class[[i]] <- var_imp_class_results[[i]] %>% arrange(desc(Importance)) %>% mutate(Position=seq(1:97)) %>% 
    rename(!! paste(names(var_imp_class_results)[i], "_rank", sep="") := Position) %>% 
    select(-c(Importance, StDev)) %>% arrange(Variable) %>% 
    mutate(Group=case_when(Variable %in% c("crystallinity_anatase", "crystallinity_rutile", "surface_area", 
                                           "primary_length", "primary_diameter", "shape_spheroid", "shape_irregular",
                                           "shape_rod", "shape_other", "coating_coated", "coating_uncoated", 
                                           "solubility_group_soluble", "solubility_group_slow", "material_group_metal",
                                           "material_group_metal.oxide") ~ "P-chem properties",
                           Variable %in% c("test_duration", "illumination_light", "water_hardness", "ph", 
                                           "temperature", "nat_org_matter_no", "nat_org_matter_yes", 
                                           "shaking_during_experiment_no", "shaking_during_experiment_yes", 
                                           "dispersion_sonicated_no", "dispersion_sonicated_yes", 
                                           "dispersion_stirred_no", "dispersion_stirred_yes",
                                           "dispersion_shaked_no", "dispersion_shaked_other", 
                                           "dispersion_vortexed_no", "dispersion_vortexed_yes", 
                                           "dispersion_homogenized_no", "dispersion_homogenized_other") ~ "Exposure conditions", 
                           Variable %in% c("MW", "Sv", "Se", "Sp", "Mv", "Me", "GD", "NssO", 
                                           "MRcons", "Vx", "VvdwZAZ", "SAscore") ~ "Molecular descriptors",
                           Variable %in% c("endpoint_mortality", "endpoint_growth.inhibition") ~ "Other",
                           .default="Species traits"))
  
  
}

## Bind everything together

var_imp_ranking_class <- bind_cols(var_imp_ranking_class) %>% rename(variable=`Variable...1`,
                                                                     group=`Group...3`) %>% select(!starts_with("Variable...")) %>% 
  select(!starts_with("Group..."))


## Relocate columns from best to worst model

var_imp_ranking_class <- var_imp_ranking_class %>% 
  select(variable, group, xgboost_rank, rf_rank, knn_rank, svm_rank, 
         nnet_rank, multinom_rank)




## Categorize the ranking into groups
## 4 categories are created based by evenly dividing the 53 variables into 4 parts

var_imp_category_class <- var_imp_ranking_class %>% 
  select(variable, group, xgboost_rank, rf_rank, knn_rank, svm_rank, 
         nnet_rank, multinom_rank) %>% 
  mutate(xgboost_rank=case_when(xgboost_rank <= 14 ~ "very important", 
                                xgboost_rank > 14 & xgboost_rank <= 28 ~ "important", 
                                xgboost_rank > 28 & xgboost_rank <= 42 ~ "less important", 
                                xgboost_rank > 42 ~ "least important"), 
         rf_rank=case_when(rf_rank <= 14 ~ "very important", 
                           rf_rank > 14 & rf_rank <= 28 ~ "important", 
                           rf_rank > 28 & rf_rank <= 42 ~ "less important", 
                           rf_rank > 42 ~ "least important"),
         knn_rank=case_when(knn_rank <= 14 ~ "very important", 
                            knn_rank > 14 & knn_rank <= 28 ~ "important", 
                            knn_rank > 28 & knn_rank <= 42 ~ "less important", 
                            knn_rank > 42 ~ "least important"),
         svm_rank=case_when(svm_rank <= 14 ~ "very important", 
                            svm_rank > 14 & svm_rank <= 28 ~ "important", 
                            svm_rank > 28 & svm_rank <= 42 ~ "less important", 
                            svm_rank > 42 ~ "least important"),
         nnet_rank=case_when(nnet_rank <= 14 ~ "very important", 
                             nnet_rank > 14 & nnet_rank <= 28 ~ "important", 
                             nnet_rank > 28 & nnet_rank <= 42 ~ "less important", 
                             nnet_rank > 42 ~ "least important"),
         multinom_rank=case_when(multinom_rank <= 14 ~ "very important", 
                                 multinom_rank > 14 & multinom_rank <= 28 ~ "important", 
                                 multinom_rank > 28 & multinom_rank <= 42 ~ "less important", 
                                 multinom_rank > 42 ~ "least important"))

## Consensus of the importance groups

var_imp_category_consensus_class <- var_imp_category_class %>% mutate(very_important = rowSums(var_imp_category_class=="very important"), 
                                                                      important = rowSums(var_imp_category_class=="important"), 
                                                                      less_important = rowSums(var_imp_category_class=="less important"),
                                                                      least_important = rowSums(var_imp_category_class=="least important")) %>% 
  select(variable, group, very_important:least_important)


var_imp_category_consensus_class %>% mutate(consensus = case_when(very_important > 4 ~ "very important", 
                                                                  important > 4 ~ "important",
                                                                  less_important > 4 ~ "less important",
                                                                  least_important > 4 ~ "least important",
                                                                  .default="mixed"))




#######################################################################################################################


## Applicability domain (based on k-nearest neighbors) ##


# Open functions to calculate applicability domain

source("C:/Users/Suri/Documents/OneDrive - Universiteit Leiden/projects/1_dmagna_acute_qsars/2_scripts/knn_applicability_domain.R")

source("C:/Users/balraadjsings/OneDrive - Universiteit Leiden/projects/1_dmagna_acute_qsars/2_scripts/knn_applicability_domain.R")


## Loop to calculate knn applicability domain for regression models

app_domain <- list()

for (i in 1:9) {
  
  # Turns data into matrix
  set.seed(12358673)
  
  ad_train_data <- models_reg %>% slice(i) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% 
    prep() %>% juice() %>% select(-endpoint_value) %>% as.matrix()
  
  set.seed(12358673)
  
  ad_test_data <- models_reg %>% slice(i) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% 
    prep() %>% bake(new_data=data_test_reg) %>% select(-endpoint_value) %>% as.matrix()
  
  # Select k based on the square root of the amount of observations in the training data
  k <- round(sqrt(nrow(ad_train_data)))
  
  # Calculate distances and limit of the training set
  train_dist <- train_distances(ad_train_data, k)
  
  # Calculate distances of the test set to the training set and assess applicability domain
  test_dist <- new_instances(ad_train_data, ad_test_data, k=k)
  
  
  
  app_domain[[models_reg %>% slice(i) %>% pull(algorithm) %>% pluck(1)]] <- list(train_data=ad_train_data, 
                                                                                 test_data=ad_test_data, 
                                                                                 k=k, 
                                                                                 dist_train=train_dist, 
                                                                                 dist_test=test_dist)
  
}


# Run function to see how many points were within or outside of applicability domain

app_domain %>% map(.f=\(x) x$dist_test %>% group_by(app_domain) %>% summarize(count=n()))

app_domain %>% map(.f=\(x) x$dist_test %>% group_by(app_domain) %>% summarize(count=n())) %>% 
  map_df(data.frame, .id = "algorithm") %>% pivot_wider(names_from=algorithm, values_from=count)

app_domain$bart$dist_train$limit_95

map_df(app_domain, \(x) data.frame(k = x$k,
                                   limit_95 = x$dist_train$limit_95,
                                   limit_99 = x$dist_train$limit_99), .id="algorithm")


## Bind all results together

app_domain_test_set <- app_domain$bart$dist_test %>% rename(bart_mean_k_distance=mean_k_distance,
                                                            bart_app_domain=app_domain) %>% 
  bind_cols(app_domain$cubist$dist_test %>% rename(cubist_mean_k_distance=mean_k_distance,
                                                   cubist_app_domain=app_domain)) %>% 
  bind_cols(app_domain$knn$dist_test %>% rename(knn_mean_k_distance=mean_k_distance,
                                                knn_app_domain=app_domain)) %>% 
  bind_cols(app_domain$linear$dist_test %>% rename(linear_mean_k_distance=mean_k_distance,
                                                   linear_app_domain=app_domain)) %>% 
  bind_cols(app_domain$mars$dist_test %>% rename(mars_mean_k_distance=mean_k_distance,
                                                 mars_app_domain=app_domain)) %>% 
  bind_cols(app_domain$nnet$dist_test %>% rename(nnet_mean_k_distance=mean_k_distance,
                                                 nnet_app_domain=app_domain)) %>% 
  bind_cols(app_domain$rf$dist_test %>% rename(rf_mean_k_distance=mean_k_distance,
                                               rf_app_domain=app_domain)) %>% 
  bind_cols(app_domain$svm$dist_test %>% rename(svm_mean_k_distance=mean_k_distance,
                                                svm_app_domain=app_domain)) %>% 
  bind_cols(app_domain$xgboost$dist_test %>% rename(xgboost_mean_k_distance=mean_k_distance,
                                                    xgboost_app_domain=app_domain))

## Calculate consensus and bind to test set

app_domain_test_set <- app_domain_test_set %>% mutate(reliable = rowSums(app_domain_test_set=="Reliable"),
                                                      caution = rowSums(app_domain_test_set=="Caution"),
                                                      unreliable = rowSums(app_domain_test_set=="Unreliable")) %>% 
  mutate(consensus = case_when(reliable > 4 ~ "reliable", 
                               caution > 4 ~ "caution",
                               unreliable > 4 ~ "unreliable",
                               .default="mixed")) %>%
  select(reliable, caution, unreliable, consensus)

quant_reg_results <- models_reg %>% filter(algorithm=="stack") %>% select(model) %>% pull(model) %>% pluck(1) %>% 
  extract_fit_engine() %>% predict(data=models_reg %>% filter(algorithm=="stack") %>% select(model) %>% 
                                     pull(model) %>% pluck(1) %>% extract_recipe() %>% bake(stack_test_reg), 
                                   quantiles=c(0.05, 0.5, 0.95), type="quantiles")

app_domain_test_set <- quant_reg_results$predictions %>% as_tibble() %>% 
  rename(.pred_lower=1, .pred=2, .pred_upper=3) %>% bind_cols(data_test_reg) %>% bind_cols(app_domain_test_set)

app_domain_test_set %>% ggplot(aes(endpoint_value, .pred, col=consensus)) + geom_point() + geom_abline() +
  geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + facet_wrap(.~consensus)

app_domain_test_set %>% ggplot(aes(endpoint_value, .pred, col=consensus)) + geom_point() + geom_abline() +
  geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + 
  geom_errorbar(aes(ymin=.pred_lower, ymax=.pred_upper)) +
  facet_wrap(.~consensus)


app_domain_test_set %>% ggplot(aes(endpoint_value, .pred, col=consensus)) + geom_point() + geom_abline() +
  geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + facet_wrap(.~material)

app_domain_test_set %>% ggplot(aes(endpoint_value, .pred, col=consensus)) + geom_point() + geom_abline() +
  geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2) + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + facet_wrap(.~species_group)


## Loop to calculate knn applicability domain for classification models

app_domain_class <- list()

for (i in 1:6) {
  
  # Turns data into matrix
  set.seed(12358673)
  
  ad_train_data <- models_class %>% slice(i) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% 
    prep() %>% juice() %>% select(-endpoint_category) %>% as.matrix()
  
  set.seed(12358673)
  
  ad_test_data <- models_class %>% slice(i) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% 
    prep() %>% bake(new_data=data_test_class) %>% select(-endpoint_category) %>% as.matrix()
  
  # Select k based on the square root of the amount of observations in the training data
  k <- round(sqrt(nrow(ad_train_data)))
  
  # Calculate distances and limit of the training set
  train_dist <- train_distances(ad_train_data, k)
  
  # Calculate distances of the test set to the training set and assess applicability domain
  test_dist <- new_instances(ad_train_data, ad_test_data, k=k)
  
  
  
  app_domain_class[[models_class %>% slice(i) %>% pull(algorithm) %>% pluck(1)]] <- list(train_data=ad_train_data, 
                                                                                         test_data=ad_test_data, 
                                                                                         k=k, 
                                                                                         dist_train=train_dist, 
                                                                                         dist_test=test_dist)
  
}


# Run function to see how many points were within or outside of applicability domain

app_domain_class %>% map(.f=\(x) x$dist_test %>% group_by(app_domain) %>% summarize(count=n()))

app_domain_class %>% map(.f=\(x) x$dist_test %>% group_by(app_domain) %>% summarize(count=n())) %>% 
  map_df(data.frame, .id = "algorithm") %>% pivot_wider(names_from=algorithm, values_from=count)

map_df(app_domain_class, \(x) data.frame(k = x$k,
                                         limit_95 = x$dist_train$limit_95,
                                         limit_99 = x$dist_train$limit_99), .id="algorithm")


## Bind all results together

app_domain_test_set_class <- app_domain_class$knn$dist_test %>% rename(knn_mean_k_distance=mean_k_distance,
                                                                       knn_app_domain=app_domain) %>% 
  bind_cols(app_domain_class$nnet$dist_test %>% rename(nnet_mean_k_distance=mean_k_distance,
                                                       nnet_app_domain=app_domain)) %>% 
  bind_cols(app_domain_class$rf$dist_test %>% rename(rf_mean_k_distance=mean_k_distance,
                                                     rf_app_domain=app_domain)) %>% 
  bind_cols(app_domain_class$svm$dist_test %>% rename(svm_mean_k_distance=mean_k_distance,
                                                      svm_app_domain=app_domain)) %>% 
  bind_cols(app_domain_class$xgboost$dist_test %>% rename(xgboost_mean_k_distance=mean_k_distance,
                                                          xgboost_app_domain=app_domain)) %>% 
  bind_cols(app_domain_class$multinom$dist_test %>% rename(multinom_mean_k_distance=mean_k_distance,
                                                           multinom_app_domain=app_domain))


## Calculate consensus and bind to test set

app_domain_test_set_class <- app_domain_test_set_class %>% mutate(reliable = rowSums(app_domain_test_set_class=="Reliable"),
                                                                  caution = rowSums(app_domain_test_set_class=="Caution"),
                                                                  unreliable = rowSums(app_domain_test_set_class=="Unreliable")) %>% 
  mutate(consensus = case_when(reliable > 3 ~ "reliable", 
                               caution > 3 ~ "caution",
                               unreliable > 3 ~ "unreliable",
                               .default="mixed")) %>%
  select(reliable, caution, unreliable, consensus) %>% 
  bind_cols(data_test_class,
            models_class %>% filter(algorithm=="xgboost") %>% select(model) %>% pull(model) %>% pluck(1) %>% 
              extract_workflow() %>% predict(data_test_class))



## Run for new data points (FIX THIS SO IT RUNS FOR EVERYTHING AT ONCE)

test <- models_reg %>% slice(1) %>% select(recipe) %>% pull(recipe) %>% pluck(1) %>% 
  prep() %>% bake(new_data=data_test_reg) %>% select(-endpoint_value) %>% as.matrix()

test_2 <- new_instances(app_domain$bart$train_data, test, k=app_domain$bart$k)

app_domain[[9]]$train_data



test_2





########################################################################################################




### Partial dependence profiles

## Check this method with DALEX package and it produces the same results

#augment.xgboost <- function(object, newdata) {
#  newdata <- as_tibble(newdata)
#  class_probs <- predict(object, newdata)
#  bind_cols(newdata, as_tibble(class_probs))
#}


## Function for partial depedence profiles of numerical variables 


fun_numeric <- function(data, variable, type, n=10) {
  
  ## Get min and max of variable to create grid
  
  min <- data %>% select(all_of(variable)) %>% min(., na.rm=T)
  max <- data %>% select(all_of(variable)) %>% max(., na.rm=T)
  
  ## Create all required predictors and prepare grid to calculate partial dependence profiles
  
  dat <- data %>% select(-c(all_of(variable))) # Select all predictors except the one of interest
  var <- tibble(!!variable := seq(min, max, length.out=n)) # Generate values for variable of interest 
  
  
  ## Combine both to create grid with all combinations
  
  grid <- expand_grid(dat, var)
  
  ## Generate predictions based on the created grid
  
  au <- augment(models_reg %>% filter(algorithm==type) %>% select(model) %>% pull(model) %>% pluck(1) %>% extract_workflow(), grid)
  
  ## Calculate average of all predictions while grouping by the variable of interest
  
  pd <- au %>%
    group_by_at(variable) %>%
    summarise(.pred = mean(.pred))
  
  
  return(list(algorithm=type, pdp=pd, raw=au))
  
}



## Function for partial depedence profiles of categorical variables


fun_categorical <- function(data, variable, type) {
  
  ## Create all required predictors and prepare grid to calculate partial dependence profiles
  
  dat <- data %>% select(-c(all_of(variable))) # Select all predictors except the one of interest
  var <- data %>% select((all_of(variable))) %>% unique() %>% na.omit() # Insert all categorical levels for variable of interest and drop NA
  
  ## Combine both to create grid with all combinations
  
  grid <- expand_grid(dat, var)
  
  ## Generate predictions based on the created grid
  
  au <- augment(models_reg %>% filter(algorithm==type) %>% select(model) %>% pull(model) %>% pluck(1) %>% extract_workflow(), grid)
  
  ## Calculate average of all predictions while grouping by the variable of interest
  
  pd <- au %>%
    group_by_at(variable) %>%
    summarise(.pred = mean(.pred))
  
  
  return(list(algorithm=type, pdp=pd, raw=au))
  
}



## Function for partial depedence profiles that combines the previous functions and automatically determines
## whether variables are numeric or categorical


fun_pdp <- function(data, variable, type=c("bart", "cubist", "knn", "linear", "mars", "nnet",
                                           "rf", "svm", "xgboost"), n=10) {
  
  ## Check for data type of variable (whether it is numeric or categorical)
  
  var_type <- data %>% select(all_of(variable)) %>% summarise_all(class) %>% pull(1)
  
  ## Calculate partial dependence profiles depending on the data type
  
  pdp <- if (var_type == "numeric") {
    fun_numeric(data, variable, type, n)
  } else {
    fun_categorical(data, variable, type)
  }
  
  
  return(pdp)
  
}



## Function to calculate partial dependence profiles for stacked model

## Uses previous functions to calculate PDP for all individual (regression and classification) models 
## and uses it as input for the stacked model. The average is then calculated like previous functions 

## Function is uses parallel processing through furrr 




fun_pdp_stack <- function(data, variable, n=10) {
  
  
  ## Create object with all model names
  
  model_names <- c("bart", "cubist", "knn", "mars", "nnet",
                   "rf", "svm", "xgboost")
  
  
  ## Calculate PDP for all individual models
  ## Parallel processing is used here
  
  set.seed(123)
  raw_pdp_results <- model_names %>% future_map(\(x) fun_pdp(data, variable, type=x, n),
                                                .options=furrr_options(packages=c("tidymodels", "tidyverse"),
                                                                       seed=TRUE)) %>% 
    set_names(nm=model_names)
  
  
  ## Extract predictions from raw_pdp_results for all models and bind rows together
  ## Parallel processing is used here
  
  stack_pdp_reg <- future_map(seq(1, length(model_names), 1), \(x) raw_pdp_results[[model_names[x]]]$raw %>% 
                                select(.pred) %>% rename(!!model_names[x] := .pred)) %>%
    bind_cols()
  
  
  ## Calculate classification probabilities based on grid created for PDP
  
  stack_pdp_class <- models_class %>% filter(algorithm=="xgboost") %>% select(model) %>% pull(model) %>% 
    pluck(1) %>% extract_workflow() %>% predict(raw_pdp_results[["xgboost"]]$raw %>% select(-endpoint_value), type="prob") %>% 
    rename(prob_very_toxic=`.pred_very toxic`,
           prob_toxic=.pred_toxic,
           prob_not_harmful=`.pred_not harmful`)
  
  
  ## Bind stack_pdp and stack_pdp_class together to create input data that can be used for stacked model
  
  stack_pdp_final <- stack_pdp_reg %>% bind_cols(stack_pdp_class)
  
  
  ## Use stacked model to predict endpoint and bind to PDP grid
  
  stack_raw_pdp <- raw_pdp_results[["xgboost"]]$raw %>% select(-c(endpoint_value, .pred)) %>% 
    bind_cols(models_reg %>% filter(algorithm=="stack") %>% select(model) %>% pull(model) %>% 
                pluck(1) %>% extract_workflow() %>% predict(stack_pdp_final))
  
  
  ## Calculate average of all predictions while grouping by the variable of interest
  
  stack_pdp <- stack_raw_pdp %>%
    group_by_at(variable) %>%
    summarise(.pred = mean(.pred))
  
  
  return(list(pdp=stack_pdp, raw=stack_raw_pdp))
  
  
}



## Function for plotting the partial dependence profiles

plot_pdp <- function(data, group) {
  
  var_type <- data$pdp[1] %>% summarise_all(class) %>% pull(1)
  var_name <- data$pdp[1] %>% colnames()
  
  if(var_type == "numeric") {
    
    if(missing(group)) {
      ggplot(data$pdp, aes(.data[[var_name]], .pred)) + geom_line() # PDP (average of all)
    } else {
      data$raw %>% group_by(pick(any_of(c(var_name, group)))) %>% 
        summarise(.pred = mean(.pred)) %>% 
        ggplot(aes(.data[[var_name]], .pred, col=.data[[group]])) + geom_line() # PDP (average per group)
    }
    
  } else {
    
    if(missing(group)) {
      ggplot(data$pdp, aes(.pred, .data[[var_name]], fill=.data[[var_name]])) + geom_col() +
        theme(legend.position = "none") # PDP (average of all)
    } else {
      data$raw %>% group_by(pick(any_of(c(var_name, group)))) %>% 
        summarise(.pred = mean(.pred)) %>% 
        ggplot(aes(.pred, .data[[var_name]], fill=.data[[var_name]])) + geom_col() + 
        facet_wrap(group) + theme(legend.position = "none")
    }
    
  }
  
}


## Function for plotting the partial dependence profiles (two graphs)


plot_pdp_multi <- function(object_1, object_2, group) {
  
  (plot_pdp(object_1, group) + theme(legend.position = "none") + labs(title=object_1$algorithm) + 
     theme(plot.title = element_text(hjust = 0.5),
           axis.title.x = element_blank(),
           axis.title.y = element_blank())) + 
    (plot_pdp(object_2, group) + theme(legend.position = "none") + labs(title=object_2$algorithm) + 
       theme(plot.title = element_text(hjust = 0.5),
             axis.title.x = element_blank(),
             axis.title.y = element_blank())) +
    plot_layout(guides = "collect") +
    plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) & 
    my_theme() & 
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) & 
    xlim(min(object_1$pdp[1]), max(object_1$pdp[1])) & ylim(-1.5, 2.5)
  
  
  
  
}






## Extract recipe

xgboost_rec <- models_reg %>% filter(algorithm=="xgboost") %>% select(model) %>% pull(model) %>% pluck(1) %>% extract_recipe()

## Extract which variables are removed by the recipe due to step_nzv and step_corr
## Also remove a few extra variables that will cause problems or are not relevant for the models

## dispersion, respiration, habitat are removed for now because they are already one hot encoded
## pdp is calculated separately for them

var_remove <- c(xgboost_rec[["steps"]][[1]][["removals"]], xgboost_rec[["steps"]][[4]][["removals"]], 
                "material", "species", "species_group", "endpoint_value",
                "dispersion_sonicated", "dispersion_stirred", "dispersion_vortexed",
                "respiration_body_cell_diffusion", "respiration_gills",
                "habitat_freshwater", "habitat_brackish", "habitat_marine") 


var_names <- colnames(data_train_reg) %>% tibble(colnames=.) %>% 
  filter(!colnames %in% var_remove)

var_names <- as_vector(var_names$colnames)



i <- c("dispersion_sonicated", "dispersion_stirred", "dispersion_shaked",
       "dispersion_vortexed", "dispersion_homogenized")

## Run PDP for individual models

plan(multisession, workers = 4)

set.seed(123)

tic()
pdp_xgboost <- future_map(seq(1:length(var_names)), \(x) fun_pdp(data=data_train_reg, variable=var_names[x], type="xgboost", n=50),
                          .options=furrr_options(packages=c("tidymodels", "tidyverse"),
                                                 seed=TRUE))
toc()



## Calculate PDP for remaining (problematic) variables
set.seed(123)
pdp_xgboost[[42]] <- fun_pdp(data=data_train_reg, variable=c("dispersion_sonicated", "dispersion_stirred", "dispersion_vortexed",
                                                             "dispersion_homogenized"), 
                             type="xgboost", n=50)

set.seed(123)
pdp_xgboost[[43]] <- fun_pdp(data=data_train_reg, variable=c("respiration_body_cell_diffusion", "respiration_gills",
                                                             "respiration_lungs", "respiration_stomata_diffusion"), 
                             type="xgboost", n=50)

set.seed(123)
pdp_xgboost[[44]] <- fun_pdp(data=data_train_reg, variable=c("habitat_freshwater", "habitat_brackish", "habitat_marine"), 
                             type="xgboost", n=50)


## Plot results

plot_pdp(pdp_xgboost[[1]], "material")


pdp_xgboost[[42]]$raw %>% 
  mutate(category=case_when(dispersion_sonicated=="no" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no" ~ "none",
                            dispersion_sonicated=="no" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="yes" & dispersion_homogenized=="no"  ~ "vortexed",
                            dispersion_sonicated=="no" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "stirred",
                            dispersion_sonicated=="no" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="yes"  ~ "homogenized",
                            dispersion_sonicated=="yes" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "sonicated",
                            dispersion_sonicated=="yes" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="yes" & dispersion_homogenized=="no"  ~ "sonicated_vortexed",
                            dispersion_sonicated=="yes" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "sonicated_stirred")) %>% 
  group_by(category, material) %>% 
  filter(category=="vortexed" | category=="stirred" | category=="homogenized" | category=="sonicated") %>%
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))


pdp_xgboost[[43]]$raw %>% 
  mutate(category=case_when(respiration_body_cell_diffusion==0 & respiration_gills==1 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "gills",
                            respiration_body_cell_diffusion==1 & respiration_gills==0 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "diffusion",
                            respiration_body_cell_diffusion==1 & respiration_gills==1 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "diffusion_gills")) %>% 
  group_by(category, material) %>% 
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))



pdp_xgboost[[44]]$raw %>% 
  mutate(category=case_when(habitat_freshwater==0 & habitat_brackish==0 & habitat_marine==0 ~ "none",
                            habitat_freshwater==0 & habitat_brackish==0 & habitat_marine==1 ~ "marine",
                            habitat_freshwater==0 & habitat_brackish==1 & habitat_marine==0 ~ "brackish",
                            habitat_freshwater==0 & habitat_brackish==1 & habitat_marine==1 ~ "brackish_marine",
                            habitat_freshwater==1 & habitat_brackish==0 & habitat_marine==0 ~ "freshwater",
                            habitat_freshwater==1 & habitat_brackish==0 & habitat_marine==1 ~ "freshwater_marine",
                            habitat_freshwater==1 & habitat_brackish==1 & habitat_marine==0 ~ "freshwater_brackish",
                            habitat_freshwater==1 & habitat_brackish==1 & habitat_marine==1 ~ "freshwater_brackish_marine")) %>% 
  group_by(category, material) %>% 
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))





## Run PDP for stacked model

plan(multisession, workers = 4)

tic()
set.seed(123)
pdp_stack <- future_map(seq(1:length(var_names)), \(x) fun_pdp_stack(data=data_train_reg, variable=var_names[x], n=50),
                        .options=furrr_options(packages=c("tidymodels", "tidyverse"),
                                               seed=TRUE)) ## Takes 3 hours
toc()





## Calculate PDP for remaining (problematic) variables

plan(multisession, workers = 4)

set.seed(123)
pdp_stack[[42]] <- fun_pdp_stack(data=data_train_reg, variable=c("dispersion_sonicated", "dispersion_stirred", "dispersion_vortexed",
                                                                 "dispersion_homogenized"), 
                                 n=50)

set.seed(123)
pdp_stack[[43]] <- fun_pdp_stack(data=data_train_reg, variable=c("respiration_body_cell_diffusion", "respiration_gills",
                                                                 "respiration_lungs", "respiration_stomata_diffusion"), 
                                 n=50)

set.seed(123)
pdp_stack[[44]] <- fun_pdp_stack(data=data_train_reg, variable=c("habitat_freshwater", "habitat_brackish", "habitat_marine"), 
                                 n=50)



## Plot results

plot_pdp(pdp_stack[[1]])
plot_pdp(pdp_stack[[1]], "material")


pdp_stack[[42]]$raw %>% 
  mutate(category=case_when(dispersion_sonicated=="no" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no" ~ "none",
                            dispersion_sonicated=="no" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="yes" & dispersion_homogenized=="no"  ~ "vortexed",
                            dispersion_sonicated=="no" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "stirred",
                            dispersion_sonicated=="no" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="yes"  ~ "homogenized",
                            dispersion_sonicated=="yes" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "sonicated",
                            dispersion_sonicated=="yes" & dispersion_stirred=="no" & 
                              dispersion_vortexed=="yes" & dispersion_homogenized=="no"  ~ "sonicated_vortexed",
                            dispersion_sonicated=="yes" & dispersion_stirred=="yes" & 
                              dispersion_vortexed=="no" & dispersion_homogenized=="no"  ~ "sonicated_stirred")) %>% 
  group_by(category, material) %>% 
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))



pdp_stack[[43]]$raw %>% 
  mutate(category=case_when(respiration_body_cell_diffusion==0 & respiration_gills==1 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "gills",
                            respiration_body_cell_diffusion==1 & respiration_gills==0 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "diffusion",
                            respiration_body_cell_diffusion==1 & respiration_gills==1 & 
                              respiration_lungs==0 & respiration_stomata_diffusion==0 ~ "diffusion_gills")) %>% 
  group_by(category, material) %>% 
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))



pdp_stack[[44]]$raw %>% 
  mutate(category=case_when(habitat_freshwater==0 & habitat_brackish==0 & habitat_marine==0 ~ "none",
                            habitat_freshwater==0 & habitat_brackish==0 & habitat_marine==1 ~ "marine",
                            habitat_freshwater==0 & habitat_brackish==1 & habitat_marine==0 ~ "brackish",
                            habitat_freshwater==0 & habitat_brackish==1 & habitat_marine==1 ~ "brackish_marine",
                            habitat_freshwater==1 & habitat_brackish==0 & habitat_marine==0 ~ "freshwater",
                            habitat_freshwater==1 & habitat_brackish==0 & habitat_marine==1 ~ "freshwater_marine",
                            habitat_freshwater==1 & habitat_brackish==1 & habitat_marine==0 ~ "freshwater_brackish",
                            habitat_freshwater==1 & habitat_brackish==1 & habitat_marine==1 ~ "freshwater_brackish_marine")) %>% 
  group_by(category, material) %>% 
  summarise(.pred = mean(.pred)) %>% 
  ggplot(aes(category, .pred, fill=material)) + geom_col() + 
  facet_wrap(.~material) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))





plot_pdp_multi(pdp_stack, pdp_xgboost, "material") & ylim(-2, 3)


grid.arrange(patchworkGrob(plot_pdp_multi(pdp_stack, pdp_xgboost, "material")), left = text_grob("Prediction", 
                                                                                                 size = 14, 
                                                                                                 face = "bold",
                                                                                                 family = "Roboto Slab",
                                                                                                 rot = 90,
                                                                                                 vjust = 1.5), 
             bottom=text_grob("Primary diameter", 
                              size = 14, 
                              face = "bold",
                              family = "Roboto Slab",
                              vjust = -0.5))








#a <- (plot_pdp(pdp_stack, "material") + theme(legend.position = "none") + labs(title="Stack") + 
#        theme(plot.title = element_text(hjust = 0.5),
#              axis.title.x = element_blank(),
#              axis.title.y = element_blank())) + 
#  (plot_pdp(pdp_xgboost, "material") + theme(legend.position = "none") + labs(title="XGBoost") + 
#     theme(plot.title = element_text(hjust = 0.5),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank())) +
#  plot_layout(guides = "collect") +
#  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5))) & 
#  my_theme() & 
#  theme(legend.position = "top",
#        legend.title = element_blank(),
#        axis.title = element_blank(),
#        axis.title.x = element_blank(),
#        axis.title.y = element_blank()) & 
#  xlim(0, 600) & ylim(-1.5, 2.5)


#grid.arrange(patchworkGrob(a), left = text_grob("Prediction", 
#                                               size = 14, 
#                                              face = "bold",
#                                             family = "Roboto Slab",
#                                                rot = 90,
#                                                vjust = 1.5), 
#             bottom=text_grob("Primary diameter", 
#                              size = 14, 
#                              face = "bold",
#                              family = "Roboto Slab",
#                              vjust = -0.5))












