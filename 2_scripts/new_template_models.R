library(tidyverse) # Collection of packages for data wrangling + plotting etc. 
library(tidymodels) # Collection of packages for machine learning
library(doParallel) # Parallel processing
library(skimr) # Inspect dataset
library(rules) # Cubist rules model


## Set working directory

setwd("...")


## Load raw data

load("./2_scripts/clean_data.RData")


## Prepare data for analysis

# Select only few materials (top 15)
# Turn EC50 data into numeric and discard all categorical data
# Filter on experiments done with pristine particles (no aging, pre-illumination, uv-radiation)
# Remove insects and plant experiments
# Filter on relevant endpoints
# Remove Poterioochromonas malhamensis because it causes problems
data_analysis <- data_final %>% 
  filter(material=="Ag" | material=="ZnO" | material=="TiO2" | material=="SiO2" | material=="CeO2" | 
           material=="CuO" | material=="Cu" | material=="Au" | material=="Fe2O3" | material=="Fe3O4" |
           material=="Al2O3" | material=="Pt" | material=="Ni" | material=="Se" | material=="NiO") %>% 
  mutate(endpoint_value=as.numeric(endpoint_value)) %>% 
  drop_na(endpoint_value) %>% 
  filter(suspension_aging=="no", pre_illumination=="no", uv_radiation=="no") %>% 
  filter(species_group!="Plant") %>% 
  filter(endpoint=="mortality" | endpoint=="growth inhibition" | endpoint=="immobilization" | 
           endpoint=="cell viability") %>% 
  filter(genus!="Poterioochromonas" & epithet!="malhamensis")

# Load solubility classes (from D. magna paper)
# Add extra data for new materials (based on literature)
solubility <- read.csv("./1_data/other/solubility.csv")

solubility <- solubility %>% bind_rows(data.frame(material=c("Al2O3", "Pt", "Ni", "Se", "NiO"), 
                                                  solubility_group=c("slow", "slow", "slow", "slow", "slow")))


# Join solubility classes to dataset
data_analysis <- data_analysis %>% left_join(solubility, by=c("material"="material"))

# Create strata for resampling
# Fill in life stage
data_analysis <- data_analysis %>% mutate(species=paste(genus, epithet, sep=" "),
                                          species_group=case_when(phylum=="Annelida" ~ "Annelida",
                                                                  phylum=="Arthropoda" ~ "Crustacea",
                                                                  phylum=="Bacillariophyta" ~ "Diatom",
                                                                  phylum=="Charophyta" | phylum=="Chlorophyta" | phylum=="Haptophyta" | phylum=="Ochrophyta" ~ "Algae",
                                                                  phylum=="Chordata" ~ "Fish",
                                                                  phylum=="Ciliophora" | phylum=="Euglenozoa" ~ "Protozoa",
                                                                  phylum=="Cnidaria" ~ "Cnidaria",
                                                                  phylum=="Cyanobacteria" ~ "Cyanobacteria",
                                                                  phylum=="Mollusca" ~ "Gastropoda", 
                                                                  phylum=="Nematoda" ~ "Nematoda",
                                                                  phylum=="Rotifera" ~ "Rotifera"),
                                          strata=paste(species_group, material, sep="_"),
                                          life_stage=case_when(species_group=="Algae" | species_group=="Cyanobacteria" | 
                                                                 species_group=="Diatom" | species_group=="Protozoa" ~ "other life stage",
                                                               str_detect(life_stage, "neonate") | str_detect(life_stage, "juvenile") | 
                                                                 str_detect(life_stage, "nauplii") | str_detect(life_stage, "larvae") | 
                                                                 str_detect(life_stage, "zoea") | str_detect(life_stage, "instar") | 
                                                                 str_detect(life_stage, "stage") | str_detect(life_stage, "phase") ~ "juvenile",
                                                               str_detect(life_stage, "embryo") | str_detect(life_stage, "egg") ~ "embryo",
                                                               str_detect(life_stage, "adult") ~ "adult",
                                                               species_group=="Crustacea" & str_detect(life_stage, "days") ~ "adult",
                                                               species=="Labeo rohita" & str_detect(life_stage, "90day old") ~ "juvenile",
                                                               species=="Oryzias latipes" & str_detect(life_stage, "month") ~ "adult",
                                                               species=="Pimephales promelas" & str_detect(life_stage, "day") ~ "juvenile",
                                                               species_group=="Crustacea" & str_detect(life_stage, "h post hatch") ~ "juvenile",
                                                               species_group=="Nematoda" & str_detect(life_stage, "4days") ~ "adult",
                                                               species=="Danio rerio" & str_detect(life_stage, "120days") ~ "adult",
                                                               species=="Barbonymus gonionotus" & str_detect(life_stage, "2 month old") ~ "juvenile"
                                          ))

# Join molecular descriptors to dataset
# Remove unnecesary variables

ochem_descr_processed <- ochem_descr_processed %>% mutate(material_group=case_when(Material=="TiO2" | Material=="ZnO" |
                                                                                     Material=="CuO" | Material=="NiO" | 
                                                                                     Material=="SiO2" | Material=="Al2O3" | 
                                                                                     Material=="Fe2O3" | Material=="Fe3O4" |
                                                                                     Material=="CeO2" ~ "metal oxide", 
                                                                                   .default="metal"))


data_analysis <- data_analysis %>% left_join(ochem_descr_processed, by=c("material"="Material")) %>% 
  select(-c(primary_method, fed_during_test, endpoint_unit, author:epithet,
            length_juvenile:feeding_strategy, SMILES))


# Change coating to binary variable
# Change shape of spherical and nearly spherical particles to spheroids and remove the unknown category
data_analysis <- data_analysis %>% mutate(coating=ifelse(coating=="uncoated", "uncoated", "coated"),
                                          shape=ifelse(shape!="unknown", shape, NA_character_),
                                          shape=ifelse(shape=="nearly spherical" | 
                                                         shape=="spherical", "spheroid", shape))


# Merge dispersion methods and turn into binary variables
data_analysis <- data_analysis %>% 
  mutate(dispersion_immediately_before=ifelse(dispersion_immediately_before=="sonicated + shaken", 
                                              "sonicated + shaked", dispersion_immediately_before)) %>% 
  separate(dispersion_immediately_before, c("dispersion_1", "dispersion_2"),
           extra = "merge", fill = "left", remove=T, sep=" \\+ ") %>% 
  separate(dispersion_before_unknown, c("dispersion_3", "dispersion_4"),
           extra = "merge", fill = "left", remove=T, sep=" \\+ ") %>% 
  mutate(dispersion_sonicated = case_when(str_detect(dispersion_1, "sonicated") ~ "yes", 
                                          str_detect(dispersion_2, "sonicated") ~ "yes",
                                          str_detect(dispersion_3, "sonicated") ~ "yes",
                                          str_detect(dispersion_4, "sonicated") ~ "yes",
                                          .default="no"),
         dispersion_stirred = case_when(str_detect(dispersion_1, "stirred") ~ "yes", 
                                        str_detect(dispersion_2, "stirred") ~ "yes",
                                        str_detect(dispersion_3, "stirred") ~ "yes",
                                        str_detect(dispersion_4, "stirred") ~ "yes",
                                        .default="no"),
         dispersion_shaked = case_when(str_detect(dispersion_1, "shaked") ~ "yes", 
                                       str_detect(dispersion_2, "shaked") ~ "yes",
                                       str_detect(dispersion_3, "shaked") ~ "yes",
                                       str_detect(dispersion_4, "shaked") ~ "yes",
                                       .default="no"),
         dispersion_vortexed = case_when(str_detect(dispersion_1, "vortexed") ~ "yes", 
                                         str_detect(dispersion_2, "vortexed") ~ "yes",
                                         str_detect(dispersion_3, "vortexed") ~ "yes",
                                         str_detect(dispersion_4, "vortexed") ~ "yes",
                                         .default="no"),
         dispersion_homogenized = case_when(str_detect(dispersion_1, "homogenized") ~ "yes", 
                                            str_detect(dispersion_2, "homogenized") ~ "yes",
                                            str_detect(dispersion_3, "homogenized") ~ "yes",
                                            str_detect(dispersion_4, "homogenized") ~ "yes",
                                            .default="no")) %>% 
  select(-c(dispersion_1:dispersion_4))

# Convert illumination conditions to proportions
# Remove one of them because they contain more or less the same information

data_analysis <- data_analysis %>% mutate(illumination_dark=illumination_dark/24,
                                          illumination_light=illumination_light/24) %>% 
  select(-illumination_dark)

# Log transform EC50 data to fix data skewness
# Remove unneccesary features for modeling
# Convert characters to factors (important for tidymodels)
data_analysis <- data_analysis %>% mutate(endpoint_value=log10(endpoint_value)) %>% 
  select(-c(test_procedure, ionic_strength)) %>%
  mutate(endpoint=case_when(endpoint=="mortality" | endpoint=="immobilization" | 
                              endpoint=="cell viability" ~ "mortality",
                            endpoint=="growth inhibition" ~ "growth inhibition")) %>% 
  mutate(across(where(is_character), as_factor))



## Remove Au, Ni, Pt material because of poor performance in models
## Remove Nematoda, Gastropoda because of poor performance in models
## Remove zeta potential because of proportion of missing data

data_analysis <- data_analysis %>% filter(!(material %in% c("Au", "Pt", "Ni"))) %>%
  filter(!(species_group %in% c("Nematoda", "Gastropoda"))) %>% 
  select(-zeta_potential)



## Inspect distribution of EC50 values
ggplot(data_analysis, aes(endpoint_value)) + geom_density()

## Inspect dataset
skim(data_analysis) %>% view()




#################################################################################################################


# Split data into 60% training and 40% test set
set.seed(123)

data_split <- initial_split(data_analysis, strata=strata,
                            prop=0.6)


data_train <- training(data_split)
data_test <- testing(data_split)

# 10-fold cross validation of training set
set.seed(345)
data_folds <- vfold_cv(data_train, v=10, repeats=1)




## Recipe for model + preprocessing steps ##

data_train <- data_train %>% select(-c(material, species_group, species, strata))
data_test <- data_test %>% select(-c(material, species_group, species, strata))



# Zero variance predictors are removed
# All categorical variables with levels that are not frequent, changed to "other"
# All categorical variables are one hot encoded (dummy variables)
# Correlated variables are removed
# Data is normalized (centered and scaled)
# Missing data is imputed with k-Nearest neighbors 

model_rec <- recipe(endpoint_value ~ ., data=data_train) %>%
  step_zv(all_predictors(), -c(habitat_freshwater, habitat_brackish, habitat_marine,
                                respiration_body_cell_diffusion, respiration_gills,
                                dispersion_stirred, dispersion_shaked, dispersion_sonicated, dispersion_vortexed,
                                dispersion_homogenized)) %>% 
  step_other(all_nominal_predictors()) %>% 
  step_normalize(all_numeric_predictors(), -c(habitat_freshwater, habitat_brackish, habitat_marine, 
                                              respiration_body_cell_diffusion, respiration_gills)) %>%
  step_corr(all_numeric_predictors()) %>%
  step_impute_knn(all_predictors(), neighbors=tune("impute")) %>%
  step_dummy(all_nominal_predictors(), one_hot = T)


################################################################################################################################


## k-Nearest neighbors ##


# Specify model and engine
knn_spec <- nearest_neighbor(neighbors=tune(), dist_power=tune()) %>% 
  set_mode("regression") %>%
  set_engine("kknn")

# Create workflow with pre-processing recipe and model specification
knn_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(knn_spec) 


## Initial hyperparameter tuning ## 

# Set hyperparameter ranges
param <- extract_parameter_set_dials(knn_workflow) %>% 
  update(neighbors = neighbors(c(5, 8)),
         impute = neighbors(c(50, 100)),
         dist_power = dist_power(c(0.1, 0.15)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Set metrics for evaluation

metrics <- metric_set(rmse, mae, rsq, ccc, huber_loss, rpiq, rpd, iic)

# Initial tuning
set.seed(456)
knn_tune <- tune_grid(knn_workflow, 
                      resamples=data_folds, 
                      grid=80,
                      control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                      metrics=metrics,
                      param_info=param)

stopCluster(cl) # Shut down cluster


# Plot hyperparameters to choose best range for final tuning
knn_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, neighbors, dist_power, impute) %>% 
  pivot_longer(neighbors:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

knn_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, neighbors, dist_power, impute) %>% 
  pivot_longer(neighbors:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")


## Show best
knn_tune %>% show_best("rsq", n=10)
knn_tune %>% show_best("rmse", n=10)
knn_tune %>% show_best("ccc", n=10)
knn_tune %>% show_best("huber_loss", n=10)


knn_best_parms <- knn_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
knn_final_rec <- finalize_recipe(model_rec, knn_best_parms)
knn_final_model <- finalize_model(knn_spec, knn_best_parms)

# Create new workflow for final model
knn_final_wf <- workflow() %>% 
  add_recipe(knn_final_rec) %>%
  add_model(knn_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
knn_final <- knn_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(knn_final) # Metrics of model (test set)


knn_tune %>% 
  collect_predictions() %>% 
  inner_join(knn_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)


# Actual vs predicted plot
collect_predictions(knn_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(knn_final) %>% mutate(offset=endpoint_value-.pred,
                                          freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))


################################################################################################################################


## Randomforest ##



# Specify model and engine + variable importance method
rf_spec <- rand_forest(trees=tune(), 
                       min_n=tune(), 
                       mtry=tune()) %>%
  set_mode("regression") %>% 
  set_engine("ranger", importance="permutation")


# Create workflow with pre-processing recipe and model specification
rf_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(rf_spec) 




## Initial hyperparameter tuning ## 

# Set hyperparameter ranges
param <- extract_parameter_set_dials(rf_workflow) %>% update(mtry = mtry(c(20, 30)),
                                                             min_n = min_n(c(1, 2)),
                                                             trees = trees(c(1000, 1500)),
                                                             impute = neighbors(c(50, 100)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
rf_tune <- tune_grid(rf_workflow, 
                     resamples=data_folds, 
                     grid=80,
                     control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                     metrics=metrics,
                     param_info=param)

stopCluster(cl) # Shut down cluster


# Plot hyperparameters to choose best range for final tuning
rf_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, mtry, trees, min_n, impute) %>% 
  pivot_longer(mtry:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

rf_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, mtry, trees, min_n, impute) %>% 
  pivot_longer(mtry:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")


## Show best
rf_tune %>% show_best("rsq", n=10)
rf_tune %>% show_best("rmse", n=10)
rf_tune %>% show_best("ccc", n=10)
rf_tune %>% show_best("huber_loss", n=10)


rf_best_parms <- rf_tune %>% select_best("rmse") # Select best parameters based on RMSE

## Final model ##


# Update hyperparameters in recipe and model 
rf_final_rec <- finalize_recipe(model_rec, rf_best_parms)
rf_final_model <- finalize_model(rf_spec, rf_best_parms)

# Create new workflow for final model
rf_final_wf <- workflow() %>% 
  add_recipe(rf_final_rec) %>%
  add_model(rf_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
rf_final <- rf_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(rf_final) # Metrics of model (test set)


rf_tune %>% collect_predictions() %>% 
  inner_join(rf_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)


# Actual vs predicted plot
collect_predictions(rf_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(rf_final) %>% mutate(offset=endpoint_value-.pred,
                                         freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))





################################################################################################################################


## Bayesian additive regression trees (BART) ##


# Specify model and engine
bart_spec <- parsnip::bart(trees=tune(),
                           prior_terminal_node_coef=tune(),
                           prior_terminal_node_expo=tune(),
                           prior_outcome_range=tune()) %>%
  set_mode("regression") %>% 
  set_engine("dbarts")


# Create workflow with pre-processing recipe and model specification
bart_workflow <- workflow() %>%
  add_recipe(model_rec) %>% 
  add_model(bart_spec) 




# Set hyperparameter ranges
param <- extract_parameter_set_dials(bart_workflow) %>% 
  update(trees = trees(c(1000, 1500)),
         prior_terminal_node_coef=prior_terminal_node_coef(c(0.7, 1), NULL),
         prior_terminal_node_expo=prior_terminal_node_expo(c(1, 2)),
         prior_outcome_range=prior_outcome_range(c(0.25, 0.75), NULL),
         impute = neighbors(c(25, 75)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
bart_tune <- tune_grid(bart_workflow, 
                       resamples=data_folds, 
                       grid=100,
                       control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                       metrics=metrics,
                       param_info = param)

stopCluster(cl) # Shut down cluster

# Plot hyperparameters to choose best range for final tuning
bart_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, prior_terminal_node_coef, prior_terminal_node_expo, prior_outcome_range, trees, impute) %>% 
  pivot_longer(prior_terminal_node_coef:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

bart_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, prior_terminal_node_coef, prior_terminal_node_expo, prior_outcome_range, trees, impute) %>% 
  pivot_longer(prior_terminal_node_coef:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

## Show best
bart_tune %>% show_best("rsq", n=10)
bart_tune %>% show_best("rmse", n=10)
bart_tune %>% show_best("ccc", n=10)
bart_tune %>% show_best("huber_loss", n=10)


bart_best_parms <- bart_tune %>% select_best("rmse") # Select best parameters based on RMSE

## Final model ##


# Update hyperparameters in recipe and model 
bart_final_rec <- finalize_recipe(model_rec, bart_best_parms)
bart_final_model <- finalize_model(bart_spec, bart_best_parms)

# Create new workflow for final model
bart_final_wf <- workflow() %>% 
  add_recipe(bart_final_rec) %>%
  add_model(bart_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
bart_final <- bart_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(bart_final) # Metrics of model (test set)


bart_tune %>% collect_predictions() %>% 
  inner_join(bart_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(bart_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(bart_final) %>% mutate(offset=endpoint_value-.pred,
                                           freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))







################################################################################################################################




## Artificial neural network ##



# Specify model and engine
nnet_spec <- mlp(hidden_units=tune(),
                 penalty=tune(),
                 epochs=tune()) %>%
  set_mode("regression") %>% 
  set_engine("nnet")


# Create workflow with pre-processing recipe and model specification
nnet_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(nnet_spec) 




## Initial hyperparameter tuning ## 


# Set hyperparameter ranges
param <- extract_parameter_set_dials(nnet_workflow) %>% update(hidden_units = hidden_units(c(10, 10)),
                                                               penalty = penalty(c(2.5, 3), NULL),
                                                               epochs = epochs(c(400, 1000)),
                                                               impute = neighbors(c(30, 80)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
nnet_tune <- tune_grid(nnet_workflow,
                       resamples=data_folds, 
                       grid=80,
                       control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                       metrics=metrics,
                       param_info = param)

stopCluster(cl) # Shut down cluster

# Plot hyperparameters to choose best range for final tuning
nnet_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, hidden_units, penalty, epochs, impute) %>% 
  pivot_longer(hidden_units:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

nnet_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, hidden_units, penalty, epochs, impute) %>% 
  pivot_longer(hidden_units:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")



## Show best
nnet_tune %>% show_best("rsq", n=10)
nnet_tune %>% show_best("rmse", n=10)
nnet_tune %>% show_best("ccc", n=10)
nnet_tune %>% show_best("huber_loss", n=10)


nnet_best_parms <- nnet_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
nnet_final_rec <- finalize_recipe(model_rec, nnet_best_parms)
nnet_final_model <- finalize_model(nnet_spec, nnet_best_parms)

# Create new workflow for final model
nnet_final_wf <- workflow() %>% 
  add_recipe(nnet_final_rec) %>%
  add_model(nnet_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
nnet_final <- nnet_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(nnet_final) # Metrics of model (test set)


nnet_tune %>% collect_predictions() %>% 
  inner_join(nnet_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(nnet_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(nnet_final) %>% mutate(offset=endpoint_value-.pred,
                                           freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))






################################################################################################################################






## Multivariate adaptive regression splines (MARS) ##




# Specify model and engine
mars_spec <- mars(num_terms=tune(),
                  prod_degree=tune()) %>%
  set_mode("regression") %>% 
  set_engine("earth")


# Create workflow with pre-processing recipe and model specification
mars_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(mars_spec) 




## Initial hyperparameter tuning ## 


# Set hyperparameter ranges
param <- extract_parameter_set_dials(mars_workflow) %>% update(num_terms=num_terms(c(30, 80)),
                                              prod_degree=prod_degree(c(1, 1)),
                                              impute = neighbors(c(50, 100)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
mars_tune <- tune_grid(mars_workflow,
                       resamples=data_folds, 
                       grid=100,
                       control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                       metrics=metrics,
                       param_info = param)

stopCluster(cl) # Shut down cluster


mars_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, num_terms, prod_degree, impute) %>% 
  pivot_longer(num_terms:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

mars_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, num_terms, prod_degree, impute) %>% 
  pivot_longer(num_terms:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

## Show best
mars_tune %>% show_best("rsq", n=10)
mars_tune %>% show_best("rmse", n=10)
mars_tune %>% show_best("ccc", n=10)
mars_tune %>% show_best("huber_loss", n=10)


mars_best_parms <- mars_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
mars_final_rec <- finalize_recipe(model_rec, mars_best_parms)
mars_final_model <- finalize_model(mars_spec, mars_best_parms)

# Create new workflow for final model
mars_final_wf <- workflow() %>% 
  add_recipe(mars_final_rec) %>%
  add_model(mars_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
mars_final <- mars_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(mars_final) # Metrics of model (test set)


mars_tune %>% collect_predictions() %>% 
  inner_join(mars_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(mars_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(mars_final) %>% mutate(offset=endpoint_value-.pred,
                                           freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))





################################################################################################################################



## Cubist rules ##



# Specify model and engine
cubist_spec <- cubist_rules(committees=tune(), neighbors=tune(), max_rules=tune()) %>%
  set_mode("regression") %>% 
  set_engine("Cubist")


# Create workflow with pre-processing recipe and model specification
cubist_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(cubist_spec) 




## Initial hyperparameter tuning ## 


# Set hyperparameter ranges
param <- extract_parameter_set_dials(cubist_workflow) %>% update(committees=committees(c(80, 100)),
                                                neighbors=neighbors(c(5, 5)),
                                                max_rules=max_rules(c(15, 25)),
                                                impute = neighbors(c(1, 25)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
cubist_tune <- tune_grid(cubist_workflow,
                         resamples=data_folds, 
                         grid=100,
                         control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                         metrics=metrics,
                         param_info = param)

stopCluster(cl) # Shut down cluster


cubist_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, committees, neighbors, max_rules, impute) %>% 
  pivot_longer(committees:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

cubist_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, committees, neighbors, max_rules, impute) %>% 
  pivot_longer(committees:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")


## Show best
cubist_tune %>% show_best("rsq", n=10)
cubist_tune %>% show_best("rmse", n=10)
cubist_tune %>% show_best("ccc", n=10)
cubist_tune %>% show_best("huber_loss", n=10)


cubist_best_parms <- cubist_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
cubist_final_rec <- finalize_recipe(model_rec, cubist_best_parms)
cubist_final_model <- finalize_model(cubist_spec, cubist_best_parms)

# Create new workflow for final model
cubist_final_wf <- workflow() %>% 
  add_recipe(cubist_final_rec) %>%
  add_model(cubist_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
cubist_final <- cubist_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(cubist_final) # Metrics of model (test set)


cubist_tune %>% collect_predictions() %>% 
  inner_join(cubist_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(cubist_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(cubist_final) %>% mutate(offset=endpoint_value-.pred,
                                             freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))




################################################################################################################################




## (Radial) Support vector machine ##




# Specify model and engine
svm_spec <- svm_rbf(cost=tune(), rbf_sigma=tune()) %>%
  set_mode("regression") %>% 
  set_engine("kernlab")


# Create workflow with pre-processing recipe and model specification
svm_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(svm_spec, formula=endpoint_value~.) 




## Initial hyperparameter tuning ## 


# Set hyperparameter range for mtry (otherwise there are errors)
param <- extract_parameter_set_dials(svm_workflow) %>% 
  update(cost=cost(c(0.95, 1), NULL),
         rbf_sigma=rbf_sigma(c(0.02, 0.04), NULL),
         impute = neighbors(c(50, 100)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
svm_tune <- tune_grid(svm_workflow,
                      resamples=data_folds, 
                      grid=80,
                      control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                      metrics=metrics,
                      param_info = param)

stopCluster(cl) # Shut down cluster

svm_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, cost, rbf_sigma, impute) %>% 
  pivot_longer(cost:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

svm_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, cost, rbf_sigma, impute) %>% 
  pivot_longer(cost:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

## Show best
svm_tune %>% show_best("rsq", n=10)
svm_tune %>% show_best("rmse", n=10)
svm_tune %>% show_best("ccc", n=10)
svm_tune %>% show_best("huber_loss", n=10)


svm_best_parms <- svm_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
svm_final_rec <- finalize_recipe(model_rec, svm_best_parms)
svm_final_model <- finalize_model(svm_spec, svm_best_parms)

# Create new workflow for final model
svm_final_wf <- workflow() %>% 
  add_recipe(svm_final_rec) %>%
  add_model(svm_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
svm_final <- svm_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(svm_final) # Metrics of model (test set)


svm_tune %>% collect_predictions() %>% 
  inner_join(svm_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(svm_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(svm_final) %>% mutate(offset=endpoint_value-.pred,
                                          freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))





################################################################################################################################



## XGBoost ##



# Specify model and engine
xgboost_spec <- boost_tree(mtry=tune(), trees=tune(), min_n=tune(), tree_depth=tune(),
                           learn_rate=tune(), loss_reduction=tune(),
                           sample_size=tune(), stop_iter=tune()) %>%
  set_mode("regression") %>% 
  set_engine("xgboost")


# Create workflow with pre-processing recipe and model specification
xgboost_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(xgboost_spec, formula=endpoint_value~.) 




## Initial hyperparameter tuning ## 


# Set hyperparameter ranges
param <- extract_parameter_set_dials(xgboost_workflow) %>% 
  update(mtry=mtry(c(25, 35)), 
         trees=trees(c(500, 1000)), 
         min_n=min_n(c(5, 10)), 
         tree_depth=tree_depth(c(30, 50)),
         learn_rate=learn_rate(c(0.01, 0.08), NULL), 
         loss_reduction=loss_reduction(c(0, 0.1), NULL),
         sample_size=sample_size(c(1, 1)), 
         stop_iter=stop_iter(c(200, 400)),
         impute = neighbors(c(1, 50)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
xgboost_tune <- tune_grid(xgboost_workflow,
                          resamples=data_folds, 
                          grid=120,
                          control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                          metrics=metrics,
                          param_info = param)

stopCluster(cl) # Shut down cluster

xgboost_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, mtry:impute) %>% 
  pivot_longer(mtry:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

xgboost_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, mtry:impute) %>% 
  pivot_longer(mtry:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

## Show best
xgboost_tune %>% show_best("rsq", n=10)
xgboost_tune %>% show_best("rmse", n=10)
xgboost_tune %>% show_best("ccc", n=10)
xgboost_tune %>% show_best("huber_loss", n=10)


xgboost_best_parms <- xgboost_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
xgboost_final_rec <- finalize_recipe(model_rec, xgboost_best_parms)
xgboost_final_model <- finalize_model(xgboost_spec, xgboost_best_parms)

# Create new workflow for final model
xgboost_final_wf <- workflow() %>% 
  add_recipe(xgboost_final_rec) %>%
  add_model(xgboost_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
xgboost_final <- xgboost_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(xgboost_final) # Metrics of model (test set)


xgboost_tune %>% collect_predictions() %>% 
  inner_join(xgboost_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(xgboost_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(xgboost_final) %>% mutate(offset=endpoint_value-.pred,
                                              freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))




################################################################################################################################


## Regularized linear model ##



# Specify model and engine + variable importance method
linear_spec <- linear_reg(penalty=tune(), mixture=tune()) %>%
  set_mode("regression") %>% 
  set_engine("glmnet")


# Create workflow with pre-processing recipe and model specification
linear_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(linear_spec, formula=endpoint_value~.) 




## Initial hyperparameter tuning ## 


# Set hyperparameter range for mtry (otherwise there are errors)
param <- extract_parameter_set_dials(linear_workflow) %>% 
  update(penalty = penalty(c(0, 0.25), NULL),
         mixture = mixture(c(0, 0.5), NULL),
         impute = neighbors(c(1, 100)))

# Parallel processing to optimize tuning (Note: can cause R to crash) 
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Initial tuning
set.seed(456)
linear_tune <- tune_grid(linear_workflow,
                         resamples=data_folds, 
                         grid=100,
                         control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                         metrics=metrics,
                         param_info = param)

stopCluster(cl) # Shut down cluster

linear_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, penalty, mixture, impute) %>% 
  pivot_longer(penalty:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

linear_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, penalty, mixture, impute) %>% 
  pivot_longer(penalty:impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

## Show best
linear_tune %>% show_best("rsq", n=10)
linear_tune %>% show_best("rmse", n=10)
linear_tune %>% show_best("ccc", n=10)
linear_tune %>% show_best("huber_loss", n=10)


linear_best_parms <- linear_tune %>% select_best("rmse") # Select best parameters based on RMSE


## Final model ##


# Update hyperparameters in recipe and model 
linear_final_rec <- finalize_recipe(model_rec, linear_best_parms)
linear_final_model <- finalize_model(linear_spec, linear_best_parms)

# Create new workflow for final model
linear_final_wf <- workflow() %>% 
  add_recipe(linear_final_rec) %>%
  add_model(linear_final_model)


# Run final model with best hyperparameters and evaluate on the test set
set.seed(678)
linear_final <- linear_final_wf %>%
  last_fit(data_split,
           metrics=metrics)


collect_metrics(linear_final) # Metrics of model (test set)



linear_tune %>% collect_predictions() %>% 
  inner_join(linear_best_parms) %>% 
  metrics(truth=endpoint_value, estimate=.pred) # Metrics of model (training set)

# Actual vs predicted plot
collect_predictions(linear_final) %>% ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)

# Amount within and outside of 0.5 residual threshold (test set)
collect_predictions(linear_final) %>% mutate(offset=endpoint_value-.pred,
                                             freq=ifelse(offset>=0.501 | offset<=-0.501, "OUT", "IN")) %>% 
  group_by(freq) %>% summarise(count=n()/nrow(data_test))















