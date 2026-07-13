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
knn_spec <- null_model() %>% 
  set_mode("regression") %>%
  set_engine("parsnip")

# Create workflow with pre-processing recipe and model specification
knn_workflow <- workflow() %>% 
  add_recipe(model_rec) %>% 
  add_model(knn_spec) 


## Initial hyperparameter tuning ## 

# Set hyperparameter ranges
param <- extract_parameter_set_dials(knn_workflow) %>% 
  update(impute = neighbors(c(50, 100)))

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
                      grid=30,
                      control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                      metrics=metrics,
                      param_info=param)

stopCluster(cl) # Shut down cluster


# Plot hyperparameters to choose best range for final tuning
knn_tune %>% collect_metrics() %>% 
  filter(.metric=="rsq") %>% 
  select(mean, impute) %>% 
  pivot_longer(impute, values_to="value",
               names_to="parameter") %>% 
  ggplot(aes(value, mean, color=parameter)) + 
  geom_point() + 
  facet_grid(~parameter, scales="free_x")

knn_tune %>% collect_metrics() %>% 
  filter(.metric=="rmse") %>% 
  select(mean, impute) %>% 
  pivot_longer(impute, values_to="value",
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
  add_model(knn_spec)


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


data_test %>% mutate(.pred=mean(endpoint_value)) %>% select(endpoint_value, .pred) %>% 
  metrics(truth=endpoint_value,
          estimate=.pred)



a <- data_test %>% group_by(material) %>% summarise(.pred=mean(endpoint_value)) 


data_test %>% left_join(a, by=c("material"="material")) %>% select(endpoint_value, .pred) %>% 
  metrics(truth=endpoint_value,
          estimate=.pred)


data_test %>% left_join(a, by=c("material"="material")) %>% select(endpoint_value, .pred) %>% 
  ggplot(aes(x=endpoint_value, y=.pred)) + geom_point() + geom_abline() + 
  xlim(-4.5, 4.5) + ylim(-4.5, 4.5) + geom_abline(intercept=0.5, linetype=2) + geom_abline(intercept=-0.5, linetype=2)










