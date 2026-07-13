library(tidyverse) # Collection of packages for data wrangling + plotting etc. 
library(tidymodels) # Collection of packages for machine learning
library(themis) # Balancing data
library(vip) # Variable importance
library(doParallel) # Parallel processing
library(doRNG) # Random seed generator during parallel processing
library(tictoc) # Timing script
library(furrr) # Parallel processing for purrr functions


## Set working directory

setwd("...")




## Create new data to test

## Load raw data

load("./2_scripts/clean_data.RData")


data_analysis <- data_final %>% 
  filter(!material %in% c("Ag", "ZnO", "TiO2", "SiO2", "CeO2", "CuO", 
                          "Cu", "Fe2O3", "Fe3O4", "Al2O3", "Se", "NiO") |
           species_group %in% c("Nematoda", "Gastropoda")) %>% 
  mutate(endpoint_value=as.numeric(endpoint_value)) %>% 
  drop_na(endpoint_value) %>% 
  filter(suspension_aging=="no", pre_illumination=="no", uv_radiation=="no") %>% 
  filter(species_group!="Plant") %>% 
  filter(genus!="Poterioochromonas" & epithet!="malhamensis") %>% 
  filter(endpoint=="mortality" | endpoint=="growth inhibition" | endpoint=="immobilization" | 
           endpoint=="cell viability")


# Load solubility classes (from D. magna paper)
# Add extra data for new materials (based on literature)
solubility <- read.csv("./1_data/other/solubility.csv")

solubility <- solubility %>% bind_rows(data.frame(material=c("Al2O3", "Pt", "Ni", "Se", "NiO"), 
                                                  solubility_group=c("slow", "slow", "slow", "slow", "slow")))


# Join solubility classes to dataset
# Fill in "slow" for all remaining materials
data_analysis <- data_analysis %>% left_join(solubility, by=c("material"="material"))

data_analysis <- data_analysis %>% mutate(solubility_group=replace_na(solubility_group, "slow"))

# Create strata for resampling
# Fill in life stage
data_analysis <- data_analysis %>% 
  mutate(species=paste(genus, epithet, sep=" "),
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
                              species=="Barbonymus gonionotus" & str_detect(life_stage, "2 month old") ~ "juvenile"))



# Join molecular descriptors to dataset
# Remove unnecesary variables

ochem_all <- read.csv("./1_data/mol_descriptors/molecular_descriptors_alvadesc_all_materials.csv")




ochem_all <- ochem_all %>% select("material", "MW", "Sv", "Se", "Sp", "Mv", "Me", "GD", "NssO", 
                                  "MRcons", "Vx", "VvdwZAZ", "SAscore", "AMW", "Mp", "Mi", "nTA")

data_analysis <- data_analysis %>% left_join(ochem_all, by=c("material"="material")) %>% 
  select(-c(primary_method, fed_during_test, endpoint_unit, author:epithet,
            length_juvenile:feeding_strategy))


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
                            endpoint=="growth inhibition" ~ "growth inhibition"))



## Select only Au, Ni, Pt material and Nematoda, Gastropoda observations
## Remove zeta potential because of proportion of missing data
data_analysis <- data_analysis %>% select(-zeta_potential) %>% 
  filter(species!="Porphyridium purpureum") %>%
  drop_na("MW") %>%
  mutate(material_group=case_when(str_detect(material, "O")==TRUE ~ "metal oxide",
                                  .default="metal")) %>%
  mutate(across(where(is_character), as_factor))


new_data <- data_analysis


###############################################################################################################################




load("./learning_curve.RData")



rf_parms <- rf_res %>% select(size:impute) %>% unique()

models <- list()

for (i in 1:nrow(rf_parms)) {
  # Split data into 60% training and 40% test set
  set.seed(123)
  
  data_split <- initial_split(data_analysis, strata=strata,
                              prop=0.6)
  
  
  set.seed(rf_parms[i, ]$set.seed)
  
  data_train <- training(data_split) %>% slice_sample(n=rf_parms[i, ]$size)
  data_test <- testing(data_split)
  
  train_full <- data_train
  test_full <- data_test
  
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
    step_impute_knn(all_predictors(), neighbors=rf_parms[i, ]$impute) %>%
    step_dummy(all_nominal_predictors(), one_hot = T)
  
  
  ################################################################################################################################
  
  
  ## Randomforest ##
  
  
  
  # Specify model and engine + variable importance method
  rf_spec <- rand_forest(trees=tune(), 
                         min_n=tune(), 
                         mtry=tune()) %>%
    set_mode("regression") %>% 
    set_engine("ranger")
  
  
  # Create workflow with pre-processing recipe and model specification
  rf_workflow <- workflow() %>% 
    add_recipe(model_rec) %>% 
    add_model(rf_spec) 
  
  
  
  
  ## Initial hyperparameter tuning ## 
  
  # Set hyperparameter ranges
  param <- extract_parameter_set_dials(rf_workflow) %>% update(mtry = mtry(c(rf_parms[i, ]$mtry, rf_parms[i, ]$mtry)),
                                                               min_n = min_n(c(rf_parms[i, ]$min_n, rf_parms[i, ]$min_n)),
                                                               trees = trees(c(rf_parms[i, ]$trees, rf_parms[i, ]$trees)))
  
  # Set metrics for evaluation
  
  metrics <- metric_set(rmse, mae, rsq, ccc, huber_loss, rpiq, rpd, iic)
  
  # Parallel processing to optimize tuning (Note: can cause R to crash) 
  ncores <- 4
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # Initial tuning
  set.seed(456)
  rf_tune <- tune_grid(rf_workflow, 
                       resamples=data_folds, 
                       grid=100,
                       control=control_grid(save_pred=T, allow_par=T, save_workflow=T),
                       metrics=metrics,
                       param_info=param)
  
  
  stopCluster(cl) # Shut down cluster
  
  
  
  
  # Plot hyperparameters to choose best range for final tuning
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
  
  
  models[[i]] <- tibble(model=list(rf_final),
         tune=list(rf_tune),
         size=rf_parms[i, ]$size,
         set.seed=rf_parms[i, ]$set.seed,
         mtry=rf_parms[i, ]$mtry,
         trees=rf_parms[i, ]$trees,
         min_n=rf_parms[i, ]$min_n,
         impute=rf_parms[i, ]$impute,
         train_data=list(train_full),
         test_data=list(test_full))
  
  
}


models <- bind_rows(models)


hierarch_metrics <- function(i) {
  
  int_val <- models %>% slice(i) %>% select(tune) %>% pull(tune) %>% pluck(1) %>% collect_metrics()
  
  int_val_materials <- models %>% slice(i) %>% select(tune) %>% pull(tune) %>% pluck(1) %>% collect_predictions() %>% 
                           select(.pred, .row, endpoint_value) %>% arrange(.row) %>% 
    bind_cols(models %>% slice(i) %>% select(train_data) %>% 
                unnest(train_data) %>% 
                select(species_group, material, species)) %>% 
    group_by(material) %>% 
    metrics(., truth=endpoint_value, estimate=.pred)
  
  int_val_species <- models %>% slice(i) %>% select(tune) %>% pull(tune) %>% pluck(1) %>% collect_predictions() %>% 
    select(.pred, .row, endpoint_value) %>% arrange(.row) %>% 
    bind_cols(models %>% slice(i) %>% select(train_data) %>% 
                unnest(train_data) %>% 
                select(species_group, material, species)) %>% 
    group_by(species_group) %>% 
    metrics(., truth=endpoint_value, estimate=.pred)
  
  ext_val <- models %>% slice(i) %>% select(model) %>% pull(model) %>% pluck(1) %>% collect_metrics()
  
  ext_val_materials <- models %>% slice(i) %>% select(model) %>% pull(model) %>% pluck(1) %>% collect_predictions() %>% 
    select(.pred, .row, endpoint_value) %>% 
    bind_cols(models %>% slice(i) %>% select(test_data) %>% 
                unnest(test_data) %>% 
                select(species_group, material, species)) %>% 
    group_by(material) %>% 
    metrics(., truth=endpoint_value, estimate=.pred)
  
  ext_val_species <- models %>% slice(i) %>% select(model) %>% pull(model) %>% pluck(1) %>% collect_predictions() %>% 
    select(.pred, .row, endpoint_value) %>% 
    bind_cols(models %>% slice(i) %>% select(test_data) %>% 
                unnest(test_data) %>% 
                select(species_group, material, species)) %>% 
    group_by(species_group) %>% 
    metrics(., truth=endpoint_value, estimate=.pred)
  
  val_new_data <- models %>% slice(i) %>% select(model) %>% pull(model) %>% pluck(1) %>% 
    extract_workflow() %>% predict(new_data) %>% 
    bind_cols(new_data %>% select(endpoint_value)) %>% metrics(truth=endpoint_value,
                                                                    estimate=.pred)
  
  
  results <- tibble(int_val=list(int_val),
                    int_val_materials=list(int_val_materials), 
                    int_val_species=list(int_val_species),
                    ext_val=list(ext_val),
                    ext_val_materials=list(ext_val_materials), 
                    ext_val_species=list(ext_val_species),
                    val_new_data=list(val_new_data))
  
  return(results)
  
}



models <- map(.x=seq(1:nrow(models)), .f=\(x) hierarch_metrics(x)) %>% 
  bind_rows() %>% 
  bind_cols(models) %>% 
  select(set.seed, size, mtry, trees, min_n, impute, 
         model, tune, train_data, test_data, int_val, int_val_materials, int_val_species,
         ext_val, ext_val_materials, ext_val_species, val_new_data)





models %>% select(int_val, impute, set.seed, size) %>% 
  unnest(int_val) %>% 
  select(set.seed, size, mtry, trees, min_n, impute, .metric, mean) %>% 
  rename(.estimate=mean) %>% 
  mutate(type="train_set") %>% 
  bind_rows(models %>% select(set.seed, size, mtry, trees, min_n, impute, ext_val) %>% 
              unnest(ext_val) %>% 
              select(-c(.estimator, .config)) %>% 
              mutate(type="test_set"), 
            models %>% select(set.seed, size, mtry, trees, min_n, impute, int_val_materials) %>% 
              unnest(int_val_materials) %>% 
              select(-c(.estimator)) %>% 
              group_by(set.seed, size, mtry, trees, min_n, impute, .metric) %>% 
              summarise(.estimate=mean(.estimate, na.rm=T)) %>%
              mutate(type="train_set_materials"),
            models %>% select(set.seed, size, mtry, trees, min_n, impute, int_val_species) %>% 
              unnest(int_val_species) %>% 
              select(-c(.estimator)) %>% 
              group_by(set.seed, size, mtry, trees, min_n, impute, .metric) %>% 
              summarise(.estimate=mean(.estimate, na.rm=T)) %>%
              mutate(type="train_set_species"),
            models %>% select(set.seed, size, mtry, trees, min_n, impute, ext_val_materials) %>% 
              unnest(ext_val_materials) %>% 
              select(-c(.estimator)) %>% 
              group_by(set.seed, size, mtry, trees, min_n, impute, .metric) %>% 
              summarise(.estimate=mean(.estimate, na.rm=T)) %>%
              mutate(type="test_set_materials"),
            models %>% select(set.seed, size, mtry, trees, min_n, impute, ext_val_species) %>% 
              unnest(ext_val_species) %>% 
              select(-c(.estimator)) %>% 
              group_by(set.seed, size, mtry, trees, min_n, impute, .metric) %>% 
              summarise(.estimate=mean(.estimate, na.rm=T)) %>%
              mutate(type="test_set_species"),
            models %>% select(set.seed, size, mtry, trees, min_n, impute, val_new_data) %>% 
              unnest(val_new_data) %>% 
              select(-c(.estimator)) %>% 
              group_by(set.seed, size, mtry, trees, min_n, impute, .metric) %>% 
              summarise(.estimate=mean(.estimate, na.rm=T)) %>%
              mutate(type="new_data")) %>% 
  group_by(size, .metric, type) %>% summarise(n=n(),
                                     std_error=(sd(.estimate)/sqrt(n)),
                                     .estimate=mean(.estimate)) %>% 
  filter(.metric=="rmse") %>%
  ggplot(aes(size, .estimate, col=type)) + 
  geom_errorbar(aes(ymin=.estimate-std_error, ymax=.estimate+std_error), col="black", linewidth=0.5,
                width=5) +
  geom_point() + geom_line()






models %>% slice(41) %>% select(model) %>% pull(model) %>% pluck(1) %>% extract_workflow() %>% 
  predict(new_data %>% filter(!material %in% c("Au", "Ag", "Ni", "ZnO", "Pt"))) %>% 
  bind_cols(new_data %>% filter(!material %in% c("Au", "Ag", "Ni", "ZnO", "Pt")) %>% select(endpoint_value, species_group, material)) %>% 
  ggplot(aes(endpoint_value, .pred)) + geom_point() + geom_abline() + 
  ylim(-4.5, 4.5) + xlim(-4.5, 4.5)




new_data %>% filter(material!="Au")




