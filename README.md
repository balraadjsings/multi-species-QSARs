# multi-species-QSARs
Code and data for "Building species trait-specific nano-QSARs: Model stacking, navigating model uncertainties and limitations, and the effect of dataset size" (https://doi.org/10.1016/j.envint.2024.108764)

### 1_data contains data files used within the scripts:
- mol_descriptors: calculated molecular descriptors 
- nanotox_database: toxicity data
- other: solubility grouping and molecular weights
- species_traits: species trait data

### 2_scripts contains R scripts to create models, calculate applicability domain and clean data:
- data_cleaning = script to clean raw data and prepare for modelling
- new_template_models = script for regression QSAR models
- new_template_classification_models = script for classification QSAR models
- new_stack = script for stacking regression and classification models into meta model
- new_model_interpretation = script for variable importance analysis, applicability domain calculation (and partial dependence profiles (not presented in paper))
- learning_curve_2 = script for generating learning curve
- null_model = script for null model

