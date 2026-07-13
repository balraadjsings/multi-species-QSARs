library(readxl) # Opening excel files
library(tidyverse)

## Set working directory

setwd("...")




## Read data and select column headers

data <- read_xlsx("./1_data/nanotox_database/nanotox_database_analysis.xlsx", sheet="data", col_names=F)

## Set headers and remove first two rows
colnames(data) <- data[2,]
data <- data[-c(1,2),]



## Select relevant columns

data <- data %>% select(2, 7:8, 10:11, 14, 19:28, 31:33, 35:37, 40, 45:46, 48, 50:55, 58:59, 62:67)

## Clean up data (remove unwanted characters and harmonize all size data)
## Size data is harmonized by taking averages from size ranges

data <- data %>% mutate(`primary_diameter (nm)`= str_replace_all(`primary_diameter (nm)`, "<", "0-") %>% 
                          str_replace_all(" \\(.*\\)", "") %>% str_remove_all("~"),
                        `primary_length (nm)`=str_remove_all(`primary_length (nm)`," \\(.*\\)") %>% 
                          ifelse(!shape_category %in% c("rod", "wire", "rod (image unclear)") & is.na(.), 0, .),
                        `surface_area (m^2/g)`=str_replace_all(`surface_area (m^2/g)`, "<", "0-") %>% 
                          str_replace_all(" \\(.*\\)", "") %>% str_remove_all("~")) %>% 
  separate(`primary_diameter (nm)`, c("diam_low", "diam_high"), extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(`primary_length (nm)`, c("length_low", "length_high"), extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(`surface_area (m^2/g)`, c("surf_low", "surf_high"), extra = "merge", fill = "left", remove=T, sep="-") %>% 
  mutate_at(c("diam_low", "diam_high", "length_low", 
              "length_high", "surf_low", "surf_high"), as.numeric) %>% 
  mutate(diameter=rowMeans(cbind(diam_low, diam_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         length=rowMeans(cbind(length_low, length_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         surface_area=rowMeans(cbind(surf_low, surf_high), na.rm=T) %>% ifelse(is.nan(.), NA, .)) %>% 
  select(-c(surf_low:diam_high))


## Create data frame with density data for calculation of specific surface area

density <- data.frame(material=unique(data$material), density=c(7.133, 10.5, 19.3, 5.606, 4.005, 
                                                                6.315, 7.55, 8.90, 7.215, 8.96, 
                                                                5.17, 6.07, 7.16, 3.6, 5.435, 
                                                                12.023, 4.86, 3.987, 2.422, 6.67,
                                                                6.51, NA, NA, NA, NA,
                                                                NA, NA, NA, NA, 8.90,
                                                                4.49, 5.25, 6.95, NA, NA,
                                                                NA, 2.70, 8.908, 4.50, 21.45,
                                                                5.68, 7.60, NA, 5.28, 5.22,
                                                                NA, 2.3446, 6.517, 19.25, 4.826, 
                                                                NA, 7.15, 7.179, 4.090, 6.02,
                                                                2.3, 7.234, 7.14, 5.010))

## Join density data with dataset

data <- data %>% left_join(density, by=c("material"="material")) %>% 
  mutate(shape_category=shape_category %>% str_replace_all(" \\(.*\\)", ""))



## Calculate specific surface area for irregular, spherical, nearly spherical and unknown shapes

data_2 <- data %>% filter(shape_category %in% c("irregular", "spherical", "nearly spherical", "unknown")) %>% 
  mutate(surface_area_calculated=6000/(density*diameter))

## Plot predicted vs observed to see how well calculated SSA matches with available data

ggplot(data_2, aes(surface_area, surface_area_calculated)) + geom_point() + geom_abline(slope=1,intercept=0) + 
  facet_grid(. ~ shape_category)

## Calculate specific surface area for rod and wire shapes (based on equivalent spherical diameter)

data_3 <- data %>% filter(shape_category %in% c("rod", "wire")) %>% 
  mutate(eq_diameter = ((3/4*diameter^2*length)^(1/3))*2,
         surface_area_calculated=6000/(density*eq_diameter))

## Plot predicted vs observed to see how well calculated SSA matches with available data

ggplot(data_3, aes(surface_area, surface_area_calculated)) + geom_point() + geom_abline(slope=1,intercept=0) + 
  facet_grid(. ~ shape_category) + xlim(0, 200) + ylim(0, 200)

## Calculate specific surface area for rod and wire shapes (based on only diameter)

data_4 <- data %>% filter(shape_category %in% c("rod", "wire")) %>% 
  mutate(eq_diameter = ((3/4*diameter^2*length)^(1/3))*2,
         surface_area_calculated=6000/(density*diameter))

## Plot predicted vs observed to see how well calculated SSA matches with available data

ggplot(data_4, aes(surface_area, surface_area_calculated)) + geom_point() + geom_abline(slope=1,intercept=0) + 
  facet_grid(. ~ shape_category) + xlim(0, 200) + ylim(0, 200) 

rm(data_2, data_3, data_4, density)


########################################################################################################

## Based on predicted vs observed results, the SSA was filled in for the shapes based on the normal SSA equation
## without considering length (no conversion of diameter to equivalent diameter for rod and wire shapes)

## Clean up the rest of the dataset and convert endpoint values to mg/l and M

data_2 <- data %>% mutate(surface_area=if_else(is.na(surface_area),
                                               case_when(shape_category %in% 
                                                           c("irregular", "spherical", "nearly spherical", 
                                                             "unknown", "rod", "wire") ~ 6000/(density*diameter)), 
                                               surface_area)) %>% 
  mutate(`crystallinity_anatase (%)`=ifelse(`crystallinity_anatase (%)`=="anatase + some rutile", 
                                            NA_character_, `crystallinity_anatase (%)`),
         `crystallinity_rutile (%)`=ifelse(`crystallinity_rutile (%)`=="anatase + some rutile", 
                                           NA_character_, `crystallinity_rutile (%)`),
         zeta_potential=case_when(str_detect(zeta_potential, "\\(.*\\)") ~ NA_character_,
                                  str_detect(zeta_potential, "^[A-Za-z]") ~ NA_character_,
                                  .default=zeta_potential) %>% str_remove_all(., "~"),
         `dispersion (immediately before)`=str_remove_all(`dispersion (immediately before)`, " \\(.*\\)"),
         `dispersion (while before experiment / unknown)`=str_remove_all(`dispersion (while before experiment / unknown)`, " \\(.*\\)"),
         fed_during_test=str_remove_all(fed_during_test, " \\(.*\\)"),
         `test_duration (h)`=str_replace_all(`test_duration (h)`, "1-4", NA_character_),
         suspension_aging=suspension_aging %>% str_replace_all("no \\(see comments\\)", "no") %>% 
           ifelse(str_detect(., "^no")==FALSE, "yes", .),
         nat_org_matter_binary=nat_org_matter_binary %>% str_remove_all(" \\(.*\\)"),
         `ionic strength (mol/l)` = `ionic strength (mol/l)` %>% 
           str_replace_all(".*\\(.*\\)", NA_character_),
         `conductivity (µS/cm)`=`conductivity (µS/cm)` %>% str_remove_all("../cm"), 
         test_ph=test_ph %>% str_replace_all(">.*", NA_character_),
         test_temperature=test_temperature %>% str_replace_all("room temperature", "20"), 
         shaken_during_experiment=shaken_during_experiment %>% str_remove_all(" \\(.*\\)"), 
         endpoint=endpoint %>% str_remove_all(" \\(.*\\)")) %>% 
  filter(!endpoint_unit %in% c("particles/ml", "µg/g")) %>% 
  mutate(unit_converted=case_when(endpoint_unit=="mg/l" ~ as.numeric(endpoint_value)*1,
                                  endpoint_unit=="µg/ml" ~ as.numeric(endpoint_value)*1,
                                  endpoint_unit=="µg/l" ~ as.numeric(endpoint_value)*0.001,
                                  endpoint_unit=="nM" ~ as.numeric(endpoint_value)*10^-9,
                                  endpoint_unit=="ppm" ~ as.numeric(endpoint_value)*1,
                                  endpoint_unit=="mM" ~ as.numeric(endpoint_value)*0.001,
                                  endpoint_unit=="µM" ~ as.numeric(endpoint_value)*10^-6,
                                  endpoint_unit=="M" ~ as.numeric(endpoint_value)*1,
                                  endpoint_unit=="mg/ml" ~ as.numeric(endpoint_value)*1000,
                                  endpoint_unit=="mol/l" ~ as.numeric(endpoint_value)*1,
                                  endpoint_unit=="g/l" ~ as.numeric(endpoint_value)*1000,
                                  endpoint_unit=="mmol/l" ~ as.numeric(endpoint_value)*0.001,
                                  endpoint_unit=="μg/l" ~ as.numeric(endpoint_value)*0.001,
                                  endpoint_unit=="ng/ml" ~ as.numeric(endpoint_value)*0.001,
                                  endpoint_unit=="ppb" ~ as.numeric(endpoint_value)*0.001),
         unit_converted_2=case_when(endpoint_unit=="mg/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                   str_remove(endpoint_value, ">") %>% 
                                                                     as.numeric(.)*1, NA_character_),
                                    endpoint_unit=="µg/ml" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                    str_remove(endpoint_value, ">") %>% 
                                                                      as.numeric(.)*1, NA_character_),
                                    endpoint_unit=="µg/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                   str_remove(endpoint_value, ">") %>% 
                                                                     as.numeric(.)*0.001, NA_character_),
                                    endpoint_unit=="nM" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                 str_remove(endpoint_value, ">") %>% 
                                                                   as.numeric(.)*10^-9, NA_character_),
                                    endpoint_unit=="ppm" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                  str_remove(endpoint_value, ">") %>% 
                                                                    as.numeric(.)*1, NA_character_),
                                    endpoint_unit=="mM" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                 str_remove(endpoint_value, ">") %>% 
                                                                   as.numeric(.)*0.001, NA_character_),
                                    endpoint_unit=="µM" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                 str_remove(endpoint_value, ">") %>% 
                                                                   as.numeric(.)*10^-6, NA_character_),
                                    endpoint_unit=="M" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                str_remove(endpoint_value, ">") %>% 
                                                                  as.numeric(.)*1, NA_character_),
                                    endpoint_unit=="mg/ml" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                    str_remove(endpoint_value, ">") %>% 
                                                                      as.numeric(.)*1000, NA_character_),
                                    endpoint_unit=="mol/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                    str_remove(endpoint_value, ">") %>% 
                                                                      as.numeric(.)*1, NA_character_),
                                    endpoint_unit=="g/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                  str_remove(endpoint_value, ">") %>% 
                                                                    as.numeric(.)*1000, NA_character_),
                                    endpoint_unit=="mmol/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                     str_remove(endpoint_value, ">") %>% 
                                                                       as.numeric(.)*0.001, NA_character_),
                                    endpoint_unit=="μg/l" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                   str_remove(endpoint_value, ">") %>% 
                                                                     as.numeric(.)*0.001, NA_character_),
                                    endpoint_unit=="ng/ml" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                    str_remove(endpoint_value, ">") %>% 
                                                                      as.numeric(.)*0.001, NA_character_),
                                    endpoint_unit=="ppb" ~ ifelse(endpoint_value %>% str_detect(., ">.*")==TRUE, 
                                                                  str_remove(endpoint_value, ">") %>% 
                                                                    as.numeric(.)*0.001, NA_character_)),
         unit_converted_3=str_c(">", unit_converted_2),
         endpoint_value=ifelse(is.na(unit_converted), unit_converted_3, unit_converted),
         endpoint_unit=ifelse(endpoint_unit %in% c("mg/l", "µg/ml", "µg/l", "ppm", "mg/ml", "g/l",
                                                   "μg/l", "ng/ml", "ppb"), "mg/l", "M")) %>% 
  select(-c("unit_converted", "unit_converted_2", "unit_converted_3"))



## Read molecular weight data and join to dataset to convert EC50 molar values to mg/l

mol_weight <- read.csv("./1_data/other/materials_mol_weight.csv")

data_2 <- data_2 %>% left_join(mol_weight, by=c("material"="materials"))

## Convert EC50 values

data_2 <- data_2 %>% mutate(unit_converted=case_when(
  endpoint_unit=="M" ~ (as.numeric(endpoint_value)*mol_weight)*1000),
  unit_converted_2=case_when(
    endpoint_unit=="M" ~ ifelse(str_detect(endpoint_value, ">.*")==TRUE, 
                                str_remove(endpoint_value, ">"), NA_character_)),
  unit_converted_3=(as.numeric(unit_converted_2)*mol_weight)*1000,
  unit_converted_3=str_c(">", unit_converted_3),
  unit_converted=ifelse(is.na(unit_converted), unit_converted_3, unit_converted),
  endpoint_value=ifelse(endpoint_unit=="M", unit_converted, endpoint_value),
  endpoint_unit=ifelse(endpoint_unit=="M", "mg/l", endpoint_unit)) %>% 
  select(-c(unit_converted, unit_converted_2, unit_converted_3, mol_weight))

data <- data_2

rm(data_2, mol_weight)

########################################################################################################


## Read species traits data

traits <- read_xlsx("./1_data/species_traits/species_traits.xlsx", sheet="clean")
traits_general <- read_xlsx("./1_data/species_traits/general_traits.xlsx", sheet="clean")

## Clean up and split certain columns

traits_2 <- traits %>% select(1:11, 13, 16:18, 25:26) %>% 
  separate_wider_delim(cols=habitat, names_sep="_", too_few = "align_start", delim = ", ") %>% 
  mutate(habitat_1=str_remove(habitat_1, " \\(.*\\)"),
         habitat_2=str_remove(habitat_2, " \\(.*\\)"),
         habitat_3=str_remove(habitat_3, " \\(.*\\)"),
         habitat_4=str_remove(habitat_4, " \\(.*\\)"),
         habitat_5=str_remove(habitat_5, " \\(.*\\)")) %>% 
  mutate(habitat_1=case_when(habitat_1=="freshwater" ~ "freshwater",
                             habitat_1=="brackish" ~ "brackish",
                             habitat_1=="marine" ~ "marine"),
         habitat_2=case_when(habitat_2=="freshwater" ~ "freshwater",
                             habitat_2=="brackish" ~ "brackish",
                             habitat_2=="marine" ~ "marine"),
         habitat_3=case_when(habitat_3=="freshwater" ~ "freshwater",
                             habitat_3=="brackish" ~ "brackish",
                             habitat_3=="marine" ~ "marine"),
         habitat_4=case_when(habitat_4=="freshwater" ~ "freshwater",
                             habitat_4=="brackish" ~ "brackish",
                             habitat_4=="marine" ~ "marine"),
         habitat_5=case_when(habitat_5=="freshwater" ~ "freshwater",
                             habitat_5=="brackish" ~ "brackish",
                             habitat_5=="marine" ~ "marine")) %>% 
  mutate(freshwater=str_count(habitat_1, "freshwater") %>% replace_na(0) + 
           str_count(habitat_2, "freshwater") %>% replace_na(0) + 
           str_count(habitat_3, "freshwater") %>% replace_na(0) + 
           str_count(habitat_4, "freshwater") %>% replace_na(0) + 
           str_count(habitat_5, "freshwater") %>% replace_na(0), 
         brackish=str_count(habitat_1, "brackish") %>% replace_na(0) + 
           str_count(habitat_2, "brackish") %>% replace_na(0) + 
           str_count(habitat_3, "brackish") %>% replace_na(0) + 
           str_count(habitat_4, "brackish") %>% replace_na(0) + 
           str_count(habitat_5, "brackish") %>% replace_na(0),
         marine=str_count(habitat_1, "marine") %>% replace_na(0) + 
           str_count(habitat_2, "marine") %>% replace_na(0) + 
           str_count(habitat_3, "marine") %>% replace_na(0) + 
           str_count(habitat_4, "marine") %>% replace_na(0) + 
           str_count(habitat_5, "marine") %>% replace_na(0)) %>% 
  select(-c(habitat_1:habitat_5)) %>% 
  mutate(`diameter egg (mm)`=`diameter egg (mm)` %>% str_replace_all("check fishbase again", NA_character_)) %>% 
  mutate(across(.cols=`diameter egg (mm)`:`mobility (type of movement)`, 
                .fns= ~ str_remove_all(., " \\(.*\\)")))

## Clean up general traits data and split certain columns

traits_3 <- traits_general %>% select(-c(13:14, 36)) %>% 
  mutate(`antenae (total)`=ifelse(`anntenae (pairs)`=="in total 4 (only 1 is paired)", "4", `anntenae (pairs)`),
         `anntenae (pairs)`=ifelse(`anntenae (pairs)`=="in total 4 (only 1 is paired)", "1", `anntenae (pairs)`),
         `Circulatory system`=case_when(str_detect(`Circulatory system`, "^closed.*") ~ "closed",
                                        str_detect(`Circulatory system`, "^open.*") ~ "open",
                                        str_detect(`Circulatory system`, "^none.*") ~ "none"),
         Respiration=case_when(Respiration=="body (diffusion)" ~ "body/cell diffusion",
                               Respiration=="none (diffusion)" ~ "body/cell diffusion",
                               Respiration=="body (diffusion) + gills" ~ "body/cell diffusion + gills",
                               Respiration=="stomata (diffusion)" ~ "stomata diffusion",
                               .default=Respiration)) %>% 
  separate_wider_delim(cols=Respiration, names_sep="_", too_few = "align_start", delim = " + ") %>% 
  mutate(`body/cell diffusion`=str_count(Respiration_1, "body/cell diffusion") %>% replace_na(0) + 
           str_count(Respiration_2, "body/cell diffusion") %>% replace_na(0),
         gills=str_count(Respiration_1, "gills") %>% replace_na(0) + 
           str_count(Respiration_2, "gills") %>% replace_na(0),
         lungs=str_count(Respiration_1, "lungs") %>% replace_na(0) + 
           str_count(Respiration_2, "lungs") %>% replace_na(0),
         `stomata diffusion`=str_count(Respiration_1, "stomata diffusion") %>% replace_na(0) + 
           str_count(Respiration_2, "stomata diffusion") %>% replace_na(0)) %>% 
  select(-c(Phylum:Family, Respiration_1, Respiration_2))

## Join both dataset

species_traits <- traits_2 %>% drop_na(Phylum) %>% left_join(traits_3, by=c("Genus"="Genus", "epithet"="epithet"))

species_traits <- species_traits %>% select(-c(`life stage`)) %>% distinct()

rm(traits, traits_2, traits_3, traits_general)


########################################################################################################

## Change species name 
data <- data %>% mutate(species=ifelse(species=="Scenedesmus subspicatus", "Desmodesmus subspicatus", species))

## Turn off scientific notation
options(scipen = 999)

## Combine traits and toxicity data (removed all Insecta because it causes issues with left joining)

#data_final <- data %>% filter(species_group!="Insecta") %>% 
  #left_join(species_traits, by=c("species"="species", "life stage"="life stage"))

data_final <- data %>% filter(species_group!="Insecta") %>% 
  left_join(species_traits, by=c("species"="species"))


## Convert columns with scientific notation to numbers

data_final <- data_final %>% mutate(`salinity`=as.numeric(`salinity (‰)`),
                                `water_hardness`=as.numeric(`water hardness (mg/l CaCO3)`),
                                `ionic_strength`=as.numeric(`ionic strength (mol/l)`),
                                `conductivity`=as.numeric(`conductivity (µS/cm)`),
                                endpoint_2=as.numeric(endpoint_value),
                                diameter_egg=as.numeric(`diameter egg (mm)`),
                                length_adult=as.numeric(`length adults (mm)`)) %>% 
  mutate(salinity=ifelse(is.na(salinity), `salinity (‰)`, salinity),
         water_hardness=ifelse(is.na(water_hardness), `water hardness (mg/l CaCO3)`, water_hardness),
         ionic_strength=ifelse(is.na(ionic_strength), `ionic strength (mol/l)`, ionic_strength),
         conductivity=ifelse(is.na(conductivity), `conductivity (µS/cm)`, conductivity),
         endpoint_value=ifelse(is.na(endpoint_2), endpoint_value, endpoint_2),
         diameter_egg=ifelse(is.na(diameter_egg), `diameter egg (mm)`, diameter_egg),
         length_adult=ifelse(is.na(length_adult), `length adults (mm)`, length_adult)) %>% 
  select(-c(`salinity (‰)`, `water hardness (mg/l CaCO3)`, `ionic strength (mol/l)`, `conductivity (µS/cm)`,
            `diameter egg (mm)`, `length adults (mm)`, endpoint_2))

## Split ranges and calculate averages

data_final <- data_final %>% 
  separate(`crystallinity_anatase (%)`, c("crystallinity_anatase_low", "crystallinity_anatase_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(`crystallinity_rutile (%)`, c("crystallinity_rutile_low", "crystallinity_rutile_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(zeta_potential, c("zeta_potential_low", "zeta_potential_high"),
           extra = "merge", fill = "left", remove=T, sep=" to ") %>% 
  separate(test_ph, c("test_ph_low", "test_ph_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(test_temperature, c("test_temperature_low", "test_temperature_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(`length juvenile (mm)`, c("length_juvenile_low", "length_juvenile_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(salinity, c("salinity_low", "salinity_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(water_hardness, c("water_hardness_low", "water_hardness_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(ionic_strength, c("ionic_strength_low", "ionic_strength_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(conductivity, c("conductivity_low", "conductivity_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(diameter_egg, c("diameter_egg_low", "diameter_egg_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  separate(length_adult, c("length_adult_low", "length_adult_high"),
           extra = "merge", fill = "left", remove=T, sep="-") %>% 
  mutate_at(c("crystallinity_anatase_low", "crystallinity_anatase_high", 
              "crystallinity_rutile_low", "crystallinity_rutile_high", 
              "zeta_potential_low", "zeta_potential_high", 
              "test_ph_low", "test_ph_high", "test_temperature_low", "test_temperature_high", 
              "length_juvenile_low", "length_juvenile_high", "salinity_low", "salinity_high", 
              "water_hardness_low", "water_hardness_high", "ionic_strength_low", "ionic_strength_high", 
              "conductivity_low", "conductivity_high", "diameter_egg_low", "diameter_egg_high", 
              "length_adult_low", "length_adult_high"), as.numeric) %>% 
  mutate(crystallinity_anatase=rowMeans(cbind(crystallinity_anatase_low, crystallinity_anatase_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         crystallinity_rutile=rowMeans(cbind(crystallinity_rutile_low, crystallinity_rutile_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         zeta_potential=rowMeans(cbind(zeta_potential_low, zeta_potential_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         test_ph=rowMeans(cbind(test_ph_low, test_ph_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         test_temperature=rowMeans(cbind(test_temperature_low, test_temperature_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         length_juvenile=rowMeans(cbind(length_juvenile_low, length_juvenile_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         salinity=rowMeans(cbind(salinity_low, salinity_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         water_hardness=rowMeans(cbind(water_hardness_low, water_hardness_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         ionic_strength=rowMeans(cbind(ionic_strength_low, ionic_strength_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         conductivity=rowMeans(cbind(conductivity_low, conductivity_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         diameter_egg=rowMeans(cbind(diameter_egg_low, diameter_egg_high), na.rm=T) %>% ifelse(is.nan(.), NA, .),
         length_adult=rowMeans(cbind(length_adult_low, length_adult_high), na.rm=T) %>% ifelse(is.nan(.), NA, .)) %>% 
  select(-c("crystallinity_anatase_low", "crystallinity_anatase_high", 
            "crystallinity_rutile_low", "crystallinity_rutile_high", 
            "zeta_potential_low", "zeta_potential_high", 
            "test_ph_low", "test_ph_high", "test_temperature_low", "test_temperature_high", 
            "length_juvenile_low", "length_juvenile_high", "salinity_low", "salinity_high", 
            "water_hardness_low", "water_hardness_high", "ionic_strength_low", "ionic_strength_high", 
            "conductivity_low", "conductivity_high", "diameter_egg_low", "diameter_egg_high", 
            "length_adult_low", "length_adult_high"))
  

########################################################################################################
  
## Final dataset (relocate and rename columns)


## Reorder columns
data_final <- data_final %>% 
  relocate(c("crystallinity_anatase", "crystallinity_rutile", "surface_area", "length", "diameter"), .after="coating") %>% 
  relocate(c("zeta_potential"), .after="primary_method") %>% 
  relocate(c("salinity", "water_hardness", "ionic_strength", "conductivity", "test_ph",
             "test_temperature"), .after="nat_org_matter_binary") %>% 
  relocate(c("length_adult", "length_juvenile", "diameter_egg"), .after="epithet")

## Rename columns
colnames(data_final) <- c("material", "shape", "coating", "crystallinity_anatase", "crystallinity_rutile",
                      "surface_area", "primary_length", "primary_diameter", "primary_method", "zeta_potential",
                      "dispersion_immediately_before", "dispersion_before_unknown", "test_procedure",
                      "species_group", "species", "life_stage", "fed_during_test", "test_duration",
                      "suspension_aging", "pre_illumination", "illumination_light", "illumination_dark",
                      "uv_radiation", "nat_org_matter", "salinity", "water_hardness", "ionic_strength",
                      "conductivity", "ph", "temperature", "shaking_during_experiment", "endpoint", 
                      "endpoint_unit", "endpoint_value", "author", "year", "title", "doi", "dataset_origin",
                      "other_comments", "density", "phylum", "sub_phylum", "class", "superorder", "order",
                      "informal_group", "family", "genus", "epithet", "length_adult", "length_juvenile", 
                      "diameter_egg", "feeding_strategy", "mobility", "habitat_freshwater", "habitat_brackish",
                      "habitat_marine", "exoskeleton", "endoskeleton", "shell", "moulting_cuticle", 
                      "carapace", "antennae_pairs", "cirri_tentacles", "cilia_tentacle_crown", "fins", "barbles",
                      "adipose_fin", "jointed_thoraic_appendages_pairs", "roots_rhizoid", "stems", "thallus", 
                      "flagella", "cilia", "alveoli", "circulatory_system", "vascular_system", "cellularity",
                      "cell_wall_silica", "prokaryote", "antennae_total", "respiration_body_cell_diffusion", 
                      "respiration_gills", "respiration_lungs", "respiration_stomata_diffusion")

data_final <- data_final %>% select(-c("dataset_origin", "density", "other_comments", "informal_group")) %>% 
  relocate(c("antennae_total"), .after="antennae_pairs") %>% 
  relocate(c("cilia", "alveoli"), .after="cilia_tentacle_crown") %>% 
  relocate(c("respiration_body_cell_diffusion", "respiration_gills", 
             "respiration_lungs", "respiration_stomata_diffusion"), .after="vascular_system")


## Convert character columns to numerical

data_final <- data_final %>% mutate(across(.cols=c("test_duration", "illumination_light",
                                               "illumination_dark"), .fns= ~ as.numeric(.x))) %>% 
  mutate(across(.cols=c("antennae_pairs", "antennae_total"), .fns= ~ as.integer(.x)))

## Fix levels in species_group and convert all characters to factors 
data_final <- data_final %>% 
  mutate(species_group=ifelse(species_group=="crustacea", "Crustacea", species_group))

## Drop species without species traits
data_final <- data_final %>% drop_na(phylum)

str(data_final)


########################################################################################################

## Molecular descriptors (OCHEM)


## Read calculated OCHEM descriptors (AlvaDesc calculator, calculated through OCHEM) --> THIS IS USED IN ANALYSIS

mol_descriptors <- read.csv("./1_data/mol_descriptors/molecular_descriptors_alvadesc.csv")
mol_descriptors_cores <- mol_descriptors %>% mutate(Material=c("TiO2", "ZnO", "Ag", "Cu", "CuO", "Au", "NiO", "SiO2", "Al2O3",
                                                               "Fe2O3", "Fe3O4", "CeO2", "Pt", "Ni", "Se"))

## Reduce variables by removing columns with NAs

mol_descriptors_cores <- mol_descriptors_cores %>% select(where(~ !any(is.na(.))))
mol_descriptors_cores <- mol_descriptors_cores %>% select(-c(CASRN:ERROR, MLOGP:ESOL, BLTF96:BLTA96)) %>% 
  select(c(SMILES, Material, MW:totalcharge, SsCH3:gmax, qpmax:LDI, Uc:PBF))

## Recipe for removing zero variance and high correlation

descr_rec <- recipe( ~ ., data=mol_descriptors_cores) %>%
  update_role(c(Material, SMILES), new_role="id") %>% 
  step_zv(all_predictors()) %>%
  step_corr(all_predictors(), threshold = 0.9, method="pearson")


## Turn recipe into data frame
ochem_descr_processed <- descr_rec %>% prep() %>% juice()



rm(descr_rec, mol_descriptors, mol_descriptors_cores)








































