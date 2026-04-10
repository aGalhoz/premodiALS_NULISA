### Libraries
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyr)
library(cowplot)
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(skimr)
library(ggExtra)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(limma)
library(showtext)
showtext_auto()  # ensures UTF-8 + font rendering
install.packages("ggnewscale")
library(ggnewscale)
library(glmnet)
library(pROC)
library(caret)
library(RColorBrewer)
library(scales)

### Directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_plasma'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_CSF'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_SERUM'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/PCA_plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/volcano_plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/signed_plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/NPQ_fluid_plate'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/ML'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/ML/high_detectable'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/ML/without_NEFL'), showWarnings = FALSE)

### Collect data
# new documentation from 12-03-2026
GeneralDocumentation <- read_delim("data input/export-2026-03-24-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_HRO_KSSGCH_MRI_NIUSASSK_267participants/GeneralDocumentation.csv", 
           delim = ";", escape_double = FALSE, trim_ws = TRUE)

GeneralDocumentation_old <- read_delim("data input/export-2025-08-27-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_HRO_KSSGCH_MRI_NIUSASSK/GeneralDocumentation.csv", 
                                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

# IDS of patients
# -> ALS patients
ALS_ID <- GeneralDocumentation %>%
  filter((ALSuncertainty == 1 | ALSFUdiagnosis == 1) & ALSFUdiagnosis %in% c(1,3,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("ALS"))

ALS_ID_old <- GeneralDocumentation_old %>%
  filter((ALSuncertainty == 1 | ALSFUdiagnosis == 1) & ALSFUdiagnosis %in% c(1,3,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("ALS"))

ALS_ID[!ALS_ID$ParticipantCode %in% ALS_ID_old$ParticipantCode,]
ALS_ID_old[!ALS_ID_old$ParticipantCode %in% ALS_ID$ParticipantCode,]

ALS_ID <- ALS_ID %>%
  filter(!ParticipantCode %in% c("DE320","TR310","TR302"))

# -> CTR 
CTR_ID <- GeneralDocumentation %>%
  filter(PGMC == 2) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("CTR"))

# -> PGMC
PGMC_ID <- GeneralDocumentation %>%
  filter(PGMC == 1) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("PGMC"))

# -> Mimic
mimic_ID <- GeneralDocumentation %>%
  filter((ALSuncertainty == 2 | ALSFUdiagnosis == 2) & ALSFUdiagnosis %in% c(2,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("mimic"))


mimic_ID_old <- GeneralDocumentation_old %>%
  filter((ALSuncertainty == 2 | ALSFUdiagnosis == 2) & ALSFUdiagnosis %in% c(2,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("mimic"))

# -> SYMP
SYMP_ID <- GeneralDocumentation %>%
  filter(PGMC == 3 & ALSuncertainty == 3 & ALSFUdiagnosis %in% c(3,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = "SYMP")

# Combine all assigned PatientIDs
assigned_IDs <- bind_rows(ALS_ID, CTR_ID, PGMC_ID, mimic_ID, SYMP_ID) %>%
  select(PatientID)

# Find remaining/unassigned patients
remaining_IDs <- GeneralDocumentation %>%
  filter(!PatientID %in% assigned_IDs$PatientID)

# Count how many are left
nrow(remaining_IDs)
# 7

NA_ID <- remaining_IDs %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = "NA")

# PGMC IDs with mutations
PGMC_mutations_ID <- GeneralDocumentation %>%
  # Select relevant columns: all MutationType* + all PreciseMutation* + PatientID etc.
  select(PatientID, ParticipantCode, PGMC, contains("MutationType"), starts_with("PreciseMutation")) %>%
  # Keep only patients in PGMC_ID
  filter(PatientID %in% PGMC_ID$PatientID) %>%
  # Pivot longer all MutationType columns
  pivot_longer(
    cols = starts_with("MutationType"),
    names_to = "mutation_column",
    values_to = "value"
  ) %>%
  # Keep only rows where mutation is present
  filter(value == 1) %>%
  # Extract mutation name from column name
  mutate(type = sub("^MutationType", "", mutation_column)) %>%
  # Rowwise to pick the correct PreciseMutation column per row
  rowwise() %>%
  mutate(
    PreciseMutationCol = if(type == "Other") "PreciseMutationOther" else paste0("PreciseMutation", type),
    # Safely pick the value if column exists, else NA
    PreciseMutation = if(PreciseMutationCol %in% names(cur_data())) {
      cur_data()[[PreciseMutationCol]]
    } else {
      NA_character_
    }
  ) %>%
  ungroup() %>%
  # Keep only relevant columns
  select(PatientID, ParticipantCode, type, PreciseMutation)

PGMC_mutations_ID <- PGMC_mutations_ID %>%
  mutate(
    mutation = if_else(
      type != "Other",
      type,
      # For "Other", extract the first word of PreciseMutation
      word(PreciseMutation, 1)
    )
  )

# ALS IDs with mutations
ALS_mutations_ID <- GeneralDocumentation %>%
  # Select relevant columns: PatientID, ParticipantCode, PGMC + all MutationType* + PreciseMutation* columns
  select(PatientID, ParticipantCode, PGMC, contains("MutationType"), starts_with("PreciseMutation")) %>%
  # Keep only ALS patients
  filter(PatientID %in% ALS_ID$PatientID) %>%
  # Pivot longer all MutationType columns
  pivot_longer(
    cols = starts_with("MutationType"),
    names_to = "mutation_column",
    values_to = "value"
  ) %>%
  # Keep only rows where mutation is present
  filter(value == 1) %>%
  # Extract mutation name from column name
  mutate(type = sub("^MutationType", "", mutation_column)) %>%
  # Rowwise to pick the correct PreciseMutation column per row
  rowwise() %>%
  mutate(
    PreciseMutationCol = if(type == "Other") "PreciseMutationOther" else paste0("PreciseMutation", type),
    # Safely pick the value if column exists, else NA
    PreciseMutation = if(PreciseMutationCol %in% names(cur_data())) {
      cur_data()[[PreciseMutationCol]]
    } else {
      NA_character_
    }
  ) %>%
  ungroup() %>%
  # Create the mutation column
  mutate(
    mutation = if_else(
      type != "Other",
      type,
      # For "Other", take first word of PreciseMutation
      word(PreciseMutation, 1)
    )
  ) %>%
  # Keep only relevant columns
  select(PatientID, ParticipantCode, type, PreciseMutation, mutation)

ALS_mutations_ID <- ALS_mutations_ID %>%
  filter(!ParticipantCode %in% c("SK307", "SK309", "SK329"))
# Remove unknown ALS mutation

# Classify DE309 as SOD1 (He is SOD1 and FIG4)
ALS_mutations_ID <- ALS_mutations_ID %>%
  filter(
    !(ParticipantCode == "DE309" & type != "SOD1")  # remove DE309 rows where type is not SOD1
  )

all_participants_IDs = do.call("rbind",
                               list(CTR_ID,
                                    ALS_ID,
                                    PGMC_ID,
                                    mimic_ID,
                                    SYMP_ID,
                                    NA_ID,
                                    PGMC_mutations_ID %>%
                                      select(PatientID,ParticipantCode,mutation) %>%
                                      rename(type = mutation)))
all_participants_IDs <- rbind(all_participants_IDs,
                              c("XH4W23T7","FR108","C9orf72"))
                              
writexl::write_xlsx(all_participants_IDs,"results/all_participants_IDs.xlsx")

# information on sample IDs
sample_ID_info <- read_excel("data input/Biospecimen Manifest Form – Banner Biomarker Program_premodiALS.xlsx", 
          sheet = "Sample Information") %>%
  select("Sample ID","Optional Informarion (patient ID)" )
# -> change NA2 info to DE107
sample_ID_info$`Optional Informarion (patient ID)` <- ifelse(sample_ID_info$`Optional Informarion (patient ID)` == "NA2",
                         "DE107",sample_ID_info$`Optional Informarion (patient ID)`)
# -> change samples of size 4 with TR116 to TR126 (writing mistake)
sample_ID_info$`Optional Informarion (patient ID)` <- ifelse(sample_ID_info$`Sample ID` %in% c("9772","9636","9626"),"TR126",
       sample_ID_info$`Optional Informarion (patient ID)`)

# data from CSF, plasma and serum fluids
protein_data <- read_excel("data input/P004_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx")

# new data from Antonia
protein_data_updated <- read_excel("data input/UPDATED_P004_NPQ.xlsx")
protein_data_old <- protein_data_updated

# new data from Marisa 04-09-2025
protein_data_updated <- read_excel("data input/Updated_P004_NPQ_09042025.xlsx")

protein_data <- protein_data_updated

# target detectability
target_detectability <- read_excel("data input/Updated_P004_NPQ_09042025.xlsx", 
           sheet = "Target Detectability")

# Patient age and sex
participants_PGMC = read_excel("data input/recruited participants.xlsx",
                               sheet = "Total PGMC") %>%
  select(Pseudonyme,sex,age)
participants_CTR = read_excel("data input/recruited participants.xlsx",
                              sheet = "Total Control") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_ALS_mimic = read_excel("data input/recruited participants.xlsx",
                                    sheet = "Total EALS-Mimics") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)

participants_Turkey = read_excel("data input/recruited participants.xlsx",
                                     sheet = "Turkey") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_Slovakia = read_excel("data input/recruited participants.xlsx",
                                 sheet = "Slovakia") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_Germany = read_excel("data input/recruited participants.xlsx",
                                   sheet = "Germany - Munich") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_Switzerland = read_excel("data input/recruited participants.xlsx",
                                   sheet = "Switzerland") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_Israel = read_excel("data input/recruited participants.xlsx",
                                   sheet = "Israel") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_France = read_excel("data input/recruited participants.xlsx",
                                   sheet = "France") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)
participants_Czech = read_excel("data input/recruited participants.xlsx",
                                   sheet = "Czech Republic") %>%
  select(Pseudonyme,Sex,`Age...7`) %>%
  rename(sex = Sex,
         age = `Age...7`)

# Sex_age_all_participants = do.call("rbind",list(participants_PGMC,
#                                                 participants_CTR,
#                                                 participants_ALS_mimic))

Sex_age_all_participants = do.call("rbind",list(participants_Czech,
                                                participants_France,
                                                participants_Germany,
                                                participants_Israel,
                                                participants_Turkey,
                                                participants_Slovakia,
                                                participants_Switzerland))

all_participants_IDs_final <- all_participants_IDs %>%
  # Join PGMC mutations
  left_join(
    PGMC_mutations_ID %>% select(PatientID, mutation) %>% rename(PGMC_mutation = mutation),
    by = "PatientID"
  ) %>%
  # Manually fix FR108
  mutate(
    PGMC_mutation = if_else(ParticipantCode == "FR108", "C9orf72", PGMC_mutation)
  ) %>%
  # Join ALS mutations
  left_join(
    ALS_mutations_ID %>% select(PatientID, mutation) %>% rename(ALS_mutation = mutation),
    by = "PatientID"
  ) %>%
  # Combine both into a single mutation column
  mutate(
    mutation = coalesce(PGMC_mutation, ALS_mutation)  # take PGMC if present, else ALS
  ) %>%
  # Optional: remove intermediate columns
  select(-PGMC_mutation, -ALS_mutation)

# Questionnaire info with ALSFRS scores
Questionnaire <- read_delim("data input/export-2026-03-20-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_HRO_KSSGCH_MRI_NIUSASSK_267participants/QuestionnaireG.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

# ECAS info
ECAS = read_delim("data input/export-2026-03-20-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_HRO_KSSGCH_MRI_NIUSASSK_267participants/Ecas.csv", 
                  delim = ";", escape_double = FALSE, trim_ws = TRUE)

# spinal/bulbar 
site_onset = read_excel("data input/all_participants_IDs.xlsx")

# disease duration
disease_duration =  read_excel("data input/all_participants_ID_visits.xlsx")

# maybe better to take the sex, age from the list of Laura - since the excel sheet from where the dates were coming was not updated
Sex_age_all_participants = site_onset %>%
  select(PatientID,sex,age) %>%
  rename(Pseudonyme = PatientID)

all_participants_IDs_final <- merge(all_participants_IDs_final,
                                    Sex_age_all_participants %>% rename(PatientID = Pseudonyme)) %>%
  filter(type %in% c("ALS","CTR","mimic","PGMC","SYMP",NA,"NA"))

writexl::write_xlsx(all_participants_IDs_final,"results/all_participants_IDs_final.xlsx")

