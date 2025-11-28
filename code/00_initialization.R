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

### Directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_plasma'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_CSF'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_SERUM'), showWarnings = FALSE)

### Collect data
# new documentation from 27-08-2025
GeneralDocumentation <- read_delim("data input/export-2025-08-27-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_HRO_KSSGCH_MRI_NIUSASSK/GeneralDocumentation.csv", 
           delim = ";", escape_double = FALSE, trim_ws = TRUE)

# IDS of patients
# -> ALS patients
ALS_ID <- GeneralDocumentation %>%
  filter((ALSuncertainty == 1 | ALSFUdiagnosis == 1) & ALSFUdiagnosis %in% c(1,3,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("ALS",68))

# -> CTR 
CTR_ID <- GeneralDocumentation %>%
  filter(PGMC == 2) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("CTR",55))

# -> PGMC
PGMC_ID <- GeneralDocumentation %>%
  filter(PGMC == 1) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("PGMC",64))

# -> Mimic
mimic_ID <- GeneralDocumentation %>%
  filter((ALSuncertainty == 2 | ALSFUdiagnosis == 2) & ALSFUdiagnosis %in% c(2,NA)) %>%
  select(PatientID,ParticipantCode) %>% 
  mutate(type = rep("mimic",23))

# PGMC IDs with mutations
PGMC_mutations_ID <- GeneralDocumentation %>%
  select(PatientID, ParticipantCode,PGMC,contains("MutationType")) %>%
  filter(PGMC==1)
PGMC_mutations_ID_tmp <- sapply(apply(PGMC_mutations_ID[,4:17],1,
                                  function(x){which(x == 1)}) %>% unique(), function(x) x + 3)
PGMC_mutations_ID <- PGMC_mutations_ID[,c(1,2,unlist(PGMC_mutations_ID_tmp))]
PGMC_mutations_ID <- PGMC_mutations_ID %>%
  mutate(type = ifelse(MutationTypeC9orf72 == 1,"C9orf72",
                       ifelse(MutationTypeSOD1 == 1,"SOD1",
                              ifelse(MutationTypeTARDBP == 1,"TARDBP",
                                     ifelse(MutationTypeFUS==1,"FUS",
                                            ifelse(MutationTypeFIG4==1,"FIG4",
                                                   ifelse(MutationTypeUBQLN2==1,"UBQLN2",
                                            ifelse(MutationTypeOther==1,"other",NA)))))))) %>%
  select(PatientID,ParticipantCode,type)

all_participants_IDs = do.call("rbind",
                               list(CTR_ID,
                                    ALS_ID,
                                    PGMC_ID,
                                    mimic_ID,
                                    PGMC_mutations_ID))
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
Sex_age_all_participants = do.call("rbind",list(participants_PGMC,
                                                participants_CTR,
                                                participants_ALS_mimic))

# Questionnaire info with ALSFRS scores
Questionnaire <- read_delim("data input/export-2025-09-04-PREMODIALS-AKDTR_BRNO_CHUFR_HMCIL_KSSGCH_MRI_NIUSASSK (ALSFRS-r)/QuestionnaireG.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)


