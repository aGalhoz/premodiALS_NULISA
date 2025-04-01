### Libraries
library(dplyr)
library(readr)
library(readxl)
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

### Directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)

### Collect data
# IDS of patients
ALS_ID <- read_excel("patients_ALS_ID_NULISA.xlsx") %>% select(ParticipantCode) %>% mutate(type = rep("ALS",64))
CTR_ID <- read_excel("patients_CTR_ID_NULISA.xlsx") %>% select(ParticipantCode) %>% mutate(type = rep("CTR",49))
PGMC_ID <- read_excel("patients_PGMC_ID_NULISA.xlsx") %>% select(ParticipantCode) %>% mutate(type = rep("PGMC",58))
mimic_ID <- read_excel("patients_mimic_ID_NULISA.xlsx") %>% select(ParticipantCode) %>% mutate(type = rep("mimic",14))

# information on sample IDs
sample_ID_info <- read_excel("Biospecimen Manifest Form – Banner Biomarker Program_premodiALS.xlsx", 
          sheet = "Sample Information") %>%
  select("Sample ID","Optional Informarion (patient ID)" )

# data from CSF, plasma and serum fluids
protein_data <- read_excel("P004_BSHRI_NULISAseq_CNSDiseasePanel_NPQCounts_2025_03_10.xlsx")

# documentation from BSIT
GeneralDocumentation <- read_delim("~/Documents/HMGU/premodiALS/smell data/export-2025-03-31-PREMODIALS-AKDTR_CHUFR_HMCIL_KSSGCH_MRI_NIUSASSK/GeneralDocumentation.csv", 
                                  delim = ";", escape_double = FALSE, trim_ws = TRUE)

