### Libraries
library(dplyr)
library(readr)
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
# raw data with smells
BSIT <- read_delim("export-2025-03-31-PREMODIALS-AKDTR_CHUFR_HMCIL_KSSGCH_MRI_NIUSASSK/BSIT.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)