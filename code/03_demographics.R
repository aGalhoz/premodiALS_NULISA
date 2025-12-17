### Demographics

# Date
date_CTR <- Sex_age_all_participants %>%
  select(Pseudonyme,age) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "CTR",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("CTR",29))
date_ALS <- Sex_age_all_participants %>%
  select(Pseudonyme,age) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "ALS",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("ALS",36))
date_mimic <- Sex_age_all_participants %>%
  select(Pseudonyme,age) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "mimic",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("mimic",8))
date_PGMC <- Sex_age_all_participants %>%
  select(Pseudonyme,age) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "PGMC",]$PatientID %>%unique())))  %>%
  mutate(patient_group = rep("PGMC",34))

skim(date_CTR)
skim(date_ALS)
skim(date_mimic)
skim(date_PGMC)

date_all <- do.call("rbind",list(date_CTR,date_ALS,date_mimic,date_PGMC))

# -> Kruskal Wallis test between samples
kruskal_test(data = date_all,age~patient_group)
dunn_test(data = date_all,age~patient_group)

## Sex analysis
Sex_CTR <- Sex_age_all_participants %>%
  select(Pseudonyme,sex) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "CTR",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("CTR",29))
Sex_ALS <- Sex_age_all_participants %>%
  select(Pseudonyme,sex) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "ALS",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("ALS",36))
Sex_mimic <- Sex_age_all_participants %>%
  select(Pseudonyme,sex) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "mimic",]$PatientID %>%unique()))) %>%
  mutate(patient_group = rep("mimic",8))
Sex_PGMC <- Sex_age_all_participants %>%
  select(Pseudonyme,sex) %>%
  filter(Pseudonyme %in% na.omit((samples_ID_type[samples_ID_type$type == "PGMC",]$PatientID %>%unique())))  %>%
  mutate(patient_group = rep("PGMC",34))
  
skim(Sex_ALS)
skim(Sex_CTR)
skim(Sex_mimic)
skim(Sex_PGMC)

Sex_all <- do.call("rbind",list(Sex_ALS,
                                Sex_CTR,
                                Sex_mimic,
                                Sex_PGMC))
Sex_all <- Sex_all %>% filter(!is.na(sex))
table(Sex_all$patient_group,Sex_all$sex)
fisher.test(Sex_all$patient_group,Sex_all$sex)
pairwise_fisher_test(table(Sex_all$patient_group,Sex_all$sex), p.adjust.method = "holm")

#  mutations for ALS
mutations_ALS <- GeneralDocumentation %>%
  select(PatientID,ParticipantCode,contains("MutationType")) %>%
  filter(PatientID %in% date_ALS$Pseudonyme) 

mutations_PGMC <- GeneralDocumentation %>%
  select(PatientID,ParticipantCode,contains("MutationType")) %>%
  filter(PatientID %in% date_PGMC$Pseudonyme) 

# new clinical table for premodiALS (10/09/2025)
clinical_table <- GeneralDocumentation %>%
  select(PatientID,ParticipantCode,ALSuncertainty,ALSFUdiagnosis,Comment,PGMC,
         DatePhenoconversion,contains("Mutation"),contains("LFU")) %>%
  left_join(Sex_age_all_participants %>% rename(PatientID = Pseudonyme)) %>%
  left_join(do.call("rbind", list(ALS_ID, CTR_ID, PGMC_ID,mimic_ID)))  %>% 
  left_join(Questionnaire %>% select(PatientID,Erfassungsdatum,contains("Gs")) %>%
              filter(!if_all(starts_with("Gs"),is.na)))
# -> reorganize mutations and reorganize order
clinical_table_extended <- clinical_table %>%
  pivot_longer(cols = starts_with("MutationType") | starts_with("PreciseMutation"),
               names_to = c(".value", "Gene"),
               names_pattern = "(MutationType|PreciseMutation)(.*)"
  ) %>% 
  distinct() %>%
  group_by(PatientID, ParticipantCode) %>%
  mutate(has_mut = any(MutationType == 1)) %>%
 # ungroup() %>%
  mutate(
    Gene = ifelse(!has_mut, NA, Gene),
    MutationType = ifelse(!has_mut, NA, MutationType),
    PreciseMutation = ifelse(!has_mut, NA, PreciseMutation)
  ) %>%
  select(-has_mut) %>%
  filter(MutationType == 1 | is.na(MutationType)) %>%
  distinct()

clinical_table_extended <- clinical_table_extended %>%
  select(PatientID,ParticipantCode,type,PGMC,ALSuncertainty,ALSFUdiagnosis,Comment,
         age,sex,Gene,PreciseMutation,
         Erfassungsdatum,contains("Gs"),DatePhenoconversion,contains("LFU")) %>% 
  mutate(LFU = ifelse(LFU == "1","yes",
                      ifelse(LFU == 2, "no",NA)),
         LFUReason = ifelse(LFUReason == "1", "Patient not willing to come",
                            ifelse(LFUReason == "2","Patient not able to come",
                                   ifelse(LFUReason == "3","Patient died",
                                          ifelse(LFUReason == "4","Reason not given",NA))))) %>%
  rename(MutationType = Gene) %>%
  mutate(type = case_when(
    ALSuncertainty == 3 & (is.na(ALSFUdiagnosis) | ALSFUdiagnosis == 3) ~ "SYMP",
    TRUE ~ type
  ))

writexl::write_xlsx(clinical_table_extended,"results/clinical_table_extended.xlsx")

# -> make slim version with ALSFRS computed for each row and summarise with format from Laura
clinical_table_slim <- clinical_table_extended %>%
  ungroup() %>%
  mutate(
    across(
      starts_with("Gs"),
      ~ recode(as.numeric(.),
               `1` = 4,
               `2` = 3,
               `3` = 2,
               `4` = 1,
               `5` = 0),
      .names = "{.col}"
    ),
    ALSFRS = rowSums(across(starts_with("Gs")), na.rm = TRUE)
  )

clinical_table_slim <- clinical_table_slim %>%
  select(PatientID,ParticipantCode,type,Comment,age,sex,MutationType,PreciseMutation,
         Erfassungsdatum,ALSFRS,DatePhenoconversion,contains("LFU")) %>%
  group_by(PatientID, ParticipantCode, type, age, sex,MutationType, PreciseMutation) %>% 
 # arrange(Erfassungsdatum, .by_group = TRUE) %>%  # make sure order is consistent
  mutate(obs = row_number()) %>%
  pivot_wider(
    names_from = obs,
    values_from = c(ALSFRS, Erfassungsdatum),
    names_sep = "_"
  ) %>%
  group_by(PatientID, ParticipantCode, type, age, sex) %>% 
  mutate(obs = row_number()) %>%
  pivot_wider(
    names_from = obs,
    values_from = c(MutationType, PreciseMutation),
    names_sep = "_"
  ) %>%
  ungroup() %>%
  select(PatientID,ParticipantCode,type,Comment,age,sex,contains("MutationType"),
         contains("PreciseMutation"), contains("Erfassungsdatum"),
         contains("ALSFRS"),DatePhenoconversion,contains("LFU"))

writexl::write_xlsx(clinical_table_slim,"results/clinical_table.xlsx")

# stats of ALSFRS for V0 of NULISA patients
clinical_table_NULISA <- clinical_table_slim %>%
  filter(PatientID %in% samples_ID_type$PatientID)

writexl::write_xlsx(clinical_table_NULISA,"results/clinical_table_NULISA.xlsx")

ALSFRS_NULISA_V0 <- clinical_table_NULISA %>%
  select(PatientID,type,ALSFRS_1) %>%
  distinct()

# -> Kruskal Wallis test between samples
kruskal_test(data = ALSFRS_NULISA_V0,ALSFRS_1~type)
dunn_test(data = ALSFRS_NULISA_V0,ALSFRS_1~type)

# stats per group
skim(ALSFRS_NULISA_V0 %>% filter(type == "CTR") %>% pull(ALSFRS_1))
skim(ALSFRS_NULISA_V0 %>% filter(type == "PGMC") %>% pull(ALSFRS_1))
skim(ALSFRS_NULISA_V0 %>% filter(type == "ALS") %>% pull(ALSFRS_1))
skim(ALSFRS_NULISA_V0 %>% filter(type == "mimic") %>% pull(ALSFRS_1))

## ECAS info
ECAS_table = clinical_table_slim %>%
  select(PatientID,ParticipantCode,type,age,sex) %>%
  left_join(ECAS %>% select(PatientID,EcasExists,Erfassungsdatum,contains("ECASALS"),
                            EcasTotalBehavior,EcasTotScore)) %>%
  mutate(EcasExists = ifelse(EcasExists == 1, "yes", ifelse(EcasExists == 2,"no",NA)),
         EcasALSSpecAllQs = ifelse(EcasALSSpecAllQs == 1, "yes", ifelse(EcasALSSpecAllQs == 2,"no",NA)),
         EcasALSnonSpecAllQs = ifelse(EcasALSnonSpecAllQs == 1, "yes", ifelse(EcasALSnonSpecAllQs == 2,"no",NA))) %>%
  group_by(PatientID,ParticipantCode,type,age,sex) %>%
  mutate(obs = row_number()) %>%
  pivot_wider(
    names_from = obs,
    values_from = c(contains("Ecas"), Erfassungsdatum),
    names_sep = "_"
  ) %>%
  select(where(~ !all(is.na(.)))) %>%
  ungroup() %>%
  select(PatientID,ParticipantCode,type,age,sex,contains("Ecas"))

writexl::write_xlsx(ECAS_table,"results/ECAS_overview.xlsx")



