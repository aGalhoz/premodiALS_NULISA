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

