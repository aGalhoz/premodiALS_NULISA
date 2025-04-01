source("00_initialization.R")

# map IDs from premodiALS with NULISA ID samples, with the status information
samples_ID_type = sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(do.call("rbind", list(ALS_ID, CTR_ID, PGMC_ID,mimic_ID))) 

# check which samples in NULISA do not have a map from premodiALS
samples_ID_NA = samples_ID_type %>%
  filter(is.na(type)) %>%
  left_join(GeneralDocumentation %>%
              select(ParticipantCode, PGMC, ALSuncertainty, ALSFUdiagnosis,LFU) %>%
              unique())
writexl::write_xlsx(samples_ID_NA,"results/samples_ID_NA.xlsx")

# get protein data with status information
protein_data_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName,SampleMatrixType,Target,UniProtID,ProteinName,NPQ) %>%
  left_join(samples_ID_type %>% select(`Sample ID`,type) %>% rename(SampleName = `Sample ID`))

# amount of unique targets
target_unique <- protein_data_IDs$Target %>% unique() #127 targets

# ANOVA on protein data per protein and type of fluid
protein_anova = protein_data_IDs %>%
  group_by(Target,SampleMatrixType) %>% 
  do(Model = aov(NPQ ~ type, data=.))
protein_anova$pvalue_anova = lapply(protein_anova$Model,function(x) {unlist(summary(x))["Pr(>F)1"]}) %>% unlist()

# Pairwise t-test between groups
protein_data_IDs$type <- factor(protein_data_IDs$type,levels = c("CTR","PGMC","ALS","mimic"))
protein_ttest = protein_data_IDs %>%
  group_by(Target,SampleMatrixType) %>%
  t_test(NPQ ~ type, p.adjust.method = "BH") %>% 
  select(-.y., -statistic, -df)

protein_pvalue = protein_ttest %>%
  left_join(protein_anova %>% select(-Model)) %>%
  arrange(pvalue_anova)
  
writexl::write_xlsx(protein_pvalue,"results/protein_pvalue.xlsx")

# Top 10 for plasma
protein_pvalue_plasma = protein_pvalue %>% 
  filter(SampleMatrixType == "PLASMA")
top10_plasma = protein_pvalue_top10_plasma$Target %>% unique()
top10_plasma = top10_plasma[1:10]

plot_plasma <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "PLASMA" & Target %in% top10_plasma & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) 
# Add statistical test p-values
protein_pvalue_plasma <- protein_pvalue_plasma %>%
  filter(Target %in% top10_plasma) %>% add_xy_position(x = "type")
plot_plasma = plot_plasma + stat_pvalue_manual(protein_pvalue_plasma, label = "p.adj.signif")

pdf("plots/top10_plasma.pdf", width = 10, height = 9)
plot_plasma
dev.off()

# Top 10 for CSF
protein_pvalue_CSF = protein_pvalue %>% 
  filter(SampleMatrixType == "CSF")
top10_CSF = protein_pvalue_CSF$Target %>% unique()
top10_CSF = top10_CSF[1:10]

plot_CSF <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "CSF" & Target %in% top10_CSF & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) 
# Add statistical test p-values
protein_pvalue_CSF <- protein_pvalue_CSF %>%
  filter(Target %in% top10_CSF) %>% add_xy_position(x = "type")
plot_CSF = plot_CSF + stat_pvalue_manual(protein_pvalue_CSF, label = "p.adj.signif")

pdf("plots/top10_CSF.pdf", width = 10, height = 9)
plot_CSF
dev.off()

# Top 10 for Serum
protein_pvalue_SERUM = protein_pvalue %>% 
  filter(SampleMatrixType == "SERUM")
top10_SERUM = protein_pvalue_SERUM$Target %>% unique()
top10_SERUM = top10_SERUM[1:10]

plot_SERUM <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "SERUM" & Target %in% top10_SERUM & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) 
# Add statistical test p-values
protein_pvalue_SERUM <- protein_pvalue_SERUM %>%
  filter(Target %in% top10_SERUM) %>% add_xy_position(x = "type")
plot_SERUM = plot_SERUM + stat_pvalue_manual(protein_pvalue_SERUM, label = "p.adj.signif")

pdf("plots/top10_SERUM.pdf", width = 10, height = 9)
plot_SERUM
dev.off()



for (variable in vector) {
  
}