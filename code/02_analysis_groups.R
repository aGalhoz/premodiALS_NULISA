source("00_initialization.R")
source("01_data_mining.R")

# ANOVA on protein data per protein and type of fluid
protein_anova = protein_data_IDs %>%
  filter(!is.na(NPQ)) %>%
  group_by(Target,SampleMatrixType) %>% 
  do(Model = aov(NPQ ~ type, data=.))
protein_anova$pvalue_anova = lapply(protein_anova$Model,function(x) {unlist(summary(x))["Pr(>F)1"]}) %>% unlist()

# Pairwise t-test between groups
protein_data_IDs$type <- factor(protein_data_IDs$type,levels = c("CTR","PGMC","ALS","mimic"))
protein_ttest = protein_data_IDs %>%
  filter(!is.na(NPQ)) %>%
  group_by(Target,SampleMatrixType) %>%
  t_test(NPQ ~ type, p.adjust.method = "BH") %>% 
  select(-.y., -statistic, -df)

protein_pvalue = protein_ttest %>%
  left_join(protein_anova %>% select(-Model)) %>%
  arrange(pvalue_anova)

# rearrange table for a better summary
protein_pvalue_summary = protein_pvalue %>%
  mutate(pvalue_PGMC_CTR = ifelse(group1 == "CTR" & group2 == "PGMC",p,NA),
         padj_PGMC_CTR = ifelse(group1 == "CTR" & group2 == "PGMC",p.adj,NA),
         padj_signif_PGMC_CTR = ifelse(group1 == "CTR" & group2 == "PGMC",p.adj.signif,NA),
         pvalue_CTR_ALS = ifelse(group1 == "CTR" & group2 == "ALS",p,NA),
         padj_CTR_ALS = ifelse(group1 == "CTR" & group2 == "ALS",p.adj,NA),
         padj_signif_CTR_ALS = ifelse(group1 == "CTR" & group2 == "ALS",p.adj.signif,NA),
         pvalue_CTR_mimic = ifelse(group1 == "CTR" & group2 == "mimic",p,NA),
         padj_CTR_mimic = ifelse(group1 == "CTR" & group2 == "mimic",p.adj,NA),
         padj_signif_CTR_mimic = ifelse(group1 == "CTR" & group2 == "mimic",p.adj.signif,NA),
         pvalue_PGMC_ALS = ifelse(group1 == "PGMC" & group2 == "ALS",p,NA),
         padj_PGMC_ALS = ifelse(group1 == "PGMC" & group2 == "ALS",p.adj,NA),
         padj_signif_PGMC_ALS = ifelse(group1 == "PGMC" & group2 == "ALS",p.adj.signif,NA),
         pvalue_PGMC_mimic = ifelse(group1 == "PGMC" & group2 == "mimic",p,NA),
         padj_PGMC_mimic = ifelse(group1 == "PGMC" & group2 == "mimic",p.adj,NA),
         padj_signif_PGMC_mimic = ifelse(group1 == "PGMC" & group2 == "mimic",p.adj.signif,NA)) %>%
  select(-c(group1,group2,n1,n2,p,p.adj,p.adj.signif)) %>%
  distinct()
protein_pvalue_summary = do.call("cbind",list(protein_pvalue_summary[,1:6] %>% distinct() %>%
                                                na.omit(),
                                              protein_pvalue_summary[,c(1:3,7:9)] %>% distinct() %>%
                                                na.omit(),
                                              protein_pvalue_summary[,c(1:3,10:12)] %>% distinct() %>%
                                                na.omit(),
                                              protein_pvalue_summary[,c(1:3,13:15)] %>% distinct() %>%
                                                na.omit(),
                                              protein_pvalue_summary[,c(1:3,16:18)] %>% distinct() %>%
                                                na.omit()))
protein_pvalue_summary = protein_pvalue_summary[,!duplicated(colnames(protein_pvalue_summary))]
protein_pvalue_summary_plasma = protein_pvalue_summary %>%
  filter(SampleMatrixType == "PLASMA")
protein_pvalue_summary_serum = protein_pvalue_summary %>%
  filter(SampleMatrixType == "SERUM")
protein_pvalue_summary_CSF = protein_pvalue_summary %>%
  filter(SampleMatrixType == "CSF")

writexl::write_xlsx(protein_pvalue,"results/protein_pvalue.xlsx")
writexl::write_xlsx(protein_pvalue_summary_plasma,"results/protein_pvalue_summary_plasma.xlsx")
writexl::write_xlsx(protein_pvalue_summary_serum,"results/protein_pvalue_summary_serum.xlsx")
writexl::write_xlsx(protein_pvalue_summary_CSF,"results/protein_pvalue_summary_CSF.xlsx")

# Top 10 for plasma
protein_pvalue_plasma = protein_pvalue %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  arrange(p.adj)
top10_plasma = protein_pvalue_plasma$Target %>% unique()
top10_plasma = top10_plasma[1:10]

plot_plasma <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "PLASMA" & Target %in% top10_plasma & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Add statistical test p-values
protein_pvalue_plasma <- protein_pvalue_plasma %>%
  filter(Target %in% top10_plasma) %>% add_xy_position(x = "type")
plot_plasma = plot_plasma + stat_pvalue_manual(protein_pvalue_plasma, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                 'ALS' = '#B2936F',
                                 'PGMC' = '#ad5291',
                                 'mimic' = '#62cda9',
                                 'other' = '#ad5291',
                                 'C9orf72' = '#55aa82',
                                 'SOD1' = '#4661b9',
                                 'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_plasma/top10_plasma.pdf", width = 10, height = 9)
plot_plasma
dev.off()

# Proteins of interest in plasma
protein_pvalue_plasma = protein_pvalue %>% 
  filter(SampleMatrixType == "PLASMA")
proteins_interest_plasma = c("NEFL","NEFH","GFAP","MAPT","pTau-181","pTau-231",
                             "pTau-217","PGF","BDNF","CRP","PARK7","APOE")

plot_plasma <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "PLASMA" & Target %in% proteins_interest_plasma & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_plasma <- protein_pvalue_plasma %>%
  filter(Target %in% proteins_interest_plasma) %>% add_xy_position(x = "type")
plot_plasma = plot_plasma + stat_pvalue_manual(protein_pvalue_plasma, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_plasma/proteins_interest_plasma.pdf", width = 10, height = 9)
plot_plasma
dev.off()

# All proteins from plasma individually
protein_pvalue_plasma = protein_pvalue %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  arrange(p.adj)
top10_plasma = protein_pvalue_plasma$Target %>% unique()
top10_plasma = top10_plasma[1:10]
proteins_plasma = c(proteins_interest_plasma,top10_plasma) %>% unique()
plots_plasma <- list()
for (i in 1:length(proteins_plasma)) {
  protein_pvalue_plasma_i <- protein_pvalue_plasma %>%
    filter(Target == proteins_plasma[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_IDs %>% filter(SampleMatrixType == "PLASMA" & Target %in% proteins_plasma[i] & !is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_plasma_i, label = "p.adj.signif",size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_plasma[[i]]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_plasma/boxplot_",proteins_plasma[i],".pdf"))
  print(plot_i)
  dev.off()
}

# Top 10 for CSF
protein_pvalue_CSF = protein_pvalue %>% 
  filter(SampleMatrixType == "CSF") %>% 
  arrange(p.adj)
top10_CSF = protein_pvalue_CSF$Target %>% unique()
top10_CSF = top10_CSF[1:10]

plot_CSF <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "CSF" & Target %in% top10_CSF & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_CSF <- protein_pvalue_CSF %>%
  filter(Target %in% top10_CSF) %>% add_xy_position(x = "type")
plot_CSF = plot_CSF + stat_pvalue_manual(protein_pvalue_CSF, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_CSF/top10_CSF.pdf", width = 10, height = 9)
plot_CSF
dev.off()

# Proteins of interest in CSF
protein_pvalue_CSF = protein_pvalue %>% 
  filter(SampleMatrixType == "CSF")
proteins_interest_CSF = c("NEFL","NEFH","IL18","TNF","IL6","HTT","CHI3L1",
                          "CCL3","UCHL1","ANXA5","CCL2","CHIT1","TARDBP","SOD1")

plot_CSF <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "CSF" & Target %in% proteins_interest_CSF & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_CSF <- protein_pvalue_CSF %>%
  filter(Target %in% proteins_interest_CSF) %>% add_xy_position(x = "type")
plot_CSF = plot_CSF + stat_pvalue_manual(protein_pvalue_CSF, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_CSF/proteins_interest_CSF.pdf", width = 10, height = 9)
plot_CSF
dev.off()

# All proteins from CSF individually
protein_pvalue_CSF = protein_pvalue %>% 
  filter(SampleMatrixType == "CSF") %>% 
  arrange(p.adj)
top10_CSF = protein_pvalue_CSF$Target %>% unique()
top10_CSF = top10_CSF[1:10]
proteins_CSF = c(proteins_interest_CSF,top10_CSF) %>% unique()
for (i in 1:length(proteins_CSF)) {
  protein_pvalue_CSF_i <- protein_pvalue_CSF %>%
    filter(Target == proteins_CSF[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_IDs %>% filter(SampleMatrixType == "CSF" & Target %in% proteins_CSF[i] & !is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_CSF_i, label = "p.adj.signif",size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_CSF[[i]]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_CSF/boxplot_",proteins_CSF[i],".pdf"))
  print(plot_i)
  dev.off()
}


# Top 10 for Serum
protein_pvalue_SERUM = protein_pvalue %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  arrange(p.adj)
top10_SERUM = protein_pvalue_SERUM$Target %>% unique()
top10_SERUM = top10_SERUM[1:10]

plot_SERUM <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "SERUM" & Target %in% top10_SERUM & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_SERUM <- protein_pvalue_SERUM %>%
  filter(Target %in% top10_SERUM) %>% add_xy_position(x = "type")
plot_SERUM = plot_SERUM + stat_pvalue_manual(protein_pvalue_SERUM, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_SERUM/top10_SERUM.pdf", width = 10, height = 9)
plot_SERUM
dev.off()

# Proteins of interest in SERUM
protein_pvalue_SERUM = protein_pvalue %>% 
  filter(SampleMatrixType == "SERUM")
proteins_interest_SERUM = c("NEFL","NEFH","pTau-181","pTau-217","pTau-231","FABP3",
                            "GFAP","MAPT","BDNF","GOT1","TAFA5","PARK7","ENO2","CRP",
                            "IL7","KLK6","IL9","FCN2","KDR","VEGFD")

plot_SERUM <- ggboxplot(
  protein_data_IDs %>% filter(SampleMatrixType == "SERUM" & Target %in% proteins_interest_SERUM & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_SERUM <- protein_pvalue_SERUM %>%
  filter(Target %in% proteins_interest_SERUM) %>% add_xy_position(x = "type")
plot_SERUM = plot_SERUM + stat_pvalue_manual(protein_pvalue_SERUM, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_SERUM/proteins_interest_SERUM.pdf", width = 10, height = 9)
plot_SERUM
dev.off()

# All proteins from SERUM individually
protein_pvalue_SERUM = protein_pvalue %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  arrange(p.adj)
top10_SERUM = protein_pvalue_SERUM$Target %>% unique()
top10_SERUM = top10_SERUM[1:10]
proteins_SERUM = c(proteins_interest_SERUM,top10_SERUM) %>% unique()
for (i in 1:length(proteins_SERUM)) {
  protein_pvalue_SERUM_i <- protein_pvalue_SERUM %>%
    filter(Target == proteins_SERUM[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_IDs %>% filter(SampleMatrixType == "SERUM" & Target %in% proteins_SERUM[i] & !is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_SERUM_i, label = "p.adj.signif", size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_SERUM[i]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_SERUM/boxplot_",proteins_SERUM[i],".pdf"))
  print(plot_i)
  dev.off()
}

##### PGMC mutations and controls

# ANOVA on protein data per protein and type of fluid
protein_anova_PGMC_CTR = protein_data_PGMC_CTR_IDs %>%
  filter(!is.na(NPQ)) %>%
  group_by(Target,SampleMatrixType) %>% 
  do(Model = aov(NPQ ~ type, data=.))
protein_anova_PGMC_CTR$pvalue_anova = lapply(protein_anova_PGMC_CTR$Model,function(x) {unlist(summary(x))["Pr(>F)1"]}) %>% unlist()

# Pairwise t-test between groups
protein_data_PGMC_CTR_IDs$type <- factor(protein_data_PGMC_CTR_IDs$type,levels = c("CTR","C9orf72","SOD1","TARDBP","FUS","other"))
protein_ttest_PGMC_CTR = protein_data_PGMC_CTR_IDs %>%
  filter(type %in% c("CTR","C9orf72","SOD1","TARDBP")) %>%
  filter(!is.na(NPQ)) %>%
  group_by(Target,SampleMatrixType) %>%
  t_test(NPQ ~ type, p.adjust.method = "BH") %>% 
  select(-.y., -statistic, -df)

protein_pvalue_PGMC_CTR = protein_ttest_PGMC_CTR %>%
  left_join(protein_anova_PGMC_CTR %>% select(-Model)) %>%
  arrange(pvalue_anova)

# rearrange table for a better summary
protein_pvalue_PGMC_CTR_summary = protein_pvalue_PGMC_CTR %>%
  mutate(pvalue_CTR_C9orf72 = ifelse(group1 == "CTR" & group2 == "C9orf72",p,NA),
         padj_CTR_C9orf72 = ifelse(group1 == "CTR" & group2 == "C9orf72",p.adj,NA),
         padj_signif_CTR_C9orf72 = ifelse(group1 == "CTR" & group2 == "C9orf72",p.adj.signif,NA),
         pvalue_CTR_SOD1 = ifelse(group1 == "CTR" & group2 == "SOD1",p,NA),
         padj_CTR_SOD1 = ifelse(group1 == "CTR" & group2 == "SOD1",p.adj,NA),
         padj_signif_CTR_SOD1 = ifelse(group1 == "CTR" & group2 == "SOD1",p.adj.signif,NA),
         pvalue_CTR_TARDBP = ifelse(group1 == "CTR" & group2 == "TARDBP",p,NA),
         padj_CTR_TARDBP = ifelse(group1 == "CTR" & group2 == "TARDBP",p.adj,NA),
         padj_signif_CTR_TARDBP = ifelse(group1 == "CTR" & group2 == "TARDBP",p.adj.signif,NA),
         pvalue_C9orf72_SOD1 = ifelse(group1 == "C9orf72" & group2 == "SOD1",p,NA),
         padj_C9orf72_SOD1 = ifelse(group1 == "C9orf72" & group2 == "SOD1",p.adj,NA),
         padj_signif_C9orf72_SOD1 = ifelse(group1 == "C9orf72" & group2 == "SOD1",p.adj.signif,NA),
         pvalue_C9orf72_TARDBP = ifelse(group1 == "C9orf72" & group2 == "TARDBP",p,NA),
         padj_C9orf72_TARDBP = ifelse(group1 == "C9orf72" & group2 == "TARDBP",p.adj,NA),
         padj_signif_C9orf72_TARDBP = ifelse(group1 == "C9orf72" & group2 == "TARDBP",p.adj.signif,NA),
         pvalue_SOD1_TARDBP = ifelse(group1 == "SOD1" & group2 == "TARDBP",p,NA),
         padj_SOD1_TARDBP = ifelse(group1 == "SOD1" & group2 == "TARDBP",p.adj,NA),
         padj_signif_SOD1_TARDBP = ifelse(group1 == "SOD1" & group2 == "TARDBP",p.adj.signif,NA)) %>%
  select(-c(group1,group2,n1,n2,p,p.adj,p.adj.signif)) %>%
  distinct()
protein_pvalue_PGMC_CTR_summary = do.call("cbind",list(protein_pvalue_PGMC_CTR_summary[,1:6] %>% distinct() %>%
                                                na.omit(),
                                                protein_pvalue_PGMC_CTR_summary[,c(1:3,7:9)] %>% distinct() %>%
                                                na.omit(),
                                                protein_pvalue_PGMC_CTR_summary[,c(1:3,10:12)] %>% distinct() %>%
                                                na.omit(),
                                                protein_pvalue_PGMC_CTR_summary[,c(1:3,13:15)] %>% distinct() %>%
                                                na.omit(),
                                                protein_pvalue_PGMC_CTR_summary[,c(1:3,16:18)] %>% distinct() %>%
                                                  na.omit(),
                                                protein_pvalue_PGMC_CTR_summary[,c(1:3,19:21)] %>% distinct() %>%
                                                  na.omit()))
protein_pvalue_PGMC_CTR_summary = protein_pvalue_PGMC_CTR_summary[,!duplicated(colnames(protein_pvalue_PGMC_CTR_summary))]
protein_pvalue_PGMC_CTR_summary_plasma = protein_pvalue_PGMC_CTR_summary %>%
  filter(SampleMatrixType == "PLASMA")
protein_pvalue_PGMC_CTR_summary_serum = protein_pvalue_PGMC_CTR_summary %>%
  filter(SampleMatrixType == "SERUM")
protein_pvalue_PGMC_CTR_summary_CSF = protein_pvalue_PGMC_CTR_summary %>%
  filter(SampleMatrixType == "CSF")

writexl::write_xlsx(protein_pvalue_PGMC_CTR,"results/protein_pvalue_PGMC_CTR.xlsx")
writexl::write_xlsx(protein_pvalue_PGMC_CTR_summary_plasma,"results/protein_pvalue_PGMC_CTR_summary_plasma.xlsx")
writexl::write_xlsx(protein_pvalue_PGMC_CTR_summary_serum,"results/protein_pvalue_PGMC_CTR_summary_serum.xlsx")
writexl::write_xlsx(protein_pvalue_PGMC_CTR_summary_CSF,"results/protein_pvalue_PGMC_CTR_summary_CSF.xlsx")

# Top 10 for plasma
protein_pvalue_plasma_PGMC_CTR_plasma = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  arrange(p.adj)
top10_plasma = protein_pvalue_plasma_PGMC_CTR_plasma$Target %>% unique()
top10_plasma = top10_plasma[1:10]

plot_plasma <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "PLASMA" & type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% top10_plasma & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_plasma_PGMC_CTR_plasma <- protein_pvalue_plasma_PGMC_CTR_plasma %>%
  filter(Target %in% top10_plasma) %>% add_xy_position(x = "type")
plot_plasma = plot_plasma + stat_pvalue_manual(protein_pvalue_plasma_PGMC_CTR_plasma, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_plasma/top10_plasma_PGMC_CTR.pdf", width = 12, height = 9)
plot_plasma
dev.off()

# Proteins of interest in PLASMA
protein_pvalue_PGMC_CTR_PLASMA = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "PLASMA")

plot_PLASMA <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "PLASMA"& type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% proteins_interest_plasma & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_PGMC_CTR_PLASMA <- protein_pvalue_PGMC_CTR_PLASMA %>%
  filter(Target %in% proteins_interest_plasma) %>% add_xy_position(x = "type")
plot_PLASMA = plot_PLASMA + stat_pvalue_manual(protein_pvalue_PGMC_CTR_PLASMA, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_plasma/proteins_interest_PLASMA_PGMC_CTR.pdf", width = 12, height = 9)
plot_PLASMA
dev.off()

# All proteins from PLASMA individually
protein_pvalue_PGMC_CTR_PLASMA = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  arrange(p.adj)
top10_plasma = protein_pvalue_plasma_PGMC_CTR_plasma$Target %>% unique()
top10_plasma = top10_plasma[1:10]
proteins_PLASMA = c(proteins_interest_plasma,top10_plasma) %>% unique()
for (i in 1:length(proteins_PLASMA)) {
  protein_pvalue_PLASMA_i <- protein_pvalue_PGMC_CTR_PLASMA %>%
    filter(Target == proteins_PLASMA[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_PGMC_CTR_IDs %>% 
      filter(SampleMatrixType == "PLASMA" & Target %in% proteins_PLASMA[i] 
             & type %in% c("CTR","C9orf72","SOD1","TARDBP") &!is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_PLASMA_i, label = "p.adj.signif",size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_PLASMA[[i]]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_plasma/boxplot_PGMC_CTR_",proteins_PLASMA[i],".pdf"))
  print(plot_i)
  dev.off()
}

# Top 10 for CSF
protein_pvalue_PGMC_CTR_CSF = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "CSF") %>% 
  arrange(p.adj)
top10_CSF = protein_pvalue_PGMC_CTR_CSF$Target %>% unique()
top10_CSF = top10_CSF[1:10]

plot_CSF <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "CSF" & type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% top10_CSF & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_PGMC_CTR_CSF <- protein_pvalue_PGMC_CTR_CSF %>%
  filter(Target %in% top10_CSF) %>% add_xy_position(x = "type")
plot_CSF = plot_CSF + stat_pvalue_manual(protein_pvalue_PGMC_CTR_CSF, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_CSF/top10_CSF_PGMC_CTR.pdf", width = 12, height = 9)
plot_CSF
dev.off()

# Proteins of interest in CSF
protein_pvalue_PGMC_CTR_CSF = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "CSF")

plot_CSF <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "CSF"& type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% proteins_interest_CSF & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_PGMC_CTR_CSF <- protein_pvalue_PGMC_CTR_CSF %>%
  filter(Target %in% proteins_interest_CSF) %>% add_xy_position(x = "type")
plot_CSF = plot_CSF + stat_pvalue_manual(protein_pvalue_PGMC_CTR_CSF, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_CSF/proteins_interest_CSF_PGMC_CTR.pdf", width = 12, height = 9)
plot_CSF
dev.off()

# All proteins from CSF individually
protein_pvalue_PGMC_CTR_CSF = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "CSF") %>% 
  arrange(p.adj)
top10_CSF = protein_pvalue_PGMC_CTR_CSF$Target %>% unique()
top10_CSF = top10_CSF[1:10]
proteins_CSF = c(proteins_interest_CSF,top10_CSF) %>% unique()
for (i in 1:length(proteins_CSF)) {
  protein_pvalue_CSF_i <- protein_pvalue_PGMC_CTR_CSF %>%
    filter(Target == proteins_CSF[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_PGMC_CTR_IDs %>% 
      filter(SampleMatrixType == "CSF" & Target %in% proteins_CSF[i] 
             & type %in% c("CTR","C9orf72","SOD1","TARDBP") &!is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_CSF_i, label = "p.adj.signif",size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_CSF[[i]]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_CSF/boxplot_PGMC_CTR_",proteins_CSF[i],".pdf"))
  print(plot_i)
  dev.off()
}

# Top 10 for Serum
protein_pvalue_PGMC_CTR_SERUM = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  arrange(p.adj)
top10_SERUM = protein_pvalue_PGMC_CTR_SERUM$Target %>% unique()
top10_SERUM = top10_SERUM[1:10]

plot_SERUM <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "SERUM" & type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% top10_SERUM & !is.na(type)), x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_PGMC_CTR_SERUM <- protein_pvalue_PGMC_CTR_SERUM %>%
  filter(Target %in% top10_SERUM) %>% add_xy_position(x = "type")
plot_SERUM = plot_SERUM + stat_pvalue_manual(protein_pvalue_PGMC_CTR_SERUM, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_SERUM/top10_SERUM_PGMC_CTR.pdf", width = 12, height = 9)
plot_SERUM
dev.off()

# Proteins of interest in SERUM
protein_pvalue_PGMC_CTR_SERUM = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "SERUM")

plot_SERUM <- ggboxplot(
  protein_data_PGMC_CTR_IDs %>% 
    filter(SampleMatrixType == "SERUM"& type %in% c("CTR","C9orf72","SOD1","TARDBP") &
             Target %in% proteins_interest_SERUM & !is.na(type)), 
  x = "type", y = "NPQ",
  fill = "type", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~Target) +
  theme(
    axis.title.x = element_blank(),
    text = element_text(size = 15),               # Base font size for everything
    axis.title = element_text(size = 18),         # Axis titles
    axis.text = element_text(size = 14),          # Axis tick labels
    strip.text = element_text(size = 16, face = "bold"),  # Facet labels
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# Add statistical test p-values
protein_pvalue_PGMC_CTR_SERUM <- protein_pvalue_PGMC_CTR_SERUM %>%
  filter(Target %in% proteins_interest_SERUM) %>% add_xy_position(x = "type")
plot_SERUM = plot_SERUM + stat_pvalue_manual(protein_pvalue_PGMC_CTR_SERUM, label = "p.adj.signif") +
  scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                'ALS' = '#B2936F',
                                'PGMC' = '#ad5291',
                                'mimic' = '#62cda9',
                                'other' = '#ad5291',
                                'C9orf72' = '#55aa82',
                                'SOD1' = '#4661b9',
                                'TARDBP' = '#B99E46')) 

pdf("plots/boxplots_SERUM/proteins_interest_SERUM_PGMC_CTR.pdf", width = 12, height = 9)
plot_SERUM
dev.off()

# All proteins from SERUM individually
protein_pvalue_PGMC_CTR_SERUM = protein_pvalue_PGMC_CTR %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  arrange(p.adj)
top10_SERUM = protein_pvalue_PGMC_CTR_SERUM$Target %>% unique()
top10_SERUM = top10_SERUM[1:10]
proteins_SERUM = c(proteins_interest_SERUM,top10_SERUM) %>% unique()
for (i in 1:length(proteins_SERUM)) {
  protein_pvalue_SERUM_i <- protein_pvalue_PGMC_CTR_SERUM %>%
    filter(Target == proteins_SERUM[i]) %>% add_xy_position(x = "type")
  plot_i = ggboxplot(
    protein_data_PGMC_CTR_IDs %>% 
      filter(SampleMatrixType == "SERUM" & Target %in% proteins_SERUM[i] 
             & type %in% c("CTR","C9orf72","SOD1","TARDBP") &!is.na(type)), 
    x = "type", y = "NPQ",
    fill = "type", palette = "npg", legend = "none",
    ggtheme = theme_pubr(border = TRUE)) +
    stat_pvalue_manual(protein_pvalue_SERUM_i, label = "p.adj.signif",size = 6) +
    scale_fill_manual(values  = c('CTR' = '#6F8EB2',  
                                  'ALS' = '#B2936F',
                                  'PGMC' = '#ad5291',
                                  'mimic' = '#62cda9',
                                  'other' = '#ad5291',
                                  'C9orf72' = '#55aa82',
                                  'SOD1' = '#4661b9',
                                  'TARDBP' = '#B99E46')) + 
    labs(title = proteins_SERUM[[i]]) +
    theme(
      axis.title.x = element_blank(),
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  pdf(paste0("plots/boxplots_SERUM/boxplot_PGMC_CTR_",proteins_SERUM[i],".pdf"))
  print(plot_i)
  dev.off()
}

## PCA with all diseased groups wiith shape of PGMC by mutation
protein_data_PCA = protein_data_IDs %>%
  left_join(protein_data_PGMC_CTR_IDs %>% rename(subtype = type)) %>%
  mutate(subtype = ifelse(subtype == "CTR" | is.na(subtype),"No subtype",
                          ifelse(subtype == "C9orf72","C9orf72",
                                 ifelse(subtype == "SOD1","SOD1",
                                        ifelse(subtype == "TARBDP","TARDBP",
                                               ifelse(subtype == "FUS","FUS",
                                                      "other"))))))

protein_data_clean <- protein_data_PCA %>%
  filter(!is.na(type))

# Step 2: Function to run PCA per SampleMatrixType
run_pca <- function(df, matrix_type) {
  df_matrix <- df %>%
    filter(SampleMatrixType == matrix_type) %>%
    select(SampleName, Target, NPQ, type, subtype)
  
  # Reshape wide
  df_wide <- df_matrix %>%
    pivot_wider(
      names_from = Target,
      values_from = NPQ,
      values_fill = 0
    )
  
  meta <- df_wide %>% select(SampleName, type, subtype)
  X <- df_wide %>% select(-SampleName, -type, -subtype)
  
  X <- X %>% select(where(~ {
    s <- sd(.x, na.rm = TRUE)
    !is.na(s) && s > 0
  }))
  
  # PCA
  pca <- prcomp(X, scale. = TRUE)
  
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta)
  
  list(scores = scores, pca = pca)
}
# Step 3: Define color palette
my_colors <- c(
  'CTR'    = '#6F8EB2',
  'ALS'    = '#B2936F',
  'PGMC'   = '#ad5291',
  'mimic'  = '#62cda9',
  'other'  = '#ad5291',
  'C9orf72'= '#55aa82',
  'SOD1'   = '#4661b9',
  'TARDBP' = '#B99E46'
)

# Step 4: Generate PCA for each matrix type
matrices <- unique(protein_data_clean$SampleMatrixType)

pca_results <- lapply(matrices, function(m) run_pca(protein_data_clean, m))
names(pca_results) <- matrices

# Step 5: Plot function
plot_pca <- function(pca_res, title) {
  pca <- pca_res$pca
  scores <- pca_res$scores %>%
    mutate(
      subtype = factor(
        subtype,
        levels = c("No subtype", "C9orf72", "SOD1", "FUS", "other")
      )
    )
  
  ggplot(scores, aes(x = PC1, y = PC2, color = type, shape = subtype)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_color_manual(values = my_colors) +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("PCA -", title),
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
    ) +
    theme(
      text = element_text(size = 16),               # Base font size for everything
      axis.title = element_text(size = 18),         # Axis titles
      axis.text = element_text(size = 16),          # Axis tick labels
      plot.title = element_text(size = 18, hjust = 0.5,face = "bold"),
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 16)
    )
}

# Step 6: Create plots for each SampleMatrixType
plots <- lapply(names(pca_results), function(m) plot_pca(pca_results[[m]], m))

# Display them one by one
pdf("plots/PCA_SERUM.pdf",width = 8,height = 6.5)
plots[[1]]
dev.off()
pdf("plots/PCA_PLASMA.pdf",width = 8,height = 6.5)
plots[[2]]
dev.off()
pdf("plots/PCA_CSF.pdf",width = 8,height = 6.5)
plots[[3]]
dev.off()


