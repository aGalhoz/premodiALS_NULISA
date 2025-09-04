# Diferential expression analyses
# -> Plasma (subgroups)
DE_subgroups_plasma = protein_data_IDs %>%
  filter(SampleMatrixType == "PLASMA") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_subgroups_plasma = DE_subgroups_plasma %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_ALS = ifelse(type == "ALS",mean_NPQ,NA),
         mean_NPQ_mimic = ifelse(type == "mimic",mean_NPQ,NA),
         mean_NPQ_PGMC = ifelse(type == "PGMC",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_subgroups_plasma = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                             list(DE_subgroups_plasma[,1:3] %>% distinct() %>% na.omit(), 
                                  DE_subgroups_plasma[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_plasma[,c(1:2,5)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_plasma[,c(1:2,6)] %>% distinct() %>% na.omit()))

DE_subgroups_plasma = DE_subgroups_plasma %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_ALS_CTR = mean_NPQ_ALS- mean_NPQ_CTR,
         log2FC_ALS_PGMC = mean_NPQ_ALS-mean_NPQ_PGMC,
         log2FC_mimic_CTR = mean_NPQ_mimic - mean_NPQ_CTR,
         log2FC_CTR_PGMC = mean_NPQ_CTR - mean_NPQ_PGMC) %>%
  left_join(protein_pvalue_summary_plasma)

# -> Plasma (PGMC subgroups)
DE_PGMC_subgroups_plasma = protein_data_CTR_PGMC_IDs %>%
  filter(SampleMatrixType == "PLASMA") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_PGMC_subgroups_plasma = DE_PGMC_subgroups_plasma %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_C9orf72 = ifelse(type == "C9orf72",mean_NPQ,NA),
         mean_NPQ_SOD1 = ifelse(type == "SOD1",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_PGMC_subgroups_plasma = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                             list(DE_PGMC_subgroups_plasma[,1:3] %>% distinct() %>% na.omit(), 
                                  DE_PGMC_subgroups_plasma[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                  DE_PGMC_subgroups_plasma[,c(1:2,5)] %>% distinct() %>% na.omit()))

DE_PGMC_subgroups_plasma = DE_PGMC_subgroups_plasma %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_C9orf72_CTR = mean_NPQ_C9orf72 - mean_NPQ_CTR,
         log2FC_SOD1_CTR = mean_NPQ_SOD1 - mean_NPQ_CTR,
         log2FC_C9orf72_SOD1 = mean_NPQ_C9orf72 - mean_NPQ_SOD1) %>%
  left_join(protein_pvalue_CTR_PGMC_summary_plasma)

# -> CSF
DE_subgroups_CSF = protein_data_IDs %>%
  filter(SampleMatrixType == "CSF") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_subgroups_CSF = DE_subgroups_CSF %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_ALS = ifelse(type == "ALS",mean_NPQ,NA),
         mean_NPQ_mimic = ifelse(type == "mimic",mean_NPQ,NA),
         mean_NPQ_PGMC = ifelse(type == "PGMC",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_subgroups_CSF = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                             list(DE_subgroups_CSF[,1:3] %>% distinct() %>% na.omit(), 
                                  DE_subgroups_CSF[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_CSF[,c(1:2,5)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_CSF[,c(1:2,6)] %>% distinct() %>% na.omit()))

DE_subgroups_CSF = DE_subgroups_CSF %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_ALS_CTR = mean_NPQ_ALS - mean_NPQ_CTR,
         log2FC_ALS_PGMC = mean_NPQ_ALS - mean_NPQ_PGMC,
         log2FC_mimic_CTR = mean_NPQ_mimic - mean_NPQ_CTR,
         log2FC_CTR_PGMC = mean_NPQ_CTR - mean_NPQ_PGMC) %>%
  left_join(protein_pvalue_summary_CSF)

# -> CSF (PGMC subgroups)
DE_PGMC_subgroups_CSF = protein_data_CTR_PGMC_IDs %>%
  filter(SampleMatrixType == "CSF") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_PGMC_subgroups_CSF = DE_PGMC_subgroups_CSF %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_C9orf72 = ifelse(type == "C9orf72",mean_NPQ,NA),
         mean_NPQ_SOD1 = ifelse(type == "SOD1",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_PGMC_subgroups_CSF = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                                  list(DE_PGMC_subgroups_CSF[,1:3] %>% distinct() %>% na.omit(), 
                                       DE_PGMC_subgroups_CSF[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                       DE_PGMC_subgroups_CSF[,c(1:2,5)] %>% distinct() %>% na.omit()))

DE_PGMC_subgroups_CSF = DE_PGMC_subgroups_CSF %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_C9orf72_CTR = mean_NPQ_C9orf72 - mean_NPQ_CTR,
         log2FC_SOD1_CTR = mean_NPQ_SOD1 - mean_NPQ_CTR,
         log2FC_C9orf72_SOD1 = mean_NPQ_C9orf72 - mean_NPQ_SOD1) %>%
  left_join(protein_pvalue_CTR_PGMC_summary_CSF)

# -> SERUM
DE_subgroups_SERUM = protein_data_IDs %>%
  filter(SampleMatrixType == "SERUM") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_subgroups_SERUM = DE_subgroups_SERUM %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_ALS = ifelse(type == "ALS",mean_NPQ,NA),
         mean_NPQ_mimic = ifelse(type == "mimic",mean_NPQ,NA),
         mean_NPQ_PGMC = ifelse(type == "PGMC",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_subgroups_SERUM = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                             list(DE_subgroups_SERUM[,1:3] %>% distinct() %>% na.omit(), 
                                  DE_subgroups_SERUM[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_SERUM[,c(1:2,5)] %>% distinct() %>% na.omit(),
                                  DE_subgroups_SERUM[,c(1:2,6)] %>% distinct() %>% na.omit()))

DE_subgroups_SERUM = DE_subgroups_SERUM %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_ALS_CTR = mean_NPQ_ALS - mean_NPQ_CTR,
         log2FC_ALS_PGMC = mean_NPQ_ALS - mean_NPQ_PGMC,
         log2FC_mimic_CTR = mean_NPQ_mimic - mean_NPQ_CTR,
         log2FC_CTR_PGMC = mean_NPQ_CTR - mean_NPQ_PGMC) %>%
  left_join(protein_pvalue_summary_serum)

# -> SERUM (PGMC subgroups)
DE_PGMC_subgroups_SERUM = protein_data_CTR_PGMC_IDs %>%
  filter(SampleMatrixType == "SERUM") %>%
  group_by(type,Target,UniProtID) %>%
  filter(!is.na(NPQ)) %>%
  summarise(mean_NPQ = mean(NPQ))

DE_PGMC_subgroups_SERUM = DE_PGMC_subgroups_SERUM %>%
  group_by(Target,UniProtID) %>%
  mutate(mean_NPQ_CTR = ifelse(type == "CTR",mean_NPQ,NA),
         mean_NPQ_C9orf72 = ifelse(type == "C9orf72",mean_NPQ,NA),
         mean_NPQ_SOD1 = ifelse(type == "SOD1",mean_NPQ,NA)) %>%
  select(-c(mean_NPQ,type))

DE_PGMC_subgroups_SERUM = Reduce(function(x, y) merge(x, y, by = c("Target","UniProtID")),
                                  list(DE_PGMC_subgroups_SERUM[,1:3] %>% distinct() %>% na.omit(), 
                                       DE_PGMC_subgroups_SERUM[,c(1:2,4)] %>% distinct() %>% na.omit(),
                                       DE_PGMC_subgroups_SERUM[,c(1:2,5)] %>% distinct() %>% na.omit()))

DE_PGMC_subgroups_SERUM = DE_PGMC_subgroups_SERUM %>%
  group_by(Target,UniProtID) %>%
  mutate(log2FC_C9orf72_CTR = mean_NPQ_C9orf72 - mean_NPQ_CTR,
         log2FC_SOD1_CTR = mean_NPQ_SOD1 - mean_NPQ_CTR,
         log2FC_C9orf72_SOD1 = mean_NPQ_C9orf72 - mean_NPQ_SOD1) %>%
  left_join(protein_pvalue_CTR_PGMC_summary_serum)

# subgroups
# ALS vs CTR (plasma)
DE_subgroups_plasma_ALS_CTR <- DE_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_CTR_ALS,padj_CTR_ALS) %>%
  mutate(log10_pval = -log10(pvalue_CTR_ALS),
         log10_padj = -log10(padj_CTR_ALS),
         DE = ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR > 0,"ALS",
                     ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR < 0,"CTR",
                            "ns")))

volcano_ALS_CTR_plasma <- volcano_plot_ALS_CTR(DE_subgroups_plasma_ALS_CTR,"log2FC_ALS_CTR",
                                       title = "Differential expression analysis ALS vs CTR in plasma")
pdf("plots/volcano_ALS_CTR_plasma.pdf")
volcano_ALS_CTR_plasma
dev.off()

# ALS vs CTR (CSF)
DE_subgroups_CSF_ALS_CTR <- DE_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_CTR_ALS,padj_CTR_ALS) %>%
  mutate(log10_pval = -log10(pvalue_CTR_ALS),
         log10_padj = -log10(padj_CTR_ALS),
         DE = ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR > 0,"ALS",
                     ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR < 0,"CTR",
                            "ns")))

volcano_ALS_CTR_CSF <- volcano_plot_ALS_CTR(DE_subgroups_CSF_ALS_CTR,"log2FC_ALS_CTR",
                                               title = "Differential expression analysis ALS vs CTR in CSF")

pdf("plots/volcano_ALS_CTR_CSF.pdf")
volcano_ALS_CTR_CSF
dev.off()

# ALS vs CTR (SERUM)
DE_subgroups_SERUM_ALS_CTR <- DE_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_ALS_CTR,pvalue_CTR_ALS,padj_CTR_ALS) %>%
  mutate(log10_pval = -log10(pvalue_CTR_ALS),
         log10_padj = -log10(padj_CTR_ALS),
         DE = ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR > 0,"ALS",
                     ifelse(padj_CTR_ALS<0.05 & log2FC_ALS_CTR < 0,"CTR",
                            "ns")))

volcano_ALS_CTR_SERUM <- volcano_plot_ALS_CTR(DE_subgroups_SERUM_ALS_CTR,"log2FC_ALS_CTR",
                                            title = "Differential expression analysis ALS vs CTR in SERUM")

pdf("plots/volcano_ALS_CTR_serum.pdf")
volcano_ALS_CTR_SERUM
dev.off()

volcano_plot_ALS_CTR <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('CTR' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'ALS'= "#3da0c2")) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    scale_size_continuous(range=c(2, 7)) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

# Mimic vs CTR
# -> plasma
DE_subgroups_plasma_mimic_CTR <- DE_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_CTR_mimic,padj_CTR_mimic) %>%
  mutate(log10_pval = -log10(pvalue_CTR_mimic),
         log10_padj = -log10(padj_CTR_mimic),
         DE = ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR > 0,"mimic",
                     ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR < 0,"CTR",
                            "ns")))

volcano_mimic_CTR_plasma <- volcano_plot_mimic_CTR(DE_subgroups_plasma_mimic_CTR,"log2FC_mimic_CTR",
                                               title = "Differential expression analysis mimic vs CTR in plasma")
pdf("plots/volcano_mimic_CTR_plasma.pdf")
volcano_mimic_CTR_plasma
dev.off()

# mimic vs CTR (CSF)
DE_subgroups_CSF_mimic_CTR <- DE_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_CTR_mimic,padj_CTR_mimic) %>%
  mutate(log10_pval = -log10(pvalue_CTR_mimic),
         log10_padj = -log10(padj_CTR_mimic),
         DE = ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR > 0,"mimic",
                     ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR < 0,"CTR",
                            "ns")))

volcano_mimic_CTR_CSF <- volcano_plot_mimic_CTR(DE_subgroups_CSF_mimic_CTR,"log2FC_mimic_CTR",
                                            title = "Differential expression analysis mimic vs CTR in CSF")

pdf("plots/volcano_mimic_CTR_CSF.pdf")
volcano_mimic_CTR_CSF
dev.off()

# mimic vs CTR (SERUM)
DE_subgroups_SERUM_mimic_CTR <- DE_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_mimic_CTR,pvalue_CTR_mimic,padj_CTR_mimic) %>%
  mutate(log10_pval = -log10(pvalue_CTR_mimic),
         log10_padj = -log10(padj_CTR_mimic),
         DE = ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR > 0,"mimic",
                     ifelse(padj_CTR_mimic<0.05 & log2FC_mimic_CTR < 0,"CTR",
                            "ns")))

volcano_mimic_CTR_SERUM <- volcano_plot_mimic_CTR(DE_subgroups_SERUM_mimic_CTR,"log2FC_mimic_CTR",
                                              title = "Differential expression analysis mimic vs CTR in SERUM")

pdf("plots/volcano_mimic_CTR_serum.pdf")
volcano_mimic_CTR_SERUM
dev.off()

volcano_plot_mimic_CTR <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('CTR' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'mimic'= "#3da0c2")) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    scale_size_continuous(range=c(2, 7)) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

# CTR vs PGMC
# -> plasma
DE_subgroups_plasma_CTR_PGMC <- DE_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_CTR_PGMC,pvalue_CTR_PGMC,padj_CTR_PGMC) %>%
  mutate(log10_pval = -log10(pvalue_CTR_PGMC),
         log10_padj = -log10(padj_CTR_PGMC),
         DE = ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC > 0,"CTR",
                     ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC < 0,"PGMC",
                            "ns")))

volcano_CTR_PGMC_plasma <- volcano_plot_CTR_PGMC(DE_subgroups_plasma_CTR_PGMC,"log2FC_CTR_PGMC",
                                               title = "Differential expression analysis CTR vs PGMC in plasma")
pdf("plots/volcano_CTR_PGMC_plasma.pdf")
volcano_CTR_PGMC_plasma
dev.off()

# CTR vs PGMC(CSF)
DE_subgroups_CSF_CTR_PGMC <- DE_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_CTR_PGMC,pvalue_CTR_PGMC,padj_CTR_PGMC) %>%
  mutate(log10_pval = -log10(pvalue_CTR_PGMC),
         log10_padj = -log10(padj_CTR_PGMC),
         DE = ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC > 0,"CTR",
                     ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC < 0,"PGMC",
                            "ns")))

volcano_CTR_PGMC_CSF <- volcano_plot_CTR_PGMC(DE_subgroups_CSF_CTR_PGMC,"log2FC_CTR_PGMC",
                                            title = "Differential expression analysis CTR vs PGMC in CSF")

pdf("plots/volcano_CTR_PGMC_CSF.pdf")
volcano_CTR_PGMC_CSF
dev.off()

# CTR vs PGMC (SERUM)
DE_subgroups_SERUM_CTR_PGMC <- DE_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_CTR_PGMC,pvalue_CTR_PGMC,padj_CTR_PGMC) %>%
  mutate(log10_pval = -log10(pvalue_CTR_PGMC),
         log10_padj = -log10(padj_CTR_PGMC),
         DE = ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC > 0,"CTR",
                     ifelse(padj_CTR_PGMC<0.05 & log2FC_CTR_PGMC < 0,"PGMC",
                            "ns")))

volcano_CTR_PGMC_SERUM <- volcano_plot_CTR_PGMC(DE_subgroups_SERUM_CTR_PGMC,"log2FC_CTR_PGMC",
                                              title = "Differential expression analysis CTR vs PGMC in SERUM")

pdf("plots/volcano_CTR_PGMC_serum.pdf")
volcano_CTR_PGMC_SERUM
dev.off()

volcano_plot_CTR_PGMC <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('PGMC' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'CTR'= "#3da0c2")) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    scale_size_continuous(range=c(2, 7)) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

# ALS vs PGMC (plasma)
DE_subgroups_plasma_ALS_PGMC <- DE_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_PGMC_ALS,padj_PGMC_ALS) %>%
  mutate(log10_pval = -log10(pvalue_PGMC_ALS),
         log10_padj = -log10(padj_PGMC_ALS),
         DE = ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC > 0,"ALS",
                     ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC < 0,"PGMC",
                            "ns")))

volcano_ALS_PGMC_plasma <- volcano_plot_ALS_PGMC(DE_subgroups_plasma_ALS_PGMC,"log2FC_ALS_PGMC",
                                               title = "Differential expression analysis ALS vs PGMC in plasma")
pdf("plots/volcano_ALS_PGMC_plasma.pdf")
volcano_ALS_PGMC_plasma
dev.off()

# ALS vs PGMC (CSF)
DE_subgroups_CSF_ALS_PGMC <- DE_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_PGMC_ALS,padj_PGMC_ALS) %>%
  mutate(log10_pval = -log10(pvalue_PGMC_ALS),
         log10_padj = -log10(padj_PGMC_ALS),
         DE = ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC > 0,"ALS",
                     ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC < 0,"PGMC",
                            "ns")))

volcano_ALS_PGMC_CSF <- volcano_plot_ALS_PGMC(DE_subgroups_CSF_ALS_PGMC,"log2FC_ALS_PGMC",
                                            title = "Differential expression analysis ALS vs PGMC in CSF")

pdf("plots/volcano_ALS_PGMC_CSF.pdf")
volcano_ALS_PGMC_CSF
dev.off()

# ALS vs PGMC (SERUM)
DE_subgroups_SERUM_ALS_PGMC <- DE_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_ALS_PGMC,pvalue_PGMC_ALS,padj_PGMC_ALS) %>%
  mutate(log10_pval = -log10(pvalue_PGMC_ALS),
         log10_padj = -log10(padj_PGMC_ALS),
         DE = ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC > 0,"ALS",
                     ifelse(padj_PGMC_ALS<0.05 & log2FC_ALS_PGMC < 0,"PGMC",
                            "ns")))

volcano_ALS_PGMC_SERUM <- volcano_plot_ALS_PGMC(DE_subgroups_SERUM_ALS_PGMC,"log2FC_ALS_PGMC",
                                              title = "Differential expression analysis ALS vs PGMC in SERUM")

pdf("plots/volcano_ALS_PGMC_serum.pdf")
volcano_ALS_PGMC_SERUM
dev.off()

volcano_plot_ALS_PGMC <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('PGMC' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'ALS'= "#3da0c2")) +
    scale_size_continuous(range=c(2, 7)) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

## C9orf72 vs SOD1
# C9orf72 vs SOD1 (PLASMA)
DE_subgroups_PLASMA_C9orf72_SOD1 <- DE_PGMC_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_C9orf72_SOD1),
         log10_padj = -log10(padj_C9orf72_SOD1),
         DE = ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 > 0,"C9orf72",
                     ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 < 0,"SOD1",
                            "ns")))

volcano_C9orf72_SOD1_PLASMA <- volcano_plot_C9orf72_SOD1(DE_subgroups_PLASMA_C9orf72_SOD1,"log2FC_C9orf72_SOD1",
                                                        title = "Differential expression analysis C9orf72 vs SOD1 in PLASMA")

pdf("plots/volcano_C9orf72_SOD1_PLASMA.pdf")
volcano_C9orf72_SOD1_PLASMA
dev.off()

# C9orf72 vs SOD1 (CSF)
DE_subgroups_CSF_C9orf72_SOD1 <- DE_PGMC_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_C9orf72_SOD1),
         log10_padj = -log10(padj_C9orf72_SOD1),
         DE = ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 > 0,"C9orf72",
                     ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 < 0,"SOD1",
                            "ns")))

volcano_C9orf72_SOD1_CSF <- volcano_plot_C9orf72_SOD1(DE_subgroups_CSF_C9orf72_SOD1,"log2FC_C9orf72_SOD1",
                                                        title = "Differential expression analysis C9orf72 vs SOD1 in CSF")

pdf("plots/volcano_C9orf72_SOD1_CSF.pdf")
volcano_C9orf72_SOD1_CSF
dev.off()

# C9orf72 vs SOD1 (SERUM)
DE_subgroups_SERUM_C9orf72_SOD1 <- DE_PGMC_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_C9orf72_SOD1,pvalue_C9orf72_SOD1,padj_C9orf72_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_C9orf72_SOD1),
         log10_padj = -log10(padj_C9orf72_SOD1),
         DE = ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 > 0,"C9orf72",
                     ifelse(padj_C9orf72_SOD1<0.05 & log2FC_C9orf72_SOD1 < 0,"SOD1",
                            "ns")))

volcano_C9orf72_SOD1_SERUM <- volcano_plot_C9orf72_SOD1(DE_subgroups_SERUM_C9orf72_SOD1,"log2FC_C9orf72_SOD1",
                                                title = "Differential expression analysis C9orf72 vs SOD1 in SERUM")

pdf("plots/volcano_C9orf72_SOD1_serum.pdf")
volcano_C9orf72_SOD1_SERUM
dev.off()

volcano_plot_C9orf72_SOD1 <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('SOD1' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'C9orf72'= "#3da0c2")) +
    scale_size_continuous(range=c(2, 7)) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

## C9orf72 vs CTR
# C9orf72 vs CTR (PLASMA)
DE_subgroups_PLASMA_C9orf72_CTR <- DE_PGMC_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_CTR_C9orf72,padj_CTR_C9orf72) %>%
  mutate(log10_pval = -log10(pvalue_CTR_C9orf72),
         log10_padj = -log10(padj_CTR_C9orf72),
         DE = ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR > 0,"C9orf72",
                     ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR < 0,"CTR",
                            "ns")))

volcano_C9orf72_CTR_PLASMA <- volcano_plot_C9orf72_CTR(DE_subgroups_PLASMA_C9orf72_CTR,"log2FC_C9orf72_CTR",
                                                         title = "Differential expression analysis C9orf72 vs CTR in PLASMA")

pdf("plots/volcano_C9orf72_CTR_PLASMA.pdf")
volcano_C9orf72_CTR_PLASMA
dev.off()

# C9orf72 vs CTR (CSF)
DE_subgroups_CSF_C9orf72_CTR <- DE_PGMC_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_CTR_C9orf72,padj_CTR_C9orf72) %>%
  mutate(log10_pval = -log10(pvalue_CTR_C9orf72),
         log10_padj = -log10(padj_CTR_C9orf72),
         DE = ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR > 0,"C9orf72",
                     ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR < 0,"CTR",
                            "ns")))

volcano_C9orf72_CTR_CSF <- volcano_plot_C9orf72_CTR(DE_subgroups_CSF_C9orf72_CTR,"log2FC_C9orf72_CTR",
                                                      title = "Differential expression analysis C9orf72 vs CTR in CSF")

pdf("plots/volcano_C9orf72_CTR_CSF.pdf")
volcano_C9orf72_CTR_CSF
dev.off()

# C9orf72 vs CTR (SERUM)
DE_subgroups_SERUM_C9orf72_CTR <- DE_PGMC_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_C9orf72_CTR,pvalue_CTR_C9orf72,padj_CTR_C9orf72) %>%
  mutate(log10_pval = -log10(pvalue_CTR_C9orf72),
         log10_padj = -log10(padj_CTR_C9orf72),
         DE = ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR > 0,"C9orf72",
                     ifelse(padj_CTR_C9orf72<0.05 & log2FC_C9orf72_CTR < 0,"CTR",
                            "ns")))

volcano_C9orf72_CTR_SERUM <- volcano_plot_C9orf72_CTR(DE_subgroups_SERUM_C9orf72_CTR,"log2FC_C9orf72_CTR",
                                                        title = "Differential expression analysis C9orf72 vs CTR in SERUM")

pdf("plots/volcano_C9orf72_CTR_serum.pdf")
volcano_C9orf72_CTR_SERUM
dev.off()

volcano_plot_C9orf72_CTR <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('CTR' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'C9orf72'= "#3da0c2")) +
    scale_size_continuous(range=c(2, 7)) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

## SOD1 vs CTR
# SOD1 vs CTR (PLASMA)
DE_subgroups_PLASMA_SOD1_CTR <- DE_PGMC_subgroups_plasma %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_CTR_SOD1,padj_CTR_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_CTR_SOD1),
         log10_padj = -log10(padj_CTR_SOD1),
         DE = ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR > 0,"SOD1",
                     ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR < 0,"CTR",
                            "ns")))

volcano_SOD1_CTR_PLASMA <- volcano_plot_SOD1_CTR(DE_subgroups_PLASMA_SOD1_CTR,"log2FC_SOD1_CTR",
                                                       title = "Differential expression analysis SOD1 vs CTR in PLASMA")

pdf("plots/volcano_SOD1_CTR_PLASMA.pdf")
volcano_SOD1_CTR_PLASMA
dev.off()

# SOD1 vs CTR (CSF)
DE_subgroups_CSF_SOD1_CTR <- DE_PGMC_subgroups_CSF %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_CTR_SOD1,padj_CTR_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_CTR_SOD1),
         log10_padj = -log10(padj_CTR_SOD1),
         DE = ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR > 0,"SOD1",
                     ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR < 0,"CTR",
                            "ns")))

volcano_SOD1_CTR_CSF <- volcano_plot_SOD1_CTR(DE_subgroups_CSF_SOD1_CTR,"log2FC_SOD1_CTR",
                                                    title = "Differential expression analysis SOD1 vs CTR in CSF")

pdf("plots/volcano_SOD1_CTR_CSF.pdf")
volcano_SOD1_CTR_CSF
dev.off()

# SOD1 vs CTR (SERUM)
DE_subgroups_SERUM_SOD1_CTR <- DE_PGMC_subgroups_SERUM %>%
  select(Target,UniProtID,log2FC_SOD1_CTR,pvalue_CTR_SOD1,padj_CTR_SOD1) %>%
  mutate(log10_pval = -log10(pvalue_CTR_SOD1),
         log10_padj = -log10(padj_CTR_SOD1),
         DE = ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR > 0,"SOD1",
                     ifelse(padj_CTR_SOD1<0.05 & log2FC_SOD1_CTR < 0,"CTR",
                            "ns")))

volcano_SOD1_CTR_SERUM <- volcano_plot_SOD1_CTR(DE_subgroups_SERUM_SOD1_CTR,"log2FC_SOD1_CTR",
                                                      title = "Differential expression analysis SOD1 vs CTR in SERUM")

pdf("plots/volcano_SOD1_CTR_serum.pdf")
volcano_SOD1_CTR_SERUM
dev.off()

volcano_plot_SOD1_CTR <- function(dataset,log2fc,title){
  ggplot(dataset,aes_string(log2fc,"log10_padj")) + 
    geom_point(aes(color=DE, size = log10_padj)) +
    geom_text_repel(data = dataset %>% filter(DE!="ns"),
                    aes(label = Target),max.overlaps = 40) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dashed",col = "darkgrey") +
    geom_vline(xintercept = 0,linetype = "dashed",col = "darkgrey") +
    scale_color_manual(values  = c('CTR' = "#C25F3D", 
                                   'ns' = 'lightgrey', 
                                   'SOD1'= "#3da0c2")) +
    scale_size_continuous(range=c(2, 7)) +
    #annotate(geom="text", x=0.75, y= (-log10(0.05)) + 0.15, label="FDR = 5%",size = 5) +
    theme_minimal() + 
    ylab(expression("-log"[10]*"(adjusted p-value)")) +
    xlab(expression("log"[2]*"(fold-change)")) + 
    ggtitle(title) +
    labs(size = expression("-log"[10]*"(adjusted p-value)"))
}

## signed p-values 
# ALS vs CTR + ALS vs PGMC
colors = c('ALS' = '#32A55E',
  'PGMC' = '#96AA9A',
  'ALS sign. both' = '#5E718B',
  'ns' = '#B4BFC5',
  'CTR' = '#CF7041',
  "CTR and PGMC" = "#E5C5BD")
signed_pvalue_ALS_CTR_PGMC <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_ALS_CTR,signed_ALS_PGMC)) +
    geom_point(aes(colour = DE_new),size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new),size = 2.25, alpha = 0.5) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    theme(plot.title = element_text(size = 15, face = "bold"))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("ALS vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("ALS vs PGMC (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+0.9,y= log10(0.05)-0.15,
             color = "grey") +
    annotate(geom = "text", label = "FDR = 5%",x=log10(0.05)-0.9,y= - log10(0.05)+0.15,
             color = "grey") +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                     aes(label=Target), size=4,  color="black",
                     nudge_y = 0.5,
                     nudge_x = 0.5, fontface="bold",  
                     min.segment.length = unit(0.1, "lines"), 
                     max.overlaps = 40) +
    scale_colour_manual(values = colors) }
# -> Plasma
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_plasma_ALS_CTR %>%
  mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
  rename(DE_ALS_CTR = DE) %>%
  select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
  DE_subgroups_plasma_ALS_PGMC %>%
    mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
    rename(DE_ALS_PGMC = DE) %>%
    select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
  ) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                         ifelse(DE_ALS_CTR == "CTR","CTR",
                                ifelse(DE_ALS_CTR == "ALS","ALS",
                                       ifelse(DE_ALS_PGMC == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

pdf("plots/signed_ALSvsCTR_ALSvsPGMC_plasma.pdf")
signed_pvalue_ALS_CTR_PGMC(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (PLASMA)")
dev.off()

# -> CSF
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_CSF_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_CSF_ALS_PGMC %>%
                                                mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
                                                rename(DE_ALS_PGMC = DE) %>%
                                                select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "ALS","ALS",
                                                     ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

pdf("plots/signed_ALSvsCTR_ALSvsPGMC_CSF.pdf")
signed_pvalue_ALS_CTR_PGMC(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (CSF)")
dev.off()

# -> Serum
signed_subgroups_ALS_CTR_and_ALS_PGMC = merge(DE_subgroups_SERUM_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_SERUM_ALS_PGMC %>%
                                                mutate(signed_ALS_PGMC = sign(log2FC_ALS_PGMC)*log10_padj) %>%
                                                rename(DE_ALS_PGMC = DE) %>%
                                                select(Target,UniProtID,signed_ALS_PGMC,DE_ALS_PGMC)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "ALS" & DE_ALS_PGMC == "ALS","ALS sign. both",
                         ifelse(DE_ALS_CTR == "CTR" & DE_ALS_PGMC == "PGMC","CTR and PGMC",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_ALS_PGMC == "ALS","ALS",
                                                     ifelse(DE_ALS_PGMC == "PGMC","PGMC","ns")))))))

pdf("plots/signed_ALSvsCTR_ALSvsPGMC_SERUM.pdf")
signed_pvalue_ALS_CTR_PGMC(signed_subgroups_ALS_CTR_and_ALS_PGMC,"ALS vs CTR + ALS vs PGMC (SERUM)")
dev.off()

# ALS vs CTR + Mimic vs CTR 
colors = c('ALS' = '#32A55E',
           'mimic' = '#96AA9A',
           'CTR sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'CTR' = '#CF7041',
           "ALS and mimic" = "#E5C5BD")
signed_pvalue_ALS_CTR_mimic <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_ALS_CTR,signed_mimic_CTR)) +
    geom_point(aes(colour = DE_new),size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new),size = 2.25, alpha = 0.5) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    theme(plot.title = element_text(size = 15, face = "bold"))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("ALS vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("mimic vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+0.9,y= log10(0.05)-0.15,
             color = "grey") +
    annotate(geom = "text", label = "FDR = 5%",x=log10(0.05)-0.9,y= - log10(0.05)+0.15,
             color = "grey") +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                    aes(label=Target), size=4,  color="black",
                    nudge_y = 0.5,
                    nudge_x = 0.5, fontface="bold",  
                    min.segment.length = unit(0.1, "lines"), 
                    max.overlaps = 40) +
    scale_colour_manual(values = colors) }
# -> Plasma
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_plasma_ALS_CTR %>%
                                                mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                rename(DE_ALS_CTR = DE) %>%
                                                select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                              DE_subgroups_plasma_mimic_CTR %>%
                                                mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                rename(DE_mimic_CTR = DE) %>%
                                                select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                         ifelse(DE_ALS_CTR == "CTR","CTR",
                                ifelse(DE_ALS_CTR == "ALS","ALS",
                                       ifelse(DE_mimic_CTR == "mimic","mimic",
                                              ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_plasma.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (PLASMA)")
dev.off()

# -> CSF
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_CSF_ALS_CTR %>%
                                                 mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                 rename(DE_ALS_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                               DE_subgroups_CSF_mimic_CTR %>%
                                                 mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                 rename(DE_mimic_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_mimic_CTR == "mimic","mimic",
                                                     ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_CSF.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (CSF)")
dev.off()

# -> Serum
signed_subgroups_ALS_CTR_and_mimic_CTR = merge(DE_subgroups_SERUM_ALS_CTR %>%
                                                 mutate(signed_ALS_CTR = sign(log2FC_ALS_CTR)*log10_padj) %>%
                                                 rename(DE_ALS_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_ALS_CTR,DE_ALS_CTR),
                                               DE_subgroups_SERUM_mimic_CTR %>%
                                                 mutate(signed_mimic_CTR = sign(log2FC_mimic_CTR)*log10_padj) %>%
                                                 rename(DE_mimic_CTR = DE) %>%
                                                 select(Target,UniProtID,signed_mimic_CTR,DE_mimic_CTR)
) %>%
  mutate(DE_new = ifelse(DE_ALS_CTR == "CTR" & DE_mimic_CTR == "CTR","CTR sign. both",
                         ifelse(DE_ALS_CTR == "ALS" & DE_mimic_CTR == "mimic","ALS and mimic",
                                ifelse(DE_ALS_CTR == "CTR","CTR",
                                       ifelse(DE_ALS_CTR == "ALS","ALS",
                                              ifelse(DE_mimic_CTR == "mimic","mimic",
                                                     ifelse(DE_mimic_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_ALSvsCTR_mimicvsCTR_SERUM.pdf")
signed_pvalue_ALS_CTR_mimic(signed_subgroups_ALS_CTR_and_mimic_CTR,"ALS vs CTR + mimic vs CTR (SERUM)")
dev.off()

# C9orf72 vs CTR + SOD1 vs CTR
colors = c('CTR' = '#32A55E',
           'C9orf72' = '#96AA9A',
           'CTR sign. both' = '#5E718B',
           'ns' = '#B4BFC5',
           'SOD1' = '#CF7041',
           "C9orf72 and SOD1" = "#E5C5BD")
signed_pvalue_C9orf72_SOD1_CTR <- function(data_to_plot,title_plot){
  ggplot(data_to_plot,aes(signed_C9orf72_CTR,signed_SOD1_CTR)) +
    geom_point(aes(colour = DE_new),size = 2.25, color = "grey", alpha = 0.5) +
    geom_point(aes(colour = DE_new),size = 2.25, alpha = 0.5) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line()) +
    ggtitle(title_plot)+
    theme(plot.title = element_text(size = 15, face = "bold"))+
    scale_x_continuous(breaks = seq(-50, 20, by = 5))+
    scale_y_continuous(breaks = seq(-30, 140, by = 10))+
    xlab(paste0("C9orf72 vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    ylab(paste0("SOD1 vs CTR (sign(LFC) x -log10(p-adj))"))+
    theme(axis.title = element_text(size = 15), axis.text=element_text(size=10))+
    geom_vline(xintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_vline(xintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",color = "grey",alpha = 0.6)+
    #geom_abline(linetype="dashed",color = "grey",alpha = 0.6) +
    annotate(geom = "text", label = "FDR = 5%",x= -log10(0.05)+0.9,y= log10(0.05)-0.15,
             color = "grey") +
    annotate(geom = "text", label = "FDR = 5%",x=log10(0.05)-0.9,y= - log10(0.05)+0.15,
             color = "grey") +
    geom_text_repel(data = subset(data_to_plot, DE_new != "ns"), 
                    aes(label=Target), size=4,  color="black",
                    nudge_y = 0.5,
                    nudge_x = 0.5, fontface="bold",  
                    min.segment.length = unit(0.1, "lines"), 
                    max.overlaps = 40) +
    scale_colour_manual(values = colors) }
# -> Plasma
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_PLASMA_C9orf72_CTR %>%
                                                mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                rename(DE_C9orf72_CTR = DE) %>%
                                                select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                              DE_subgroups_PLASMA_SOD1_CTR %>%
                                                mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                rename(DE_SOD1_CTR = DE) %>%
                                                select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_plasma.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (PLASMA)")
dev.off()

# -> CSF
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_CSF_C9orf72_CTR %>%
                                                    mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                    rename(DE_C9orf72_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                                  DE_subgroups_CSF_SOD1_CTR %>%
                                                    mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                    rename(DE_SOD1_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_CSF.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (CSF)")
dev.off()

# -> SERUM
signed_subgroups_C9orf72_CTR_and_SOD1_CTR = merge(DE_subgroups_SERUM_C9orf72_CTR %>%
                                                    mutate(signed_C9orf72_CTR = sign(log2FC_C9orf72_CTR)*log10_padj) %>%
                                                    rename(DE_C9orf72_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_C9orf72_CTR,DE_C9orf72_CTR),
                                                  DE_subgroups_SERUM_SOD1_CTR %>%
                                                    mutate(signed_SOD1_CTR = sign(log2FC_SOD1_CTR)*log10_padj) %>%
                                                    rename(DE_SOD1_CTR = DE) %>%
                                                    select(Target,UniProtID,signed_SOD1_CTR,DE_SOD1_CTR)
) %>%
  mutate(DE_new = ifelse(DE_C9orf72_CTR == "CTR" & DE_SOD1_CTR == "CTR","CTR sign. both",
                         ifelse(DE_C9orf72_CTR == "C9orf72" & DE_SOD1_CTR == "SOD1","C9orf72 and SOD1",
                                ifelse(DE_C9orf72_CTR == "CTR","CTR",
                                       ifelse(DE_C9orf72_CTR == "C9orf72","C9orf72",
                                              ifelse(DE_SOD1_CTR == "SOD1","SOD1",
                                                     ifelse(DE_SOD1_CTR == "CTR","CTR","ns")))))))

pdf("plots/signed_C9orf72vsCTR_SOD1vsCTR_SERUM.pdf")
signed_pvalue_C9orf72_SOD1_CTR(signed_subgroups_C9orf72_CTR_and_SOD1_CTR,"C9orf72 vs CTR + SOD1 vs CTR (SERUM)")
dev.off()


