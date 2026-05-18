# Differential expression analyses of the interleukin's group all together

###############################################
### Helper Functions
###############################################

# 1. Get Interleukins
extract_il_ligands <- function(df) {
  df %>% 
    filter(str_detect(Target, "^IL[0-9]")) %>% 
    filter(!str_detect(Target, "R(L)?[0-9]*$|BP|RB|RA|RN|ST|RAP"))
}

# 2. Plot Interleukins trends in each fluid
plot_il_trends <- function(data, stats, labels_df, comp_name, title_suffix) {
  
  plot_df <- data %>% filter(Comparison == comp_name)
  plot_stats <- stats %>% filter(Comparison == comp_name)
  plot_labels <- labels_df %>% filter(Comparison == comp_name)
  # define seed 
  jitter_pos <- position_jitter(width = 0.1, seed = 42)
  
  ggplot(plot_df, aes(x = Fluid, y = LFC, fill = Fluid)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(position = jitter_pos, alpha = 0.7, size = 2) +
    geom_text_repel(
      data = plot_labels,
      aes(label = Target),
      position = jitter_pos,
      size = 3.5,
      fontface = "italic",
      box.padding = 0.2,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = 'grey50',
      max.overlaps = Inf) +
    geom_text(data = plot_stats, aes(x = Fluid, y = max_lfc + 0.3, label = p_label), 
              fontface = "bold", size = 4) +
    scale_fill_brewer(palette = "Pastel1") +
    theme_classic(base_size = 14) +
    labs(title = paste("Interleukin family distribution for", title_suffix),
         subtitle = "Signed Wilcoxon Test",
         y = "Log2 Fold Change", x = "") +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
}

plot_il_trends_sex <- function(data, stats, labels_df, comp_name, title_suffix) {
  
  plot_df <- data %>% filter(Comparison == comp_name)
  plot_stats <- stats %>% filter(Comparison == comp_name)
  plot_labels <- labels_df %>% filter(Comparison == comp_name)
  
  plot_df$Sex <- factor(plot_df$Sex, levels = c("all", "female", "male"))
  plot_stats$Sex <- factor(plot_stats$Sex, levels = c("all", "female", "male"))
  plot_labels$Sex <- factor(plot_labels$Sex, levels = c("all", "female", "male"))
  
  jitter_pos <- position_jitter(width = 0.1, seed = 42)
  
  ggplot(plot_df, aes(x = Fluid, y = LFC, fill = Fluid)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(position = jitter_pos, alpha = 0.7, size = 2) +
    geom_text_repel(
      data = plot_labels,
      aes(label = Target),
      position = jitter_pos,
      size = 3.2,
      fontface = "italic",
      box.padding = 0.2,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = 'grey50',
      max.overlaps = Inf
    ) +
    geom_text(
      data = plot_stats,
      aes(x = Fluid, y = max_lfc + 0.3, label = p_label),
      fontface = "bold",
      size = 3.8
    ) +
    scale_fill_brewer(palette = "Pastel1") +
    facet_wrap(~Sex, nrow = 1) +  
    theme_classic(base_size = 14) +
    labs(
      title = paste("Interleukin family distribution for", title_suffix),
      subtitle = "All and gender stratification",
      y = "Log2 Fold Change",
      x = ""
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),  
      panel.spacing = unit(1.2, "lines")
    )
}

plot_il_trends_sex_CSF <- function(data, stats, labels_df, comp_name, title_suffix) {
  
  plot_df <- data %>% filter(Comparison == comp_name)
  plot_stats <- stats %>% filter(Comparison == comp_name)
  plot_labels <- labels_df %>% filter(Comparison == comp_name)
  
  # clean labels
  plot_df$Sex <- factor(plot_df$Sex, levels = c("all", "female", "male"))
  plot_stats$Sex <- factor(plot_stats$Sex, levels = c("all", "female", "male"))
  plot_labels$Sex <- factor(plot_labels$Sex, levels = c("all", "female", "male"))
  
  # jitter
  jitter_pos <- position_jitterdodge(
    jitter.width = 0.15,
    dodge.width = 0.8,
    seed = 42
  )
  
  ggplot(plot_df, aes(x = Fluid, y = LFC, fill = Sex)) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    
    geom_violin(
      aes(group = interaction(Fluid, Sex),alpha = Sex),
      position = position_dodge(width = 0.8),
      color = NA,
      width = 0.6
    ) +
    
    geom_boxplot(
      aes(group = interaction(Fluid, Sex)),
      position = position_dodge(width = 0.8),
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.7
    ) +
    
    geom_jitter(
      aes(color = Sex),
      position = jitter_pos,
      alpha = 0.8,
      size = 2,
      show.legend = FALSE
    ) +
    
    geom_text_repel(
      data = plot_labels,
      aes(label = Target, color = Sex),
      position = jitter_pos,
      size = 3,
      fontface = "italic",
      box.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = plot_stats,
      aes(x = Fluid, y = max_lfc + 0.35, label = p_label, group = Sex),
      position = position_dodge(width = 0.8),
      size = 3.5,
      fontface = "bold"
    ) +
    
    scale_alpha_manual(values = c(
      "all" = 0.2,
      "female" = 0.35,
      "male" = 0.35
    )) +
    
    scale_fill_manual(values = c(
      "all" = "#7f7f7f",
      "female" = "#dd8452",
      "male" = "#4c72b0"
    )) +
    
    scale_color_manual(values = c(
      "all" = "#7f7f7f",
      "female" = "#dd8452",
      "male" = "#4c72b0"
    )) +
    
    theme_classic(base_size = 14) +
    
    labs(
      title = title_suffix,
      y = "Log2 Fold Change",
      x = "",
      fill = "Group"
    ) +
    
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(size=10),
      axis.text.y = element_text(size = 10)) + 
    guides(alpha = "none") 
}

# 3. Get top n drivers based on LFC
get_top_drivers <- function(data, comp,n_top = 10) {
  data %>%
    filter(Comparison == comp) %>%
    group_by(Fluid) %>%
    slice_max(order_by = abs(LFC), n = n_top) %>%
    arrange(Fluid, desc(LFC)) %>%
    select(Fluid, Target, LFC)
}

###############################################################################
# Run pipeline
###############################################################################

## 1. Enrichment of interleukins in ALS and PGMCs

# 1.1 Check log2FC trends in the fluids for ALS vs CTR and PGMC vs CTR
combined_ils <- list()

for(t in names(DE_results)) {
  df_il <- DE_results[[t]]$standard_adj %>%
    extract_il_ligands() %>%
    mutate(Fluid = t) %>%
    select(Target, log2FC_ALS_CTR, log2FC_PGMC_CTR, Fluid)
  
  combined_ils[[t]] <- df_il
}

il_long <- bind_rows(combined_ils) %>%
  pivot_longer(cols = c(log2FC_ALS_CTR, log2FC_PGMC_CTR), 
               names_to = "Comparison", 
               values_to = "LFC") %>%
  mutate(Comparison = gsub("log2FC_", "", Comparison))

stats_summary <- il_long %>%
  group_by(Fluid, Comparison) %>%
  summarise(
    p_val = wilcox.test(LFC, mu = 0)$p.value,
    median_lfc = median(LFC),
    max_lfc = max(LFC),
    .groups = 'drop'
  ) %>%
  mutate(p_label = ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 4))))

top_drivers_als <- get_top_drivers(il_long, "ALS_CTR",8)
top_drivers_pgmc <- get_top_drivers(il_long, "PGMC_CTR",8)

# Combine top driver dataframes 
all_top_drivers <- bind_rows(
  top_drivers_als %>% mutate(Comparison = "ALS_CTR"),
  top_drivers_pgmc %>% mutate(Comparison = "PGMC_CTR")
)

# Create Plots with labels
p_als <- plot_il_trends(il_long, stats_summary, all_top_drivers, "ALS_CTR", "ALS vs Controls")
p_pgmc <- plot_il_trends(il_long, stats_summary, all_top_drivers, "PGMC_CTR", "PGMC vs Controls")

# Save plots
pdf("plots/Interleukin_Collective_Analysis.pdf", width = 8, height = 6)
print(p_als)
print(p_pgmc)
dev.off()

# 1.2 Check log2FC trends in CSF for ALS vs CTR and PGMC vs CTR, and stratified by gender 
combined_ils <- list()

for(t in names(DE_results)) {
  df_il <- DE_results[[t]]$standard_adj %>%
    extract_il_ligands() %>%
    mutate(Fluid = t) %>%
    select(Target, log2FC_ALS_CTR, log2FC_PGMC_CTR, Fluid)
  
  df_il_female = summarise_mean_NPQ(results_ALL[[t]]$data_adjusted %>%
                                      filter(sex == "F"), 
                                    matrix_type = t, adjusted = TRUE) %>%
    add_log2fc_standard() %>%
    extract_il_ligands() %>%
    mutate(Fluid = t) %>%
    select(Target, log2FC_ALS_CTR, log2FC_PGMC_CTR, Fluid) %>%
    rename(log2FC_ALS_CTR_female = log2FC_ALS_CTR,
           log2FC_PGMC_CTR_female = log2FC_PGMC_CTR)
  
  df_il_male = summarise_mean_NPQ(results_ALL[[t]]$data_adjusted %>%
                                      filter(sex == "M"), 
                                    matrix_type = t, adjusted = TRUE) %>%
    add_log2fc_standard() %>%
    extract_il_ligands() %>%
    mutate(Fluid = t) %>%
    select(Target, log2FC_ALS_CTR, log2FC_PGMC_CTR, Fluid) %>%
    rename(log2FC_ALS_CTR_male = log2FC_ALS_CTR,
           log2FC_PGMC_CTR_male = log2FC_PGMC_CTR)
  
  
  combined_ils[[t]] <- df_il %>%
    left_join(df_il_female) %>%
    left_join(df_il_male)
}

il_long <- bind_rows(combined_ils) %>%
  pivot_longer(
    cols = starts_with("log2FC"),
    names_to = c("Comparison", "Sex"),
    names_pattern = "log2FC_(ALS_CTR|PGMC_CTR)(?:_(female|male))?",
    values_to = "LFC") %>%
  mutate(
    Sex = ifelse(Sex == "", "all", Sex))

stats_summary <- il_long %>%
  group_by(Fluid, Comparison, Sex) %>%
  summarise(
    p_val = wilcox.test(LFC, mu = 0)$p.value,
    median_lfc = median(LFC),
    max_lfc = max(LFC),
    .groups = 'drop'
  ) %>%
  mutate(
    p_label = ifelse(p_val < 0.001, "p < 0.001",
                     paste0("p = ", round(p_val, 4))))

top_drivers_als <- il_long %>%
  filter(Comparison == "ALS_CTR") %>%
  group_by(Fluid, Sex) %>%
  slice_max(abs(LFC), n = 8) %>%
  ungroup()

top_drivers_pgmc <- il_long %>%
  filter(Comparison == "PGMC_CTR") %>%
  group_by(Fluid, Sex) %>%
  slice_max(abs(LFC), n = 8) %>%
  ungroup()

all_top_drivers <- bind_rows(top_drivers_als, top_drivers_pgmc)

# Create Plots with labels
p_als <- plot_il_trends_sex(il_long,stats_summary,all_top_drivers,"ALS_CTR","ALS vs Controls")

p_pgmc <- plot_il_trends_sex(il_long,stats_summary,all_top_drivers,"PGMC_CTR","PGMC vs Controls")

# Save plots
pdf("plots/Interleukin_Collective_Analysis_sex.pdf", width = 12, height = 6)
print(p_als)
print(p_pgmc)
dev.off()

# 1.3 Log2FC trends only for CSF
il_long_CSF <- il_long %>% filter(Fluid == "CSF")
stats_summary_CSF <- stats_summary %>% filter(Fluid == "CSF")
all_top_drivers_CSF <- all_top_drivers %>% filter(Fluid == "CSF")

p_A <- plot_il_trends_sex_CSF(il_long_CSF,stats_summary_CSF,all_top_drivers_CSF,
  "ALS_CTR","ALS vs Controls")

p_B <- plot_il_trends_sex_CSF(il_long_CSF,stats_summary_CSF,all_top_drivers_CSF,
  "PGMC_CTR","PGMC vs Controls")

# Save plots
pdf("plots/Interleukin_Collective_Analysis_sex_CSF.pdf", width = 6, height = 6)
print(p_A)
print(p_B)
dev.off()


## 2. Make heatmap of IL33 with the ALSFRSR sections
## First approach: based on the correlations between IL33 and NPQ_adjusted data of CSF
CSF_data <- results_ALL[["CSF"]]$data_adjusted

CSF_data_ALSFRS <- inner_join(CSF_data, ALSFRS_NULISA_V0) %>%
  filter(Target == "IL33") %>%
  filter(type == "ALS") %>%
  filter(!is.na(ALSFRS_1) & !is.na(NPQ_adj)) %>%
  left_join(disease_duration %>%
              ungroup() %>%
              mutate(
                across(
                  starts_with("Gs"),
                  ~ recode(as.numeric(.),
                           `0` = 4,
                           `1` = 3,
                           `2` = 2,
                           `3` = 1,
                           `4` = 0),
                  .names = "{.col}")) %>%
              select(ParticipantCode,PatientID,contains("Gs"),Visit) %>%
              filter(Visit == "V0")) %>%
  mutate(ALSFRS_bulbar = rowSums(select(., Gs1:Gs3), na.rm = TRUE),
         ALSFRS_spinal = rowSums(select(., Gs4:Gs9), na.rm = TRUE),
         ALSFRS_respiratory = rowSums(select(., Gs10:Gs12), na.rm = TRUE))

cor_test_fun <- function(x, y) {
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  c(cor = test$estimate, p = test$p.value)
}

# for each ALSFRS-R question
results_Qs <- lapply(paste0("Gs", 1:12), function(g) {
  out <- cor_test_fun(CSF_data_ALSFRS[[g]], CSF_data_ALSFRS$NPQ_adj)
  data.frame(
    Question = g,
    Correlation = out["cor.rho"],
    p_value = out["p"],
    Type = "Question"
  )
})

# for each category of ALSFRS-R question
results_domains <- lapply(c("ALSFRS_bulbar","ALSFRS_spinal","ALSFRS_respiratory"), 
                          function(g) {
  out <- cor_test_fun(CSF_data_ALSFRS[[g]], CSF_data_ALSFRS$NPQ_adj)
  data.frame(
    Question = g,
    Correlation = out["cor.rho"],
    p_value = out["p"],
    Type = "Domain"
  )
})

cor_df <- bind_rows(results_Qs,results_domains)

cor_df <- cor_df %>%
  mutate(
    order = case_when(
      Question == "ALSFRS_bulbar" ~ 13,
      Question == "ALSFRS_spinal" ~ 14,
      Question == "ALSFRS_respiratory" ~ 15,
      TRUE ~ as.numeric(gsub("Gs", "", Question))
    ),
    Label = case_when(
      Question == "ALSFRS_bulbar" ~ "Bulbar related (Q1–3)",
      Question == "ALSFRS_spinal" ~ "Spinal related (Q4–9)",
      Question == "ALSFRS_respiratory" ~ "Respiratory related (Q10–12)",
      Question == "Gs1" ~ "Speech (Q1)",
      Question == "Gs2" ~ "Salivation (Q2)",
      Question == "Gs3" ~ "Swallowing (Q3)",
      Question == "Gs4" ~ "Handwriting (Q4)",
      Question == "Gs5" ~ "Cutting food and handling utensils (Q5)",
      Question == "Gs6" ~ "Dressing and hygiene (Q6)",
      Question == "Gs7" ~ "Turning in bed and adjusting bed clothes (Q7)",
      Question == "Gs8" ~ "Walking (Q8)",
      Question == "Gs9" ~ "Climbing stairs (Q9)",
      Question == "Gs10" ~ "Dyspnea and shortness of breath (Q10)",
      Question == "Gs11" ~ "Sleep disturbance due to breathing (Q11)",
      Question == "Gs12" ~ "Mechanical ventilation (Q12)"
    )
  )

cor_df <- cor_df %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

cor_df$order <- ifelse(cor_df$Type == "Domain",
                       cor_df$order + 2,  
                       cor_df$order)


pdf("plots/IL33_heatmap_version1.pdf", width = 8, height = 6)
ggplot(cor_df, aes(x = "IL33", y = reorder(Label, -order), fill = Correlation)) +
  geom_tile(width = 0.6, height = 0.8) +
  geom_text(aes(label = case_when(
    is.na(Correlation) ~ "NA",
    sig == "" ~ sprintf("%.2f", Correlation),
    TRUE ~ sprintf("%.2f (%s)", Correlation, sig))),
  size = 3.5) +
  scale_fill_gradient2(
    low = "#3B6FB6",   # muted blue
    mid = "white",
    high = "#B63B3B",  # muted red
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = expression("Spearman "*rho)) +
  facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "",
    y = "",
    title = "Correlation of IL33 with ALSFRS-R categories") +
  theme_minimal(base_size = 12) +
  theme( panel.grid = element_blank(),
         strip.text.y = element_text(size = 11, face = "bold"),
         strip.background = element_blank(),
         axis.text.y = element_text(size = 10, color = "black"),
         axis.text.x = element_blank(),
         legend.position = "right",
         legend.title = element_text(size = 10),
         legend.text = element_text(size = 9),
         plot.title = element_text(size = 13, face = "bold", hjust = 0),
         panel.spacing = unit(0.8, "lines"),
         plot.margin = margin(10, 15, 10, 10)) 
dev.off()

## Second approach: put the z-score of the ALSFRS-R and the expression of IL33 separately
df_std <- CSF_data_ALSFRS %>%
  select(ParticipantCode, NPQ_adj, starts_with("Gs")) %>%
  mutate(across(c(contains("Gs")), ~ scale(.)[,1]),
         IL33 = scale(NPQ_adj)[,1]) 

# Reshape for plotting
df_long <- df_std %>%
  pivot_longer(cols = -ParticipantCode, names_to = "Variable", values_to = "Value") %>%
  filter(Variable != "NPQ_adj") %>%
  mutate(
    Domain = case_when(
      Variable %in% paste0("Gs",1:3) ~ "Bulbar",
      Variable %in% paste0("Gs",4:9) ~ "Spinal",
      Variable %in% paste0("Gs",10:12) ~ "Respiratory",
      Variable == "IL33" ~ "IL33"),
    order = case_when(
      Variable == "IL33" ~ 13,
      TRUE ~ as.numeric(gsub("Gs", "", Variable))),
    Variable = case_when(
      Variable == "Gs1" ~ "Speech (Q1)",
      Variable == "Gs2" ~ "Salivation (Q2)",
      Variable == "Gs3" ~ "Swallowing (Q3)",
      Variable == "Gs4" ~ "Handwriting (Q4)",
      Variable == "Gs5" ~ "Cutting food and handling utensils (Q5)",
      Variable == "Gs6" ~ "Dressing and hygiene (Q6)",
      Variable == "Gs7" ~ "Turning in bed and adjusting bed clothes (Q7)",
      Variable == "Gs8" ~ "Walking (Q8)",
      Variable == "Gs9" ~ "Climbing stairs (Q9)",
      Variable == "Gs10" ~ "Dyspnea and shortness of breath (Q10)",
      Variable == "Gs11" ~ "Sleep disturbance due to breathing (Q11)",
      Variable == "Gs12" ~ "Mechanical ventilation (Q12)",
      Variable == "IL33" ~ "IL33"),
    Type = ifelse(Variable == "IL33", "IL33", "ALSFRS")
    )

df_long$Domain <- factor(df_long$Domain,levels = c("Bulbar",
                                                   "Spinal","Respiratory",
                                                   "IL33"))

# plot
pdf("plots/IL33_heatmap_version2.pdf", width = 8, height = 6)
ggplot() +
  # ALSFRS layer
  geom_tile(
    data = subset(df_long, Type=="ALSFRS"),
    aes(x = ParticipantCode, y = reorder(Variable,-order), fill = Value)) +
  scale_fill_gradient2(
    name = "Z-score",
    low = "#1a1de5", mid = "white", high = "#e31c3c",
    na.value = "grey85",
    midpoint = 0, limits = c(-4,1)) +
  ggnewscale::new_scale_fill() +
  # IL33 layer
  geom_tile(
    data = subset(df_long, Variable == "IL33"),
    aes(x = ParticipantCode, y = Variable, fill = Value)
  ) +
  scale_fill_gradient2(
    name = "Scaled Exp.",
    low = "#1EB980", mid = "white", high = "#6A1B9A",
    midpoint = 0, limits = c(-2,2)
  ) +
  facet_grid(Domain ~ ., scales="free_y", space="free_y") +
  labs(
    x = "Patients",
    y = "",
    title = "Patient-level relation between ALSFRS-R and IL33"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(face="bold"),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0)
  )
dev.off()

## Third approach: put the original ALSFRS-R and the expression of IL33 separately
df_std_original <- CSF_data_ALSFRS %>%
  select(ParticipantCode, NPQ_adj, starts_with("Gs")) %>%
  mutate(IL33 = scale(NPQ_adj)[,1]) %>%
  arrange(IL33)

# Reshape for plotting
df_long_original <- df_std_original %>%
  pivot_longer(cols = -ParticipantCode, names_to = "Variable", values_to = "Value") %>%
  filter(Variable != "NPQ_adj") %>%
  mutate(
    Domain = case_when(
      Variable %in% paste0("Gs",1:3) ~ "Bulbar",
      Variable %in% paste0("Gs",4:9) ~ "Spinal",
      Variable %in% paste0("Gs",10:12) ~ "Respiratory",
      Variable == "IL33" ~ "IL33"),
    order = case_when(
      Variable == "IL33" ~ 13,
      TRUE ~ as.numeric(gsub("Gs", "", Variable))),
    Variable = case_when(
      Variable == "Gs1" ~ "Speech (Q1)",
      Variable == "Gs2" ~ "Salivation (Q2)",
      Variable == "Gs3" ~ "Swallowing (Q3)",
      Variable == "Gs4" ~ "Handwriting (Q4)",
      Variable == "Gs5" ~ "Cutting food and handling utensils (Q5)",
      Variable == "Gs6" ~ "Dressing and hygiene (Q6)",
      Variable == "Gs7" ~ "Turning in bed and adjusting bed clothes (Q7)",
      Variable == "Gs8" ~ "Walking (Q8)",
      Variable == "Gs9" ~ "Climbing stairs (Q9)",
      Variable == "Gs10" ~ "Dyspnea and shortness of breath (Q10)",
      Variable == "Gs11" ~ "Sleep disturbance due to breathing (Q11)",
      Variable == "Gs12" ~ "Mechanical ventilation (Q12)",
      Variable == "IL33" ~ "IL33"),
    Type = ifelse(Variable == "IL33", "IL33", "ALSFRS")
  )

df_long_original$Domain <- factor(df_long_original$Domain,levels = c("Bulbar",
                                                   "Spinal","Respiratory",
                                                   "IL33"))

df_long_original$ParticipantCode <- factor(df_long_original$ParticipantCode,
                                           levels = df_std_original$ParticipantCode)

# plot
pdf("plots/IL33_heatmap_version3.pdf", width = 8, height = 6)
ggplot() +
  # ALSFRS layer
  geom_tile(
    data = subset(df_long_original, Type=="ALSFRS"),
    aes(x = ParticipantCode, y = reorder(Variable,-order), fill = Value)) +
  scale_fill_gradient2(
    name = "ALSFRS-R Score",
    low = "#0c1cf3", mid = "#d1d1f0", high = "#f4f5f9",
    na.value = "grey85",
    midpoint = 2, limits = c(0,4)) +
  ggnewscale::new_scale_fill() +
  # IL33 layer
  geom_tile(
    data = subset(df_long_original, Variable == "IL33"),
    aes(x = ParticipantCode, y = Variable, fill = Value)
  ) +
  scale_fill_gradient2(
    name = "Scaled Exp.",
    low = "#1EB980", mid = "white", high = "#6A1B9A",
    midpoint = 0, limits = c(-2,2)
  ) +
  facet_grid(Domain ~ ., scales="free_y", space="free_y") +
  labs(
    x = "Patients",
    y = "",
    title = "Patient-level relation between ALSFRS-R and IL33"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text.y = element_text(face="bold"),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0)
  )
dev.off()
