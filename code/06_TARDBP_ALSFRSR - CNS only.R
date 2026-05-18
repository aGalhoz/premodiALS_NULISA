### ALSFRS-R stratified by categories for TARDBP

## 2. Make heatmap of TARDBP with the ALSFRSR sections
## First approach: based on the correlations between TARDBP and NPQ_adjusted data of SERUM
SERUM_data <- results_ALL[["SERUM"]]$data_adjusted

SERUM_data_ALSFRS <- inner_join(SERUM_data, ALSFRS_NULISA_V0) %>%
  filter(Target == "TARDBP") %>%
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
  out <- cor_test_fun(SERUM_data_ALSFRS[[g]], SERUM_data_ALSFRS$NPQ_adj)
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
                            out <- cor_test_fun(SERUM_data_ALSFRS[[g]], SERUM_data_ALSFRS$NPQ_adj)
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


pdf("plots/TARDBP_heatmap_version1.pdf", width = 8, height = 6)
ggplot(cor_df, aes(x = "TARDBP", y = reorder(Label, -order), fill = Correlation)) +
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
    title = "Correlation of TARDBP with ALSFRS-R categories") +
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

## Second approach: put the z-score of the ALSFRS-R and the expression of TARDBP separately
df_std <- SERUM_data_ALSFRS %>%
  select(ParticipantCode, NPQ_adj, starts_with("Gs")) %>%
  mutate(across(c(contains("Gs")), ~ scale(.)[,1]),
         TARDBP = scale(NPQ_adj)[,1]) 

# Reshape for plotting
df_long <- df_std %>%
  pivot_longer(cols = -ParticipantCode, names_to = "Variable", values_to = "Value") %>%
  filter(Variable != "NPQ_adj") %>%
  mutate(
    Domain = case_when(
      Variable %in% paste0("Gs",1:3) ~ "Bulbar",
      Variable %in% paste0("Gs",4:9) ~ "Spinal",
      Variable %in% paste0("Gs",10:12) ~ "Respiratory",
      Variable == "TARDBP" ~ "TARDBP"),
    order = case_when(
      Variable == "TARDBP" ~ 13,
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
      Variable == "TARDBP" ~ "TARDBP"),
    Type = ifelse(Variable == "TARDBP", "TARDBP", "ALSFRS")
  )

df_long$Domain <- factor(df_long$Domain,levels = c("Bulbar",
                                                   "Spinal","Respiratory",
                                                   "TARDBP"))

# plot
pdf("plots/TARDBP_heatmap_version2.pdf", width = 8, height = 6)
ggplot() +
  # ALSFRS layer
  geom_tile(
    data = subset(df_long, Type=="ALSFRS"),
    aes(x = ParticipantCode, y = reorder(Variable,-order), fill = Value)) +
  scale_fill_gradient2(
    name = "Z-score",
    low = "#1a1de5", mid = "white", high = "#e31c3c",
    na.value = "grey85",
    midpoint = 0, limits = c(-5,2)) +
  ggnewscale::new_scale_fill() +
  # TARDBP layer
  geom_tile(
    data = subset(df_long, Variable == "TARDBP"),
    aes(x = ParticipantCode, y = Variable, fill = Value)
  ) +
  scale_fill_gradient2(
    name = "Scaled Exp.",
    low = "#1EB980", mid = "white", high = "#6A1B9A",
    midpoint = 0, limits = c(-2,4)
  ) +
  facet_grid(Domain ~ ., scales="free_y", space="free_y") +
  labs(
    x = "Patients",
    y = "",
    title = "Patient-level relation between ALSFRS-R and TARDBP"
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

## Third approach: put the original ALSFRS-R and the expression of TARDBP separately
df_std_original <- SERUM_data_ALSFRS %>%
  select(ParticipantCode, NPQ_adj, starts_with("Gs")) %>%
  mutate(TARDBP = scale(NPQ_adj)[,1]) %>%
  arrange(TARDBP)

# Reshape for plotting
df_long_original <- df_std_original %>%
  pivot_longer(cols = -ParticipantCode, names_to = "Variable", values_to = "Value") %>%
  filter(Variable != "NPQ_adj") %>%
  mutate(
    Domain = case_when(
      Variable %in% paste0("Gs",1:3) ~ "Bulbar",
      Variable %in% paste0("Gs",4:9) ~ "Spinal",
      Variable %in% paste0("Gs",10:12) ~ "Respiratory",
      Variable == "TARDBP" ~ "TARDBP"),
    order = case_when(
      Variable == "TARDBP" ~ 13,
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
      Variable == "TARDBP" ~ "TARDBP"),
    Type = ifelse(Variable == "TARDBP", "TARDBP", "ALSFRS")
  )

df_long_original$Domain <- factor(df_long_original$Domain,levels = c("Bulbar",
                                                                     "Spinal","Respiratory",
                                                                     "TARDBP"))

df_long_original$ParticipantCode <- factor(df_long_original$ParticipantCode,
                                           levels = df_std_original$ParticipantCode)

# plot
pdf("plots/TARDBP_heatmap_version3.pdf", width = 8, height = 6)
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
  # TARDBP layer
  geom_tile(
    data = subset(df_long_original, Variable == "TARDBP"),
    aes(x = ParticipantCode, y = Variable, fill = Value)
  ) +
  scale_fill_gradient2(
    name = "Scaled Exp.",
    low = "#1EB980", mid = "white", high = "#6A1B9A",
    midpoint = 0, limits = c(-2,4)
  ) +
  facet_grid(Domain ~ ., scales="free_y", space="free_y") +
  labs(
    x = "Patients",
    y = "",
    title = "Patient-level relation between ALSFRS-R and TARDBP"
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
