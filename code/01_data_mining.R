source("00_initialization.R")

###############################################
### Helper Functions
###############################################

## 1. Assign group labels (ALS / CTR / PGMC / mimic)
assign_group_labels <- function(sample_df, ...) {
  bind_rows(...) %>%
    rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
    left_join(sample_df, by = "ParticipantCode")
}

## 2. Collect sample counts per fluid
count_samples_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(type)) %>%
    distinct(SampleName, type) %>%
    group_by(type) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}

count_subtypes_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid, !is.na(subtype)) %>%
    distinct(SampleName, subtype) %>%
    group_by(subtype) %>%
    summarise(nr_samples = n(), .groups = "drop") %>%
    mutate(biofluid = fluid)
}


## 3. Get mean per fluid
get_mean_per_fluid <- function(df, fluid) {
  df %>%
    filter(SampleMatrixType == fluid) %>%
    distinct(Target, UniProtID, NPQ) %>%
    group_by(Target, UniProtID) %>%
    summarise(NPQ = mean(NPQ), .groups = "drop") %>%
    rename_with(~ paste0("NPQ_", fluid), "NPQ")
}

## 4. Correlation plot generator
corr_plot <- function(data, x, y, title) {
  
  cor_test <- cor.test(data[[x]], data[[y]], use = "complete.obs")
  r_val <- round(cor_test$estimate, 2)
  p_val <- ifelse(cor_test$p.value < 0.001,
                  formatC(cor_test$p.value, format = "e", digits = 2),
                  round(cor_test$p.value, 3))
  
  cor_text <- paste0("italic(R) == ", r_val,
                     " * ',' ~ italic(p) == ", p_val)
  
  p <- ggscatter(data, x = x, y = y,
                 add = "reg.line",
                 add.params = list(color = "blue",fill = "lightblue"),
                 conf.int = TRUE) +
    annotate("text",
             x = min(data[[x]], na.rm = TRUE),
             y = max(data[[y]], na.rm = TRUE),
             hjust = 0,
             label = cor_text,
             parse = TRUE, size = 5) +
    labs(title = title, x = x, y = y) +
    theme( plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
           axis.title = element_text(size = 16), 
           axis.text = element_text(size = 14), 
           legend.title = element_text(size = 14), 
           legend.text = element_text(size = 14))
  
  ggMarginal(p, type = "histogram", color = "lightblue", fill = "lightblue")
}

make_corr <- function(df1, df2, label1, label2, filename) {
  merged <- inner_join(df1, df2, by = c("Target", "UniProtID"))
  p <- corr_plot(merged,
                 paste0("NPQ_", label1),
                 paste0("NPQ_", label2),
                 paste("Correlation", label1, "vs", label2))
  
  pdf(filename)
  print(p)
  dev.off()
  p
}

## 5. Detectability by fluid and plate
plot_detectability <- function(df, fluid, by_plate = FALSE) {
  
  df_f <- df %>% filter(SampleMatrixType == fluid)
  
  if (by_plate) {
    df_f <- df_f %>%
      group_by(PlateID, Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  } else {
    df_f <- df_f %>%
      group_by(Target) %>%
      mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
             mean_pct   = mean(TargetDetectability, na.rm = TRUE)) %>%
      ungroup()
  }
  
  df_f <- df_f %>%
    mutate(Target = factor(Target,
                           levels = unique(Target[order(mean_value, decreasing = TRUE)])
    ))
  
  p <- ggplot(df_f,
              aes(x = Target, y = TargetDetectability_value,
                  fill = mean_pct)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_gradient(low = "blue", high = "orange") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Target Detectability -", fluid),
      x = "Target",
      y = "Detectability (NPQ - LOD)",
      fill = "Detectability (%)"
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  if (by_plate) p <- p + facet_wrap(~PlateID, ncol = 1)
  
  p
}

###############################################
### 1. Assign Group Labels
###############################################

samples_ID_type <- sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(bind_rows(ALS_ID, CTR_ID, PGMC_ID, mimic_ID), by = "ParticipantCode")

samples_PGMC_CTR_ID_type = sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(bind_rows(CTR_ID,PGMC_mutations_ID), by = "ParticipantCode")

###############################################
### 2. Missing ID mapping
###############################################

samples_ID_NA <- samples_ID_type %>%
  filter(is.na(type)) %>%
  left_join(
    GeneralDocumentation %>%
      select(ParticipantCode, PGMC, ALSuncertainty, ALSFUdiagnosis, LFU) %>%
      distinct(),
    by = "ParticipantCode"
  )

writexl::write_xlsx(samples_ID_NA, "results/samples_ID_NA.xlsx")

###############################################
### 3. Merge protein data with sample type
###############################################

protein_data_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName, SampleMatrixType, Target, UniProtID, ProteinName, NPQ) %>%
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`), by = "SampleName")

# get count table of APOE
table_APOE = protein_data_IDs %>%
  left_join(samples_PGMC_CTR_ID_type %>% rename(subtype = type,
                                                SampleName = `Sample ID`)) %>%
  mutate(
    subtype = ifelse(subtype == "CTR" | is.na(subtype), "No subtype",
                     ifelse(subtype %in% c("FUS","UBQLN2","FIG4","other"),"others",subtype))) %>%
  filter(Target == "APOE4") %>%
  mutate(APOE_status =  ifelse(NPQ > 10, "Carrier",ifelse(NPQ <= 10, "Non-carrier",NA)))

carrier_table <- table_APOE %>%
  group_by(SampleMatrixType, type, subtype, APOE_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = APOE_status,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = Carrier + `Non-carrier`
  )
writexl::write_xlsx(carrier_table,"results/APOE4_carriers.xlsx")

###############################################
### 4. Sample Counts Across Fluids
###############################################

fluids <- c("SERUM", "PLASMA", "CSF")

sample_counts <- bind_rows(
  lapply(fluids, function(f) count_samples_per_fluid(protein_data_IDs, f))
)

writexl::write_xlsx(sample_counts, "results/samples_biofluid_overview.xlsx")

###############################################
### 5. Correlation Between Fluids
###############################################

df_SERUM  <- get_mean_per_fluid(protein_data, "SERUM")
df_PLASMA <- get_mean_per_fluid(protein_data, "PLASMA")
df_CSF    <- get_mean_per_fluid(protein_data, "CSF")

p1 = make_corr(df_PLASMA, df_SERUM,  "PLASMA", "SERUM", "plots/correlation_plasma_serum.pdf")
p2 = make_corr(df_SERUM,  df_CSF,    "SERUM",  "CSF",   "plots/correlation_serum_CSF.pdf")
p3 = make_corr(df_PLASMA, df_CSF,    "PLASMA", "CSF",   "plots/correlation_plasma_CSF.pdf")

# -> Correlation plots all together 
combined <- plot_grid(p2, p1, p3, nrow = 1, align = "h") 
pdf("plots/correlation_combined.pdf", width = 18, height = 6) 
print(combined) 
dev.off()

###############################################
### 6. NPQ Distributions by Plate
###############################################

unique_targets <- unique(protein_data$Target)

for (target in unique_targets) {
  
  df_t <- protein_data %>%
    filter(Target == target) %>%
    mutate(SampleMatrixType = factor(SampleMatrixType,
                                     levels = c("PLASMA", "SERUM", "CSF", "CONTROL")))
  
  p <- ggboxplot(df_t,
                 x = "SampleMatrixType", y = "NPQ",
                 add = "jitter", color = "PlateID") +
    theme_bw()
  
  pdf(paste0("plots/NPQ_fluid_plate/NPQ_", target, ".pdf"))
  print(p)
  dev.off()
}

###############################################
### 7. Target Detectability per Fluid & Plate
###############################################

td <- protein_data  %>%
  select(PlateID, SampleMatrixType, Target, NPQ) %>%
  left_join(target_detectability %>% rename(Target = TargetName), by = c("Target","PlateID")) %>%
  mutate(TargetDetectability_value = NPQ - TargetLOD_NPQ,
         TargetDetectability =
           as.numeric(TargetDetectability) * 100) %>%
  mutate(TargetDetectability_value =
           ifelse(Target %in% c("APOE", "CRP"),
                  abs(TargetDetectability_value),
                  TargetDetectability_value))

fluids <- c("CSF", "PLASMA", "SERUM")

# Save plots
for (fluid in fluids) {
  p1 <- plot_detectability(td, fluid, by_plate = FALSE)
  p2 <- plot_detectability(td, fluid, by_plate = TRUE)
  
  Cairo::CairoPDF(paste0("plots/TargetDetectability_", fluid, ".pdf"),
                  width = 28, height = 6)
  print(p1)
  dev.off()
  
  Cairo::CairoPDF(paste0("plots/TargetDetectability_byPlate_", fluid, ".pdf"),
                  width = 28, height = 8)
  print(p2)
  dev.off()
}