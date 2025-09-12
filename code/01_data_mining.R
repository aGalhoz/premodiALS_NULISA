source("00_initialization.R")

## Data mining of all proteomics data by groupps

##### Diseased groups CTR, ALS, PGMC and mimic
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
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`))

# how many samples per biofluid
nr_samples_serum <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  pull(SampleName) %>% unique() %>% length()

samples_serum <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "SERUM") %>% 
  select(SampleName,type) %>% 
  distinct() %>%
  group_by(type) %>% 
  left_join(samples_ID_type %>% rename(SampleName = `Sample ID`)) %>%
  summarise(nr_samples = n())

nr_samples_plasma <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  pull(SampleName) %>% unique() %>% length()

samples_plasma <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "PLASMA") %>% 
  select(SampleName,type) %>% 
  distinct() %>%
  group_by(type) %>%
  summarise(nr_samples = n())

nr_samples_CSF <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "CSF") %>% 
  pull(SampleName) %>% unique() %>% length()

samples_CSF <- protein_data_IDs[!is.na(protein_data_IDs$type),] %>% 
  distinct() %>% na.omit() %>% 
  filter(SampleMatrixType == "CSF") %>% 
  select(SampleName,type) %>% 
  distinct() %>%
  group_by(type) %>%
  summarise(nr_samples = n())


# amount of unique targets
target_unique <- protein_data_IDs$Target %>% unique() #127 targets

##### PGMC mutations and controls
samples_PGMC_CTR_ID_type = sample_ID_info %>%
  rename(ParticipantCode = `Optional Informarion (patient ID)`) %>%
  left_join(do.call("rbind", list(CTR_ID, PGMC_mutations_ID))) 

# get protein data with status information
protein_data_PGMC_CTR_IDs <- protein_data %>%
  filter(SampleType == "Sample") %>%
  select(SampleName,SampleMatrixType,Target,UniProtID,ProteinName,NPQ) %>%
  left_join(samples_PGMC_CTR_ID_type %>% select(`Sample ID`,type) %>% rename(SampleName = `Sample ID`))

##### Correlation plots between fluids
protein_data_IDs_SERUM <- protein_data %>%
  filter(SampleMatrixType == "SERUM") %>%
  select(NPQ,Target,UniProtID) %>%
  distinct() %>%
  group_by(Target,UniProtID) %>%
  summarise(NPQ_Serum = mean(NPQ)) %>%
  distinct()

protein_data_IDs_CSF <- protein_data %>%
  filter(SampleMatrixType == "CSF") %>%
  select(NPQ,Target,UniProtID) %>%
  distinct() %>%
  group_by(Target,UniProtID) %>%
  summarise(NPQ_CSF = mean(NPQ)) %>%
  distinct()

protein_data_IDs_PLASMA <- protein_data %>%
  filter(SampleMatrixType == "PLASMA") %>%
  select(NPQ,Target,UniProtID) %>%
  distinct() %>%
  group_by(Target,UniProtID) %>%
  summarise(NPQ_PLASMA = mean(NPQ)) %>%
  distinct()

# function to create correlation plot
corr_plot <- function(data, x_col, y_col, title, reg_color, hist_color){
  
  # Compute Pearson correlation
  cor_test <- cor.test(data[[x_col]], data[[y_col]], use = "complete.obs")
  r_val <- round(cor_test$estimate, 2)
  p_val <- ifelse(cor_test$p.value < 0.001,
                  formatC(cor_test$p.value, format = "e", digits = 2),
                  round(cor_test$p.value, 3))
  cor_text <- paste0("italic(R) == ", r_val, " * ',' ~ italic(p) == ", p_val)
  
  base_plot <- ggscatter(
    data, x = x_col, y = y_col,
    add = "reg.line",
    add.params = list(color = reg_color, fill = hist_color),
    conf.int = TRUE
  ) +
    annotate(geom = "text",
      x = min(data[[x_col]], na.rm = TRUE), 
          y = max(data[[y_col]], na.rm = TRUE) + 1, 
          label = cor_text,
      size = 5, hjust = 0, parse = TRUE
    ) +
    labs(title = title, x = x_col, y = y_col) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    ) 
    
  
  # Add marginal histogram
  ggMarginal(base_plot, type = "histogram", color = hist_color, fill = hist_color)
}

# -> Serum vs CSF
protein_data_IDS_SERUM_CSF = merge(protein_data_IDs_SERUM,
                                   protein_data_IDs_CSF,
                                   by = c("Target","UniProtID"))
p1 <- corr_plot(protein_data_IDS_SERUM_CSF, 
                       "NPQ_Serum", "NPQ_CSF", "Correlation Serum vs CSF", 
                       "blue", "lightblue")

pdf("plots/correlation_serum_CSF.pdf")
p1
dev.off()

# -> Plasma vs Serum
protein_data_IDS_PLASMA_SERUM = merge(protein_data_IDs_PLASMA,
                                      protein_data_IDs_SERUM,
                                      by = c("Target","UniProtID"))
p2 <- corr_plot(protein_data_IDS_PLASMA_SERUM, 
                       "NPQ_PLASMA", "NPQ_Serum", "Correlation Plasma vs Serum", 
                       "blue", "lightblue")
pdf("plots/correlation_plasma_serum.pdf")
p2
dev.off()

# -> Plasma vs CSF
protein_data_IDS_PLASMA_CSF = merge(protein_data_IDs_PLASMA,
                                      protein_data_IDs_CSF,
                                      by = c("Target","UniProtID"))
p3 <- corr_plot(protein_data_IDS_PLASMA_CSF, 
                       "NPQ_PLASMA", "NPQ_CSF", "Correlation Plasma vs CSF", 
                       "blue", "lightblue")
pdf("plots/correlation_plasma_CSF.pdf")
p3
dev.off()

# -> Correlation plots all together
combined <- plot_grid(p2, p1, p3, nrow = 1, align = "h")

pdf("plots/correlation_combined.pdf", width = 18, height = 6)
print(combined)
dev.off()

# check differences in NPQ values depending on the plate
unique_targets <- protein_data$Target %>% unique()

for (i in 1:length(unique_targets)) {
  target = unique_targets[i]
  protein_data_target <- protein_data %>%
    filter(Target == target)
  protein_data_target$SampleMatrixType <- factor(protein_data_target$SampleMatrixType,
                                                 levels = c("PLASMA","SERUM","CSF","CONTROL"))
  p <- ggboxplot(protein_data_target,x = "SampleMatrixType",y = "NPQ",
                 add = "jitter",color = "PlateID")
  pdf(paste0("plots/NPQ_fluid_plate/NPQ_",unique_targets[i],".pdf"))
  print(p)
  dev.off() 
}

# check differences in target detectability per target & plate
target_detectability_plate_fluid <- protein_data %>% 
  select(PlateID,SampleMatrixType,Target,NPQ) %>%
  left_join(target_detectability %>%
              rename(Target = TargetName)) %>%
  mutate(TargetDetectability_value = NPQ - TargetLOD_NPQ,
         TargetDetectability = as.numeric(TargetDetectability) * 100)

target_detectability_plate_fluid$TargetDetectability_value <- ifelse(target_detectability_plate_fluid$Target %in% c("APOE","CRP"),
                                                                     abs(target_detectability_plate_fluid$TargetDetectability_value),
                                                                     target_detectability_plate_fluid$TargetDetectability_value)

fluids = c("CSF","PLASMA","SERUM")

for (i in 1:length(fluids)) {
  fluid = fluids[i]
  target_detectability_plate_fluid_i <- target_detectability_plate_fluid %>%
    filter(SampleMatrixType == fluid) %>%
    group_by(Target) %>%
    mutate(mean_value = mean(TargetDetectability_value, na.rm = T),
           mean_value_percentage = mean(TargetDetectability,na.rm = T)) %>%
    ungroup() %>%
    mutate(Target = factor(Target, 
                           levels = unique(Target[order(mean_value,decreasing = T)])))
  p <- ggplot(target_detectability_plate_fluid_i,
              aes(x = Target, y = TargetDetectability_value, 
                  fill = mean_value_percentage, group = Target)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_gradient(low = "blue", high = "orange") +
    labs(
      title = paste("Target Detectability -", fluid),
      x = "Target",
      y = "Detectability (NPQ - LOD)",
      fill = "Detectability (%)",
      color = "Detectability (%)"
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Arial Unicode MS"),
      axis.text.x = element_text(angle = 60, hjust = 1),
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15)
    )
  Cairo::CairoPDF(paste0("plots/TargetDetectability",fluid,".pdf"),width = 28,height = 6)
  print(p)
  dev.off() 
}

# stratified by plate
for (i in 1:length(fluids)) {
  fluid <- fluids[i]
  
  # Filter data for this fluid
  target_detectability_plate_fluid_i <- target_detectability_plate_fluid %>%
    filter(SampleMatrixType == fluid) %>%
    group_by(PlateID, Target) %>%
    mutate(mean_value = mean(TargetDetectability_value, na.rm = TRUE),
           mean_value_percentage = mean(TargetDetectability,na.rm = T)) %>%
    ungroup() %>%
    mutate(Target = factor(Target, 
                           levels = unique(Target[order(mean_value, decreasing = TRUE)]))) 
  
  # Plot stratified by PlateID
  p <- ggplot(target_detectability_plate_fluid_i,
              aes(x = Target, y = TargetDetectability_value, 
                  fill = mean_value_percentage)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~PlateID, ncol = 1) +  # stratify by plate
    scale_fill_gradient(low = "blue", high = "orange") +
    labs(
      title = paste("Target Detectability -", fluid),
      x = "Target",
      y = "Detectability (NPQ - LOD)",
      fill = "Detectability (%)"
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Arial Unicode MS"),
      axis.text.x = element_text(angle = 60, hjust = 1),
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15)
    )
  
  Cairo::CairoPDF(paste0("plots/TargetDetectability_byPlate_", fluid, ".pdf"), width = 28,height = 8)
  print(p)
  dev.off()
}
