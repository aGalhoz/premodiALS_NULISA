# Analyses of relevant biomarkers with clinical variables

###############################################################################
# Run pipeline
###############################################################################

# 1. Define the interesting protein lists
targets_cns <- c("NEFH", "NEFL", "MAPT", "GDNF", "FABP3", "TAFA5", "pTau-181",
                 "pTau-231", "pTau-217","HTT","TARDBP")
targets_immune <- c("TAFA5", "CX3CL1", "CEACAM5", "VEGFD", "IL7", "IL27", "IL33")

current_targets <- targets_cns

# 2. Loop over fluids to check site of onset stratification
comparisons <- list( c("spinal", "bulbar"))

for (fluid_name in names(results_ALL)) {
  
  message(paste("Processing fluid:", fluid_name))
  
  fluid_data <- results_ALL[[fluid_name]]$data_adjusted
  
  plot_data <- inner_join(fluid_data, site_onset) %>%
    mutate(SiteDiseaseOnset = tolower(trimws(as.character(SiteDiseaseOnset)))) %>% 
    filter(Target %in% current_targets) %>%
    filter(SiteDiseaseOnset %in% c("spinal", "bulbar")) %>% 
    filter(!is.na(NPQ_adj))
  
  if (nrow(plot_data) == 0) {
    message(paste("No matching 'spinal' or 'bulbar' data found in", fluid_name))
    next
  }
  
  # Plot
  p <- ggplot(plot_data, aes(x = SiteDiseaseOnset, y = NPQ_adj, fill = SiteDiseaseOnset)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
    facet_wrap(~Target, scales = "free", ncol = 3) +
    stat_compare_means(comparisons = comparisons, 
                       method = "wilcox.test", 
                       label = "p.format",      
                       size = 3.5,              
                       label.y.npc = 0.9) +        
    theme_classic(base_size = 12) +
    scale_fill_manual(values = c("spinal" = "#3498db", "bulbar" = "#e74c3c")) +
    labs(
      title = paste(fluid_name, ": Bulbar vs Spinal Segregation"),
      x = "Site of Onset",
      y = "Adjusted NPQ Expression"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  # save plots
  file_name <- paste0("plots/correlations/Bulbar_Spinal_", fluid_name, ".pdf")
  
  ggsave(file_name, 
         plot = p, 
         width = 10, 
         height = 14, 
         units = "in", 
         device = cairo_pdf)
  
  message(paste("Saved:", file_name))
}

# 3. Loop over fluids to check correlation with ALSFRS-R
for (fluid_name in names(results_ALL)) {
  
  message(paste("Processing fluid:", fluid_name))
  
  fluid_data <- results_ALL[[fluid_name]]$data_adjusted
  
  # Merge with clinical data 
  plot_data <- inner_join(fluid_data, ALSFRS_NULISA_V0) %>%
    filter(Target %in% current_targets) %>%
    filter(type == "ALS") %>%
    filter(!is.na(ALSFRS_1) & !is.na(NPQ_adj))
  
  if (nrow(plot_data) == 0) {
    message(paste("No data for correlation in", fluid_name))
    next
  }
  
  # correlation Plot
  p <- ggplot(plot_data, aes(x = ALSFRS_1, y = NPQ_adj)) +
    geom_point(color = "grey40", alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", color = "#001f3f", fill = "#001f3f", 
                alpha = 0.1, linewidth = 1) +
    facet_wrap(~Target, scales = "free_y", ncol = 3) +
    stat_cor(method = "spearman", 
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top", size = 3.2) + 
    theme_bw(base_size = 12) + 
    labs(
      title = paste(fluid_name, ": ALSFRS-R Correlation"),
      x = "ALSFRS-R Total Score",
      y = "Adjusted NPQ Expression"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey95") 
    )
  
  # save plots
  if(!dir.exists("plots/correlations")) dir.create("plots/correlations", recursive = TRUE)
  file_name <- paste0("plots/correlations/ALSFRS_Corr_", fluid_name, ".pdf")
  
  ggsave(file_name, 
         plot = p, 
         width = 12, 
         height = 12, 
         units = "in", 
         device = cairo_pdf)
  
  message(paste("Saved:", file_name))
}

# 4. Loop over fluids to check correlation with disease duration
for (fluid_name in names(results_ALL)) {
  
  message(paste("Processing fluid:", fluid_name))
  
  fluid_data <- results_ALL[[fluid_name]]$data_adjusted
  
  # Merge with clinical data 
  plot_data <- inner_join(fluid_data, disease_duration) %>%
    filter(Target %in% current_targets) %>%
    filter(type == "ALS") %>%
    filter(Visit == "V0") %>%
    filter(!is.na(`Disease duration`) & !is.na(NPQ_adj))
  
  if (nrow(plot_data) == 0) {
    message(paste("No data for correlation in", fluid_name))
    next
  }
  
  # correlation Plot
  p <- ggplot(plot_data, aes(x = `Disease duration`, y = NPQ_adj)) +
    geom_point(color = "grey40", alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", color = "#001f3f", fill = "#001f3f", 
                alpha = 0.1, linewidth = 1) +
    facet_wrap(~Target, scales = "free_y", ncol = 3) +
    stat_cor(method = "spearman", 
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top", size = 3.2) + 
    theme_bw(base_size = 12) + 
    labs(
      title = paste(fluid_name, ": Disease Duration Correlation"),
      x = "Disease Duration (months)",
      y = "Adjusted NPQ Expression"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(color = "grey95") 
    )
  
  # save plots
  if(!dir.exists("plots/correlations")) dir.create("plots/correlations", recursive = TRUE)
  file_name <- paste0("plots/correlations/Disease_Duration_Corr_", fluid_name, ".pdf")
  
  ggsave(file_name, 
         plot = p, 
         width = 12, 
         height = 12, 
         units = "in", 
         device = cairo_pdf)
  
  message(paste("Saved:", file_name))
}

# 5. Loop over fluids to check correlation with age
group_modes <- list(
  "ALS" = c("ALS"),
  "Alldiseased" = c("ALS", "CTR", "PGMC", "mimic")
)

for (mode_name in names(group_modes)) {
  
  groups_to_keep <- group_modes[[mode_name]]
  
  for (fluid_name in names(results_ALL)) {
    
    message(paste("Processing:", fluid_name, "-", mode_name))
    
    fluid_data <- results_ALL[[fluid_name]]$data
    
    plot_data <- inner_join(fluid_data, Sex_age_all_participants) %>%  
      filter(SampleMatrixType == fluid_name) %>%
      filter(Target %in% current_targets) %>%
      filter(type %in% groups_to_keep) %>%
      filter(!is.na(age) & !is.na(NPQ)) %>%
      select(age,NPQ,ParticipantCode,PatientID,Target,SampleName) %>%
      distinct()
    
    if (nrow(plot_data) == 0) next
    
    print(plot_data)
    
    p <- ggplot(plot_data, aes(x = age, y = NPQ)) +
      geom_point(color = "grey40", alpha = 0.5, size = 1.5) +
      geom_smooth(method = "lm", color = "#001f3f", fill = "#001f3f",
                  alpha = 0.1, linewidth = 1) +
      facet_wrap(~Target, scales = "free_y", ncol = 3) +
      stat_cor(method = "spearman",
               aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
               label.x.npc = "left", label.y.npc = "top", size = 3.2) +
      theme_bw(base_size = 12) +
      labs(
        title = paste(fluid_name, ": Age Correlation for", ifelse(mode_name == "Alldiseased",
                                                                  "all groups",
                                                                  mode_name)),
        x = "Age",
        y = "NPQ Expression"
      ) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major.y = element_line(color = "grey95") 
      )
    
    file_name <- paste0("plots/correlations/Age_Corr_", mode_name, "_", fluid_name, ".pdf")
    
    ggsave(file_name, plot = p, width = 12, height = 12, device = cairo_pdf)
    
    message("Saved: ", file_name)
  }
}

# 6. Loop over fluids to check stratification by gender
comparisons <- list(c("male", "female"))

for (mode_name in names(group_modes)) {
  
  groups_to_keep <- group_modes[[mode_name]]
  
  for (fluid_name in names(results_ALL)) {
    
    message(paste("Processing:", fluid_name, "-", mode_name))
    
    fluid_data <- results_ALL[[fluid_name]]$data
    
    plot_data <- inner_join(fluid_data, Sex_age_all_participants) %>%
      filter(SampleMatrixType == fluid_name) %>%
      filter(Target %in% current_targets) %>%
      filter(type %in% groups_to_keep) %>%
      filter(sex %in% c("M", "F")) %>%
      filter(!is.na(NPQ)) %>%
      mutate(Sex = recode(sex, "M"="male", "F"="female")) %>%
      select(Target,Sex,NPQ,ParticipantCode,PatientID,SampleName) %>%
      distinct()
    
    if (nrow(plot_data) == 0) next
    
    p <- ggplot(plot_data, aes(x = Sex, y = NPQ, fill = Sex)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
      facet_wrap(~Target, scales = "free_y", ncol = 3) +
      stat_compare_means(comparisons = comparisons,
                         method = "wilcox.test",
                         label = "p.format",
                         size = 3.5,
                         label.y.npc = 0.9) +
      theme_classic(base_size = 12) +
      scale_fill_manual(values = c("male" = "#4c72b0", "female" = "#dd8452")) +
      labs(
        title = paste(fluid_name, ": Sex Stratification for", ifelse(mode_name == "Alldiseased",
                                                                     "all groups",
                                                                     mode_name)),
        x = "Sex",
        y = "NPQ Expression"
      ) +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
      )
    
    file_name <- paste0("plots/correlations/Sex_", mode_name, "_", fluid_name, ".pdf")
    
    ggsave(file_name, plot = p, width = 10, height = 14, device = cairo_pdf)
    
    message("Saved: ", file_name)
  }
}