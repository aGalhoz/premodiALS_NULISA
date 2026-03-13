# Analyses of relevant biomarkers with clinical variables

###############################################################################
# Run pipeline
###############################################################################

# 1. Define the interesting protein lists
targets_cns <- c("NEFH", "NEFL", "MAPT", "GDNF", "FABP3", "TAFA5", "pTau181", "pTau231", "pTau217")
targets_immune <- c("TAFA5", "CX3CL1", "CEACAM5", "VEGFD", "IL7", "IL27", "IL33")

current_targets <- targets_immune

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
  if(!dir.exists("plots/site_onset")) dir.create("plots/site_onset", recursive = TRUE)
  file_name <- paste0("plots/site_onset/Bulbar_Spinal_", fluid_name, ".pdf")
  
  ggsave(file_name, 
         plot = p, 
         width = 10, 
         height = 10, 
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
      panel.grid.major.y = element_line(color = "grey95") # Very subtle horizontal guide
    )
  
  # save plots
  if(!dir.exists("plots/correlations")) dir.create("plots/correlations", recursive = TRUE)
  file_name <- paste0("plots/correlations/ALSFRS_Corr_", fluid_name, ".pdf")
  
  ggsave(file_name, 
         plot = p, 
         width = 12, 
         height = 10, 
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
         height = 10, 
         units = "in", 
         device = cairo_pdf)
  
  message(paste("Saved:", file_name))
}

