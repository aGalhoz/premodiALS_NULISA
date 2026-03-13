# Differential expression analyses of the interleukin's group all together

###############################################
### Helper Functions
###############################################

# 1. Get Interleukins
extract_il_ligands <- function(df) {
  df %>% 
    filter(str_detect(Target, "^IL[0-9]")) %>% 
    filter(!str_detect(Target, "R[0-9]|BP|RB|RA|RN|ST|RAP"))
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

# 1. Check log2FC trends in the fluids for ALS vs CTR and PGMC vs CTR
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