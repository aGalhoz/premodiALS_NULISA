## ================================================================================
## PGMC selection & protein signature overlap: Venn diagrams + summary heatmap


## ------------------------------------------------------------------
## HELPER FUNCTIONS
## ------------------------------------------------------------------

## PGMC patients classified as "at risk" (ALS_risk_score > 0) 
get_selected_pgmc <- function(model_result) {
  model_result$results %>%
    left_join(protein_data_IDs %>% select(SampleName, ParticipantCode, type), by = "SampleName") %>%
    distinct() %>%
    filter(type == "PGMC", ALS_risk_score > 0) %>%
    pull(ParticipantCode)
}

## Proteins with nonzero coefficient in the locked model 
get_selected_proteins <- function(model_result, covariate_cols = c("age", "sex_male"),
                                  normalize_panel_suffix = TRUE) {
  coef_df <- model_result$coefficients
  
  value_col <- setdiff(names(coef_df), "feature")
  if (length(value_col) != 1) {
    stop("Expected exactly one non-'feature' column in model_result$coefficients, found: ",
         paste(value_col, collapse = ", "), ". Inspect names(model_result$coefficients) and adjust.")
  }
  
  proteins <- coef_df %>%
    filter(feature != "(Intercept)", !feature %in% covariate_cols, .data[[value_col]] != 0) %>%
    pull(feature)
  
  if (normalize_panel_suffix) {
    proteins <- gsub("_CNS$|_IMMUNE$", "", proteins)
  }
  proteins
}

## Venn diagram helpers
shorten_labels <- function(names_vec) {
  token_lists <- strsplit(names_vec, "_")
  max_len <- max(sapply(token_lists, length))
  token_lists <- lapply(token_lists, function(t) { length(t) <- max_len; t })
  token_matrix <- do.call(rbind, token_lists)
  
  varies <- apply(token_matrix, 2, function(col) length(unique(col)) > 1)
  if (!any(varies)) return(names_vec)  # nothing varies -- fall back to full names
  
  apply(token_matrix[, varies, drop = FALSE], 1, function(row) paste(na.omit(row), collapse = "_"))
}

## Compute every Venn region 
compute_venn_regions <- function(sets) {
  set_names <- names(sets)
  n <- length(set_names)
  region_list <- list()
  
  for (k in n:1) {   # largest combinations first (most informative) down to singles
    combos <- combn(set_names, k, simplify = FALSE)
    for (combo in combos) {
      in_these <- Reduce(intersect, sets[combo])
      other_sets <- setdiff(set_names, combo)
      exclusive <- if (length(other_sets) > 0) {
        setdiff(in_these, unique(unlist(sets[other_sets])))
      } else in_these
      region_list[[paste(combo, collapse = " & ")]] <- exclusive
    }
  }
  region_list
}

save_venn_regions <- function(sets, output_path) {
  regions <- compute_venn_regions(sets)
  lines <- unlist(lapply(names(regions), function(nm) {
    c(paste0("=== ", nm, " (n=", length(regions[[nm]]), ") ==="),
      if (length(regions[[nm]]) > 0) paste(sort(regions[[nm]]), collapse = ", ") else "(none)",
      "")
  }))
  writeLines(lines, output_path)
}

make_venn_patients <- function(model_names, title, plot_path = NULL) {
  sets <- lapply(model_registry[model_names], get_selected_pgmc)
  names(sets) <- shorten_labels(names(model_registry[model_names]))
  
  p <- ggVennDiagram(sets, label = "count") +
    scale_fill_gradient(low = "grey95", high = "#AD5291") +
    labs(title = title) +
    theme(legend.position = "none",
          plot.margin = margin(t = 1, r = 3, b = 1, l = 3, unit = "cm"))
  
  if (!is.null(plot_path)) {
    ggsave(plot_path, p, width = 9, height = 7)
    save_venn_regions(sets, gsub("\\.pdf$", "_ids.txt", plot_path))
  }
  p
}

make_venn_proteins <- function(model_names, title, plot_path = NULL, normalize_panel_suffix = TRUE) {
  sets <- lapply(model_registry[model_names], get_selected_proteins,
                 normalize_panel_suffix = normalize_panel_suffix)
  names(sets) <- shorten_labels(names(model_registry[model_names]))
  
  p <- ggVennDiagram(sets, label = "count") +
    scale_fill_gradient(low = "grey95", high = "#6F8EB2") +
    labs(title = title) +
    theme(legend.position = "none",
          plot.margin = margin(t = 1, r = 3, b = 1, l = 3, unit = "cm"))
  
  if (!is.null(plot_path)) {
    ggsave(plot_path, p, width = 9, height = 7)
    save_venn_regions(sets, gsub("\\.pdf$", "_proteins.txt", plot_path))
  }
  p
}


# selection of PGMCs
selection_long <- purrr::imap_dfr(model_registry, function(model_result, model_name) {
  selected <- get_selected_pgmc(model_result)
  tibble(ParticipantCode = all_pgmc_ids,
         model = model_name,
         selected = as.integer(all_pgmc_ids %in% selected))
})

selection_wide <- selection_long %>%
  pivot_wider(names_from = model, values_from = selected) %>%
  column_to_rownames("ParticipantCode")

# get data basd on fluid
model_fluid <- sapply(names(model_registry), function(nm) {
  if (grepl("SERUM", nm))  "SERUM"
  else if (grepl("PLASMA", nm)) "PLASMA"
  else if (grepl("CSF", nm))    "CSF"
  else NA_character_
})

## get visit data
get_visit_value <- function(colname, visit_label) {
  lookup <- clinical_data_extra_heatmap %>%
    filter(Visit == visit_label) %>%
    distinct(ParticipantCode, .data[[colname]])
  
  if (any(duplicated(lookup$ParticipantCode))) {
    warning("Multiple values of '", colname, "' at Visit ", visit_label,
            " for at least one patient -- check clinical_data_extra_heatmap for duplicate rows.")
  }
  deframe(lookup)
}

to_row_vec <- function(lookup) {
  setNames(lookup[rownames(selection_matrix)], rownames(selection_matrix))
}

## ------------------------------------------------------------------
## WORKFLOW
## ------------------------------------------------------------------

# Get models' outputs
model_registry <- list(
  Lasso_SERUM_both   = lasso_PGMC_serum_signature_both,
  Lasso_PLASMA_both  = lasso_PGMC_plasma_signature_both,
  EN_SERUM_both      = EN_PGMC_serum_signature_both,
  EN_PLASMA_both     = EN_PGMC_plasma_signature_both,
  EN_CSF_both        = EN_PGMC_CSF_signature_both,
  
  Lasso_SERUM_CNS    = lasso_PGMC_serum_signature_CNS,
  Lasso_PLASMA_CNS   = lasso_PGMC_plasma_signature_CNS,
  Lasso_CSF_CNS      = EN_PGMC_CSF_signature_CNS,
  EN_SERUM_CNS       = EN_PGMC_serum_signature_CNS,
  EN_PLASMA_CNS      = EN_PGMC_plasma_signature_CNS,
  EN_CSF_CNS         = EN_PGMC_CSF_signature_CNS,
  
  Lasso_SERUM_IMMUNE  = lasso_PGMC_serum_signature_IMMUNE,
  Lasso_PLASMA_IMMUNE = lasso_PGMC_plasma_signature_IMMUNE,
  EN_SERUM_IMMUNE     = EN_PGMC_serum_signature_IMMUNE,
  EN_PLASMA_IMMUNE    = EN_PGMC_plasma_signature_IMMUNE
)

## Venns per fluid across panel modes (both / CNS / Immune) separately for Lasso and Elastic Net

## Lasso: SERUM and PLASMA only
make_venn_patients(c("Lasso_SERUM_both","Lasso_SERUM_CNS","Lasso_SERUM_IMMUNE"),
                   "PGMC selected as at-risk: Lasso, SERUM",
                   "plots/venn/patients_lasso_SERUM.pdf")
make_venn_proteins(c("Lasso_SERUM_both","Lasso_SERUM_CNS","Lasso_SERUM_IMMUNE"),
                   "Signature overlap: Lasso, SERUM",
                   "plots/venn/proteins_lasso_SERUM.pdf")

make_venn_patients(c("Lasso_PLASMA_both","Lasso_PLASMA_CNS","Lasso_PLASMA_IMMUNE"),
                   "PGMC selected as at-risk: Lasso, PLASMA",
                   "plots/venn/patients_lasso_PLASMA.pdf")
make_venn_proteins(c("Lasso_PLASMA_both","Lasso_PLASMA_CNS","Lasso_PLASMA_IMMUNE"),
                   "Signature overlap: Lasso, PLASMA",
                   "plots/venn/proteins_lasso_PLASMA.pdf")

## Elastic Net: SERUM, PLASMA, CSF
for (fl in c("SERUM", "PLASMA", "CSF")) {
  models <- c(paste0("EN_", fl, "_both"), paste0("EN_", fl, "_CNS"), paste0("EN_", fl, "_IMMUNE"))
  make_venn_patients(models, paste0("PGMC selected as at-risk: Elastic Net, ", fl),
                     paste0("plots/venn/patients_EN_", fl, ".pdf"))
  make_venn_proteins(models, paste0("Signature overlap: Elastic Net, ", fl),
                     paste0("plots/venn/proteins_EN_", fl, ".pdf"))
}

## Venn set across fluids, within one panel mode 
for (pm in c("both", "CNS", "IMMUNE")) {
  models <- c(paste0("EN_SERUM_", pm), paste0("EN_PLASMA_", pm), paste0("EN_CSF_", pm))
  make_venn_patients(models, paste0("PGMC selected as at-risk across fluids (Elastic Net, ", pm, ")"),
                     paste0("plots/venn/patients_EN_across_fluids_", pm, ".pdf"))
  make_venn_proteins(models, paste0("Signature overlap across fluids (Elastic Net, ", pm, ")"),
                     paste0("plots/venn/proteins_EN_across_fluids_", pm, ".pdf"),
                     normalize_panel_suffix = FALSE)  
}

## Lasso across fluids is only a 2-set comparison (SERUM vs PLASMA, no CSF)
for (pm in c("both", "CNS", "IMMUNE")) {
  models <- c(paste0("Lasso_SERUM_", pm), paste0("Lasso_PLASMA_", pm))
  make_venn_patients(models, paste0("PGMC selected as at-risk, Lasso SERUM vs PLASMA (", pm, ")"),
                     paste0("plots/venn/patients_Lasso_across_fluids_", pm, ".pdf"))
  make_venn_proteins(models, paste0("Signature overlap, Lasso SERUM vs PLASMA (", pm, ")"),
                     paste0("plots/venn/proteins_Lasso_across_fluids_", pm, ".pdf"),
                     normalize_panel_suffix = FALSE)
}


## ================================================================================
## Summary heatmap for all fluids

clinical_data_extra_heatmap = read_excel("data input/all_participants_ID_visits.xlsx")

# gt PGMCs and data of fluids
all_pgmc_ids <- protein_data_IDs %>% filter(type == "PGMC") %>% pull(ParticipantCode) %>% unique()
selection_matrix <- as.matrix(selection_wide[, names(model_registry)])  
column_split_fluid <- factor(model_fluid[colnames(selection_matrix)],
                             levels = c("SERUM", "PLASMA", "CSF"))

visits <- c("V0", "V1", "V2", "V3")

# participant label
mutation_lookup <- clinical_data_extra_heatmap %>%
  filter(!is.na(mutation)) %>%
  distinct(ParticipantCode, mutation)
if (any(duplicated(mutation_lookup$ParticipantCode))) {
  mutation_lookup <- mutation_lookup %>% distinct(ParticipantCode, .keep_all = TRUE)
}
mutation_vec <- deframe(mutation_lookup)[rownames(selection_matrix)]
row_labels_vec <- ifelse(is.na(mutation_vec),
                         rownames(selection_matrix),
                         paste0(rownames(selection_matrix), " (", mutation_vec, ")"))

# age and motor signs
age_v0_vec <- to_row_vec(get_visit_value("age", "V0"))

motor_mat <- sapply(visits, function(v) to_row_vec(get_visit_value("TotalMotorSigns", v)))
rownames(motor_mat) <- rownames(selection_matrix)
colnames(motor_mat) <- visits

# event markers 
visit_index <- function(v) match(v, visits)

pgmc_visits <- clinical_data_extra_heatmap %>%
  filter(ParticipantCode %in% rownames(selection_matrix)) %>%
  mutate(vi = visit_index(Visit))

# PHENOCONVERSION
pheno_df <- pgmc_visits %>%
  filter(type == "PGMC", DiseaseOnset == "ALS", `Disease duration` > 0) %>%
  group_by(ParticipantCode) %>%
  slice_min(`Disease duration`, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ParticipantCode, vi)

# MOTOR SIGNS POSITIVE: first visit with motor signs > 0
motor_pos_df <- pgmc_visits %>%
  filter(!is.na(TotalMotorSigns), TotalMotorSigns > 0, !is.na(vi)) %>%
  group_by(ParticipantCode) %>%
  slice_min(vi, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ParticipantCode, vi) %>%
  anti_join(pheno_df, by = c("ParticipantCode", "vi"))

# MOTOR SIGNS NEGATIVE: motor signs > 0 at a visit and 0 at the next visit
motor_neg_df <- pgmc_visits %>%
  filter(!is.na(vi)) %>%
  arrange(ParticipantCode, vi) %>%
  group_by(ParticipantCode) %>%
  mutate(next_signs = lead(TotalMotorSigns), next_vi = lead(vi)) %>%
  ungroup() %>%
  filter(!is.na(TotalMotorSigns), TotalMotorSigns > 0,
         next_vi == vi + 1, !is.na(next_signs), next_signs == 0) %>%
  select(ParticipantCode, vi)

make_keys <- function(df) paste(df$ParticipantCode, df$vi)
pheno_keys <- make_keys(pheno_df)
pos_keys   <- make_keys(motor_pos_df)
neg_keys   <- make_keys(motor_neg_df)

# DEATH
stopifnot(exists("disease_duration"))
death_ids <- disease_duration %>%
  filter(type == "PGMC", !is.na(LFUDeathDate)) %>%
  distinct(ParticipantCode) %>%
  pull(ParticipantCode)
death_df <- pgmc_visits %>%
  filter(ParticipantCode %in% death_ids, !is.na(vi)) %>%
  group_by(ParticipantCode) %>%
  slice_max(vi, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ParticipantCode, vi)
death_keys <- make_keys(death_df)

# colours
mutation_levels <- unique(na.omit(mutation_vec))
mutation_levels <- c(intersect(c("C9orf72", "SOD1", "TARDBP", "FUS"), mutation_levels),
                     setdiff(mutation_levels, c("C9orf72", "SOD1", "TARDBP", "FUS", "NEK1")),
                     intersect("NEK1", mutation_levels))
mutation_palette <- c("C9orf72" = "#55AA82", "SOD1" = "#4661B9", "TARDBP" = "#B99E46",
                      "FUS" = "#B96546", "NEK1" = "#8E6BBF")
mutation_col <- mutation_palette[mutation_levels]
names(mutation_col) <- mutation_levels
extra <- is.na(mutation_col)
if (any(extra)) mutation_col[extra] <- hcl.colors(sum(extra), "Pastel 1")
binary_col <- c(`0` = "grey95", `1` = "#AD5291")
age_col    <- colorRamp2(range(age_v0_vec, na.rm = TRUE), c("grey95", "#7B4F9E"))
motor_col  <- colorRamp2(c(0, 8), c("white", "red"))

# legend for the markers
mk <- function(pch, fill) {
  function(x, y, w, h) grid.points(x, y, pch = pch, size = unit(2.5, "mm"),
                                   gp = gpar(col = "black", fill = fill, lwd = 0.8))
}
# cross for death (vertical bar with a shorter bar near the top)
draw_cross <- function(x, y, size = unit(2.8, "mm")) {
  gp <- gpar(col = "black", lwd = 1.6, lineend = "butt")
  grid.segments(x0 = x, x1 = x, y0 = y - 0.5 * size, y1 = y + 0.5 * size, gp = gp)
  grid.segments(x0 = x - 0.3 * size, x1 = x + 0.3 * size,
                y0 = y + 0.15 * size, y1 = y + 0.15 * size, gp = gp)
}

lgd_markers <- Legend(
  title = "Events",
  labels = c("Phenoconversion", "Motor signs positive", "Motor signs negative", "Death"),
  graphics = list(mk(21, "black"), mk(24, "white"), mk(25, "white"),
                  function(x, y, w, h) draw_cross(x, y)),
  background = "white"
)

# heatmap builder
n_models <- ncol(selection_matrix)             # models tried
freq_models <- rowSums(selection_matrix)       # in how many of the models each PGMC was selected

make_pgmc_summary_heatmap <- function(fluids, output_path, width = 17, height = 10, show_freq = FALSE) {
  
  cols_keep <- colnames(selection_matrix)[column_split_fluid %in% fluids]
  mat <- selection_matrix[, cols_keep, drop = FALSE]
  multi <- length(fluids) > 1
  
  # no clustering: PGMC selected most often first
  ord <- order(-rowSums(mat), -freq_models[rownames(mat)], rownames(mat))
  mat <- mat[ord, , drop = FALSE]
  ids <- rownames(mat)
  motor_ord <- motor_mat[ids, , drop = FALSE]
  freq_ord <- freq_models[ids]
  
  ht_main <- Heatmap(
    mat,
    name = "At-risk\n(score > 0)",
    col = binary_col,
    column_split = if (multi) droplevels(column_split_fluid[column_split_fluid %in% fluids]) else NULL,
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_labels = ids,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title = if (multi) as.character(fluids) else paste0("PGMC risk classification in ", fluids),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    # extra column: in how many of the models the PGMC was selected
    right_annotation = if (show_freq) rowAnnotation(
      freq = anno_barplot(freq_ord, ylim = c(0, n_models), width = unit(2, "cm"),
                          gp = gpar(fill = "#AD5291", col = NA),
                          axis_param = list(at = pretty(c(0, n_models), 3), gp = gpar(fontsize = 7))),
      n = anno_text(paste0(freq_ord, "/", n_models), gp = gpar(fontsize = 7)),
      annotation_label = c(paste0("Selected in models")),
      show_annotation_name = c(TRUE, FALSE),
      annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
    ) else NULL,
    heatmap_legend_param = list(
      title = "At-risk\n(score > 0)",
      at = c("0", "1"),
      labels = c("Not selected", "Selected")
    )
  )
  
  age_anno <- rowAnnotation(
    Mutation = factor(mutation_vec[ids], levels = mutation_levels),
    Age = age_v0_vec[ids],
    col = list(Mutation = mutation_col, Age = age_col),
    na_col = "white",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size = unit(0.35, "cm")
  )
  
  ht_motor <- Heatmap(
    motor_ord,
    name = "Total motor signs",
    col = motor_col,
    na_col = "grey90",
    left_annotation = age_anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_labels = ids,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = "Motor signs per visit",
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    width = unit(ncol(motor_ord) * 6, "mm"),
    heatmap_legend_param = list(at = 0:8),
    layer_fun = function(j, i, x, y, w, h, fill) {
      # tiles first, then markers on top so neighbouring tiles do not cover them
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.5))
      keys <- paste(rownames(motor_ord)[i], j)
      
      sel <- keys %in% pos_keys                       # motor signs positive: left edge of the visit
      if (any(sel)) grid.points(x[sel] - w[sel] / 2, y[sel], pch = 24, size = unit(1.8, "mm"),
                                gp = gpar(fill = "white", col = "black", lwd = 0.8))
      sel <- keys %in% neg_keys                       # motor signs negative: boundary to the next visit
      if (any(sel)) grid.points(x[sel] + w[sel] / 2, y[sel], pch = 25, size = unit(1.8, "mm"),
                                gp = gpar(fill = "white", col = "black", lwd = 0.8))
      sel <- keys %in% pheno_keys                     # phenoconversion: left edge of the visit
      if (any(sel)) grid.points(x[sel] - w[sel] / 2, y[sel], pch = 21, size = unit(2.2, "mm"),
                                gp = gpar(fill = "black", col = "black", lwd = 0.5))
      sel <- keys %in% death_keys                     # death: right edge of the last visit
      if (any(sel)) draw_cross(x[sel] + w[sel] / 2, y[sel])
    }
  )
  
  pdf(output_path, width = width, height = height)
  draw(ht_main + ht_motor,
       annotation_legend_list = list(lgd_markers),
       merge_legends = TRUE)
  dev.off()
  
  invisible(ht_main + ht_motor)
}

# all fluids together and per fluid
ht_all    <- make_pgmc_summary_heatmap(c("SERUM", "PLASMA", "CSF"),
                                       "plots/venn/PGMC_summary_heatmap_all_fluids.pdf", width = 17,
                                       show_freq = TRUE)
ht_serum  <- make_pgmc_summary_heatmap("SERUM",  "plots/venn/PGMC_summary_heatmap_SERUM.pdf",  width = 10)
ht_plasma <- make_pgmc_summary_heatmap("PLASMA", "plots/venn/PGMC_summary_heatmap_PLASMA.pdf", width = 10)
ht_csf    <- make_pgmc_summary_heatmap("CSF",    "plots/venn/PGMC_summary_heatmap_CSF.pdf",    width = 10)
