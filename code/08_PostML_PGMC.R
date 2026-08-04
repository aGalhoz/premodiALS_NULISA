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
  proteins <- model_result$coefficients %>%
    filter(feature != "(Intercept)", !feature %in% covariate_cols, s1 != 0) %>%
    pull(feature)
  
  if (normalize_panel_suffix) {
    proteins <- gsub("_CNS$|_IMMUNE$", "", proteins)
  }
  proteins
}

## Venn diagram helpers
make_venn_patients <- function(model_names, title, plot_path = NULL) {
  sets <- lapply(model_registry[model_names], get_selected_pgmc)
  names(sets) <- names(model_registry[model_names])
  
  p <- ggVennDiagram(sets, label = "count") +
    scale_fill_gradient(low = "grey95", high = "#AD5291") +
    labs(title = title) +
    theme(legend.position = "none")
  
  if (!is.null(plot_path)) ggsave(plot_path, p, width = 7, height = 6)
  p
}

make_venn_proteins <- function(model_names, title, plot_path = NULL, normalize_panel_suffix = TRUE) {
  sets <- lapply(model_registry[model_names], get_selected_proteins,
                 normalize_panel_suffix = normalize_panel_suffix)
  names(sets) <- names(model_registry[model_names])
  
  p <- ggVennDiagram(sets, label = "count") +
    scale_fill_gradient(low = "grey95", high = "#6F8EB2") +
    labs(title = title) +
    theme(legend.position = "none")
  
  if (!is.null(plot_path)) ggsave(plot_path, p, width = 7, height = 6)
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

## heatmap by fluid
make_fluid_summary_heatmap <- function(fluid_name, output_path) {
  
  cols_this_fluid <- colnames(selection_matrix)[column_split_fluid == fluid_name]
  mat_fluid <- selection_matrix[, cols_this_fluid, drop = FALSE]
  
  row_anno_fluid <- rowAnnotation(
    Converted              = converted_vec,
    MotorSigns_V0          = motor_v0_vec,
    MotorSigns_V1          = motor_v1_vec,
    MotorSigns_V2          = motor_v2_vec,
    TimeToPhenoconv        = phenoconversion_vec,
    NEFL                   = nefl_by_fluid[[fluid_name]],
    Urine_Neopterin_V0     = urine_neopterin_vec,
    Urine_p75ECD_V0        = urine_p75_vec,
    Urine_Neopterin_Delta  = urine_neopterin_delta_vec,
    Urine_p75ECD_Delta     = urine_p75_delta_vec,
    col = list(
      Converted             = c(`0` = "grey95", `1` = "#B2242A"),
      MotorSigns_V0         = motor_col,
      MotorSigns_V1         = motor_col,
      MotorSigns_V2         = motor_col,
      TimeToPhenoconv       = phenoconv_col,
      NEFL                  = nefl_col,
      Urine_Neopterin_V0    = urine_neo_col,
      Urine_p75ECD_VO       = urine_p75_col,
      Urine_Neopterin_Delta = urine_neo_delta_col,
      Urine_p75ECD_Delta    = urine_p75_delta_col
    ),
    na_col = "white",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size = unit(0.35, "cm")
  )
  
  ht_fluid <- Heatmap(
    mat_fluid,
    name = "At-risk\n(score > 0)",
    col = binary_col,
    right_annotation = row_anno_fluid,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    column_title = paste0("PGMC risk classification in ", fluid_name),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title = "At-risk\n(score > 0)",
      at = c("0", "1"),
      labels = c("Not selected", "Selected")
    )
  )
  
  pdf(output_path, width = 10, height = 10)
  draw(ht_fluid)
  dev.off()
  
  ht_fluid
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
  EN_SERUM_CNS       = EN_PGMC_serum_signature_CNS,
  EN_PLASMA_CNS      = EN_PGMC_plasma_signature_CNS,
  EN_CSF_CNS         = EN_PGMC_CSF_signature_CNS,
  
  Lasso_SERUM_IMMUNE  = lasso_PGMC_serum_signature_IMMUNE,
  Lasso_PLASMA_IMMUNE = lasso_PGMC_plasma_signature_IMMUNE,
  EN_SERUM_IMMUNE     = EN_PGMC_serum_signature_IMMUNE,
  EN_PLASMA_IMMUNE    = EN_PGMC_plasma_signature_IMMUNE,
  EN_CSF_IMMUNE       = EN_PGMC_CSF_signature_IMMUNE
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
                     normalize_panel_suffix = FALSE)  # same panel mode -> no suffix mismatch to resolve
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

clinical_data_extra_heatmap = read_excel("data input/all_participants_ID_visits_heatmap.xlsx")

# gt PGMCs and data of fluids
all_pgmc_ids <- protein_data_IDs %>% filter(type == "PGMC") %>% pull(ParticipantCode) %>% unique()
selection_matrix <- as.matrix(selection_wide[, names(model_registry)])  
column_split_fluid <- factor(model_fluid[colnames(selection_matrix)],
                             levels = c("SERUM", "PLASMA", "CSF"))

## add info on converters
converted_ids <- c("DE102", "TR119", "TR122","TR112")
converted_vec <- as.integer(rownames(selection_matrix) %in% converted_ids)
names(converted_vec) <- rownames(selection_matrix)

## Motor signs at all visits
motor_v0_vec <- to_row_vec(get_visit_value("TotalMotorSigns", "V0"))
motor_v1_vec <- to_row_vec(get_visit_value("TotalMotorSigns", "V1"))
motor_v2_vec <- to_row_vec(get_visit_value("TotalMotorSigns", "V2"))

##  NEFL at V0
nefl_serum_vec  <- to_row_vec(get_visit_value("Serum NEFL",  "V0"))
nefl_plasma_vec <- to_row_vec(get_visit_value("Plasma NEFL", "V0"))
nefl_csf_vec    <- to_row_vec(get_visit_value("CSF NEFL",    "V0"))

## urine at V0
urine_neopterin_vec <- to_row_vec(get_visit_value("umol neopterin/ mol creatinine", "V0"))
urine_p75_vec       <- to_row_vec(get_visit_value("ng p75ECD/mg creatinine",       "V0"))

## urine markers, delta between V0 and V1 
neopterin_v0 <- get_visit_value("umol neopterin/ mol creatinine", "V0")
neopterin_v1 <- get_visit_value("umol neopterin/ mol creatinine", "V1")
common_neo   <- intersect(names(neopterin_v0), names(neopterin_v1))
neopterin_delta_lookup <- setNames(neopterin_v1[common_neo] - neopterin_v0[common_neo], common_neo)
urine_neopterin_delta_vec <- to_row_vec(neopterin_delta_lookup)

p75_v0     <- get_visit_value("ng p75ECD/mg creatinine", "V0")
p75_v1     <- get_visit_value("ng p75ECD/mg creatinine", "V1")
common_p75 <- intersect(names(p75_v0), names(p75_v1))
p75_delta_lookup <- setNames(p75_v1[common_p75] - p75_v0[common_p75], common_p75)
urine_p75_delta_vec <- to_row_vec(p75_delta_lookup)

## get time to phenoconversion
phenoconversion_lookup <- clinical_data_extra_heatmap %>%
  filter(ParticipantCode %in% converted_ids, !is.na(`Disease duration`)) %>%
  distinct(ParticipantCode, `Disease duration`)

if (any(duplicated(phenoconversion_lookup$ParticipantCode))) {
  warning("Multiple distinct non-missing 'Disease duration' values found for at least one ",
          "converted patient across visits -- check clinical_data_extra_heatmap. ",
          "Taking the first one found per patient; verify this is correct.")
  phenoconversion_lookup <- phenoconversion_lookup %>% distinct(ParticipantCode, .keep_all = TRUE)
}

missing_conversion <- setdiff(converted_ids, phenoconversion_lookup$ParticipantCode)
if (length(missing_conversion) > 0) {
  warning("No non-missing 'Disease duration' value found for: ",
          paste(missing_conversion, collapse = ", "), " -- check these patients' records.")
}

phenoconversion_map <- phenoconversion_lookup %>%
  mutate(time_to_phenoconversion = -`Disease duration`) %>%
  select(ParticipantCode, time_to_phenoconversion) %>%
  deframe()

phenoconversion_vec <- setNames(rep(NA_real_, nrow(selection_matrix)), rownames(selection_matrix))
phenoconversion_vec[names(phenoconversion_map)] <- phenoconversion_map


## combined heatmap with results of all fluids
binary_col     <- c(`0` = "grey95", `1` = "#AD5291")
motor_col      <- colorRamp2(range(c(motor_v0_vec, motor_v1_vec, motor_v2_vec), na.rm = TRUE), c("grey95", "#2A6DB2"))
phenoconv_col  <- colorRamp2(range(phenoconversion_vec, na.rm = TRUE), c("#FEE8C8", "#B2242A")) 
nefl_col       <- colorRamp2(range(c(nefl_serum_vec, nefl_plasma_vec, nefl_csf_vec), na.rm = TRUE),
                             c("grey95", "#21918C"))
urine_neo_col  <- colorRamp2(range(urine_neopterin_vec, na.rm = TRUE), c("grey95", "#B2936F"))
urine_p75_col  <- colorRamp2(range(urine_p75_vec, na.rm = TRUE), c("grey95", "#B2936F"))         

neo_delta_range <- max(abs(urine_neopterin_delta_vec), na.rm = TRUE)
p75_delta_range <- max(abs(urine_p75_delta_vec), na.rm = TRUE)
urine_neo_delta_col <- colorRamp2(c(-neo_delta_range, 0, neo_delta_range), c("#2166AC", "white", "#B2182B"))
urine_p75_delta_col <- colorRamp2(c(-p75_delta_range, 0, p75_delta_range), c("#2166AC", "white", "#B2182B"))

row_anno <- rowAnnotation(
  Converted              = converted_vec,
  MotorSigns_V0          = motor_v0_vec,
  MotorSigns_V1          = motor_v1_vec,
  MotorSigns_V2          = motor_v2_vec,
  TimeToPhenoconv        = phenoconversion_vec,
  NEFL_Serum             = nefl_serum_vec,
  NEFL_Plasma            = nefl_plasma_vec,
  NEFL_CSF               = nefl_csf_vec,
  Urine_Neopterin_V0     = urine_neopterin_vec,
  Urine_p75ECD_V0        = urine_p75_vec,
  Urine_Neopterin_Delta  = urine_neopterin_delta_vec,
  Urine_p75ECD_Delta     = urine_p75_delta_vec,
  col = list(
    Converted             = c(`0` = "grey95", `1` = "#B2242A"),
    MotorSigns_V0         = motor_col,
    MotorSigns_V1         = motor_col,
    MotorSigns_V2         = motor_col,
    TimeToPhenoconv       = phenoconv_col,
    NEFL_Serum            = nefl_col,
    NEFL_Plasma           = nefl_col,
    NEFL_CSF              = nefl_col,
    Urine_Neopterin_V0    = urine_neo_col,
    Urine_p75ECD_VO       = urine_p75_col,
    Urine_Neopterin_Delta = urine_neo_delta_col,
    Urine_p75ECD_Delta    = urine_p75_delta_col
  ),
  na_col = "white",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = unit(0.35, "cm")
)

ht <- Heatmap(
  selection_matrix,
  name = "At-risk\n(score > 0)",
  col = binary_col,
  right_annotation = row_anno,
  column_split = column_split_fluid,
  cluster_column_slices = FALSE,   
  cluster_rows = TRUE,
  cluster_columns = FALSE,         
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_title = c("SERUM", "PLASMA", "CSF"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "At-risk\n(score > 0)",
    at = c("0", "1"),
    labels = c("Not selected", "Selected")
  )
)

pdf("plots/venn/PGMC_summary_heatmap_all_fluids.pdf", width = 17, height = 10)
draw(ht)
dev.off()

# same heatmap but by fluid 
nefl_by_fluid <- list(
  SERUM  = nefl_serum_vec,
  PLASMA = nefl_plasma_vec,
  CSF    = nefl_csf_vec
)

ht_serum  <- make_fluid_summary_heatmap("SERUM",  "plots/venn/PGMC_summary_heatmap_SERUM.pdf")
ht_plasma <- make_fluid_summary_heatmap("PLASMA", "plots/venn/PGMC_summary_heatmap_PLASMA.pdf")
ht_csf    <- make_fluid_summary_heatmap("CSF",    "plots/venn/PGMC_summary_heatmap_CSF.pdf")

