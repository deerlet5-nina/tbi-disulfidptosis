suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(data.table)
  library(edgeR)
  library(fgsea)
  library(limma)
  library(org.Hs.eg.db)
})

options(stringsAsFactors = FALSE)
set.seed(20260710)

required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Missing environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

optional_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value) || !file.exists(value)) return(NA_character_)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

output_dir <- required_env("NCR_OUTPUT_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
module_file <- required_env("NCR_MODULES")
modules <- fread(module_file, sep = "\t", data.table = FALSE)
modules$gene_symbol <- toupper(modules$gene_symbol)
pathways <- split(modules$gene_symbol, modules$module)
module_class <- setNames(modules$module_class, modules$module)
module_class <- module_class[!duplicated(names(module_class))]

all_gene_stats <- list()
all_module_effects <- list()
all_pathway_results <- list()
all_sample_scores <- list()
all_interactions <- list()
dataset_inventory <- list()

append_gene_stats <- function(x) all_gene_stats[[length(all_gene_stats) + 1L]] <<- x
append_module_effects <- function(x) all_module_effects[[length(all_module_effects) + 1L]] <<- x
append_pathway_results <- function(x) all_pathway_results[[length(all_pathway_results) + 1L]] <<- x
append_sample_scores <- function(x) all_sample_scores[[length(all_sample_scores) + 1L]] <<- x
append_interactions <- function(x) all_interactions[[length(all_interactions) + 1L]] <<- x
append_inventory <- function(x) dataset_inventory[[length(dataset_inventory) + 1L]] <<- x

bh <- function(p) p.adjust(p, method = "BH")

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("///.*$", "", x)
  x <- sub(";.*$", "", x)
  x <- sub(",.*$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

collapse_counts <- function(counts, symbols) {
  symbols <- clean_symbol(symbols)
  keep <- !is.na(symbols)
  counts <- as.matrix(counts[keep, , drop = FALSE])
  storage.mode(counts) <- "numeric"
  rowsum(counts, group = symbols[keep], reorder = FALSE)
}

map_human_ensembl_counts <- function(frame, id_column) {
  ids <- sub("\\..*$", "", as.character(frame[[id_column]]))
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(ids),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  symbols <- unname(mapped[ids])
  collapse_counts(frame[, setdiff(names(frame), id_column), drop = FALSE], symbols)
}

standard_gene_table <- function(tt, dataset, contrast, species, tissue, region, time, model, statistic_type = "moderated_t") {
  data.frame(
    dataset = dataset,
    contrast = contrast,
    species = species,
    tissue = tissue,
    region = region,
    time = time,
    model = model,
    gene_symbol = toupper(rownames(tt)),
    logFC = tt$logFC,
    statistic = tt$t,
    statistic_type = statistic_type,
    p_value = tt$P.Value,
    FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

run_voom <- function(counts, design, coefficient) {
  stopifnot(identical(colnames(counts), rownames(design)))
  dge <- DGEList(counts = round(counts))
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  tt <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  list(table = tt, expression = v$E, fit = fit)
}

run_limma_expression <- function(expression, design, coefficient) {
  stopifnot(identical(colnames(expression), rownames(design)))
  fit <- eBayes(lmFit(expression, design), robust = TRUE)
  tt <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  list(table = tt, expression = expression, fit = fit)
}

module_scores <- function(expression) {
  rownames(expression) <- toupper(rownames(expression))
  expression <- expression[!duplicated(rownames(expression)), , drop = FALSE]
  z <- t(scale(t(expression)))
  z[!is.finite(z)] <- NA_real_
  scores <- lapply(names(pathways), function(module) {
    present <- intersect(pathways[[module]], rownames(z))
    if (length(present) == 0L) return(rep(NA_real_, ncol(z)))
    colMeans(z[present, , drop = FALSE], na.rm = TRUE)
  })
  scores <- as.data.frame(scores, check.names = FALSE)
  names(scores) <- names(pathways)
  rownames(scores) <- colnames(expression)
  scores
}

hedges_g <- function(case, control) {
  case <- case[is.finite(case)]
  control <- control[is.finite(control)]
  n1 <- length(case)
  n0 <- length(control)
  if (n1 < 2L || n0 < 2L) return(NA_real_)
  pooled <- sqrt(((n1 - 1) * var(case) + (n0 - 1) * var(control)) / (n1 + n0 - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  d <- (mean(case) - mean(control)) / pooled
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  correction * d
}

two_group_module_effects <- function(scores, meta, group_column, case, control, dataset, contrast, species, tissue, region, time, model) {
  rows <- lapply(names(scores), function(module) {
    case_values <- scores[meta[[group_column]] == case, module]
    control_values <- scores[meta[[group_column]] == control, module]
    test <- tryCatch(t.test(case_values, control_values), error = function(e) NULL)
    data.frame(
      dataset = dataset,
      contrast = contrast,
      species = species,
      tissue = tissue,
      region = region,
      time = time,
      model = model,
      module = module,
      module_class = unname(module_class[module]),
      n_case = sum(is.finite(case_values)),
      n_control = sum(is.finite(control_values)),
      mean_case = mean(case_values, na.rm = TRUE),
      mean_control = mean(control_values, na.rm = TRUE),
      effect = hedges_g(case_values, control_values),
      effect_type = "Hedges_g",
      p_value = if (is.null(test)) NA_real_ else test$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- bh(out$p_value)
  out
}

linear_module_effects <- function(scores, meta, formula_rhs, coefficient, dataset, contrast, species, tissue, region, time, model) {
  rows <- lapply(names(scores), function(module) {
    dat <- cbind(meta, score = scores[meta$sample_id, module])
    fit <- lm(as.formula(paste("score ~", formula_rhs)), data = dat)
    co <- summary(fit)$coefficients
    if (!coefficient %in% rownames(co)) {
      beta <- p <- NA_real_
    } else {
      beta <- co[coefficient, "Estimate"]
      p <- co[coefficient, "Pr(>|t|)"]
    }
    data.frame(
      dataset = dataset,
      contrast = contrast,
      species = species,
      tissue = tissue,
      region = region,
      time = time,
      model = model,
      module = module,
      module_class = unname(module_class[module]),
      n_case = NA_integer_,
      n_control = NA_integer_,
      mean_case = NA_real_,
      mean_control = NA_real_,
      effect = beta,
      effect_type = "adjusted_beta",
      p_value = p,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$FDR <- bh(out$p_value)
  out
}

fgsea_from_gene_table <- function(gene_table) {
  stats <- gene_table$statistic
  names(stats) <- gene_table$gene_symbol
  ok <- is.finite(stats) & !is.na(names(stats)) & nzchar(names(stats))
  stats <- stats[ok]
  if (length(stats) < 100L) return(data.frame())
  split_stats <- split(stats, names(stats))
  stats <- vapply(split_stats, function(v) v[which.max(abs(v))], numeric(1))
  stats <- stats + rank(names(stats), ties.method = "first") * 1e-12
  result <- suppressWarnings(fgseaMultilevel(
    pathways = pathways,
    stats = sort(stats, decreasing = TRUE),
    minSize = 3,
    maxSize = 500,
    eps = 0
  ))
  if (nrow(result) == 0L) return(data.frame())
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
  result$dataset <- gene_table$dataset[1]
  result$contrast <- gene_table$contrast[1]
  result$species <- gene_table$species[1]
  result$tissue <- gene_table$tissue[1]
  result$region <- gene_table$region[1]
  result$time <- gene_table$time[1]
  result$model <- gene_table$model[1]
  result$module_class <- unname(module_class[result$pathway])
  result
}

record_analysis <- function(gene_table, module_table = NULL, scores = NULL, score_meta = NULL) {
  append_gene_stats(gene_table)
  enrichment <- fgsea_from_gene_table(gene_table)
  if (nrow(enrichment)) append_pathway_results(enrichment)
  if (!is.null(module_table)) append_module_effects(module_table)
  if (!is.null(scores) && !is.null(score_meta)) {
    long <- cbind(sample_id = rownames(scores), scores)
    long <- melt(as.data.table(long), id.vars = "sample_id", variable.name = "module", value.name = "score")
    long <- merge(long, as.data.table(score_meta), by = "sample_id", all.x = TRUE)
    long$analysis_contrast <- gene_table$contrast[1]
    append_sample_scores(as.data.frame(long))
  }
}

message("GSE193407: adjusted chronic CTE stage analysis")
g193_counts_frame <- fread(required_env("NCR_GSE193_COUNTS"), data.table = FALSE)
g193_meta <- fread(required_env("NCR_GSE193_META"), data.table = FALSE)
g193_counts <- map_human_ensembl_counts(g193_counts_frame, "Name")
rownames(g193_meta) <- g193_meta$Core_ID
common <- intersect(colnames(g193_counts), rownames(g193_meta))
g193_meta <- g193_meta[common, , drop = FALSE]
g193_counts <- g193_counts[, common, drop = FALSE]
complete <- complete.cases(g193_meta[, c("AGE", "RIN", "CTE_STAGE")])
g193_meta <- g193_meta[complete, , drop = FALSE]
g193_counts <- g193_counts[, rownames(g193_meta), drop = FALSE]
g193_meta$sample_id <- rownames(g193_meta)
g193_meta$stage <- as.numeric(g193_meta$CTE_STAGE)
g193_meta$age_z <- as.numeric(scale(g193_meta$AGE))
g193_meta$rin_z <- as.numeric(scale(g193_meta$RIN))
design <- model.matrix(~ age_z + rin_z + stage, data = g193_meta)
rownames(design) <- g193_meta$sample_id
g193_stage <- run_voom(g193_counts, design, "stage")
g193_stage_stats <- standard_gene_table(
  g193_stage$table, "GSE193407", "continuous_CTE_stage", "human", "postmortem brain", "BA9 cortex", "chronic course",
  "limma-voom; stage adjusted for age and RIN"
)
g193_scores <- module_scores(g193_stage$expression)
g193_module <- linear_module_effects(
  g193_scores, g193_meta, "age_z + rin_z + stage", "stage", "GSE193407", "continuous_CTE_stage", "human",
  "postmortem brain", "BA9 cortex", "chronic course", "module stage trend adjusted for age and RIN"
)
score_meta <- g193_meta[, c("sample_id", "stage", "AGE", "RIN")]
score_meta$dataset <- "GSE193407"
score_meta$condition <- ifelse(score_meta$stage == 0, "Stage0", paste0("Stage", score_meta$stage))
record_analysis(g193_stage_stats, g193_module, g193_scores, score_meta)

late_keep <- g193_meta$stage %in% c(0, 3, 4)
late_meta <- g193_meta[late_keep, , drop = FALSE]
late_meta$group <- factor(ifelse(late_meta$stage == 0, "Stage0", "Late"), levels = c("Stage0", "Late"))
late_counts <- g193_counts[, late_meta$sample_id, drop = FALSE]
late_design <- model.matrix(~ age_z + rin_z + group, data = late_meta)
rownames(late_design) <- late_meta$sample_id
g193_late <- run_voom(late_counts, late_design, "groupLate")
g193_late_stats <- standard_gene_table(
  g193_late$table, "GSE193407", "late_CTE_stage3_4_vs_stage0", "human", "postmortem brain", "BA9 cortex", "late chronic",
  "limma-voom; late-stage contrast adjusted for age and RIN"
)
late_scores <- module_scores(g193_late$expression)
late_module <- two_group_module_effects(
  late_scores, late_meta, "group", "Late", "Stage0", "GSE193407", "late_CTE_stage3_4_vs_stage0", "human",
  "postmortem brain", "BA9 cortex", "late chronic", "donor-level module contrast"
)
late_score_meta <- late_meta[, c("sample_id", "stage", "AGE", "RIN", "group")]
late_score_meta$dataset <- "GSE193407"
late_score_meta$condition <- late_score_meta$group
record_analysis(g193_late_stats, late_module, late_scores, late_score_meta)

coupling_outcomes <- c("Regulatory_bridge", "Cytoskeletal_endpoint", "WRC_actin_branching", "Focal_adhesion_tension", "Actomyosin_endpoint")
coupling_rows <- lapply(coupling_outcomes, function(outcome) {
  dat <- data.frame(
    outcome = as.numeric(scale(g193_scores[, outcome])),
    transport = as.numeric(scale(g193_scores[, "Transporter_entry"])),
    stage = as.numeric(scale(g193_meta$stage)),
    age = g193_meta$age_z,
    rin = g193_meta$rin_z
  )
  fit <- lm(outcome ~ transport + stage + age + rin, data = dat)
  interaction_fit <- lm(outcome ~ transport * stage + age + rin, data = dat)
  data.frame(
    dataset = "GSE193407",
    outcome_module = outcome,
    predictor_module = "Transporter_entry",
    adjusted_beta = coef(summary(fit))["transport", "Estimate"],
    p_value = coef(summary(fit))["transport", "Pr(>|t|)"],
    stage_interaction_beta = coef(summary(interaction_fit))["transport:stage", "Estimate"],
    stage_interaction_p = coef(summary(interaction_fit))["transport:stage", "Pr(>|t|)"],
    n = nrow(dat)
  )
})
g193_coupling <- do.call(rbind, coupling_rows)
g193_coupling$FDR <- bh(g193_coupling$p_value)
g193_coupling$stage_interaction_FDR <- bh(g193_coupling$stage_interaction_p)
fwrite(g193_coupling, file.path(output_dir, "GSE193407_transport_actin_coupling.csv"))
append_inventory(data.frame(dataset = "GSE193407", species = "human", design = "CTE stages 0-4", samples = nrow(g193_meta), primary_use = "adjusted chronic stage trend and transporter-actin coupling"))
rm(g193_counts_frame, g193_counts, g193_stage, g193_late); gc()

message("GSE319253: donor-collapsed external CTE replication")
g319_frame <- fread(required_env("NCR_GSE319_COUNTS"), data.table = FALSE)
g319_counts <- map_human_ensembl_counts(g319_frame, "gene_id")
replicate_names <- colnames(g319_counts)
donor_names <- sub("_REP[0-9]+$", "", replicate_names)
g319_counts <- t(rowsum(t(g319_counts), group = donor_names, reorder = FALSE))
g319_meta <- data.frame(sample_id = colnames(g319_counts))
g319_meta$condition <- factor(ifelse(grepl("^CTE", g319_meta$sample_id), "CTE", "Control"), levels = c("Control", "CTE"))
rownames(g319_meta) <- g319_meta$sample_id
design <- model.matrix(~ condition, data = g319_meta)
rownames(design) <- g319_meta$sample_id
g319 <- run_voom(g319_counts, design, "conditionCTE")
g319_stats <- standard_gene_table(g319$table, "GSE319253", "CTE_vs_control_donor_collapsed", "human", "postmortem brain", "superior frontal cortex", "chronic CTE", "limma-voom after summing technical replicates within donor")
g319_scores <- module_scores(g319$expression)
g319_module <- two_group_module_effects(g319_scores, g319_meta, "condition", "CTE", "Control", "GSE319253", "CTE_vs_control_donor_collapsed", "human", "postmortem brain", "superior frontal cortex", "chronic CTE", "donor-level module contrast")
g319_score_meta <- g319_meta
g319_score_meta$dataset <- "GSE319253"
record_analysis(g319_stats, g319_module, g319_scores, g319_score_meta)
append_inventory(data.frame(dataset = "GSE319253", species = "human", design = "6 CTE and 6 control donors; two technical replicates summed", samples = nrow(g319_meta), primary_use = "external chronic CTE replication"))
rm(g319_frame, g319_counts, g319); gc()

message("GSE209552: acute severe human TBI bulk RNA analysis")
g209_path <- required_env("NCR_GSE209_COUNTS")
g209_frame <- suppressWarnings(fread(g209_path, data.table = FALSE))
if (ncol(g209_frame) == 14L && names(g209_frame)[1] == "V1") {
  names(g209_frame)[1] <- "row_symbol"
} else if (!"row_symbol" %in% names(g209_frame)) {
  stop("Unexpected GSE209552 count-table width")
}
annotation_columns <- c("row_symbol", "Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_columns <- setdiff(names(g209_frame), annotation_columns)
g209_counts <- collapse_counts(g209_frame[, sample_columns, drop = FALSE], g209_frame$row_symbol)
g209_meta <- data.frame(sample_id = colnames(g209_counts))
g209_meta$condition <- factor(ifelse(grepl("ctrl", g209_meta$sample_id, ignore.case = TRUE), "Control", "TBI"), levels = c("Control", "TBI"))
rownames(g209_meta) <- g209_meta$sample_id
design <- model.matrix(~ condition, data = g209_meta)
rownames(design) <- g209_meta$sample_id
g209 <- run_voom(g209_counts, design, "conditionTBI")
g209_stats <- standard_gene_table(g209$table, "GSE209552", "acute_severe_TBI_vs_control", "human", "resected brain", "mixed lesion tissue", "4h_to_8d", "limma-voom; small acute surgical cohort")
g209_scores <- module_scores(g209$expression)
g209_module <- two_group_module_effects(g209_scores, g209_meta, "condition", "TBI", "Control", "GSE209552", "acute_severe_TBI_vs_control", "human", "resected brain", "mixed lesion tissue", "4h_to_8d", "sample-level module contrast")
g209_score_meta <- g209_meta
g209_score_meta$dataset <- "GSE209552"
record_analysis(g209_stats, g209_module, g209_scores, g209_score_meta)
append_inventory(data.frame(dataset = "GSE209552", species = "human", design = "4 acute severe TBI and 3 postmortem control bulk samples", samples = nrow(g209_meta), primary_use = "acute directionality; snRNA analyzed separately"))
rm(g209_frame, g209_counts, g209); gc()

message("GSE163415: treatment-adjusted mouse CCI time-region analysis")
g163_3 <- fread(required_env("NCR_GSE163_3D_COUNTS"), data.table = FALSE)
g163_29 <- fread(required_env("NCR_GSE163_29D_COUNTS"), data.table = FALSE)
g163_frame <- merge(g163_3, g163_29, by = "Gene", all = TRUE)
g163_frame[is.na(g163_frame)] <- 0
g163_counts <- collapse_counts(g163_frame[, setdiff(names(g163_frame), "Gene"), drop = FALSE], g163_frame$Gene)
g163_meta <- data.frame(sample_id = colnames(g163_counts))
parts <- strsplit(g163_meta$sample_id, "-", fixed = TRUE)
g163_meta$time <- vapply(parts, `[`, character(1), 1)
g163_meta$region <- vapply(parts, `[`, character(1), 2)
g163_meta$treatment <- factor(vapply(parts, `[`, character(1), 3), levels = c("Veh", "Drug"))
g163_meta$condition <- factor(vapply(parts, `[`, character(1), 4), levels = c("NoTBI", "TBI"))

for (time_value in unique(g163_meta$time)) {
  for (region_value in unique(g163_meta$region)) {
    keep <- g163_meta$time == time_value & g163_meta$region == region_value
    meta <- g163_meta[keep, , drop = FALSE]
    counts <- g163_counts[, meta$sample_id, drop = FALSE]
    rownames(meta) <- meta$sample_id
    design <- model.matrix(~ treatment + condition, data = meta)
    rownames(design) <- meta$sample_id
    fit <- run_voom(counts, design, "conditionTBI")
    contrast <- paste(time_value, region_value, "TBI_vs_NoTBI_adjusted_for_treatment", sep = "_")
    stats <- standard_gene_table(fit$table, "GSE163415", contrast, "mouse", "brain", region_value, time_value, "limma-voom additive injury effect adjusted for candesartan stratum")
    scores <- module_scores(fit$expression)
    module_effect <- two_group_module_effects(scores, meta, "condition", "TBI", "NoTBI", "GSE163415", contrast, "mouse", "brain", region_value, time_value, "module contrast across treatment strata")
    score_meta <- meta
    score_meta$dataset <- "GSE163415"
    record_analysis(stats, module_effect, scores, score_meta)

    interaction_design <- model.matrix(~ treatment * condition, data = meta)
    rownames(interaction_design) <- meta$sample_id
    interaction_fit <- run_voom(counts, interaction_design, "treatmentDrug:conditionTBI")
    interaction_stats <- standard_gene_table(interaction_fit$table, "GSE163415", paste(time_value, region_value, "treatment_by_injury", sep = "_"), "mouse", "brain", region_value, time_value, "limma-voom candesartan-by-injury interaction")
    append_interactions(interaction_stats[interaction_stats$gene_symbol %in% unique(modules$gene_symbol), , drop = FALSE])
    interaction_scores <- module_scores(interaction_fit$expression)
    interaction_module <- linear_module_effects(
      interaction_scores, meta, "treatment * condition", "treatmentDrug:conditionTBI",
      "GSE163415", paste(time_value, region_value, "treatment_by_injury", sep = "_"),
      "mouse", "brain", region_value, time_value, "module-level candesartan-by-injury interaction"
    )
    append_interactions(interaction_module)
  }
}
append_inventory(data.frame(dataset = "GSE163415", species = "mouse", design = "CCI by 3/29 DPI, three regions and vehicle/candesartan strata", samples = nrow(g163_meta), primary_use = "time-region replication with treatment adjustment"))
rm(g163_3, g163_29, g163_frame, g163_counts); gc()

parse_sample_titles <- function(series_path) {
  lines <- readLines(gzfile(series_path), warn = FALSE)
  title_line <- lines[grepl("^!Sample_title", lines)][1]
  if (is.na(title_line)) stop("No !Sample_title in ", series_path)
  fields <- strsplit(title_line, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", fields)
}

message("GSE298240: chronic repetitive closed-head injury with sex and treatment adjustment")
g298_frame <- fread(required_env("NCR_GSE298_COUNTS"), data.table = FALSE)
g298_sample_columns <- setdiff(names(g298_frame), "Gene_name")
g298_counts <- collapse_counts(g298_frame[, g298_sample_columns, drop = FALSE], g298_frame$Gene_name)
g298_titles <- parse_sample_titles(required_env("NCR_GSE298_SERIES"))
if (length(g298_titles) != length(g298_sample_columns)) stop("GSE298240 title/sample mismatch")
g298_meta <- data.frame(sample_id = g298_sample_columns, title = g298_titles)
g298_meta$sex <- factor(sub("_.*$", "", g298_meta$title), levels = c("F", "M"))
g298_meta$injury <- factor(ifelse(grepl("5xCHI", g298_meta$title, ignore.case = TRUE), "Injury", "Sham"), levels = c("Sham", "Injury"))
g298_meta$treatment <- factor(ifelse(grepl("\\+ SB", g298_meta$title, ignore.case = TRUE), "SB", "Vehicle"), levels = c("Vehicle", "SB"))
rownames(g298_meta) <- g298_meta$sample_id
design <- model.matrix(~ sex + treatment + injury, data = g298_meta)
rownames(design) <- g298_meta$sample_id
g298 <- run_voom(g298_counts, design, "injuryInjury")
g298_stats <- standard_gene_table(g298$table, "GSE298240", "chronic_5xCHI_vs_sham_adjusted", "mouse", "brain", "somatosensory cortex", "chronic", "limma-voom injury effect adjusted for sex and p38-inhibitor SB239063 stratum")
g298_scores <- module_scores(g298$expression)
g298_module <- linear_module_effects(g298_scores, g298_meta, "sex + treatment + injury", "injuryInjury", "GSE298240", "chronic_5xCHI_vs_sham_adjusted", "mouse", "brain", "somatosensory cortex", "chronic", "injury effect adjusted for sex and p38-inhibitor SB239063 stratum")
g298_score_meta <- g298_meta
g298_score_meta$dataset <- "GSE298240"
g298_score_meta$condition <- g298_score_meta$injury
record_analysis(g298_stats, g298_module, g298_scores, g298_score_meta)

for (spec in list(
  list(formula = "treatment + sex * injury", coefficient = "sexM:injuryInjury", label = "sex_by_injury"),
  list(formula = "sex + treatment * injury", coefficient = "treatmentSB:injuryInjury", label = "treatment_by_injury")
)) {
  effects <- linear_module_effects(g298_scores, g298_meta, spec$formula, spec$coefficient, "GSE298240", spec$label, "mouse", "brain", "somatosensory cortex", "chronic", spec$label)
  append_interactions(effects)
}
append_inventory(data.frame(dataset = "GSE298240", species = "mouse", design = "chronic 5xCHI/sham by sex and p38-inhibitor SB239063/vehicle", samples = nrow(g298_meta), primary_use = "chronic repetitive-injury replication and heterogeneity"))
rm(g298_frame, g298_counts, g298); gc()

message("GSE111452: rat FPI trajectory from 24 hours to 12 months")
for (platform in c("GPL22740", "GPL15084")) {
  expr_path <- required_env(paste0("NCR_GSE111_", platform, "_EXPR"))
  meta_path <- required_env(paste0("NCR_GSE111_", platform, "_META"))
  expression <- as.data.frame(fread(expr_path), check.names = FALSE)
  rownames(expression) <- expression[[1]]
  expression <- as.matrix(expression[, -1, drop = FALSE])
  storage.mode(expression) <- "numeric"
  expression <- normalizeBetweenArrays(expression, method = "quantile")
  meta_all <- fread(meta_path, data.table = FALSE)
  rownames(meta_all) <- meta_all$sample_id
  for (region_value in unique(meta_all$region)) {
    for (time_value in unique(meta_all$time)) {
      meta <- meta_all[meta_all$region == region_value & meta_all$time == time_value & meta_all$condition %in% c("Sham", "TBI"), , drop = FALSE]
      if (sum(meta$condition == "TBI") < 2L || sum(meta$condition == "Sham") < 2L) next
      expr <- expression[, meta$sample_id, drop = FALSE]
      meta$condition <- factor(meta$condition, levels = c("Sham", "TBI"))
      design <- model.matrix(~ condition, data = meta)
      rownames(design) <- meta$sample_id
      fit <- run_limma_expression(expr, design, "conditionTBI")
      contrast <- paste(platform, region_value, time_value, "TBI_vs_sham", sep = "_")
      stats <- standard_gene_table(fit$table, "GSE111452", contrast, "rat", "brain", region_value, time_value, "limma empirical-Bayes microarray contrast")
      scores <- module_scores(fit$expression)
      module_effect <- two_group_module_effects(scores, meta, "condition", "TBI", "Sham", "GSE111452", contrast, "rat", "brain", region_value, time_value, "sample-level module contrast")
      score_meta <- meta
      score_meta$dataset <- "GSE111452"
      record_analysis(stats, module_effect, scores, score_meta)
    }
  }
  append_inventory(data.frame(dataset = paste0("GSE111452_", platform), species = "rat", design = "FPI cortex/hippocampus by 24h-12mo", samples = nrow(meta_all), primary_use = "long-horizon regional trajectory"))
  rm(expression, meta_all); gc()
}

g223_expr_path <- optional_env("NCR_GSE223_EXPR")
g223_meta_path <- optional_env("NCR_GSE223_META")
if (!is.na(g223_expr_path) && !is.na(g223_meta_path)) {
  message("GSE223245: exploratory peripheral severity trend")
  expression <- as.data.frame(fread(g223_expr_path), check.names = FALSE)
  rownames(expression) <- expression[[1]]
  expression <- as.matrix(expression[, -1, drop = FALSE])
  storage.mode(expression) <- "numeric"
  meta <- fread(g223_meta_path, data.table = FALSE)
  rownames(meta) <- meta$sample_id
  expression <- expression[, meta$sample_id, drop = FALSE]
  scores <- module_scores(expression)
  module_rows <- lapply(names(scores), function(module) {
    test <- suppressWarnings(cor.test(scores[, module], meta$severity, method = "spearman", exact = FALSE))
    data.frame(
      dataset = "GSE223245", contrast = "ordinal_peripheral_severity", species = "human", tissue = "whole blood/PBMC",
      region = "peripheral", time = "acute", model = "Spearman severity trend", module = module,
      module_class = unname(module_class[module]), n_case = sum(meta$severity > 0), n_control = sum(meta$severity == 0),
      mean_case = mean(scores[meta$severity > 0, module]), mean_control = mean(scores[meta$severity == 0, module]),
      effect = unname(test$estimate), effect_type = "Spearman_r", p_value = test$p.value
    )
  })
  module_table <- do.call(rbind, module_rows)
  module_table$FDR <- bh(module_table$p_value)
  gene_rows <- lapply(seq_len(nrow(expression)), function(i) {
    test <- suppressWarnings(cor.test(expression[i, ], meta$severity, method = "spearman", exact = FALSE))
    c(statistic = unname(test$estimate), p_value = test$p.value)
  })
  gene_mat <- do.call(rbind, gene_rows)
  gene_table <- data.frame(
    dataset = "GSE223245", contrast = "ordinal_peripheral_severity", species = "human", tissue = "whole blood/PBMC",
    region = "peripheral", time = "acute", model = "Spearman severity trend", gene_symbol = toupper(rownames(expression)),
    logFC = NA_real_, statistic = gene_mat[, "statistic"], statistic_type = "Spearman_r", p_value = gene_mat[, "p_value"]
  )
  gene_table$FDR <- bh(gene_table$p_value)
  score_meta <- meta
  score_meta$dataset <- "GSE223245"
  score_meta$condition <- score_meta$group
  record_analysis(gene_table, module_table, scores, score_meta)
  append_inventory(data.frame(dataset = "GSE223245", species = "human", design = "control/mild/moderate/severe peripheral samples", samples = nrow(meta), primary_use = "exploratory severity context only"))
}

g104_path <- optional_env("NCR_GSE104_STATS")
if (!is.na(g104_path)) {
  message("GSE104687: importing donor-level regional statistics")
  g104 <- fread(g104_path, data.table = FALSE)
  for (comparison_value in unique(g104$comparison)) {
    sub <- g104[g104$comparison == comparison_value, , drop = FALSE]
    gene_table <- data.frame(
      dataset = "GSE104687", contrast = comparison_value, species = "human", tissue = "postmortem brain",
      region = sub$region[1], time = "remote history", model = "donor-level regional Welch analysis from preserved source table",
      gene_symbol = toupper(sub$gene_symbol), logFC = sub$logFC_TBI_vs_Control, statistic = sub$t_statistic,
      statistic_type = "Welch_t", p_value = sub$p_value, FDR = sub$adj_p_value
    )
    record_analysis(gene_table)
  }
  append_inventory(data.frame(dataset = "GSE104687", species = "human", design = "same-donor remote multi-region postmortem comparison", samples = max(g104$n_donor_TBI + g104$n_donor_Control, na.rm = TRUE), primary_use = "remote-regional boundary check"))
}

gene_stats <- rbindlist(all_gene_stats, fill = TRUE)
module_effects <- rbindlist(all_module_effects, fill = TRUE)
pathway_results <- rbindlist(all_pathway_results, fill = TRUE)
sample_scores <- rbindlist(all_sample_scores, fill = TRUE)
interactions <- rbindlist(all_interactions, fill = TRUE)
inventory <- unique(rbindlist(dataset_inventory, fill = TRUE))

core_genes <- unique(modules$gene_symbol[modules$module == "Eight_gene_panel"])
core_stats <- gene_stats[gene_symbol %in% core_genes]

fwrite(gene_stats, file.path(output_dir, "Table_S_all_bulk_gene_statistics.csv.gz"))
fwrite(core_stats, file.path(output_dir, "Table_core_8gene_effects.csv"))
fwrite(module_effects, file.path(output_dir, "Table_module_effects.csv"))
fwrite(pathway_results, file.path(output_dir, "Table_ranked_module_enrichment.csv"))
fwrite(sample_scores, file.path(output_dir, "Table_sample_module_scores.csv.gz"))
fwrite(interactions, file.path(output_dir, "Table_sex_treatment_interactions.csv"))
fwrite(inventory, file.path(output_dir, "Table_dataset_inventory.csv"))
fwrite(modules, file.path(output_dir, "Table_mechanism_module_definitions.csv"))

capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
message("Expanded bulk analysis completed: ", output_dir)
