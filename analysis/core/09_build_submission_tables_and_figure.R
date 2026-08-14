#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) tolower(args[[1]]) else "final"
stopifnot(phase %in% c("preview", "final"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
point_dir <- file.path(root, "03_outputs", "04_formal_point")
boot_dir <- file.path(root, "03_outputs", paste0("06_bootstrap_summary_", phase))
gate_dir <- file.path(root, "03_outputs", "01_gate0")
sensitivity_dir <- file.path(root, "03_outputs", "07_sensitivity_point")
threshold_dir <- file.path(root, "03_outputs", "09_threshold_sensitivity")
table_dir <- file.path(root, "03_outputs", paste0("08_submission_tables_", phase))
figure_dir <- file.path(root, if (phase == "final") "04_figures" else "04_figures_preview")
panel_dir <- file.path(figure_dir, "panels")
source_dir <- file.path(figure_dir, "source_data")
qc_dir <- file.path(figure_dir, "qc")
log_dir <- file.path(root, "06_logs")
for (d in c(table_dir, figure_dir, panel_dir, source_dir, qc_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required <- c(
  file.path(boot_dir, "bootstrap_percentile_intervals.csv"),
  file.path(boot_dir, "mechanism_promotion_by_cohort.csv"),
  file.path(boot_dir, "bootstrap_promotion_decision.csv"),
  file.path(point_dir, "shapley_decomposition_point.csv"),
  file.path(point_dir, "model_coefficients_model_based.csv"),
  file.path(gate_dir, "gate0_event_counts_by_cohort.csv")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) stop("Missing required inputs: ", paste(missing_required, collapse = "; "))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
exposure_order <- c("education", "wealth")
exposure_label <- c(education = "Education", wealth = "Wealth")

ci <- fread(file.path(boot_dir, "bootstrap_percentile_intervals.csv"))
promotion <- fread(file.path(boot_dir, "mechanism_promotion_by_cohort.csv"))
decision <- fread(file.path(boot_dir, "bootstrap_promotion_decision.csv"))
shapley <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[estimand == "dfle"]
coefficients <- fread(file.path(point_dir, "model_coefficients_model_based.csv"))
events <- fread(file.path(gate_dir, "gate0_event_counts_by_cohort.csv"))

primary_absolute <- ci[metric %in% c("gap_dfle", "gap_tle", "gap_disabled"), .(
  cohort, exposure, estimand = sub("^gap_", "", metric), point_estimate,
  bootstrap_valid_replicates = valid_replicates, ci_low, ci_high
)]
primary_absolute[, `:=`(
  .cohort_order = match(cohort, cohort_order),
  .exposure_order = match(exposure, exposure_order)
)]
setorder(primary_absolute, .cohort_order, .exposure_order, estimand)
primary_absolute[, c(".cohort_order", ".exposure_order") := NULL]

mechanism <- promotion[, .(
  cohort, exposure, dfle_gap, recovery_relapse_contribution_years,
  recovery_relapse_percent, rr_ci_low, rr_ci_high,
  rr_percent_ci_low, rr_percent_ci_high, rr_ci_excludes_zero,
  point_threshold_met, validity_gate_pass, mechanism_trigger_cohort
)]

block_ci <- ci[grepl("^contribution_", metric), .(
  cohort, exposure, block = sub("^contribution_", "", metric),
  point_estimate, bootstrap_valid_replicates = valid_replicates, ci_low, ci_high
)]
model_hr <- coefficients[term == "seslow", .(
  cohort, exposure, process, low_vs_high_hazard_ratio = hazard_ratio,
  model_based_ci_low = hr_ci_low_model, model_based_ci_high = hr_ci_high_model
)]

fwrite(events, file.path(table_dir, "TableS1_gate0_event_counts.csv"))
fwrite(primary_absolute, file.path(table_dir, "Table1_absolute_life_year_gaps.csv"))
fwrite(mechanism, file.path(table_dir, "Table2_recovery_relapse_mechanism.csv"))
fwrite(block_ci, file.path(table_dir, "TableS2_shapley_contributions.csv"))
fwrite(model_hr, file.path(table_dir, "TableS3_transition_hazard_ratios.csv"))
fwrite(decision, file.path(table_dir, "TableS4_binding_promotion_decision.csv"))
if (file.exists(file.path(sensitivity_dir, "sensitivity_direction_and_magnitude.csv"))) {
  invisible(file.copy(
    file.path(sensitivity_dir, "sensitivity_direction_and_magnitude.csv"),
    file.path(table_dir, "TableS5_sensitivity_direction_and_magnitude.csv"), overwrite = TRUE
  ))
}
if (file.exists(file.path(threshold_dir, "threshold_gate0_event_counts.csv"))) {
  invisible(file.copy(
    file.path(threshold_dir, "threshold_gate0_event_counts.csv"),
    file.path(table_dir, "TableS6_threshold_gate0_event_counts.csv"), overwrite = TRUE
  ))
}
threshold_parts <- list()
for (threshold_name in c("adl_only", "at_least_two", "permissive_partial")) {
  threshold_path <- file.path(threshold_dir, threshold_name, "point_promotion_screen.csv")
  if (file.exists(threshold_path)) {
    z <- fread(threshold_path)
    z[, `:=`(threshold = threshold_name, model_status = "COMPLETED")]
  } else {
    z <- CJ(cohort = cohort_order, exposure = exposure_order)
    z[, `:=`(
      recovery_relapse_contribution_years = NA_real_, dfle_gap = NA_real_,
      recovery_relapse_percent = NA_real_, point_threshold_met = FALSE,
      threshold = threshold_name,
      model_status = "STOP: full-pipeline nonconvergence"
    )]
  }
  threshold_parts[[length(threshold_parts) + 1L]] <- z
}
threshold_combined <- rbindlist(threshold_parts, fill = TRUE)
threshold_combined[, `:=`(
  .threshold_order = match(threshold, c("adl_only", "at_least_two", "permissive_partial")),
  .cohort_order = match(cohort, cohort_order),
  .exposure_order = match(exposure, exposure_order)
)]
setorder(threshold_combined, .threshold_order, .cohort_order, .exposure_order)
threshold_combined[, c(".threshold_order", ".cohort_order", ".exposure_order") := NULL]
fwrite(threshold_combined, file.path(table_dir, "TableS7_threshold_mechanism_sensitivity.csv"))

panel_b <- primary_absolute[estimand == "dfle", .(
  cohort, exposure, point = point_estimate, ci_low, ci_high,
  bootstrap_valid_replicates
)]
panel_c <- mechanism[, .(
  cohort, exposure, point = recovery_relapse_contribution_years,
  ci_low = rr_ci_low, ci_high = rr_ci_high,
  percent = recovery_relapse_percent, percent_ci_low = rr_percent_ci_low,
  percent_ci_high = rr_percent_ci_high, mechanism_trigger_cohort
)]
panel_d <- block_ci[, .(cohort, exposure, block, point = point_estimate, ci_low, ci_high)]

fwrite(panel_b, file.path(source_dir, "Figure1B_dfle_gap.csv"))
fwrite(panel_c, file.path(source_dir, "Figure1C_recovery_relapse_contribution.csv"))
fwrite(panel_d, file.path(source_dir, "Figure1D_shapley_blocks.csv"))
fwrite(events, file.path(source_dir, "Figure1A_gate0_context.csv"))

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

line_svg <- function(x1, y1, x2, y2, class = "axis", extra = "") {
  sprintf('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" class="%s" %s/>', x1, y1, x2, y2, class, extra)
}
rect_svg <- function(x, y, w, h, fill, stroke = "none", rx = 0, class = "") {
  sprintf('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="%.2f" fill="%s" stroke="%s" class="%s"/>', x, y, w, h, rx, fill, stroke, class)
}
circle_svg <- function(x, y, r, fill, stroke = "white", sw = 1) {
  sprintf('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" stroke="%s" stroke-width="%.2f"/>', x, y, r, fill, stroke, sw)
}
text_svg <- function(x, y, label, class = "label", anchor = "start", extra = "") {
  sprintf('<text x="%.2f" y="%.2f" class="%s" text-anchor="%s" %s>%s</text>', x, y, class, anchor, extra, xml_escape(label))
}
polygon_svg <- function(points, fill) sprintf('<polygon points="%s" fill="%s"/>', points, fill)
arrow_svg <- function(x1, y1, x2, y2, color = "#424A53", head_length = 6, head_half_width = 3.5) {
  dx <- x2 - x1
  dy <- y2 - y1
  distance <- sqrt(dx^2 + dy^2)
  ux <- dx / distance
  uy <- dy / distance
  base_x <- x2 - head_length * ux
  base_y <- y2 - head_length * uy
  perp_x <- -uy * head_half_width
  perp_y <- ux * head_half_width
  points <- sprintf(
    "%.2f,%.2f %.2f,%.2f %.2f,%.2f",
    x2, y2, base_x + perp_x, base_y + perp_y, base_x - perp_x, base_y - perp_y
  )
  c(line_svg(x1, y1, x2, y2, "axis", sprintf('stroke="%s"', color)), polygon_svg(points, color))
}

panel_groups <- function(letter, axis, data, labels, legend, annotations, transform = NULL) {
  tr <- if (is.null(transform)) "" else sprintf(' transform="%s"', transform)
  c(
    sprintf('<g id="panel_%s"%s>', letter, tr),
    sprintf('<g id="panel_%s_axis">', letter), axis, "</g>",
    sprintf('<g id="panel_%s_data">', letter), data, "</g>",
    sprintf('<g id="panel_%s_labels">', letter), labels, "</g>",
    sprintf('<g id="panel_%s_legend">', letter), legend, "</g>",
    sprintf('<g id="panel_%s_annotations">', letter), annotations, "</g>",
    "</g>"
  )
}

svg_document <- function(width, height, content) {
  c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', width, height, width, height),
    '<style><![CDATA[
      text { font-family: Arial, Helvetica, sans-serif; fill: #20252B; }
      .panel-letter { font-size: 18px; font-weight: 700; }
      .title { font-size: 13px; font-weight: 700; }
      .subtitle { font-size: 9px; fill: #58616B; }
      .label { font-size: 9px; }
      .small { font-size: 7.5px; }
      .axis { stroke: #424A53; stroke-width: 1; fill: none; }
      .grid { stroke: #D9DEE3; stroke-width: 0.7; fill: none; }
      .zero { stroke: #20252B; stroke-width: 1.1; fill: none; }
      .ci { stroke-width: 1.8; fill: none; stroke-linecap: round; }
    ]]></style>',
    rect_svg(0, 0, width, height, "#FFFFFF"), content, "</svg>"
  )
}

write_svg <- function(path, width, height, content) {
  writeLines(svg_document(width, height, content), path, useBytes = TRUE)
}

make_panel_a <- function(transform = NULL) {
  axis <- c(rect_svg(5, 5, 340, 288, "none", "#D9DEE3", 2), line_svg(15, 46, 335, 46, "grid"))
  node_x <- c(I0 = 25, D1 = 125, R1 = 225, D2 = 225, DEAD = 125)
  node_y <- c(I0 = 105, D1 = 105, R1 = 75, D2 = 150, DEAD = 225)
  node_fill <- c(I0 = "#D9F0EA", D1 = "#FDE3D0", R1 = "#E7E1F4", D2 = "#F6CFD8", DEAD = "#D8DCE1")
  data <- character()
  labels <- character()
  for (st in names(node_x)) {
    data <- c(data, rect_svg(node_x[[st]], node_y[[st]], 78, 34, node_fill[[st]], "#56606B", 7))
  }
  arrows <- list(
    c(103, 122, 125, 122),
    c(203, 112, 225, 96),
    c(254, 109, 254, 150),
    c(274, 150, 274, 109),
    c(164, 139, 164, 225),
    c(245, 109, 184, 225),
    c(264, 184, 190, 232),
    c(103, 130, 145, 225)
  )
  for (a in arrows) {
    data <- c(data, arrow_svg(a[1], a[2], a[3], a[4]))
  }
  node_labels <- c(I0 = "Independent", D1 = "First difficulty", R1 = "Recovered", D2 = "Relapsed", DEAD = "Death")
  for (st in names(node_x)) labels <- c(labels, text_svg(node_x[[st]] + 39, node_y[[st]] + 21, node_labels[[st]], "label", "middle"))
  labels <- c(labels,
    text_svg(113, 111, "onset", "small", "middle"),
    text_svg(216, 93, "recovery", "small", "end"),
    text_svg(244, 137, "relapse", "small", "end"),
    text_svg(282, 137, "recovery", "small"),
    text_svg(171, 193, "death from any living state", "small")
  )
  legend <- c(
    rect_svg(22, 262, 12, 8, "#D9F0EA", "#56606B", 1), text_svg(39, 270, "disability-free", "small"),
    rect_svg(116, 262, 12, 8, "#FDE3D0", "#56606B", 1), text_svg(133, 270, "with difficulty", "small"),
    rect_svg(227, 262, 12, 8, "#E7E1F4", "#56606B", 1), text_svg(244, 270, "recovery history", "small")
  )
  annotations <- c(
    text_svg(14, 28, "A", "panel-letter"),
    text_svg(42, 27, "Frozen five-state process", "title"),
    text_svg(42, 40, "Any difficulty in common 5 ADL + 4 IADL; attrition is not death", "subtitle")
  )
  panel_groups("A", axis, data, labels, legend, annotations, transform)
}

forest_panel <- function(letter, z, title, subtitle, x_title, colors, transform = NULL) {
  z <- copy(z)
  z[, cohort_order_n := match(cohort, cohort_order)]
  z[, exposure_order_n := match(exposure, exposure_order)]
  setorder(z, cohort_order_n, exposure_order_n)
  finite_range <- range(c(z$ci_low, z$ci_high, 0), finite = TRUE)
  pad <- max(0.25, diff(finite_range) * 0.10)
  x_min <- floor((finite_range[1] - pad) * 2) / 2
  x_max <- ceiling((finite_range[2] + pad) * 2) / 2
  plot_left <- 82
  plot_right <- 332
  plot_top <- 65
  plot_bottom <- 245
  xmap <- function(x) plot_left + (x - x_min) / (x_max - x_min) * (plot_right - plot_left)
  cohort_y <- setNames(c(82, 124, 166, 208), cohort_order)
  offset <- c(education = -6, wealth = 6)
  tick_step <- if ((x_max - x_min) <= 4) 1 else 2
  ticks <- seq(ceiling(x_min / tick_step) * tick_step, floor(x_max / tick_step) * tick_step, by = tick_step)
  axis <- c(rect_svg(5, 5, 340, 288, "none", "#D9DEE3", 2))
  labels <- character()
  for (tick in ticks) {
    axis <- c(axis, line_svg(xmap(tick), plot_top, xmap(tick), plot_bottom, if (tick == 0) "zero" else "grid"))
    labels <- c(labels, text_svg(xmap(tick), 260, format(tick, trim = TRUE), "small", "middle"))
  }
  for (cohort_name in cohort_order) {
    axis <- c(axis, line_svg(plot_left, cohort_y[[cohort_name]] + 18, plot_right, cohort_y[[cohort_name]] + 18, "grid"))
    labels <- c(labels, text_svg(74, cohort_y[[cohort_name]] + 3, cohort_name, "label", "end"))
  }
  axis <- c(axis, line_svg(plot_left, plot_bottom, plot_right, plot_bottom, "axis"))
  labels <- c(labels, text_svg((plot_left + plot_right) / 2, 279, x_title, "small", "middle"))
  data <- character()
  for (i in seq_len(nrow(z))) {
    yy <- cohort_y[[z$cohort[[i]]]] + offset[[z$exposure[[i]]]]
    color <- colors[[z$exposure[[i]]]]
    data <- c(data,
      line_svg(xmap(z$ci_low[[i]]), yy, xmap(z$ci_high[[i]]), yy, "ci", sprintf('stroke="%s"', color)),
      line_svg(xmap(z$ci_low[[i]]), yy - 3, xmap(z$ci_low[[i]]), yy + 3, "ci", sprintf('stroke="%s"', color)),
      line_svg(xmap(z$ci_high[[i]]), yy - 3, xmap(z$ci_high[[i]]), yy + 3, "ci", sprintf('stroke="%s"', color)),
      circle_svg(xmap(z$point[[i]]), yy, 4, color)
    )
  }
  legend <- c(
    circle_svg(205, 49, 4, colors[["education"]]), text_svg(214, 52, "Education", "small"),
    circle_svg(275, 49, 4, colors[["wealth"]]), text_svg(284, 52, "Wealth", "small")
  )
  annotations <- c(
    text_svg(14, 28, letter, "panel-letter"),
    text_svg(42, 27, title, "title"),
    text_svg(42, 40, subtitle, "subtitle")
  )
  panel_groups(letter, axis, data, labels, legend, annotations, transform)
}

make_panel_d <- function(transform = NULL) {
  z <- copy(panel_d)
  z[, row_key := paste(cohort, ifelse(exposure == "education", "Edu", "Wealth"), sep = " · ")]
  row_order <- as.vector(t(outer(cohort_order, c("Edu", "Wealth"), paste, sep = " · ")))
  z[, row_n := match(row_key, row_order)]
  setorder(z, row_n)
  block_order <- c("onset", "recovery", "relapse", "post_disability_mortality", "pre_disability_mortality")
  block_colors <- c(
    onset = "#1B9E77", recovery = "#D95F02", relapse = "#7570B3",
    post_disability_mortality = "#E7298A", pre_disability_mortality = "#66A61E"
  )
  x_min <- min(-6.5, floor(min(z[, sum(point[point < 0]), by = row_key]$V1) * 2) / 2)
  x_max <- max(1, ceiling(max(z[, sum(point[point > 0]), by = row_key]$V1) * 2) / 2)
  plot_left <- 100
  plot_right <- 334
  plot_top <- 62
  plot_bottom <- 236
  xmap <- function(x) plot_left + (x - x_min) / (x_max - x_min) * (plot_right - plot_left)
  row_y <- setNames(seq(77, 217, length.out = length(row_order)), row_order)
  ticks <- seq(ceiling(x_min), floor(x_max), by = 1)
  axis <- c(rect_svg(5, 5, 340, 288, "none", "#D9DEE3", 2))
  labels <- character()
  for (tick in ticks) {
    axis <- c(axis, line_svg(xmap(tick), plot_top, xmap(tick), plot_bottom, if (tick == 0) "zero" else "grid"))
    labels <- c(labels, text_svg(xmap(tick), 249, format(tick, trim = TRUE), "small", "middle"))
  }
  axis <- c(axis, line_svg(plot_left, plot_bottom, plot_right, plot_bottom, "axis"))
  data <- character()
  for (rk in row_order) {
    zz <- z[row_key == rk]
    neg_cursor <- 0
    pos_cursor <- 0
    for (block_name in block_order) {
      value <- zz[block == block_name, point]
      if (!length(value)) next
      if (value < 0) {
        x0 <- neg_cursor + value
        x1 <- neg_cursor
        neg_cursor <- x0
      } else {
        x0 <- pos_cursor
        x1 <- pos_cursor + value
        pos_cursor <- x1
      }
      data <- c(data, rect_svg(xmap(x0), row_y[[rk]] - 5, max(0.5, xmap(x1) - xmap(x0)), 10, block_colors[[block_name]], "white", 0))
    }
    labels <- c(labels, text_svg(94, row_y[[rk]] + 3, rk, "small", "end"))
  }
  labels <- c(labels, text_svg((plot_left + plot_right) / 2, 263, "Contribution to low–high DFLE gap (years)", "small", "middle"))
  legend_labels <- c("Onset", "Recovery", "Relapse", "Post-disability death", "Pre-disability death")
  legend <- character()
  legend_x <- c(16, 78, 148, 210, 16)
  legend_y <- c(273, 273, 273, 273, 286)
  for (i in seq_along(block_order)) {
    legend <- c(legend,
      rect_svg(legend_x[[i]], legend_y[[i]] - 7, 9, 7, block_colors[[block_order[[i]]]], "none"),
      text_svg(legend_x[[i]] + 13, legend_y[[i]], legend_labels[[i]], "small")
    )
  }
  annotations <- c(
    text_svg(14, 28, "D", "panel-letter"),
    text_svg(42, 27, "What generates the DFLE divide?", "title"),
    text_svg(42, 40, "Five-block Shapley decomposition; negative values widen the deficit", "subtitle")
  )
  panel_groups("D", axis, data, labels, legend, annotations, transform)
}

colors <- c(education = "#0072B2", wealth = "#D55E00")
panel_a_content <- make_panel_a()
panel_b_content <- forest_panel(
  "B", panel_b, "Remaining disability-free life expectancy",
  "Low minus high SES at age 60; 95% person-bootstrap intervals",
  "DFLE gap (years)", colors
)
panel_c_content <- forest_panel(
  "C", panel_c, "Recovery and relapse contribution",
  "Combined contribution to the low–high DFLE gap; 95% intervals",
  "Recovery + relapse contribution (years)", colors
)
panel_d_content <- make_panel_d()

write_svg(file.path(panel_dir, "Figure1A_state_model.svg"), 350, 298, panel_a_content)
write_svg(file.path(panel_dir, "Figure1B_dfle_gap.svg"), 350, 298, panel_b_content)
write_svg(file.path(panel_dir, "Figure1C_recovery_relapse.svg"), 350, 298, panel_c_content)
write_svg(file.path(panel_dir, "Figure1D_shapley_decomposition.svg"), 350, 298, panel_d_content)

composite <- c(
  make_panel_a("translate(10,10)"),
  forest_panel(
    "B", panel_b, "Remaining disability-free life expectancy",
    "Low minus high SES at age 60; 95% person-bootstrap intervals",
    "DFLE gap (years)", colors, "translate(370,10)"
  ),
  forest_panel(
    "C", panel_c, "Recovery and relapse contribution",
    "Combined contribution to the low–high DFLE gap; 95% intervals",
    "Recovery + relapse contribution (years)", colors, "translate(10,318)"
  ),
  make_panel_d("translate(370,318)")
)
write_svg(file.path(figure_dir, "Figure1_recovery_divide.svg"), 730, 626, composite)

manifest <- list(
  figure = "Figure 1",
  title = "The recovery divide",
  phase = phase,
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  source_tables = list.files(source_dir, full.names = FALSE),
  svg_files = c("Figure1_recovery_divide.svg", file.path("panels", list.files(panel_dir, pattern = "[.]svg$"))),
  contract = "FROZEN_ANALYSIS_CONTRACT_v1.0 with amendments v1.0.1-v1.0.3 and sparse-origin v1.1/bootstrap v1.2",
  rendering = "native SVG primitives and live Arial text; no raster objects"
)
write_json(manifest, file.path(figure_dir, "figure_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
caption <- c(
  "# Figure 1 caption",
  "",
  "**Figure 1 | The recovery divide in disability-free life expectancy.**",
  "",
  "**a,** Frozen interval-observed five-state framework. Independence before any observed difficulty (I0) transitions to first difficulty (D1), recovery after difficulty (R1), relapse after recovery (D2), and verified death; non-death attrition is censoring and never assigned to death. **b,** Modelled low-minus-high socioeconomic differences in remaining disability-free life expectancy at age 60 by cohort and education or wealth. **c,** Combined Shapley contribution of recovery and relapse to the low-minus-high disability-free life-expectancy gap. Points in b and c are frozen point estimates and horizontal lines are 95% person-cluster bootstrap percentile intervals; each bootstrap re-estimates wealth cut points, all transition models, transition matrices, life tables, and the decomposition. **d,** Five-block Shapley decomposition of the modelled low-minus-high gap into disability onset, recovery, relapse, post-disability mortality, and pre-disability mortality. Negative contributions widen the deficit among lower-SES adults. Shapley components are descriptive modelled replacements, not intervention effects or causal mediation. CHARLS, China Health and Retirement Longitudinal Study; ELSA, English Longitudinal Study of Ageing; HRS, Health and Retirement Study; MHAS, Mexican Health and Aging Study; DFLE, disability-free life expectancy."
)
writeLines(caption, file.path(figure_dir, "Figure1_caption.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_submission_figure_", phase, ".txt")))

cat("Submission tables and SVG figure completed\n")
cat("phase=", phase, "\n", sep = "")
cat("figure_dir=", figure_dir, "\n", sep = "")
cat("table_dir=", table_dir, "\n", sep = "")
