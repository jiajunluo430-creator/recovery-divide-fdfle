#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) tolower(args[[1L]]) else "preview"
stopifnot(phase %in% c("preview", "final"))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
point_dir <- file.path(root, "03_outputs", "11_exploratory_upgrade_point")
boot_dir <- file.path(root, "03_outputs", paste0("13_exploratory_bootstrap_summary_", phase))
support_dir <- file.path(root, "03_outputs", "10_exploratory_upgrade_support")
sensitivity_dir <- file.path(root, "03_outputs", "16_exploratory_sensitivity_summary")
table_dir <- file.path(root, "03_outputs", paste0("17_exploratory_submission_tables_", phase))
figure_dir <- file.path(root, if (phase == "final") "04_figures_v2" else "04_figures_v2_preview")
panel_dir <- file.path(figure_dir, "panels")
source_dir <- file.path(figure_dir, "source_data")
qc_dir <- file.path(figure_dir, "qc")
log_dir <- file.path(root, "06_logs")
for (directory in c(table_dir, figure_dir, panel_dir, source_dir, qc_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required <- c(
  file.path(boot_dir, "low_education_wealth_promotion_by_cohort.csv"),
  file.path(boot_dir, "recovery_phase_promotion_by_cohort.csv"),
  file.path(point_dir, "shapley_decomposition_point.csv"),
  file.path(support_dir, "joint_ses_process_support.csv"),
  file.path(sensitivity_dir, "low_education_wealth_sensitivity_robustness.csv")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) stop("Missing figure inputs: ", paste(missing_required, collapse = "; "))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
block_order <- c("onset", "recovery", "relapse", "post_disability_mortality", "pre_disability_mortality")
block_colors <- c(
  onset = "#4D8C80",
  recovery = "#D46A1F",
  relapse = "#7663A8",
  post_disability_mortality = "#B65A78",
  pre_disability_mortality = "#7C8B4A"
)

set_cohort_order <- function(z, extra = character()) {
  z[, .cohort_order := match(cohort, cohort_order)]
  setorderv(z, c(".cohort_order", extra))
  z[, .cohort_order := NULL]
  invisible(z)
}

lowedu <- fread(file.path(boot_dir, "low_education_wealth_promotion_by_cohort.csv"))
set_cohort_order(lowedu)
panel_b <- lowedu[, .(
  cohort,
  point = dfle_gap,
  ci_low = ci_low_gap_dfle,
  ci_high = ci_high_gap_dfle,
  valid_replicates = valid_replicates_gap_dfle,
  promoted = cohort_promotion_gate
)]

panel_c <- fread(file.path(point_dir, "shapley_decomposition_point.csv"))[
  module == "wealth_within_low_education" & estimand == "dfle",
  .(
    cohort, block,
    point = contribution_years,
    total_gap = total_low_minus_high_gap,
    percent = contribution_percent
  )
]
panel_c[, block := factor(block, levels = block_order)]
set_cohort_order(panel_c, "block")

panel_d <- fread(file.path(boot_dir, "recovery_phase_promotion_by_cohort.csv"))[
  stratum %in% c("early_recovery", "sustained_recovery"),
  .(
    cohort,
    phase = stratum,
    point = hazard_ratio,
    ci_low = bootstrap_hr_ci_low,
    ci_high = bootstrap_hr_ci_high,
    valid_replicates,
    interval_support = bootstrap_interval_support
  )
]
panel_d[, phase := factor(phase, levels = c("early_recovery", "sustained_recovery"))]
set_cohort_order(panel_d, "phase")

support <- fread(file.path(support_dir, "joint_ses_process_support.csv"))[
  joint_ses4 %in% c("low_low", "lowedu_highwealth") & process %in% c("recovery", "relapse")
]
support_wide <- dcast(
  support,
  cohort + joint_ses4 ~ process,
  value.var = "events",
  fill = 0
)
panel_a <- support_wide[, .(
  minimum_recovery_events = min(recovery),
  minimum_relapse_events = min(relapse)
), by = cohort]
set_cohort_order(panel_a)

sensitivity <- fread(file.path(sensitivity_dir, "low_education_wealth_sensitivity_robustness.csv"))
set_cohort_order(sensitivity)

fwrite(panel_a, file.path(source_dir, "Figure2A_low_education_support.csv"))
fwrite(panel_b, file.path(source_dir, "Figure2B_low_education_dfle_gap.csv"))
fwrite(panel_c, file.path(source_dir, "Figure2C_low_education_shapley_blocks.csv"))
fwrite(panel_d, file.path(source_dir, "Figure2D_recovery_phase_relapse_hr.csv"))
fwrite(sensitivity, file.path(source_dir, "Figure2_sensitivity_context.csv"))

fwrite(lowedu, file.path(table_dir, "Table3_low_education_wealth_mechanism.csv"))
fwrite(panel_c, file.path(table_dir, "TableS8_low_education_shapley_blocks.csv"))
fwrite(panel_d, file.path(table_dir, "TableS9_recovery_phase_relapse_hr.csv"))
fwrite(sensitivity, file.path(table_dir, "TableS10_low_education_sensitivity_robustness.csv"))

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

line_svg <- function(id, x1, y1, x2, y2, class = "axis", extra = "") {
  sprintf('<line id="%s" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" class="%s" %s/>', id, x1, y1, x2, y2, class, extra)
}

rect_svg <- function(id, x, y, width, height, fill, stroke = "none", rx = 0, class = "") {
  sprintf('<rect id="%s" x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="%.2f" fill="%s" stroke="%s" class="%s"/>', id, x, y, width, height, rx, fill, stroke, class)
}

circle_svg <- function(id, x, y, radius, fill, stroke = "white", stroke_width = 1) {
  sprintf('<circle id="%s" cx="%.2f" cy="%.2f" r="%.2f" fill="%s" stroke="%s" stroke-width="%.2f"/>', id, x, y, radius, fill, stroke, stroke_width)
}

text_svg <- function(id, x, y, label, class = "label", anchor = "start", extra = "") {
  sprintf('<text id="%s" x="%.2f" y="%.2f" class="%s" text-anchor="%s" %s>%s</text>', id, x, y, class, anchor, extra, xml_escape(label))
}

polygon_svg <- function(id, points, fill) {
  sprintf('<polygon id="%s" points="%s" fill="%s"/>', id, points, fill)
}

arrow_svg <- function(id, x1, y1, x2, y2, color = "#424A53", head_length = 6, head_half_width = 3.5) {
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
  c(
    line_svg(paste0(id, "_shaft"), x1, y1, x2, y2, "axis", sprintf('stroke="%s"', color)),
    polygon_svg(paste0(id, "_head"), points, color)
  )
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
      .subtitle { font-size: 8.5px; fill: #58616B; }
      .label { font-size: 9px; }
      .small { font-size: 7.5px; }
      .tiny { font-size: 6.8px; }
      .axis { stroke: #424A53; stroke-width: 1; fill: none; }
      .grid { stroke: #D9DEE3; stroke-width: 0.7; fill: none; }
      .zero { stroke: #20252B; stroke-width: 1.1; fill: none; }
      .ci { stroke-width: 1.8; fill: none; stroke-linecap: round; }
    ]]></style>',
    content,
    "</svg>"
  )
}

write_svg <- function(path, width, height, content) {
  writeLines(svg_document(width, height, content), path, useBytes = TRUE)
}

make_panel_a <- function(transform = NULL) {
  axis <- c(
    rect_svg("A_border", 5, 5, 340, 288, "none", "#D9DEE3", 2),
    line_svg("A_header_rule", 15, 47, 335, 47, "grid")
  )
  data <- c(
    rect_svg("A_low_education", 92, 61, 166, 30, "#E8EDF2", "#59636E", 6),
    rect_svg("A_high_wealth", 32, 111, 120, 38, "#DCEFE9", "#4D8C80", 6),
    rect_svg("A_low_wealth", 198, 111, 120, 38, "#FBE1D2", "#D46A1F", 6),
    arrow_svg("A_split_high", 150, 91, 92, 111, "#59636E"),
    arrow_svg("A_split_low", 200, 91, 258, 111, "#59636E")
  )
  process_x <- c(onset = 28, recovery = 105, relapse = 182, death = 259)
  process_fill <- c(onset = "#E8EDF2", recovery = "#F7D7C1", relapse = "#E7E0F1", death = "#E8EDF2")
  for (process_name in names(process_x)) {
    data <- c(data, rect_svg(
      paste0("A_process_", process_name), process_x[[process_name]], 173, 64, 25,
      process_fill[[process_name]], "#737C85", 4
    ))
  }
  labels <- c(
    text_svg("A_lowedu_text", 175, 80, "Frozen low-education stratum", "label", "middle"),
    text_svg("A_highwealth_text", 92, 135, "High wealth", "label", "middle"),
    text_svg("A_lowwealth_text", 258, 135, "Low wealth", "label", "middle"),
    text_svg("A_contrast_text", 175, 159, "reported contrast: low minus high", "small", "middle")
  )
  for (process_name in names(process_x)) {
    labels <- c(labels, text_svg(
      paste0("A_process_label_", process_name), process_x[[process_name]] + 32, 190,
      tools::toTitleCase(gsub("_", " ", process_name)), "small", "middle"
    ))
  }
  cohort_y <- c(220, 236, 252, 268)
  for (i in seq_along(cohort_order)) {
    z <- panel_a[cohort == cohort_order[[i]]]
    labels <- c(labels,
      text_svg(paste0("A_support_cohort_", i), 28, cohort_y[[i]], cohort_order[[i]], "small"),
      text_svg(
        paste0("A_support_counts_", i), 91, cohort_y[[i]],
        paste0("min recovery ", z$minimum_recovery_events, "  |  min relapse ", z$minimum_relapse_events),
        "small"
      )
    )
  }
  legend <- c(
    rect_svg("A_legend_recovery", 237, 211, 9, 7, block_colors[["recovery"]]),
    text_svg("A_legend_recovery_text", 250, 218, "recovery", "tiny"),
    rect_svg("A_legend_relapse", 292, 211, 9, 7, block_colors[["relapse"]]),
    text_svg("A_legend_relapse_text", 305, 218, "relapse", "tiny")
  )
  annotations <- c(
    text_svg("A_letter", 14, 28, "A", "panel-letter"),
    text_svg("A_title", 42, 27, "Wealth divide within low education", "title"),
    text_svg("A_subtitle", 42, 40, "Support was audited before outcome modelling", "subtitle")
  )
  panel_groups("A", axis, data, labels, legend, annotations, transform)
}

forest_panel <- function(letter, z, title, subtitle, x_title, color, zero = 0, log_scale = FALSE, transform = NULL) {
  z <- copy(z)
  set_cohort_order(z, if ("phase" %in% names(z)) "phase" else character())
  values <- c(z$ci_low, z$ci_high, zero)
  if (log_scale) {
    values <- log(values)
    x_min <- floor(min(values, na.rm = TRUE) * 5) / 5
    x_max <- ceiling(max(values, na.rm = TRUE) * 5) / 5
    xmap <- function(x) 90 + (log(x) - x_min) / (x_max - x_min) * 240
    ticks <- c(0.5, 0.75, 1, 1.5, 2)
    ticks <- ticks[log(ticks) >= x_min & log(ticks) <= x_max]
  } else {
    finite_range <- range(values, finite = TRUE)
    pad <- max(0.25, diff(finite_range) * 0.10)
    x_min <- floor((finite_range[[1]] - pad) * 2) / 2
    x_max <- ceiling((finite_range[[2]] + pad) * 2) / 2
    xmap <- function(x) 90 + (x - x_min) / (x_max - x_min) * 240
    tick_step <- if ((x_max - x_min) <= 5) 1 else 2
    ticks <- seq(ceiling(x_min / tick_step) * tick_step, floor(x_max / tick_step) * tick_step, by = tick_step)
  }
  plot_top <- 65
  plot_bottom <- 244
  cohort_y <- setNames(c(84, 126, 168, 210), cohort_order)
  offsets <- if ("phase" %in% names(z)) c(early_recovery = -6, sustained_recovery = 6) else c(all = 0)
  axis <- c(rect_svg(paste0(letter, "_border"), 5, 5, 340, 288, "none", "#D9DEE3", 2))
  labels <- character()
  for (i in seq_along(ticks)) {
    tick <- ticks[[i]]
    axis <- c(axis, line_svg(
      paste0(letter, "_grid_", i), xmap(tick), plot_top, xmap(tick), plot_bottom,
      if (abs(tick - zero) < 1e-10) "zero" else "grid"
    ))
    labels <- c(labels, text_svg(
      paste0(letter, "_tick_", i), xmap(tick), 258,
      if (log_scale) format(tick, trim = TRUE) else format(tick, trim = TRUE), "small", "middle"
    ))
  }
  for (i in seq_along(cohort_order)) {
    cohort_name <- cohort_order[[i]]
    axis <- c(axis, line_svg(
      paste0(letter, "_row_", i), 90, cohort_y[[cohort_name]] + 18, 330, cohort_y[[cohort_name]] + 18, "grid"
    ))
    labels <- c(labels, text_svg(
      paste0(letter, "_cohort_", i), 82, cohort_y[[cohort_name]] + 3, cohort_name, "label", "end"
    ))
  }
  axis <- c(axis, line_svg(paste0(letter, "_xaxis"), 90, plot_bottom, 330, plot_bottom, "axis"))
  labels <- c(labels, text_svg(paste0(letter, "_xtitle"), 210, 278, x_title, "small", "middle"))
  data <- character()
  for (i in seq_len(nrow(z))) {
    stratum_name <- if ("phase" %in% names(z)) as.character(z$phase[[i]]) else "all"
    yy <- cohort_y[[z$cohort[[i]]]] + offsets[[stratum_name]]
    point_color <- if (length(color) > 1L) color[[stratum_name]] else color[[1L]]
    data <- c(data,
      line_svg(paste0(letter, "_ci_", i), xmap(z$ci_low[[i]]), yy, xmap(z$ci_high[[i]]), yy, "ci", sprintf('stroke="%s"', point_color)),
      line_svg(paste0(letter, "_cap_l_", i), xmap(z$ci_low[[i]]), yy - 3, xmap(z$ci_low[[i]]), yy + 3, "ci", sprintf('stroke="%s"', point_color)),
      line_svg(paste0(letter, "_cap_h_", i), xmap(z$ci_high[[i]]), yy - 3, xmap(z$ci_high[[i]]), yy + 3, "ci", sprintf('stroke="%s"', point_color)),
      circle_svg(paste0(letter, "_point_", i), xmap(z$point[[i]]), yy, 4, point_color)
    )
  }
  legend <- character()
  if ("phase" %in% names(z)) {
    legend <- c(
      circle_svg(paste0(letter, "_legend_early"), 198, 51, 4, color[["early_recovery"]]),
      text_svg(paste0(letter, "_legend_early_text"), 207, 54, "Early", "small"),
      circle_svg(paste0(letter, "_legend_sustained"), 260, 51, 4, color[["sustained_recovery"]]),
      text_svg(paste0(letter, "_legend_sustained_text"), 269, 54, "Sustained", "small")
    )
  }
  annotations <- c(
    text_svg(paste0(letter, "_letter"), 14, 28, letter, "panel-letter"),
    text_svg(paste0(letter, "_title"), 42, 27, title, "title"),
    text_svg(paste0(letter, "_subtitle"), 42, 40, subtitle, "subtitle")
  )
  panel_groups(letter, axis, data, labels, legend, annotations, transform)
}

make_panel_c <- function(transform = NULL) {
  z <- copy(panel_c)
  x_min <- min(-6.5, floor(min(z[, sum(point[point < 0]), by = cohort]$V1) * 2) / 2)
  x_max <- max(1.5, ceiling(max(z[, sum(point[point > 0]), by = cohort]$V1) * 2) / 2)
  xmap <- function(x) 95 + (x - x_min) / (x_max - x_min) * 235
  cohort_y <- setNames(c(84, 126, 168, 210), cohort_order)
  ticks <- seq(ceiling(x_min), floor(x_max), by = 1)
  axis <- c(rect_svg("C_border", 5, 5, 340, 288, "none", "#D9DEE3", 2))
  labels <- character()
  for (i in seq_along(ticks)) {
    tick <- ticks[[i]]
    axis <- c(axis, line_svg(
      paste0("C_grid_", i), xmap(tick), 64, xmap(tick), 238,
      if (tick == 0) "zero" else "grid"
    ))
    labels <- c(labels, text_svg(paste0("C_tick_", i), xmap(tick), 253, tick, "small", "middle"))
  }
  data <- character()
  for (i in seq_along(cohort_order)) {
    cohort_name <- cohort_order[[i]]
    zz <- z[cohort == cohort_name]
    neg_cursor <- 0
    pos_cursor <- 0
    for (j in seq_along(block_order)) {
      block_name <- block_order[[j]]
      value <- zz[as.character(block) == block_name, point]
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
      data <- c(data, rect_svg(
        paste0("C_", cohort_name, "_", block_name),
        xmap(x0), cohort_y[[cohort_name]] - 7, max(0.5, xmap(x1) - xmap(x0)), 14,
        block_colors[[block_name]], "white"
      ))
    }
    labels <- c(labels,
      text_svg(paste0("C_cohort_", i), 87, cohort_y[[cohort_name]] + 3, cohort_name, "label", "end"),
      text_svg(
        paste0("C_gap_", i), 333, cohort_y[[cohort_name]] + 3,
        paste0("gap ", sprintf("%.1f", unique(zz$total_gap)), " y"), "tiny", "end"
      )
    )
  }
  axis <- c(axis, line_svg("C_xaxis", 95, 238, 330, 238, "axis"))
  labels <- c(labels, text_svg("C_xtitle", 212, 270, "Contribution to low–high DFLE gap (years)", "small", "middle"))
  legend_labels <- c("Onset", "Recovery", "Relapse", "Post-disability death", "Pre-disability death")
  legend <- character()
  legend_x <- c(17, 76, 145, 210, 17)
  legend_y <- c(279, 279, 279, 279, 291)
  for (i in seq_along(block_order)) {
    legend <- c(legend,
      rect_svg(paste0("C_legend_box_", i), legend_x[[i]], legend_y[[i]] - 7, 9, 7, block_colors[[block_order[[i]]]]),
      text_svg(paste0("C_legend_text_", i), legend_x[[i]] + 12, legend_y[[i]], legend_labels[[i]], "tiny")
    )
  }
  annotations <- c(
    text_svg("C_letter", 14, 28, "C", "panel-letter"),
    text_svg("C_title", 42, 27, "Which process generates the residual divide?", "title"),
    text_svg("C_subtitle", 42, 40, "Five-block Shapley decomposition within low education", "subtitle")
  )
  panel_groups("C", axis, data, labels, legend, annotations, transform)
}

panel_a_content <- make_panel_a()
panel_b_content <- forest_panel(
  "B", panel_b,
  "Disability-free years within low education",
  paste0("Low minus high wealth; ", if (phase == "final") "final" else "preview", " 95% person-bootstrap intervals"),
  "DFLE gap (years)", c(all = "#D46A1F")
)
panel_c_content <- make_panel_c()
phase_colors <- c(early_recovery = "#3D7EA6", sustained_recovery = "#7663A8")
panel_d_content <- forest_panel(
  "D", panel_d,
  "Relapse after observed recovery",
  paste0("Low/high wealth HR (>1: more relapse in low wealth); ", phase, " intervals"),
  "Low-versus-high relapse hazard ratio", phase_colors, zero = 1, log_scale = TRUE
)

write_svg(file.path(panel_dir, "Figure2A_low_education_contrast.svg"), 350, 298, panel_a_content)
write_svg(file.path(panel_dir, "Figure2B_low_education_dfle.svg"), 350, 298, panel_b_content)
write_svg(file.path(panel_dir, "Figure2C_low_education_decomposition.svg"), 350, 298, panel_c_content)
write_svg(file.path(panel_dir, "Figure2D_recovery_durability.svg"), 350, 298, panel_d_content)

composite <- c(
  make_panel_a("translate(10,10)"),
  forest_panel(
    "B", panel_b,
    "Disability-free years within low education",
    paste0("Low minus high wealth; ", if (phase == "final") "final" else "preview", " 95% person-bootstrap intervals"),
    "DFLE gap (years)", c(all = "#D46A1F"), transform = "translate(370,10)"
  ),
  make_panel_c("translate(10,318)"),
  forest_panel(
    "D", panel_d,
    "Relapse after observed recovery",
    paste0("Low/high wealth HR (>1: more relapse in low wealth); ", phase, " intervals"),
    "Low-versus-high relapse hazard ratio", phase_colors, zero = 1, log_scale = TRUE,
    transform = "translate(370,318)"
  )
)
write_svg(file.path(figure_dir, "Figure2_low_education_recovery_divide.svg"), 730, 626, composite)

manifest <- list(
  figure = "Figure 2",
  title = "The wealth recovery divide within low education",
  phase = phase,
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  source_tables = list.files(source_dir, full.names = FALSE),
  svg_files = c(
    "Figure2_low_education_recovery_divide.svg",
    file.path("panels", list.files(panel_dir, pattern = "[.]svg$"))
  ),
  contract = "Exploratory upgrade charter v2.0 with addenda v2.1-v2.2",
  rendering = "native SVG primitives with live Arial text; no image elements"
)
write_json(manifest, file.path(figure_dir, "figure_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

caption <- c(
  "# Figure 2 caption",
  "",
  "**Figure 2 | The wealth recovery divide among adults with low education.**",
  "",
  paste0(
    "**a,** The versioned C2 contrast compares frozen low versus high entry-wave wealth within the harmonised low-education group; values show the smaller recovery and relapse event counts across the two decisive wealth cells. ",
    "**b,** Modelled low-minus-high wealth differences in remaining disability-free life expectancy at age 60. ",
    "**c,** Five-block Shapley decomposition of the low-minus-high difference into onset, recovery, relapse, post-disability mortality, and pre-disability mortality. Negative contributions widen the low-wealth deficit. ",
    "**d,** Low-versus-high wealth hazard ratios for relapse immediately after an observed recovery and after recovery persisted for another observed wave. Points and horizontal lines in b and d are point estimates and ",
    if (phase == "final") "final 95% person-bootstrap percentile intervals." else "100-replicate preview 95% person-bootstrap percentile intervals; these are not final inference.",
    " Shapley components are descriptive modelled replacements, not intervention effects or causal mediation. All panels use the unchanged v1 functional state history. CHARLS, China Health and Retirement Longitudinal Study; ELSA, English Longitudinal Study of Ageing; HRS, Health and Retirement Study; MHAS, Mexican Health and Aging Study; DFLE, disability-free life expectancy; HR, hazard ratio."
  )
)
writeLines(caption, file.path(figure_dir, "Figure2_caption.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, paste0("sessionInfo_exploratory_figure_", phase, ".txt")))

cat("Exploratory submission tables and Figure 2 SVG completed\n")
cat("phase=", phase, "\n", sep = "")
cat("figure_dir=", figure_dir, "\n", sep = "")
