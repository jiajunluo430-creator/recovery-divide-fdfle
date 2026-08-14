#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
revision_root <- file.path(root, "09_manuscript_bmc_medicine_v1", "revision_round_1")
outputs <- file.path(revision_root, "03_outputs")
tables <- file.path(revision_root, "04_tables_v1_1")
figure_dir_override <- Sys.getenv("D4_BMC_FIGURE_DIR", unset = "")
figure_dir <- if (nzchar(figure_dir_override)) {
  normalizePath(figure_dir_override, winslash = "/", mustWork = FALSE)
} else {
  file.path(revision_root, "05_figures_v1_1")
}
panel_dir <- file.path(figure_dir, "panels")
source_dir <- file.path(figure_dir, "source_data")
qc_dir <- file.path(figure_dir, "qc")
log_dir <- if (nzchar(figure_dir_override)) file.path(figure_dir, "logs") else file.path(revision_root, "06_logs")
for (d in c(figure_dir, panel_dir, source_dir, qc_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
country_labels <- c(CHARLS = "China", ELSA = "England", HRS = "United States", MHAS = "Mexico")

colors <- list(
  education = "#2C6DA4",
  wealth = "#D56C1E",
  high = "#3178A8",
  low = "#C34E4E",
  difficulty_high = "#78A6C4",
  difficulty_low = "#D88B85",
  initial = "#6B7D8D",
  onset = "#4E9085",
  recovery = "#D56C1E",
  relapse = "#6F5AA8",
  post = "#B55B7A",
  pre = "#718146",
  broad = "#2678A8",
  adl = "#7B4FA3",
  two = "#D47A1F",
  death = "#747B82",
  censor = "#999FA5"
)

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

svg_text <- function(id, x, y, label, class = "label", anchor = "start", extra = "") {
  sprintf('<text id="%s" x="%.2f" y="%.2f" class="%s" text-anchor="%s" %s>%s</text>',
          id, x, y, class, anchor, extra, xml_escape(label))
}

svg_line <- function(id, x1, y1, x2, y2, class = "axis", extra = "") {
  sprintf('<line id="%s" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" class="%s" %s/>',
          id, x1, y1, x2, y2, class, extra)
}

svg_rect <- function(id, x, y, width, height, fill = "none", stroke = "none", rx = 0, class = "", extra = "") {
  sprintf('<rect id="%s" x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="%.2f" fill="%s" stroke="%s" class="%s" %s/>',
          id, x, y, width, height, rx, fill, stroke, class, extra)
}

svg_circle <- function(id, x, y, radius, fill, stroke = "white", stroke_width = 1.2) {
  sprintf('<circle id="%s" cx="%.2f" cy="%.2f" r="%.2f" fill="%s" stroke="%s" stroke-width="%.2f"/>',
          id, x, y, radius, fill, stroke, stroke_width)
}

svg_polygon <- function(id, points, fill, stroke = "none") {
  sprintf('<polygon id="%s" points="%s" fill="%s" stroke="%s"/>', id, points, fill, stroke)
}

svg_arrow <- function(id, x1, y1, x2, y2, color = "#525960", dashed = FALSE,
                      head_length = 9, head_half_width = 5) {
  dx <- x2 - x1
  dy <- y2 - y1
  distance <- sqrt(dx^2 + dy^2)
  ux <- dx / distance
  uy <- dy / distance
  base_x <- x2 - head_length * ux
  base_y <- y2 - head_length * uy
  perp_x <- -uy * head_half_width
  perp_y <- ux * head_half_width
  points <- sprintf("%.2f,%.2f %.2f,%.2f %.2f,%.2f",
                    x2, y2, base_x + perp_x, base_y + perp_y, base_x - perp_x, base_y - perp_y)
  dash <- if (dashed) 'stroke-dasharray="6 4"' else ""
  c(
    svg_line(paste0(id, "_shaft"), x1, y1, x2, y2, "transition", sprintf('stroke="%s" %s', color, dash)),
    svg_polygon(paste0(id, "_head"), points, color)
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

svg_document <- function(width, height, content, reference_width = width, final_width_mm = 170) {
  final_width_pt <- final_width_mm * 72 / 25.4
  style_scale <- reference_width / final_width_pt
  physical_width_mm <- width / reference_width * final_width_mm
  physical_height_mm <- height / reference_width * final_width_mm
  style <- sprintf(
    '<style><![CDATA[
      text { font-family: Arial, Helvetica, sans-serif; fill: #202428; }
      .panel-label { font-size: %.3fpx; font-weight: 700; }
      .facet-title { font-size: %.3fpx; font-weight: 700; }
      .label { font-size: %.3fpx; }
      .small { font-size: %.3fpx; }
      .tiny { font-size: %.3fpx; }
      .micro { font-size: %.3fpx; }
      .state-main { font-size: %.3fpx; font-weight: 700; }
      .state-sub { font-size: %.3fpx; fill: #525960; }
      .axis { stroke: #525960; stroke-width: %.3f; fill: none; }
      .grid { stroke: #DCE0E3; stroke-width: %.3f; fill: none; }
      .separator { stroke: #C8CDD1; stroke-width: %.3f; fill: none; }
      .zero { stroke: #252A2E; stroke-width: %.3f; fill: none; }
      .ci { stroke-width: %.3f; fill: none; stroke-linecap: round; }
      .transition { stroke-width: %.3f; fill: none; }
    ]]></style>',
    11.5 * style_scale,
    9.2 * style_scale,
    8.2 * style_scale,
    7.8 * style_scale,
    7.3 * style_scale,
    7.0 * style_scale,
    8.5 * style_scale,
    7.0 * style_scale,
    0.80 * style_scale,
    0.45 * style_scale,
    0.50 * style_scale,
    0.90 * style_scale,
    1.20 * style_scale,
    1.10 * style_scale
  )
  c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" width="%.3fmm" height="%.3fmm" ',
                   'viewBox="0 0 %d %d" data-final-width-mm="%.1f">'),
            physical_width_mm, physical_height_mm, width, height, final_width_mm),
    style,
    content,
    "</svg>"
  )
}

write_svg <- function(path, width, height, content, reference_width = width, final_width_mm = 170) {
  writeLines(
    svg_document(width, height, content, reference_width = reference_width, final_width_mm = final_width_mm),
    path,
    useBytes = TRUE
  )
}

facet_strip <- function(letter, title, width) {
  title_lines <- strsplit(title, "|", fixed = TRUE)[[1L]]
  has_letter <- nzchar(letter)
  title_x <- if (has_letter) 68 else 12
  strip_height <- if (length(title_lines) <= 2L) 72 else if (length(title_lines) == 3L) 112 else 144
  output <- c(svg_rect(paste0(if (has_letter) letter else "single", "_facet_strip"),
                       0, 0, width, strip_height, "#F0F1F2"))
  if (has_letter) {
    output <- c(output, svg_text(paste0(letter, "_panel_label"), 12, strip_height / 2 + 8,
                                 paste0("(", tolower(letter), ")"), "panel-label"))
  }
  if (length(title_lines) == 1L) {
    title_y <- strip_height / 2 + 7
  } else {
    title_y <- seq(25, strip_height - 25, length.out = length(title_lines))
  }
  for (i in seq_along(title_lines)) {
    output <- c(output, svg_text(
      paste0(if (has_letter) letter else "single", "_facet_title_", i),
      title_x, title_y[[i]], title_lines[[i]], "facet-title"
    ))
  }
  output
}

format_tick <- function(x) {
  ifelse(abs(x - round(x)) < 1e-8, as.character(round(x)), formatC(x, format = "f", digits = 1))
}

grouped_forest_panel <- function(letter, z, title, width, height, x_min, x_max, ticks,
                                 x_title, series_levels, series_colors, series_labels,
                                 log_scale = FALSE, left = 112, right = width - 18,
                                 top = 92, bottom = NULL, legend_rows = 1L,
                                 legend_title = NULL) {
  z <- copy(z)
  legend_title_offset <- if (is.null(legend_title)) 0 else 20
  if (is.null(bottom)) bottom <- height - 145 - 26 * (legend_rows - 1L) - legend_title_offset
  z[, cohort_order__ := match(cohort, cohort_order)]
  z[, series_order__ := match(series, series_levels)]
  setorder(z, cohort_order__, series_order__)
  xmap <- if (log_scale) {
    function(x) left + (log(x) - log(x_min)) / (log(x_max) - log(x_min)) * (right - left)
  } else {
    function(x) left + (x - x_min) / (x_max - x_min) * (right - left)
  }
  centers <- seq(top + 20, bottom - 20, length.out = length(cohort_order))
  names(centers) <- cohort_order
  offsets <- if (length(series_levels) == 2L) c(-7, 7) else seq(-12, 12, length.out = length(series_levels))
  names(offsets) <- series_levels
  axis <- labels <- data <- legend <- annotations <- character()
  annotations <- c(annotations, facet_strip(letter, title, width))
  for (i in seq_along(ticks)) {
    tick <- ticks[[i]]
    cls <- if ((log_scale && abs(tick - 1) < 1e-10) || (!log_scale && abs(tick) < 1e-10)) "zero" else "grid"
    axis <- c(axis, svg_line(paste0(letter, "_grid_", i), xmap(tick), top, xmap(tick), bottom, cls))
    labels <- c(labels, svg_text(paste0(letter, "_tick_", i), xmap(tick), bottom + 20, format_tick(tick), "small", "middle"))
  }
  axis <- c(axis, svg_line(paste0(letter, "_xaxis"), left, bottom, right, bottom, "axis"))
  for (i in seq_along(cohort_order)) {
    cc <- cohort_order[[i]]
    labels <- c(labels, svg_text(paste0(letter, "_cohort_", cc), left - 10, centers[[cc]] + 4, cc, "label", "end"))
    if (i < length(cohort_order)) {
      sep <- mean(c(centers[[cc]], centers[[cohort_order[[i + 1L]]]]))
      axis <- c(axis, svg_line(paste0(letter, "_sep_", i), left, sep, right, sep, "separator"))
    }
  }
  for (i in seq_len(nrow(z))) {
    yy <- centers[[z$cohort[[i]]]] + offsets[[z$series[[i]]]]
    col <- series_colors[[z$series[[i]]]]
    data <- c(data,
              svg_line(paste0(letter, "_ci_", i), xmap(z$ci_low[[i]]), yy, xmap(z$ci_high[[i]]), yy,
                       "ci", sprintf('stroke="%s"', col)),
              svg_line(paste0(letter, "_cap_l_", i), xmap(z$ci_low[[i]]), yy - 3.5, xmap(z$ci_low[[i]]), yy + 3.5,
                       "ci", sprintf('stroke="%s"', col)),
              svg_line(paste0(letter, "_cap_h_", i), xmap(z$ci_high[[i]]), yy - 3.5, xmap(z$ci_high[[i]]), yy + 3.5,
                       "ci", sprintf('stroke="%s"', col)),
              svg_circle(paste0(letter, "_point_", i), xmap(z$point[[i]]), yy, 4.5, col))
  }
  labels <- c(labels, svg_text(paste0(letter, "_xtitle"), (left + right) / 2,
                               height - 90 - 26 * (legend_rows - 1L) - legend_title_offset,
                               x_title, "small", "middle"))
  if (!is.null(legend_title)) {
    legend <- c(legend, svg_text(paste0(letter, "_legend_title"), left,
                                 height - 68 - 24 * (legend_rows - 1L),
                                 legend_title, "micro", extra = 'font-weight="700"'))
  }
  legend_cols <- ceiling(length(series_levels) / legend_rows)
  legend_x <- seq(left, left + (right - left) * (legend_cols - 1L) / legend_cols,
                  length.out = legend_cols)
  for (i in seq_along(series_levels)) {
    ss <- series_levels[[i]]
    legend_col <- (i - 1L) %% legend_cols + 1L
    legend_row <- (i - 1L) %/% legend_cols + 1L
    legend_y <- height - 42 - (legend_rows - legend_row) * 24
    legend <- c(legend,
                svg_circle(paste0(letter, "_legend_point_", i), legend_x[[legend_col]], legend_y, 4.5, series_colors[[ss]]),
                svg_text(paste0(letter, "_legend_text_", i), legend_x[[legend_col]] + 11,
                         legend_y + 4, series_labels[[ss]], "tiny"))
  }
  list(axis = axis, data = data, labels = labels, legend = legend, annotations = annotations)
}

long_forest_panel <- function(letter, z, title, width, height, x_min, x_max, ticks, x_title,
                              series_levels, series_colors, series_labels, log_scale = FALSE,
                              left = 168, right = width - 18, top = 88, bottom = NULL,
                              show_letter = TRUE, legend_rows = 1L,
                              show_cohort_column = FALSE, legend_title = NULL) {
  z <- copy(z)
  legend_title_offset <- if (is.null(legend_title)) 0 else 20
  if (is.null(bottom)) bottom <- height - 164 - 26 * (legend_rows - 1L) - legend_title_offset
  z[, row_order__ := seq_len(.N)]
  xmap <- if (log_scale) {
    function(x) left + (log(x) - log(x_min)) / (log(x_max) - log(x_min)) * (right - left)
  } else {
    function(x) left + (x - x_min) / (x_max - x_min) * (right - left)
  }
  yy <- seq(top + 12, bottom - 14, length.out = nrow(z))
  axis <- labels <- data <- legend <- annotations <- character()
  annotations <- c(annotations, facet_strip(if (show_letter) letter else "", title, width))
  if (show_cohort_column) {
    cohort_values <- unique(as.character(z$cohort))
    for (j in seq_along(cohort_values)) {
      cc <- cohort_values[[j]]
      idx <- which(as.character(z$cohort) == cc)
      if (j %% 2L == 0L) {
        band_top <- if (min(idx) == 1L) top else mean(c(yy[[min(idx) - 1L]], yy[[min(idx)]]))
        band_bottom <- if (max(idx) == nrow(z)) bottom else mean(c(yy[[max(idx)]], yy[[max(idx) + 1L]]))
        axis <- c(axis, svg_rect(paste0(letter, "_cohort_band_", cc), 0, band_top,
                                 width, band_bottom - band_top, "#F7F8F9"))
      }
      group_y <- mean(yy[idx])
      labels <- c(labels, svg_text(
        paste0(letter, "_cohort_group_", cc), 28, group_y, cc, "small", "middle",
        extra = sprintf('font-weight="700" transform="rotate(-90 28 %.2f)"', group_y)
      ))
    }
    axis <- c(axis, svg_line(paste0(letter, "_cohort_column_rule"), 60, top, 60, bottom, "separator"))
  }
  for (i in seq_along(ticks)) {
    tick <- ticks[[i]]
    cls <- if ((log_scale && abs(tick - 1) < 1e-10) || (!log_scale && abs(tick) < 1e-10)) "zero" else "grid"
    axis <- c(axis, svg_line(paste0(letter, "_grid_", i), xmap(tick), top, xmap(tick), bottom, cls))
    labels <- c(labels, svg_text(paste0(letter, "_tick_", i), xmap(tick), bottom + 17, format_tick(tick), "micro", "middle"))
  }
  axis <- c(axis, svg_line(paste0(letter, "_xaxis"), left, bottom, right, bottom, "axis"))
  for (i in seq_len(nrow(z))) {
    if (i > 1L && z$cohort[[i]] != z$cohort[[i - 1L]]) {
      axis <- c(axis, svg_line(paste0(letter, "_group_sep_", i), 8, mean(c(yy[[i - 1L]], yy[[i]])), right,
                               mean(c(yy[[i - 1L]], yy[[i]])), "separator"))
    }
    row_text <- if (show_cohort_column) z$row_label[[i]] else paste0(z$cohort[[i]], "  ", z$row_label[[i]])
    labels <- c(labels, svg_text(paste0(letter, "_row_", i), left - 8, yy[[i]] + 3,
                                row_text, "micro", "end"))
    col <- series_colors[[z$series[[i]]]]
    data <- c(data,
              svg_line(paste0(letter, "_ci_", i), xmap(z$ci_low[[i]]), yy[[i]], xmap(z$ci_high[[i]]), yy[[i]],
                       "ci", sprintf('stroke="%s"', col)),
              svg_line(paste0(letter, "_cap_l_", i), xmap(z$ci_low[[i]]), yy[[i]] - 3, xmap(z$ci_low[[i]]), yy[[i]] + 3,
                       "ci", sprintf('stroke="%s"', col)),
              svg_line(paste0(letter, "_cap_h_", i), xmap(z$ci_high[[i]]), yy[[i]] - 3, xmap(z$ci_high[[i]]), yy[[i]] + 3,
                       "ci", sprintf('stroke="%s"', col)),
              svg_circle(paste0(letter, "_point_", i), xmap(z$point[[i]]), yy[[i]], 3.6, col))
  }
  labels <- c(labels, svg_text(paste0(letter, "_xtitle"), (left + right) / 2,
                               height - 110 - 26 * (legend_rows - 1L) - legend_title_offset,
                               x_title, "micro", "middle"))
  if (!is.null(legend_title)) {
    legend <- c(legend, svg_text(paste0(letter, "_legend_title"), 10,
                                 height - 68 - 24 * (legend_rows - 1L),
                                 legend_title, "micro", extra = 'font-weight="700"'))
  }
  legend_cols <- ceiling(length(series_levels) / legend_rows)
  legend_x <- seq(10, width - 150, length.out = legend_cols)
  for (i in seq_along(series_levels)) {
    ss <- series_levels[[i]]
    legend_col <- (i - 1L) %% legend_cols + 1L
    legend_row <- (i - 1L) %/% legend_cols + 1L
    legend_y <- height - 50 - (legend_rows - legend_row) * 24
    legend <- c(legend, svg_rect(paste0(letter, "_legend_box_", i), legend_x[[legend_col]], legend_y, 9, 7,
                                 series_colors[[ss]]),
                svg_text(paste0(letter, "_legend_text_", i), legend_x[[legend_col]] + 13, legend_y + 7,
                         series_labels[[ss]], "micro"))
  }
  list(axis = axis, data = data, labels = labels, legend = legend, annotations = annotations)
}

# ------------------------------
# Locked aggregate inputs
# ------------------------------

required_files <- c(
  table1 = file.path(tables, "source_data", "Table1_numeric_source.csv"),
  education = file.path(outputs, "04_primary_education_inference", "bootstrap_percentile_ci_summary.csv"),
  wealth = file.path(outputs, "04_primary_wealth_inference", "bootstrap_percentile_ci_summary.csv"),
  focused = file.path(outputs, "04_lowedu_inference", "bootstrap_percentile_ci_summary.csv"),
  adl = file.path(outputs, "04_lowedu_inference_adl_only", "bootstrap_percentile_ci_summary.csv"),
  two = file.path(outputs, "04_lowedu_inference_at_least_two", "bootstrap_percentile_ci_summary.csv"),
  continuous = file.path(outputs, "15_consolidated_evidence", "continuous_time_process_intensity_ratios.csv")
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required)) {
  stop("Missing final figure input(s): ", paste(names(missing_required), missing_required, collapse = " | "))
}

table1 <- fread(required_files[["table1"]])
education <- fread(required_files[["education"]])
wealth <- fread(required_files[["wealth"]])
focused <- fread(required_files[["focused"]])
adl <- fread(required_files[["adl"]])
two <- fread(required_files[["two"]])
continuous <- fread(required_files[["continuous"]])

# ------------------------------
# Figure 1: flow plus state model
# ------------------------------

fig1_flow <- melt(
  table1[, .(
    cohort,
    `Eligible entrants` = eligible_entrants,
    `Primary education` = primary_education_people,
    `Primary wealth` = primary_wealth_people,
    `Focused secondary` = focused_lowedu_wealth_people
  )],
  id.vars = "cohort", variable.name = "stage", value.name = "people"
)
fig1_flow[, cohort := factor(cohort, levels = cohort_order)]
fig1_flow[, stage := factor(stage, levels = c("Eligible entrants", "Primary education", "Primary wealth", "Focused secondary"))]
fwrite(fig1_flow, file.path(source_dir, "Figure1A_cohort_flow.csv"))

flow_colors <- c(
  `Eligible entrants` = "#BFC7CD",
  `Primary education` = colors$education,
  `Primary wealth` = colors$wealth,
  `Focused secondary` = "#6F5AA8"
)

fig1a_axis <- fig1a_data <- fig1a_labels <- fig1a_legend <- fig1a_annotations <- character()
fig1a_annotations <- c(fig1a_annotations, facet_strip("A", "Cohort entry and analysis samples", 560))
x_left <- 112
x_right <- 540
x_max <- 35000
xmap_flow <- function(x) x_left + x / x_max * (x_right - x_left)
ticks_flow <- c(0, 10000, 20000, 30000)
for (i in seq_along(ticks_flow)) {
  fig1a_axis <- c(fig1a_axis,
                  svg_line(paste0("A_grid_", i), xmap_flow(ticks_flow[[i]]), 90,
                           xmap_flow(ticks_flow[[i]]), 365, if (i == 1L) "zero" else "grid"))
  fig1a_labels <- c(fig1a_labels,
                    svg_text(paste0("A_tick_", i), xmap_flow(ticks_flow[[i]]), 386,
                             formatC(ticks_flow[[i]] / 1000, format = "f", digits = 0), "small", "middle"))
}
fig1a_axis <- c(fig1a_axis, svg_line("A_xaxis", x_left, 365, x_right, 365, "axis"))
cohort_y <- setNames(c(120, 185, 250, 315), cohort_order)
bar_offsets <- setNames(c(-15, -5, 5, 15), levels(fig1_flow$stage))
for (cc in cohort_order) {
  fig1a_labels <- c(fig1a_labels, svg_text(paste0("A_cohort_", cc), x_left - 10, cohort_y[[cc]] + 4,
                                            cc, "label", "end"))
}
for (i in seq_len(nrow(fig1_flow))) {
  cc <- as.character(fig1_flow$cohort[[i]])
  ss <- as.character(fig1_flow$stage[[i]])
  yy <- cohort_y[[cc]] + bar_offsets[[ss]]
  xx <- xmap_flow(fig1_flow$people[[i]])
  fig1a_data <- c(fig1a_data,
                  svg_rect(paste0("A_bar_", i), x_left, yy - 3.4, max(0.5, xx - x_left), 6.8,
                           flow_colors[[ss]], "none", 1.5))
  if (ss == "Eligible entrants") {
    fig1a_labels <- c(fig1a_labels, svg_text(paste0("A_n_", cc), min(xx + 6, 528), yy + 3,
                                              format(fig1_flow$people[[i]], big.mark = ","), "micro"))
  }
}
fig1a_labels <- c(fig1a_labels, svg_text("A_xtitle", (x_left + x_right) / 2, 408,
                                          "Participants (thousands)", "small", "middle"))
legend_x <- c(20, 300, 20, 300)
legend_y <- c(425, 425, 452, 452)
for (i in seq_along(levels(fig1_flow$stage))) {
  ss <- levels(fig1_flow$stage)[[i]]
  fig1a_legend <- c(fig1a_legend,
                    svg_rect(paste0("A_legend_box_", i), legend_x[[i]], legend_y[[i]], 10, 8, flow_colors[[ss]]),
                    svg_text(paste0("A_legend_text_", i), legend_x[[i]] + 14, legend_y[[i]] + 8, ss, "micro"))
}
panel1a <- panel_groups("A", fig1a_axis, fig1a_data, fig1a_labels, fig1a_legend, fig1a_annotations)
write_svg(file.path(panel_dir, "Figure1A_cohort_flow.svg"), 560, 470, panel1a, reference_width = 1340)

fig1b_axis <- character()
fig1b_data <- character()
fig1b_labels <- character()
fig1b_legend <- character()
fig1b_annotations <- facet_strip("B", "History-aware states, death, and censoring", 760)
state_x <- c(I0 = 28, D1 = 190, R1 = 370, D2 = 550)
state_fill <- c(I0 = "#E6F0EE", D1 = "#F4E2D4", R1 = "#F4E2D4", D2 = "#E9E3F3")
for (ss in names(state_x)) {
  fig1b_data <- c(fig1b_data, svg_rect(paste0("B_state_", ss), state_x[[ss]], 100, 142, 98,
                                        state_fill[[ss]], "#525960", 8))
}
fig1b_data <- c(fig1b_data,
                svg_rect("B_death", 292, 292, 142, 70, "#E2E5E7", "#525960", 8),
                svg_rect("B_censor", 530, 295, 190, 68, "#F4F5F6", "#8D949A", 6,
                         extra = 'stroke-dasharray="6 4"'))
fig1b_labels <- c(fig1b_labels,
                  svg_text("B_I0_main", 99, 130, "Independent", "state-main", "middle"),
                  svg_text("B_I0_sub_1", 99, 156, "before observed", "state-sub", "middle"),
                  svg_text("B_I0_sub_2", 99, 178, "difficulty", "state-sub", "middle"),
                  svg_text("B_D1_main", 261, 137, "Difficulty", "state-main", "middle"),
                  svg_text("B_D1_sub", 261, 166, "D1", "state-sub", "middle"),
                  svg_text("B_R1_main", 441, 137, "Recovered", "state-main", "middle"),
                  svg_text("B_R1_sub", 441, 166, "R1", "state-sub", "middle"),
                  svg_text("B_D2_main", 621, 137, "Relapsed", "state-main", "middle"),
                  svg_text("B_D2_sub", 621, 166, "D2", "state-sub", "middle"),
                  svg_text("B_death_text", 363, 333, "Verified death", "state-main", "middle"),
                  svg_text("B_censor_main", 625, 316, "Non-death loss", "state-sub", "middle"),
                  svg_text("B_censor_mid", 625, 338, "or unknown status", "state-sub", "middle"),
                  svg_text("B_censor_sub", 625, 358, "censored — never death", "state-sub", "middle"))
fig1b_data <- c(fig1b_data,
                svg_arrow("B_onset", 170, 146, 190, 146, colors$onset),
                svg_arrow("B_recovery_first", 332, 146, 370, 146, colors$recovery),
                svg_arrow("B_relapse", 512, 132, 550, 132, colors$relapse),
                svg_arrow("B_recovery_repeat", 550, 167, 512, 167, colors$recovery),
                svg_arrow("B_death_I0", 99, 198, 314, 292, colors$death),
                svg_arrow("B_death_D1", 261, 198, 340, 292, colors$death),
                svg_arrow("B_death_R1", 441, 198, 386, 292, colors$death),
                svg_arrow("B_death_D2", 621, 198, 414, 292, colors$death),
                svg_arrow("B_censor_arrow", 495, 228, 558, 295, colors$censor, dashed = TRUE))
fig1b_annotations <- c(fig1b_annotations,
                       svg_text("B_predeath_note", 170, 247, "pre-difficulty death block", "micro", "middle"),
                       svg_text("B_postdeath_note", 430, 261, "shared post-difficulty death block", "micro", "middle"),
                       svg_text("B_recovery_note", 441, 218, "D1/D2 share recovery contrast", "micro", "middle"))
legend_specs <- list(
  c("Onset", colors$onset), c("Recovery", colors$recovery), c("Relapse", colors$relapse),
  c("Mortality", colors$death), c("Censoring", colors$censor)
)
legend_x2 <- c(25, 155, 300, 435, 580)
for (i in seq_along(legend_specs)) {
  fig1b_legend <- c(fig1b_legend,
                    svg_line(paste0("B_legend_line_", i), legend_x2[[i]], 438, legend_x2[[i]] + 25, 438,
                             "transition", sprintf('stroke="%s"', legend_specs[[i]][[2]])),
                    svg_text(paste0("B_legend_text_", i), legend_x2[[i]] + 31, 442,
                             legend_specs[[i]][[1]], "micro"))
}
panel1b <- panel_groups("B", fig1b_axis, fig1b_data, fig1b_labels, fig1b_legend, fig1b_annotations)
write_svg(file.path(panel_dir, "Figure1B_state_framework.svg"), 760, 470, panel1b, reference_width = 1340)

figure1 <- c(
  panel_groups("A", fig1a_axis, fig1a_data, fig1a_labels, fig1a_legend, fig1a_annotations,
               transform = "translate(0,0)"),
  panel_groups("B", fig1b_axis, fig1b_data, fig1b_labels, fig1b_legend, fig1b_annotations,
               transform = "translate(580,0)")
)
write_svg(file.path(figure_dir, "Figure1_flow_and_multistate_framework.svg"), 1340, 470, figure1)

# ------------------------------
# Figure 2: primary gaps and focused absolute levels
# ------------------------------

fig2a <- rbindlist(list(
  education[metric == "population_fdfle_gap", .(
    cohort, series = "education", point = point_estimate, ci_low, ci_high
  )],
  wealth[metric == "population_fdfle_gap", .(
    cohort, series = "wealth", point = point_estimate, ci_low, ci_high
  )]
))
fig2a[, cohort := factor(cohort, levels = cohort_order)]
fwrite(fig2a, file.path(source_dir, "Figure2A_primary_education_wealth_fdfle_gaps.csv"))

fig2b_specs <- data.table(
  metric = c("population_high_fdfle", "population_low_fdfle",
             "population_high_difficulty_years", "population_low_difficulty_years"),
  series = c("high_fdfle", "low_fdfle", "high_difficulty", "low_difficulty")
)
fig2b <- focused[fig2b_specs, on = "metric", nomatch = 0,
                 .(cohort, series = i.series, point = point_estimate, ci_low, ci_high)]
fig2b[, cohort := factor(cohort, levels = cohort_order)]
fwrite(fig2b, file.path(source_dir, "Figure2B_focused_absolute_life_years.csv"))

p2a <- grouped_forest_panel(
  "A", fig2a, "Primary population FDFLE differences", 650, 480,
  x_min = -11, x_max = 1, ticks = c(-10, -5, 0),
  x_title = "Low minus high socioeconomic group (years)",
  series_levels = c("education", "wealth"),
  series_colors = c(education = colors$education, wealth = colors$wealth),
  series_labels = c(education = "Education", wealth = "Wealth"),
  left = 115, right = 628
)
panel2a <- panel_groups("A", p2a$axis, p2a$data, p2a$labels, p2a$legend, p2a$annotations)
write_svg(file.path(panel_dir, "Figure2A_primary_fdfle_gaps.svg"), 650, 480, panel2a, reference_width = 1320)

p2b <- grouped_forest_panel(
  "B", fig2b, "Focused secondary absolute years|within lower education", 650, 480,
  x_min = 0, x_max = 22, ticks = c(0, 5, 10, 15, 20),
  x_title = "Truncated remaining years from age 60 through 100",
  series_levels = c("high_fdfle", "low_fdfle", "high_difficulty", "low_difficulty"),
  series_colors = c(
    high_fdfle = colors$high, low_fdfle = colors$low,
    high_difficulty = colors$difficulty_high, low_difficulty = colors$difficulty_low
  ),
  series_labels = c(
    high_fdfle = "High wealth FDFLE", low_fdfle = "Low wealth FDFLE",
    high_difficulty = "High wealth difficulty", low_difficulty = "Low wealth difficulty"
  ),
  left = 115, right = 628, legend_rows = 2
)
panel2b <- panel_groups("B", p2b$axis, p2b$data, p2b$labels, p2b$legend, p2b$annotations)
write_svg(file.path(panel_dir, "Figure2B_focused_absolute_years.svg"), 650, 480, panel2b, reference_width = 1320)

figure2 <- c(
  panel_groups("A", p2a$axis, p2a$data, p2a$labels, p2a$legend, p2a$annotations,
               transform = "translate(0,0)"),
  panel_groups("B", p2b$axis, p2b$data, p2b$labels, p2b$legend, p2b$annotations,
               transform = "translate(670,0)")
)
write_svg(file.path(figure_dir, "Figure2_primary_gaps_and_absolute_years.svg"), 1320, 480, figure2)

# ------------------------------
# Figure 3: primary post-onset summary, decomposition, definition sensitivity
# ------------------------------

fig3a <- rbindlist(list(
  education[metric == "population_recovery_relapse", .(
    cohort, series = "education", point = point_estimate, ci_low, ci_high
  )],
  wealth[metric == "population_recovery_relapse", .(
    cohort, series = "wealth", point = point_estimate, ci_low, ci_high
  )]
))
fig3a[, cohort_order__ := match(cohort, cohort_order)]
fig3a[, series_order__ := match(series, c("education", "wealth"))]
setorder(fig3a, cohort_order__, series_order__)
fig3a[, c("cohort_order__", "series_order__") := NULL]
fwrite(fig3a, file.path(source_dir, "Figure3A_primary_post_onset_contributions.csv"))

component_map <- data.table(
  metric = c(
    "population_contribution_initial_state",
    "population_contribution_onset",
    "population_contribution_recovery",
    "population_contribution_relapse",
    "population_contribution_post_difficulty_mortality",
    "population_contribution_pre_difficulty_mortality"
  ),
  series = c("initial", "onset", "recovery", "relapse", "post", "pre"),
  row_label = c("Initial state", "Onset", "Recovery", "Relapse", "Post-diff. death", "Pre-diff. death"),
  component_order = 1:6
)
fig3b <- focused[component_map, on = "metric", nomatch = 0,
                 .(cohort, series = i.series, row_label = i.row_label,
                   component_order = i.component_order, point = point_estimate, ci_low, ci_high)]
fig3b[, cohort_order__ := match(cohort, cohort_order)]
setorder(fig3b, cohort_order__, component_order)
fig3b[, cohort_order__ := NULL]
fwrite(fig3b, file.path(source_dir, "Figure3B_six_block_decomposition_with_ci.csv"))

fig3c <- rbindlist(list(
  focused[metric == "population_recovery_relapse", .(
    cohort, series = "broad", row_label = "Broad (any)", point = point_estimate, ci_low, ci_high
  )],
  adl[metric == "population_recovery_relapse", .(
    cohort, series = "adl", row_label = "ADL only", point = point_estimate, ci_low, ci_high
  )],
  two[metric == "population_recovery_relapse", .(
    cohort, series = "two", row_label = "At least two", point = point_estimate, ci_low, ci_high
  )]
))
fig3c[, cohort_order__ := match(cohort, cohort_order)]
fig3c[, series_order__ := match(series, c("broad", "adl", "two"))]
setorder(fig3c, cohort_order__, series_order__)
fig3c[, c("cohort_order__", "series_order__") := NULL]
fwrite(fig3c, file.path(source_dir, "Figure3C_functional_definition_sensitivity.csv"))

process_map <- data.table(
  process = c("onset", "recovery", "relapse", "death_post", "death_pre"),
  series = c("onset", "recovery", "relapse", "post", "pre"),
  row_label = c("Onset", "Recovery", "Relapse", "Post-diff. death", "Pre-diff. death"),
  process_order = 1:5
)
figs1 <- continuous[process_map, on = "process", nomatch = 0,
                    .(cohort, series = i.series, row_label = i.row_label,
                      process_order = i.process_order, point = intensity_ratio,
                      ci_low, ci_high)]
figs1[, cohort_order__ := match(cohort, cohort_order)]
setorder(figs1, cohort_order__, process_order)
figs1[, cohort_order__ := NULL]
fwrite(figs1, file.path(source_dir, "AdditionalFile2_continuous_time_intensity_ratios.csv"))

p3a <- grouped_forest_panel(
  "A", fig3a, "Primary contributions", 440, 900,
  x_min = -5.5, x_max = 1.6, ticks = c(-5, -3, -1, 0, 1),
  x_title = "Recovery + relapse contribution (years)",
  series_levels = c("education", "wealth"),
  series_colors = c(education = colors$education, wealth = colors$wealth),
  series_labels = c(education = "Education", wealth = "Wealth"),
  left = 92, right = 425, top = 88,
  legend_title = "Socioeconomic indicator"
)
panel3a <- panel_groups("A", p3a$axis, p3a$data, p3a$labels, p3a$legend, p3a$annotations)
write_svg(file.path(panel_dir, "Figure3A_primary_post_onset.svg"), 440, 900, panel3a, reference_width = 1480)

p3b <- long_forest_panel(
  "B", fig3b, "Six-component decomposition", 560, 900,
  x_min = -3.2, x_max = 1.6, ticks = c(-3, -2, -1, 0, 1),
  x_title = "Contribution to FDFLE gap (years)",
  series_levels = c("initial", "onset", "recovery", "relapse", "post", "pre"),
  series_colors = c(
    initial = colors$initial, onset = colors$onset, recovery = colors$recovery,
    relapse = colors$relapse, post = colors$post, pre = colors$pre
  ),
  series_labels = c(
    initial = "Initial", onset = "Onset", recovery = "Recovery",
    relapse = "Relapse", post = "Post death", pre = "Pre death"
  ),
  left = 170, right = 542, top = 88, legend_rows = 2,
  show_cohort_column = TRUE, legend_title = "Decomposition component"
)
panel3b <- panel_groups("B", p3b$axis, p3b$data, p3b$labels, p3b$legend, p3b$annotations)
write_svg(file.path(panel_dir, "Figure3B_six_block_decomposition.svg"), 560, 900, panel3b, reference_width = 1480)

p3c <- long_forest_panel(
  "C", fig3c, "By outcome definition", 440, 900,
  x_min = -3.2, x_max = 2.2, ticks = c(-3, -2, -1, 0, 1, 2),
  x_title = "Recovery + relapse (years)",
  series_levels = c("broad", "adl", "two"),
  series_colors = c(broad = colors$broad, adl = colors$adl, two = colors$two),
  series_labels = c(broad = "Broad", adl = "ADL-only", two = "≥2 difficulties"),
  left = 155, right = 425, top = 88,
  show_cohort_column = TRUE, legend_title = "Outcome definition"
)
panel3c <- panel_groups("C", p3c$axis, p3c$data, p3c$labels, p3c$legend, p3c$annotations)
write_svg(file.path(panel_dir, "Figure3C_definition_sensitivity.svg"), 440, 900, panel3c, reference_width = 1480)

ps1 <- long_forest_panel(
  "A", figs1, "Continuous-time process contrasts|Focused secondary wealth contrast within lower education", 520, 720,
  x_min = 0.45, x_max = 3.0, ticks = c(0.5, 1, 2, 3),
  x_title = "Low-versus-high wealth intensity ratio (log scale)",
  series_levels = c("onset", "recovery", "relapse", "post", "pre"),
  series_colors = c(
    onset = colors$onset, recovery = colors$recovery, relapse = colors$relapse,
    post = colors$post, pre = colors$pre
  ),
  series_labels = c(
    onset = "Onset", recovery = "Recovery", relapse = "Relapse",
    post = "Post death", pre = "Pre death"
  ),
  log_scale = TRUE,
  left = 190,
  show_letter = FALSE
)
panels1 <- panel_groups("A", ps1$axis, ps1$data, ps1$labels, ps1$legend, ps1$annotations)
write_svg(file.path(figure_dir, "Additional_file_2_continuous_time_processes.svg"), 520, 720, panels1)

figure3 <- c(
  panel_groups("A", p3a$axis, p3a$data, p3a$labels, p3a$legend, p3a$annotations,
               transform = "translate(0,0)"),
  panel_groups("B", p3b$axis, p3b$data, p3b$labels, p3b$legend, p3b$annotations,
               transform = "translate(460,0)"),
  panel_groups("C", p3c$axis, p3c$data, p3c$labels, p3c$legend, p3c$annotations,
               transform = "translate(1040,0)")
)
write_svg(file.path(figure_dir, "Figure3_decomposition_and_robustness.svg"), 1480, 900, figure3)

figure_legends <- c(
  "# Main figure titles and legends",
  "",
  "## Figure 1. Cohort flow and history-aware multistate framework",
  "",
  "Panel A shows cohort entry and final participant counts for primary education, primary wealth, and the focused secondary low-education wealth analysis. Panel B shows the observed-history states and fitted transition blocks. D1 and D2 share the recovery socioeconomic contrast; D1, R1, and D2 share the post-difficulty mortality contrast. Verified death is absorbing. Non-death loss and unknown vital or observation status are censoring branches and were never coded as death.",
  "",
  "## Figure 2. Primary socioeconomic differences and focused absolute life years",
  "",
  "Panel A shows low-minus-high education and wealth differences in population-initialised functional-difficulty-free life expectancy across the four cohorts. Panel B shows focused secondary high- and low-wealth functional-difficulty-free years and years with functional difficulty within lower education. Points and intervals are household-cluster percentile-bootstrap estimates and 95% confidence intervals. Negative differences indicate fewer functional-difficulty-free years in the lower socioeconomic group. Absolute quantities are truncated remaining years from age 60 through age 100.",
  "",
  "## Figure 3. Post-onset contributions, decomposition, and functional-definition sensitivity",
  "",
  "Panel A shows primary recovery-plus-relapse contributions for all eight education and wealth contrasts. Panel B shows the focused secondary six-component decomposition within low education. Panel C shows the focused secondary recovery-plus-relapse contribution by outcome definition. Points and intervals are household-cluster bootstrap estimates and 95% confidence intervals. Negative contributions widen the lower-group FDFLE deficit; Shapley components are descriptive accounting quantities.",
  "",
  "## Additional file 2. Continuous-time process contrasts",
  "",
  "Focused secondary low-versus-high wealth process intensity ratios within low education from continuous-time panel Markov models are shown on a logarithmic scale. Intervals are model based."
)
writeLines(figure_legends, file.path(figure_dir, "MAIN_FIGURE_TITLES_AND_LEGENDS.md"), useBytes = TRUE)

manifest <- list(
  package = "Direction4 BMC Medicine first-submission vector figures",
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  inputs = as.list(required_files),
  assembled_svgs = basename(list.files(figure_dir, pattern = "^Figure[123].*[.]svg$", full.names = TRUE)),
  panel_svgs = basename(list.files(panel_dir, pattern = "[.]svg$", full.names = TRUE)),
  source_tables = basename(list.files(source_dir, pattern = "[.]csv$", full.names = TRUE)),
  svg_contract = "live Arial text; semantic panel/axis/data/labels/legend/annotations groups; vector primitives; no image elements",
  participant_level_data_opened = FALSE
)
write_json(manifest, file.path(figure_dir, "FIGURE_ASSET_MANIFEST_v1.1.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_56_v1_1_vector_figures.txt"), useBytes = TRUE)

cat("v1.1 SVG figures completed in ", figure_dir, "\n", sep = "")
