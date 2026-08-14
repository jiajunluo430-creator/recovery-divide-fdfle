#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240, warn = 1)
if (nzchar(.d4_rlib <- Sys.getenv("RECOVERY_DIVIDE_RLIB", unset = ""))) .libPaths(c(.d4_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("RECOVERY_DIVIDE_ROOT", unset = normalizePath(".", winslash = "/", mustWork = FALSE))
final_dir <- file.path(root, "03_outputs", "06_bootstrap_summary_final")
point_dir <- file.path(root, "03_outputs", "04_formal_point")
sensitivity_dir <- file.path(root, "03_outputs", "07_sensitivity_point")
threshold_root <- file.path(root, "03_outputs", "09_threshold_sensitivity")
report_dir <- file.path(root, "05_reports")
log_dir <- file.path(root, "06_logs")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(final_dir, "bootstrap_validity_qc.csv"),
  file.path(final_dir, "bootstrap_percentile_intervals.csv"),
  file.path(final_dir, "mechanism_promotion_by_cohort.csv"),
  file.path(final_dir, "pairwise_country_heterogeneity.csv"),
  file.path(final_dir, "bootstrap_promotion_decision.csv"),
  file.path(sensitivity_dir, "sensitivity_direction_and_magnitude.csv"),
  file.path(sensitivity_dir, "sensitivity_rate_support_qc.csv")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) stop("Missing required final inputs: ", paste(missing_required, collapse = "; "))

cohort_order <- c("CHARLS", "ELSA", "HRS", "MHAS")
exposure_order <- c("education", "wealth")
exposure_cn <- c(education = "教育", wealth = "财富")

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}
ci_text <- function(point, low, high, digits = 2) {
  paste0(fmt(point, digits), " [", fmt(low, digits), ", ", fmt(high, digits), "]")
}
md_table <- function(z) {
  z <- as.data.frame(z, stringsAsFactors = FALSE)
  z[] <- lapply(z, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "NA"
    gsub("\\|", "\\\\|", x)
  })
  header <- paste0("| ", paste(names(z), collapse = " | "), " |")
  separator <- paste0("|", paste(rep("---", ncol(z)), collapse = "|"), "|")
  rows <- apply(z, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, separator, rows)
}

validity <- fread(file.path(final_dir, "bootstrap_validity_qc.csv"))
ci <- fread(file.path(final_dir, "bootstrap_percentile_intervals.csv"))
promotion <- fread(file.path(final_dir, "mechanism_promotion_by_cohort.csv"))
heterogeneity <- fread(file.path(final_dir, "pairwise_country_heterogeneity.csv"))
decision <- fread(file.path(final_dir, "bootstrap_promotion_decision.csv"))

gap_ci <- ci[metric == "gap_dfle", .(
  cohort, exposure, gap_point = point_estimate, gap_ci_low = ci_low, gap_ci_high = ci_high,
  gap_valid = valid_replicates
)]
main <- merge(promotion, gap_ci, by = c("cohort", "exposure"), all.x = TRUE)
main[, `:=`(.cohort_order = match(cohort, cohort_order), .exposure_order = match(exposure, exposure_order))]
setorder(main, .cohort_order, .exposure_order)

validity_display <- validity[, .(
  队列 = cohort,
  SES = unname(exposure_cn[exposure]),
  有效次数 = valid_replicates,
  失败次数 = failed_replicates,
  门槛 = ifelse(validity_gate_pass, "PASS", "FAIL"),
  失败原因 = ifelse(is.na(failure_messages) | !nzchar(failure_messages), "—", failure_messages)
)]
main_display <- main[, .(
  队列 = cohort,
  SES = unname(exposure_cn[exposure]),
  `DFLE差值 low-high, y [95% CI]` = ci_text(gap_point, gap_ci_low, gap_ci_high),
  `恢复+复发贡献, y [95% CI]` = ci_text(rr_point, rr_ci_low, rr_ci_high),
  `贡献比例, % [95% CI]` = ci_text(rr_percent_point, rr_percent_ci_low, rr_percent_ci_high, 1),
  `预冻结机制触发` = ifelse(mechanism_trigger_cohort, "是", "否")
)]
hetero_trigger <- heterogeneity[heterogeneity_trigger_pair == TRUE]
hetero_display <- hetero_trigger[, .(
  SES = unname(exposure_cn[exposure]),
  过程 = block,
  对比 = paste0(cohort_a, "−", cohort_b),
  `贡献差, y [95% CI]` = ci_text(point_difference_a_minus_b, ci_low, ci_high)
)]

report5 <- c(
  "# 正式 person-bootstrap 结果与绑定决策",
  "",
  "状态：500 次/队列的最终 person-cluster bootstrap；财富切点、五个过程模型、转移矩阵、生命表和 Shapley 分解均在每次重抽样中重新估计。",
  "",
  "## 绑定决策",
  "",
  if (isTRUE(decision$top_journal_promotion_gate[[1]]))
    "**GO：预冻结的顶刊升级门槛通过。** 该结论不是由单个 SES 系数显著性触发，而是由 recovery+relapse 的绝对年数、贡献比例、同向重复和/或跨国机制差异共同触发。"
  else
    "**PIVOT：预冻结的顶刊升级门槛未通过。** 保留全部真实结果，但不得把常规 onset gradient 包装成 recovery-divide 顶刊主线。",
  "",
  paste0("- 全部 cohort×SES 有效次数门槛：", ifelse(decision$all_cohort_exposure_validity_gates_pass[[1]], "通过", "未通过"), "。"),
  paste0("- 至少两个队列同向 recovery+relapse 机制门槛：", ifelse(decision$mechanism_same_direction_two_cohorts_met[[1]], "通过", "未通过"), "。"),
  paste0("- ≥0.75 年且区间排除 0 的跨国过程差异门槛：", ifelse(decision$country_heterogeneity_trigger_met[[1]], "通过", "未通过"), "。"),
  "",
  "## Bootstrap 有效性",
  "",
  md_table(validity_display),
  "",
  "`invalid standardized hazard` 表示该重抽样在固定模型下产生非有限/无效标准化 hazard；按合同计为失败，不做模型修补。最终门槛要求每个 cohort×SES 至少 450 次有效。",
  "",
  "## 主估计",
  "",
  md_table(main_display),
  "",
  "负值表示低 SES 相比高 SES 少有 disability-free years；负的过程贡献表示该过程扩大低 SES 的 DFLE 缺口。百分比区间在总缺口接近 0 的重抽样中可变得很宽，因此顶刊判断同时要求绝对年数。",
  "",
  "## 通过门槛的跨国机制差异",
  "",
  if (nrow(hetero_display)) md_table(hetero_display) else "无跨国过程差异满足预冻结阈值。",
  "",
  "## 解释边界",
  "",
  "- 这是观察性、描述性 transition-contribution decomposition，不是把教育或财富设为可干预处理的因果分解。",
  "- onset 仍是多数缺口的最大单一来源；主线是 recovery/relapse 提供了不能被 onset-only 故事解释的额外绝对年数。",
  "- MHAS 财富是预先保留的透明阴性/异质性结果，不得通过重分组或改阈值消除。",
  "- 绑定区间使用 percentile bootstrap；模型自身的 Wald 区间仅作诊断。"
)
writeLines(report5, file.path(report_dir, "05_bootstrap_final_results.md"), useBytes = TRUE)

sensitivity <- fread(file.path(sensitivity_dir, "sensitivity_direction_and_magnitude.csv"))
rate_qc <- fread(file.path(sensitivity_dir, "sensitivity_rate_support_qc.csv"))
sensitivity_order <- unique(sensitivity$sensitivity)
sensitivity[, `:=`(.sensitivity_order = match(sensitivity, sensitivity_order),
  .cohort_order = match(cohort, cohort_order), .exposure_order = match(exposure, exposure_order))]
setorder(sensitivity, .sensitivity_order, .cohort_order, .exposure_order)
sensitivity_display <- sensitivity[sensitivity != "primary_reproduction", .(
  敏感性 = sensitivity,
  队列 = cohort,
  SES = unname(exposure_cn[exposure]),
  `DFLE差值, y` = fmt(dfle_gap_low_minus_high),
  `恢复+复发, y` = fmt(recovery_relapse_contribution),
  `贡献比例, %` = fmt(recovery_relapse_percent, 1),
  `DFLE方向一致` = ifelse(dfle_gap_direction_concordant, "是", "否"),
  `机制方向一致` = ifelse(rr_direction_concordant, "是", "否")
)]
rate_display <- rate_qc[sensitivity %in% c("exclude_pandemic_crossing", "pre_pandemic_end_before_2020") &
  (hazards_gt_10 > 0 | hazards_gt_100 > 0), .(
    敏感性 = sensitivity, 队列 = cohort, SES = unname(exposure_cn[exposure]), 过程 = process,
    `最大annual hazard` = format(max_annual_hazard, scientific = TRUE, digits = 3),
    `hazard>10个数` = hazards_gt_10, `hazard>100个数` = hazards_gt_100
  )]

threshold_names <- c("adl_only", "at_least_two", "permissive_partial")
threshold_parts <- list()
for (threshold_name in threshold_names) {
  path <- file.path(threshold_root, threshold_name, "point_promotion_screen.csv")
  if (file.exists(path)) {
    z <- fread(path)
    z[, `:=`(threshold = threshold_name, model_status = "COMPLETED")]
  } else {
    z <- CJ(cohort = cohort_order, exposure = exposure_order)
    z[, `:=`(
      recovery_relapse_contribution_years = NA_real_, dfle_gap = NA_real_,
      recovery_relapse_percent = NA_real_, point_threshold_met = FALSE,
      threshold = threshold_name,
      model_status = "STOP: reproducible full-pipeline nonconvergence"
    )]
  }
  threshold_parts[[length(threshold_parts) + 1L]] <- z
}
threshold <- rbindlist(threshold_parts)
threshold[, `:=`(.threshold_order = match(threshold, threshold_names),
  .cohort_order = match(cohort, cohort_order), .exposure_order = match(exposure, exposure_order))]
setorder(threshold, .threshold_order, .cohort_order, .exposure_order)
threshold_display <- threshold[, .(
  阈值 = threshold, 队列 = cohort, SES = unname(exposure_cn[exposure]),
  `DFLE差值, y` = fmt(dfle_gap),
  `恢复+复发, y` = fmt(recovery_relapse_contribution_years),
  `贡献比例, %` = fmt(recovery_relapse_percent, 1),
  `点估计门槛` = ifelse(model_status == "COMPLETED", ifelse(point_threshold_met, "通过", "未通过"), "不可估"),
  `模型状态` = model_status
)]

report6 <- c(
  "# 预设敏感性分析",
  "",
  "## 一次只改变一个实现条件",
  "",
  md_table(sensitivity_display),
  "",
  "解释：`unweighted`、`interval_1_to_3_years` 与 `exclude_explicit_proxy` 保持主状态定义不变。CHARLS 早期波次 proxy 为结构性不可得，因此 `exclude_explicit_proxy` 只排除明确标记的 proxy，不能等同于完整 self-only。",
  "",
  "## 疫情相关敏感性的数值支持",
  "",
  if (nrow(rate_display)) md_table(rate_display) else "两项疫情敏感性未发现 annual hazard >10。",
  "",
  "删除关键跨期区间后，部分年龄×波次×状态组合失去支持；若出现极端标准化 hazard，则相应生命表结果标记为 numerically non-informative，保留但不作为反驳或支持主机制的证据。模型未被重新调参。",
  "",
  "## 三个预冻结功能阈值",
  "",
  md_table(threshold_display),
  "",
  "三套阈值均从 individual-wave item data 重新构建 I0/D1/R1/D2 历史和死亡风险集。ADL-only 检验更严重的基本活动受限，≥2 difficulties 检验更高负担，permissive partial-item 允许任一已观察困难判为 difficulty。只报告完成全部 40 个模型与生命表的阈值；permissive partial-item 在 MHAS education death_post 的完整流程中重复不收敛，按合同 STOP，未以孤立模型或类别合并救援。主论文的不确定性仍以冻结 primary state 的 500 次 bootstrap 为绑定依据。",
  "",
  "## 稳健性结论",
  "",
  "- 不加权、1–3 年间隔和排除明确 proxy 是主结论的核心实现敏感性；方向性失败必须逐项保留。",
  "- ADL-only 若明显减小 recovery/relapse 比例，说明 recovery divide 更适用于 broad functional difficulty，而不能外推为 severe ADL disability 的同等机制。",
  "- MHAS 财富持续作为异质/阴性结果，不以任何阈值选择来救。",
  "- 疫情限制若因支持断裂产生极端 hazard，则结论是该敏感性不可识别，而不是 recovery/relapse 为 0。"
)
writeLines(report6, file.path(report_dir, "06_sensitivity_results.md"), useBytes = TRUE)

report7 <- c(
  "# 投稿故事与经核验期刊 shortlist",
  "",
  "核验日期：2026-08-10。JIF 取期刊/出版方当前页面；SCIE 状态按 Clarivate Master Journal List 精确 ISSN 检索。",
  "",
  "## 推荐投稿顺序",
  "",
  "| 顺序 | 期刊 | 当前 JIF | Clarivate | 决策 |",
  "|---:|---|---:|---|---|",
  "| 1 | The Lancet Healthy Longevity | 14.6 | SCIE | scope-first 主投：最贴合 healthy longevity 与 DFLE transition mechanism |",
  "| 2 | Nature Aging | 25.0 (2025) | SCIE | conceptual stretch：需强调 recovery/relapse 作为 ageing inequality architecture，而非数据库比较 |",
  "| 3 | The Lancet Public Health | 25.2 | SCIE + SSCI | policy/inequality highest stretch：需把跨国制度差异与 population relevance 做强；描述性边界降低其当前适配度 |",
  "| 4 | Nature Communications | 18.1 (2025) | SCIE | broad backup：方法闭合和四国异质性充分时可执行 |",
  "",
  "`npj Aging` 当前 JIF 13.0 但仅 ESCI，按硬合同排除；`JAMA Network Open` 当前 JIF 9.7，亦排除。",
  "",
  "## 主投稿命题",
  "",
  "**The recovery divide: socioeconomic inequalities in disability-free life expectancy arise not only from onset, but from unequal return to independence and relapse after recovery.**",
  "",
  "建议正文结构：",
  "",
  "1. 四国 common-contract Gate-0 与 five-state framework；",
  "2. 教育/财富 low−high DFLE 绝对年数差；",
  "3. onset、recovery、relapse、post-disability death、pre-disability death 的闭合 Shapley 分解；",
  "4. 哪些国家/SES 的 recovery+relapse 贡献经 bootstrap 支持；",
  "5. MHAS 财富的弱/反向机制作为制度异质性而非异常值；",
  "6. 阈值、权重、proxy、间隔与疫情支持敏感性。",
  "",
  "## 与既有研究的真实增量",
  "",
  "既有工作已经覆盖 HRS–ELSA 的 SES–DFLE 总差、HRS onset/recovery/mortality 趋势、England 教育差异，以及 MHAS/US–Mexico multistate transitions。不能宣称首次发现 SES 与 disability 或首次使用 multistate life table。可防守增量是：在同一九项功能合同下比较中国、英格兰、美国、墨西哥，并把 low−high DFLE gap 以绝对年数闭合到 onset、recovery、relapse、失能后死亡和失能前死亡。",
  "",
  "## 声明边界",
  "",
  "避免使用 causal effect、mediation 或 policy effect。可使用 contribution、accounted for、consistent with 和 cross-country mechanism heterogeneity。若政策 counterfactual 模块后续升级，必须作为描述性 transition equalisation scenario，而非教育/财富干预效应。",
  "",
  "官方核验链接见 `05_reports/02_live_collision_and_journal_audit.md`。"
)
writeLines(report7, file.path(report_dir, "07_submission_story_and_journal_shortlist.md"), useBytes = TRUE)

report8 <- c(
  "# 可复现执行与投稿候选包交接",
  "",
  "## 运行顺序",
  "",
  "1. `00_audit_source_schema.R`",
  "2. `02_audit_elsa_wave11.R`、`02c_audit_charls2020_schema.R`、`02d_audit_charls2020_adapter.R`",
  "3. `01_build_gate0_panel.R`",
  "4. `03_build_formal_risksets.R`",
  "5. `04_fit_formal_models_and_decompose.R`",
  "6. `05_life_tables_shapley_from_checkpoint.R`",
  "7. `06_bootstrap_cohort.R COHORT 500 final`（四队列可并行）",
  "8. `07_summarize_bootstrap.R final`",
  "9. `08_sensitivity_point_models.R`",
  "10. `08b_build_threshold_sensitivity_risksets.R` 与 `08c_run_threshold_sensitivity_models.ps1`",
  "11. `09_build_submission_tables_and_figure.R final`",
  "12. `10_export_figure_with_illustrator.jsx`、`10b_run_figure_qc.ps1` 与 Illustrator open QC，完成 Figure 1。",
  "13. `11_write_final_reports.R`，锁定 v1 主分析报告。",
  "14. 运行 v2 支持与点估计模块：`13_exploratory_upgrade_support.R`、`14_exploratory_upgrade_point_models.R`、`14a_recheck_age_interaction_stop.R`、`18_audit_recovery_exhaustion_support.R`。",
  "15. `14b_launch_exploratory_sensitivities.ps1`，再由 `19_summarize_exploratory_sensitivities.R` 汇总。",
  "16. `15b_launch_exploratory_preview.ps1` 通过预览门槛后，运行 `15c_launch_exploratory_final.ps1`；`15d`/`15e` 仅用于缺失 shard 的恢复，不得重复已完成工作。",
  "17. `16_summarize_exploratory_bootstrap.R final` 只汇总一次；随后运行 `22_write_exploratory_final_report_and_tables.R`。",
  "18. `21_build_exploratory_figure.R`、`21b_export_exploratory_figure_with_illustrator.jsx`、`21c_run_exploratory_figure_qc.ps1` 与 Illustrator open QC，完成 Figure 2。",
  "19. `26_write_binding_manuscript.py`、`23_build_manuscript_table1.R`、`24_build_submission_documents.py` 生成终稿源文件、Table 1、正文和投稿信。",
  "20. `24b_run_submission_document_qa.ps1` 后逐页视觉复核，并运行 `27_finalize_strobe_from_render.py`。",
  "21. 所有 final gate 通过后，运行 `25_build_v2_submission_candidate.ps1`。",
  "",
  "非交互式 R 固定为 `Rscript`。每个主步骤均在 `06_logs` 保存 timestamp log 或 session information。",
  "",
  "## 发布边界",
  "",
  "- `02_derived/*.rds` 和阈值内部风险集含个体级数据，不得进入投稿/公开包。",
  "- 原始 `.dta/.sav/.sas7bdat` 不复制到候选包。",
  "- 可发布内容限于冻结合同、变量/状态字典、聚合 Gate 表、模型/生命表聚合结果、代码、报告和 vector figures。",
  "- `90_invalidated` 保留纠错前内部输出供审计，但不进入投稿候选包。",
  "- v2 候选包还排除 preview、test、原始 bootstrap shard、模型 checkpoint 与 participant-level risk set；只纳入通过 QC 的 final 汇总。",
  "",
  "## 最终绑定产物",
  "",
  "- v2 结果报告：`05_reports/15_exploratory_bootstrap_final.md`。",
  "- v2 主表：`03_outputs/17_exploratory_submission_tables_final`。",
  "- Figure 1/2：`04_figures`、`04_figures_v2`，均有 raster-free SVG/PDF 和 Illustrator 结构/open QC。",
  "- 正文与投稿信：`08_manuscript_tlhl_v2/submission_documents`。",
  "- Appendix Methods 与 STROBE：`08_manuscript_tlhl_v2/supplement`。",
  "- 非披露发布候选：`07_submission_candidate_v2_20260810` 及同名 ZIP；`PACKAGE_QA.txt` 与 `SHA256_MANIFEST.csv` 为最终封装证据。",
  "",
  "## 外部文件需求",
  "",
  "当前四队列 primary analysis 不需要用户追加下载。ELSA 后期死亡完整性不足已通过冻结的 death-supported origin-wave 范围处理，而不是把 EOL 子样本当作完整死亡登记。若未来扩展 ELSA late-wave mortality，需获取覆盖 wave 7–11 全体样本、可区分死亡与非死亡失访的官方 mortality status/date release；在此之前不得扩展死亡风险窗。"
)
writeLines(report8, file.path(report_dir, "08_reproducibility_handoff.md"), useBytes = TRUE)

completed_threshold <- threshold[model_status == "COMPLETED"]
adl_trigger_n <- completed_threshold[threshold == "adl_only" & point_threshold_met == TRUE, .N]
two_trigger_n <- completed_threshold[threshold == "at_least_two" & point_threshold_met == TRUE, .N]
qualifying_primary_n <- main[mechanism_trigger_cohort == TRUE, .N]
overall_label <- if (isTRUE(decision$top_journal_promotion_gate[[1]])) "GO" else "PIVOT"
executive <- c(
  "# Direction 4 最终执行决策",
  "",
  paste0("## **", overall_label, "**"),
  "",
  if (overall_label == "GO")
    paste0("Primary any-difficulty state 的绑定 bootstrap 升级门槛通过，共 ", qualifying_primary_n,
      " 个 cohort×SES 组合达到预冻结 recovery+relapse 年数、比例和区间标准。项目可进入投稿候选阶段。")
  else
    "Primary bootstrap 升级门槛未通过；保留 onset/DFLE 结果，但 recovery-divide 顶刊主线需降级。",
  "",
  "## 必须同时保留的 claim boundary",
  "",
  paste0("- ADL-only 点估计门槛通过 ", adl_trigger_n, " 个组合；≥2 difficulties 通过 ", two_trigger_n, " 个组合。"),
  "- 因而主张应明确限定为 common 5 ADL + 4 IADL 中的 **any reported difficulty**；不得声称相同机制同等适用于 severe/multiple disability。",
  "- permissive partial-item threshold 在完整正式流程中重复不收敛，按 STOP 保留，不能用孤立拟合替代。",
  "- 疫情区间限制产生极端 age-period mortality hazards，标记为 numerically non-informative，既不支持也不反驳主线。",
  "",
  "## 投稿动作",
  "",
  if (overall_label == "GO")
    "首选 The Lancet Healthy Longevity；The Lancet Public Health 与 Nature Aging 保留为 stretch，前提是摘要和讨论把严重度阈值特异性当作真实结果，而不是隐藏。"
  else
    "停止 stretch framing，转向仍满足 JIF>10/SCIE 且接受描述性 multistate inequality 的期刊；不得用新阈值或新模型救结果。"
)
writeLines(executive, file.path(report_dir, "09_executive_decision.md"), useBytes = TRUE)

edu <- main[exposure == "education"]
wealth <- main[exposure == "wealth"]
qualifying_labels <- main[mechanism_trigger_cohort == TRUE, paste0(cohort, " ", exposure)]
abstract <- c(
  "# Structured abstract candidate",
  "",
  "## Background",
  "",
  "Socioeconomic inequalities in disability-free life expectancy are established, but it is unclear whether these inequalities arise mainly from disability onset or also from unequal recovery, relapse after recovery, and mortality after disability. We compared these processes across China, England, the United States, and Mexico.",
  "",
  "## Methods",
  "",
  "We harmonised longitudinal data from CHARLS, ELSA, HRS, and MHAS for adults aged 60 years or older. The primary functional state was any difficulty in a common set of five activities of daily living and four instrumental activities. We modelled interval-observed transitions among independence, first difficulty, recovery, relapse, and verified death using survey-weighted complementary log-log models. Age-specific transition matrices yielded remaining disability-free life expectancy at age 60. We decomposed low-versus-high education and wealth gaps into onset, recovery, relapse, post-disability mortality, and pre-disability mortality using Shapley values. Uncertainty was estimated with 500 person-cluster bootstrap replicates that re-estimated wealth cut points and the full analysis.",
  "",
  "## Findings",
  "",
  paste0(
    "Low education was associated with ", fmt(abs(max(edu$gap_point)), 2), "–", fmt(abs(min(edu$gap_point)), 2),
    " fewer disability-free years across cohorts. Recovery and relapse jointly contributed ",
    fmt(min(abs(edu$rr_point)), 2), "–", fmt(max(abs(edu$rr_point)), 2),
    " years (", fmt(min(edu$rr_percent_point), 1), "%–", fmt(max(edu$rr_percent_point), 1), "%) of the education gap. ",
    "The corresponding wealth contribution was adverse in CHARLS, ELSA, and HRS but not MHAS. ",
    "The frozen recovery-relapse uncertainty criterion was met for ", paste(qualifying_labels, collapse = ", "), "."
  ),
  "",
  "## Interpretation",
  "",
  "Socioeconomic differences in later-life functional health were not explained by disability onset alone. Unequal recovery and relapse accounted for a material share of broad difficulty-free life-expectancy gaps in several settings, while Mexico showed a distinct wealth pattern. The mechanism attenuated under stricter disability thresholds, indicating that the recovery divide should be interpreted as a feature of broad functional difficulty rather than assumed to apply equally to multiple or severe disability. These estimates are descriptive transition contributions and not causal effects of education or wealth."
)
writeLines(abstract, file.path(report_dir, "10_structured_abstract_candidate.md"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_final_report_writer.txt"))
cat("Final reports written to ", report_dir, "\n", sep = "")
