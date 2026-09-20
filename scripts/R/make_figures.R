# make_figures.R
# -----------------------------------------------------------------------
# Publication figures from the combined supplement CSVs (aggregate_supplement).
# Each figure is written as vector PDF and 300 dpi PNG, and the exact data
# plotted is written next to it as CSV (<out>/data/<figure>.csv).
# Legends and axes use full words, no abbreviations.
#
# Usage: Rscript make_figures.R <figure_data_dir> <output_dir>
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript make_figures.R <figure_data_dir> <output_dir>")
data_dir <- args[1]; out_dir <- args[2]
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)

rd <- function(name) {
  f <- file.path(data_dir, paste0(name, ".csv"))
  if (!file.exists(f)) stop("missing input: ", f)
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- labels: full words -----------------------------------------------
metric_label <- c(
  RAP_mean = "Mean right atrial pressure", PAP_mean = "Mean pulmonary arterial pressure",
  PAP_max = "Systolic pulmonary arterial pressure", PAP_min = "Diastolic pulmonary arterial pressure",
  SAP_mean = "Mean systemic arterial pressure", SAP_max = "Systolic systemic arterial pressure",
  SAP_min = "Diastolic systemic arterial pressure", QpQs = "Pulmonary-to-systemic flow ratio",
  CO_Lmin = "Systemic flow (cardiac output)", Q_shunt_Lmin = "Shunt flow",
  Qp_Lmin = "Pulmonary flow", SVR = "Systemic vascular resistance", PVR = "Pulmonary vascular resistance",
  LVEDV = "Left ventricular end-diastolic volume", LVESV = "Left ventricular end-systolic volume",
  RVEDV = "Right ventricular end-diastolic volume", RVESV = "Right ventricular end-systolic volume",
  LVEF = "Left ventricular ejection fraction", RVEF = "Right ventricular ejection fraction",
  LVP_max = "Left ventricular systolic pressure", RVP_max = "Right ventricular systolic pressure",
  LAP_mean = "Mean left atrial pressure", LV_RV_systolic_gradient = "Left-to-right ventricular systolic gradient")
param_label <- c(
  "R.SAR" = "Systemic arterial resistance", "R.SC" = "Systemic capillary resistance",
  "R.SVEN" = "Systemic venous resistance", "R.PAR" = "Pulmonary arterial resistance",
  "R.PCOX" = "Pulmonary capillary resistance", "R.PVEN" = "Pulmonary venous resistance",
  "C.SAR" = "Systemic arterial compliance", "C.SVEN" = "Systemic venous compliance",
  "C.PAR" = "Pulmonary arterial compliance", "C.PVEN" = "Pulmonary venous compliance",
  "E.LV.EA" = "Left ventricular active elastance", "E.LV.EB" = "Left ventricular passive elastance",
  "E.RV.EA" = "Right ventricular active elastance", "E.RV.EB" = "Right ventricular passive elastance",
  "E.LA.EA" = "Left atrial active elastance", "E.RA.EA" = "Right atrial active elastance",
  "V0.LV" = "Left ventricular unstressed volume", "V0.RV" = "Right ventricular unstressed volume",
  "V0.LA" = "Left atrial unstressed volume", "V0.RA" = "Right atrial unstressed volume",
  "vsd.Cd" = "Shunt discharge coefficient",
  "group.R_sys_scale" = "Systemic resistance scale", "group.R_pul_scale" = "Pulmonary resistance scale")
role_label <- c(fitted = "Calibration target", evaluation_only = "Compared, not fitted",
                predicted_evaluated = "Predicted, compared with measurement")
lab <- function(x, dict) ifelse(x %in% names(dict), dict[x], x)
# Protocol variable names spelled out (no abbreviations in figures).
protocol_label <- function(x) {
  r <- c("^PA " = "Pulmonary arterial ", "^RA mean pressure" = "Mean right atrial pressure",
         ", RFA\\b" = ", femoral artery catheter", ", DAo\\b" = ", descending aorta catheter",
         ", NIBP\\b" = ", non-invasive cuff", ", form$" = ", protocol form (cuff)",
         "\\(RFA, protocol formula\\)" = "(femoral artery catheter, protocol formula)",
         "^LVEDV" = "Left ventricular end-diastolic volume", "^LVESV" = "Left ventricular end-systolic volume",
         "^RVEDV" = "Right ventricular end-diastolic volume", "^RVESV" = "Right ventricular end-systolic volume",
         "^LVEF" = "Left ventricular ejection fraction", "^RVEF" = "Right ventricular ejection fraction",
         "^LV systolic" = "Left ventricular systolic", "^RV systolic" = "Right ventricular systolic",
         "^LV-RV systolic gradient" = "Left-to-right ventricular systolic gradient",
         "Pulmonary flow Qp" = "Pulmonary flow", "Systemic flow Qs" = "Systemic flow",
         "^Qp/Qs$" = "Pulmonary-to-systemic flow ratio", "^SVR" = "Systemic vascular resistance",
         "pre-form row" = "pre-closure form row")
  for (k in names(r)) x <- gsub(k, r[[k]], x, perl = TRUE)
  x
}
governed <- c("RAP_mean", "PAP_mean", "PAP_max", "PAP_min", "SAP_mean", "SAP_max", "SAP_min", "QpQs", "CO_Lmin")

seed_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
theme_paper <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey92"))

save_fig <- function(p, name, data, w = 18, h = 12) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "cm")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  write.csv(data, file.path(out_dir, "data", paste0(name, ".csv")), row.names = FALSE)
  cat("  written:", name, "\n")
}
run_fig <- function(name, fn) {
  tryCatch(fn(), error = function(e) cat("  FAILED:", name, "-", conditionMessage(e), "\n"))
}

# ---- F01/F02: tornado plots --------------------------------------------
tornado <- function(d, title, name, with_range) {
  d <- d %>% filter(Output %in% governed) %>%
    mutate(OutputLabel = lab(Output, metric_label), ParameterLabel = lab(Parameter, param_label))
  top <- d %>% group_by(Output) %>% slice_max(ST, n = 8, with_ties = FALSE) %>% ungroup() %>%
    mutate(ParameterLabel = tidytext_reorder(ParameterLabel, ST, OutputLabel))
  p <- ggplot(top, aes(x = ST, y = ParameterLabel)) +
    geom_col(aes(fill = "Total-order index"), width = 0.7) +
    geom_point(aes(x = S1, shape = "First-order index"), size = 1.2) +
    geom_vline(xintercept = 0.10, linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~OutputLabel, scales = "free_y", ncol = 3, labeller = label_wrap_gen(30)) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = "grey60", name = NULL) + scale_shape_manual(values = 16, name = NULL) +
    labs(x = "Sobol sensitivity index (dashed line: screening threshold 0.10)", y = NULL) +
    theme_paper
  if (with_range) p <- p + geom_errorbar(aes(xmin = ST_min, xmax = ST_max), width = 0.3, linewidth = 0.3, orientation = "y")
  save_fig(p, name, top, w = 19, h = 20)
}
tidytext_reorder <- function(x, by, within) {
  new <- paste(x, within, sep = "___")
  stats::reorder(new, by)
}
run_fig("F01_gsa_tornado_screening", function() {
  g <- rd("gsa_indices_initial_all"); g <- g[g$Seed == g$Seed[1], ]
  tornado(g, "Sensitivity before calibration (screening)", "F01_gsa_tornado_screening", FALSE)
})
run_fig("F02_gsa_tornado_calibrated", function() {
  g <- rd("gsa_indices_final_summary") %>%
    transmute(Output, Parameter, ST = ST_mean, S1 = S1_mean, ST_min, ST_max)
  tornado(g, "Sensitivity at the calibrated state (mean and range across runs)", "F02_gsa_tornado_calibrated", TRUE)
})

# ---- F03: calibration effect (before vs after) -------------------------
run_fig("F03_calibration_before_after", function() {
  d <- rd("protocol_comparison_pre_all") %>% filter(Role == "fitted", ModelMetric %in% governed) %>%
    distinct(Seed, ModelMetric, .keep_all = TRUE) %>%
    select(Seed, ModelMetric, Before = ErrorPctBaseline, After = ErrorPct) %>%
    pivot_longer(c(Before, After), names_to = "Stage", values_to = "ErrorPct") %>%
    mutate(AbsErrorPct = abs(ErrorPct), Label = lab(ModelMetric, metric_label),
           Stage = factor(Stage, c("Before", "After"), c("Before calibration (scaled)", "After calibration")))
  order <- d %>% filter(grepl("^After", Stage)) %>% group_by(Label) %>%
    summarise(o = mean(AbsErrorPct), .groups = "drop") %>% arrange(o) %>% pull(Label)
  d$Label <- factor(d$Label, order)
  s <- d %>% group_by(Label, Stage) %>% summarise(Mean = mean(AbsErrorPct), .groups = "drop")
  p <- ggplot(d, aes(x = AbsErrorPct, y = Label, colour = Stage)) +
    geom_col(data = s, aes(x = Mean, y = Label, fill = Stage), inherit.aes = FALSE,
             position = position_dodge(0.8), width = 0.7, alpha = 0.35) +
    geom_point(position = position_dodge(0.8), size = 1.3) +
    geom_vline(xintercept = 10, linetype = "dashed", linewidth = 0.3) +
    scale_colour_manual(values = c("grey40", "#0072B2"), name = NULL) +
    scale_fill_manual(values = c("grey40", "#0072B2"), name = NULL) +
    labs(x = "Absolute error relative to measurement (%) (dashed line: 10% acceptance)", y = NULL) + theme_paper
  save_fig(p, "F03_calibration_before_after", d)
})

# ---- F04/F05: model vs protocol, pre and post --------------------------
error_panel <- function(d, name, title) {
  d <- d %>% filter(!is.na(ErrorPct)) %>%
    mutate(Label = protocol_label(Variable),
           RoleLabel = lab(Role, role_label), Seed = factor(Seed))
  p <- ggplot(d, aes(x = ErrorPct, y = reorder(Label, abs(ErrorPct)), colour = Seed, shape = RoleLabel)) +
    annotate("rect", xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf, alpha = 0.12) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    # The runs can coincide exactly (the accepted candidate is deterministic),
    # so separate them vertically; otherwise one run hides the others.
    geom_point(size = 1.6, alpha = 0.9,
               position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = seed_palette, name = "Run (multi-start seed)") +
    scale_shape_manual(values = c(16, 1, 17), name = NULL) +
    labs(x = "Model minus measurement (% of measurement); shaded band: within 10%", y = NULL) +
    theme_paper + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    theme(legend.box = "vertical")
  save_fig(p, name, d, w = 19, h = 14)
}
run_fig("F04_pre_model_vs_protocol", function() {
  error_panel(rd("protocol_comparison_pre_all"), "F04_pre_model_vs_protocol",
              "Pre-closure: calibrated model against every comparable protocol measurement")
})
run_fig("F05_post_prediction_vs_protocol", function() {
  error_panel(rd("protocol_comparison_post_all"), "F05_post_prediction_vs_protocol",
              "Post-closure: prediction from pre-closure parameters against protocol measurements")
})

# ---- F06: pressure waveforms --------------------------------------------
run_fig("F06_pressure_waveforms", function() {
  ts <- bind_rows(
    rd("timeseries_pre_calibrated_all") %>% mutate(State = "Before closure (calibrated)"),
    rd("timeseries_post_prediction_all") %>% mutate(State = "After closure (predicted)"))
  d <- ts %>% select(Seed, State, phase, P_SAR_mmHg, P_PAR_mmHg, P_RA_mmHg, P_LV_mmHg, P_RV_mmHg) %>%
    pivot_longer(starts_with("P_"), names_to = "Compartment", values_to = "Pressure") %>%
    mutate(Compartment = recode(Compartment, P_SAR_mmHg = "Systemic artery", P_PAR_mmHg = "Pulmonary artery",
                                P_RA_mmHg = "Right atrium", P_LV_mmHg = "Left ventricle", P_RV_mmHg = "Right ventricle"),
           Seed = factor(Seed), State = factor(State, c("Before closure (calibrated)", "After closure (predicted)")))
  p <- ggplot(d, aes(phase, Pressure, colour = Seed, linetype = State)) +
    geom_line(linewidth = 0.5) + facet_wrap(~Compartment, scales = "free_y", ncol = 3, labeller = label_wrap_gen(30)) +
    scale_colour_manual(values = seed_palette, name = "Run (multi-start seed)") +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    labs(x = "Fraction of the cardiac cycle", y = "Pressure (mmHg)") + theme_paper +
    guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) + theme(legend.box = "vertical")
  save_fig(p, "F06_pressure_waveforms", d, w = 19, h = 12)
})

# ---- F07: pressure-volume loops -----------------------------------------
run_fig("F07_pressure_volume_loops", function() {
  ts <- bind_rows(
    rd("timeseries_pre_baseline_all") %>% mutate(State = "Before closure, before calibration"),
    rd("timeseries_pre_calibrated_all") %>% mutate(State = "Before closure (calibrated)"),
    rd("timeseries_post_prediction_all") %>% mutate(State = "After closure (predicted)"))
  d <- bind_rows(
    ts %>% transmute(Seed, State, Ventricle = "Left ventricle", Volume = state_V_LV, Pressure = P_LV_mmHg, phase),
    ts %>% transmute(Seed, State, Ventricle = "Right ventricle", Volume = state_V_RV, Pressure = P_RV_mmHg, phase)) %>%
    arrange(Seed, State, Ventricle, phase) %>%
    mutate(Seed = factor(Seed),
           State = factor(State, c("Before closure, before calibration", "Before closure (calibrated)", "After closure (predicted)")))
  measured <- rd("protocol_comparison_post_all") %>% filter(ModelMetric %in% c("LVEDV", "LVESV", "RVEDV", "RVESV")) %>%
    distinct(ModelMetric, ComparedValue) %>%
    mutate(Ventricle = ifelse(grepl("^LV", ModelMetric), "Left ventricle", "Right ventricle"),
           Kind = ifelse(grepl("EDV", ModelMetric), "Measured end-diastolic volume", "Measured end-systolic volume"))
  p <- ggplot(d, aes(Volume, Pressure, colour = State, group = interaction(Seed, State))) +
    geom_path(linewidth = 0.45) +
    geom_vline(data = measured, aes(xintercept = ComparedValue, linetype = Kind), colour = "black", linewidth = 0.3) +
    facet_wrap(~Ventricle, scales = "free") +
    scale_colour_manual(values = c("grey60", "#0072B2", "#D55E00"), name = NULL) +
    scale_linetype_manual(values = c("dashed", "dotted"), name = NULL) +
    labs(x = "Volume (mL)", y = "Pressure (mmHg)") + theme_paper +
    guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) + theme(legend.box = "vertical")
  save_fig(p, "F07_pressure_volume_loops", d, w = 18, h = 10)
})

# ---- F08: shunt flow --------------------------------------------------------
run_fig("F08_shunt_flow", function() {
  d <- rd("timeseries_pre_calibrated_all") %>% transmute(Seed = factor(Seed), phase, Flow = Q_VSD_mLs)
  p <- ggplot(d, aes(phase, Flow, colour = Seed)) + geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(linewidth = 0.5) + scale_colour_manual(values = seed_palette, name = "Run (multi-start seed)") +
    labs(x = "Fraction of the cardiac cycle", y = "Shunt flow (mL/s), positive from left to right ventricle") + theme_paper
  save_fig(p, "F08_shunt_flow", d, w = 12, h = 8)
})

# ---- F09: calibrated parameters across runs ---------------------------------
run_fig("F09_parameters_across_runs", function() {
  d <- rd("parameters_calibration_all") %>%
    mutate(Label = lab(Parameter, param_label), Seed = factor(Seed),
           Status = ifelse(Calibrated == 1, "Calibrated (selected by sensitivity analysis)", "Fixed"),
           LowerRatio = lb / baseline_scaled, UpperRatio = ub / baseline_scaled)
  bounds <- d %>% distinct(Label, LowerRatio, UpperRatio)
  p <- ggplot(d, aes(y = Label)) +
    geom_errorbar(data = bounds, aes(xmin = LowerRatio, xmax = UpperRatio), width = 0, colour = "grey75", linewidth = 2, orientation = "y") +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_point(aes(x = RatioToScaled, colour = Seed, shape = Status), size = 1.7) +
    scale_x_log10() + scale_colour_manual(values = seed_palette, name = "Run (multi-start seed)") +
    scale_shape_manual(values = c(16, 1), name = NULL) +
    labs(x = "Value relative to the size-scaled value (log scale); grey bar: calibration bounds", y = NULL) +
    theme_paper + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) + theme(legend.box = "vertical")
  save_fig(p, "F09_parameters_across_runs", d, w = 18, h = 11)
})

# ---- F10/F11: identifiability ---------------------------------------------
run_fig("F10_identifiability_condition_number", function() {
  d <- rd("identifiability_reduction_all") %>% mutate(Seed = factor(Seed))
  p <- ggplot(d, aes(SubsetSize, ConditionNumber, colour = Seed)) + geom_line() + geom_point(size = 1.4) +
    scale_y_log10() + scale_x_continuous(breaks = sort(unique(d$SubsetSize))) +
    scale_colour_manual(values = seed_palette, name = "Run (multi-start seed)") +
    labs(x = "Number of parameters in the best-conditioned subset", y = "Condition number of the sensitivity matrix (log scale)") +
    theme_paper
  save_fig(p, "F10_identifiability_condition_number", d, w = 12, h = 9)
})
run_fig("F11_identifiability_correlation", function() {
  pr <- rd("identifiability_pairs_all")
  d <- bind_rows(pr, pr %>% rename(ParameterA = ParameterB, ParameterB = ParameterA)) %>%
    mutate(A = lab(ParameterA, param_label), B = lab(ParameterB, param_label), Seed = paste("Run, seed", Seed))
  p <- ggplot(d, aes(A, B, fill = Correlation)) + geom_tile() +
    geom_tile(data = filter(d, abs(Correlation) >= 0.9), colour = "black", linewidth = 0.5, fill = NA) +
    facet_wrap(~Seed) + scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", limits = c(-1, 1),
                                            name = "Correlation between sensitivity columns") +
    labs(x = NULL, y = NULL) + theme_paper + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  save_fig(p, "F11_identifiability_correlation", d, w = 22, h = 12)
})

# ---- F12: multi-start ---------------------------------------------------------
run_fig("F12_multistart", function() {
  d <- rd("multistart_all") %>% mutate(Seed = factor(Seed), Chosen = ifelse(Selected == 1, "Selected start", "Other start"))
  p <- ggplot(d, aes(Seed, PrimaryRMSE, shape = Chosen)) + geom_point(size = 2, position = position_jitter(width = 0.08, seed = 1)) +
    scale_shape_manual(values = c(1, 16), name = NULL) +
    labs(x = "Run (multi-start seed)", y = "Relative root-mean-square error over governed targets") + theme_paper
  save_fig(p, "F12_multistart", d, w = 10, h = 9)
})

# ---- F13: surrogate quality -----------------------------------------------------
run_fig("F13_gsa_surrogate_quality", function() {
  d <- rd("gsa_surrogate_quality_all") %>% filter(Output %in% c(governed, "LVESV", "LVEF", "LVEDV", "RVEDV", "RVESV")) %>%
    group_by(Stage, Output) %>% summarise(LOO = mean(ModifiedLOO), LOOmin = min(ModifiedLOO), LOOmax = max(ModifiedLOO), .groups = "drop") %>%
    mutate(Label = lab(Output, metric_label),
           Stage = recode(Stage, before_calibration = "Before calibration", after_calibration = "After calibration"))
  p <- ggplot(d, aes(LOO, reorder(Label, LOO), colour = Stage)) +
    geom_pointrange(aes(xmin = LOOmin, xmax = LOOmax), size = 0.25, position = position_dodge(0.5)) +
    scale_colour_manual(values = c("grey40", "#0072B2"), name = NULL) +
    labs(x = "Leave-one-out error of the polynomial chaos surrogate", y = NULL) + theme_paper
  save_fig(p, "F13_gsa_surrogate_quality", d, w = 14, h = 10)
})

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
cat("Figures done:", out_dir, "\n")
