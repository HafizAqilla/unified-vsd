# make_manuscript_figures.R
# -----------------------------------------------------------------------
# Main-text figures for the P1 manuscript, laid out after Lechuga et al.
# (2026) and Owashi et al. (2020): sensitivity bar chart + heat map,
# simulated waveforms against the measured values, parity plots, and the
# pressure-volume loops with a loop-metrics table. The three multi-start
# seeds give identical results, so seed 101 is plotted.
#
# Each figure is written as vector PDF and 300 dpi PNG, with the plotted
# data next to it (<out>/data/<figure>.csv).
#
# Usage: Rscript make_manuscript_figures.R <figure_data_dir> <output_dir> [noise_summary_dir]
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript make_manuscript_figures.R <figure_data_dir> <output_dir> [noise_summary_dir]")
data_dir <- args[1]; out_dir <- args[2]
noise_dir <- if (length(args) >= 3) args[3] else ""
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
SEED <- 101

rd <- function(name) {
  f <- file.path(data_dir, paste0(name, ".csv"))
  if (!file.exists(f)) stop("missing input: ", f)
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Seed" %in% names(d)) d <- d[d$Seed == SEED, ]
  d
}

col_pre  <- "#0072B2"   # calibrated pre-closure
col_post <- "#D55E00"   # predicted post-closure
col_base <- "grey55"    # scaled baseline
col_meas <- "black"

theme_paper <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        legend.position = "bottom", plot.tag = element_text(face = "bold"))

save_fig <- function(p, name, data, w = 18, h = 12) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "cm")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  write.csv(data, file.path(out_dir, "data", paste0(name, ".csv")), row.names = FALSE)
  cat("  written:", name, "\n")
}
run_fig <- function(name, fn) tryCatch(fn(), error = function(e) cat("  FAILED:", name, "-", conditionMessage(e), "\n"))

governed <- c("RAP_mean", "PAP_mean", "SAP_mean", "CO_Lmin", "QpQs", "PAP_min", "PAP_max", "SAP_min", "SAP_max")
target_short <- c(RAP_mean = "Mean RA
pressure", PAP_mean = "Mean PA
pressure", SAP_mean = "Mean arterial
pressure",
                  CO_Lmin = "Qs", QpQs = "Qp/Qs", PAP_min = "Diastolic PA
pressure",
                  PAP_max = "Systolic PA
pressure", SAP_min = "Diastolic RFA
pressure",
                  SAP_max = "Systolic RFA
pressure")
cand_label <- c("group.R_sys_scale" = "Systemic resistance scale", "group.R_pul_scale" = "Pulmonary resistance scale",
                "R.SVEN" = "Systemic venous resistance", "C.SAR" = "Systemic arterial compliance",
                "C.PAR" = "Pulmonary arterial compliance", "E.RV.EB" = "RV passive elastance",
                "E.RV.EA" = "RV active elastance", "E.LV.EB" = "LV passive elastance",
                "E.LV.EA" = "LV active elastance", "E.RA.EA" = "RA active elastance",
                "E.LA.EA" = "LA active elastance", "V0.RV" = "RV unstressed volume",
                "V0.LV" = "LV unstressed volume", "vsd.Cd" = "Shunt discharge coefficient")
group_of <- c(R.SAR = "group.R_sys_scale", R.SC = "group.R_sys_scale",
              R.PAR = "group.R_pul_scale", R.PCOX = "group.R_pul_scale", R.PVEN = "group.R_pul_scale")

# ---- M1: sensitivity screen ----------------------------------------------
run_fig("M1_gsa_screening", function() {
  g <- rd("gsa_indices_initial_all") %>% filter(Output %in% governed) %>%
    mutate(Candidate = ifelse(Parameter %in% names(group_of), group_of[Parameter], Parameter)) %>%
    filter(Candidate %in% names(cand_label)) %>%
    group_by(Candidate, Output) %>% summarise(ST = max(ST), .groups = "drop")
  mx <- g %>% group_by(Candidate) %>% summarise(MaxST = max(ST), .groups = "drop") %>%
    mutate(Retained = MaxST >= 0.10, Label = cand_label[Candidate]) %>% arrange(desc(MaxST))
  lvl <- rev(mx$Label)                       # largest at the top
  # (A) tornado-style: sorted horizontal bars, threshold as a vertical line
  pa <- ggplot(mx, aes(y = factor(Label, levels = lvl), x = MaxST, fill = Retained)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0.10, linetype = "dashed", linewidth = 0.35) +
    geom_text(aes(x = pmax(MaxST, 0.10), label = sprintf("%.2f", MaxST)), hjust = -0.2, size = 2.3) +
    scale_fill_manual(values = c(`TRUE` = col_pre, `FALSE` = "grey70"),
                      labels = c(`TRUE` = "Retained and calibrated", `FALSE` = "Held at scaled value"),
                      breaks = c("TRUE", "FALSE"), name = NULL) +
    scale_x_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2), expand = c(0, 0)) +
    labs(y = NULL, x = expression("Largest total-order index, " * S[T])) +
    theme_paper +
    theme(legend.position = c(0.62, 0.12), legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 7), legend.background = element_blank())
  # (B) heat map, same row order as (A)
  hm <- g %>% mutate(Label = factor(cand_label[Candidate], levels = lvl),
                     Target = factor(gsub("
", " ", target_short[Output]), levels = gsub("
", " ", target_short[governed])))
  pb <- ggplot(hm, aes(x = Target, y = Label, fill = ST)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = filter(hm, ST >= 0.10), colour = "black", fill = NA, linewidth = 0.45) +
    geom_text(aes(label = ifelse(ST >= 0.05, sprintf("%.2f", ST), "")), size = 2.2,
              colour = ifelse(hm$ST > 0.45, "white", "black")) +
    scale_fill_gradient(low = "white", high = col_pre, limits = c(0, 0.75), name = expression(S[T])) +
    labs(x = NULL, y = NULL) + theme_paper +
    theme(legend.position = "right", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), axis.text.x = element_text(size = 6.8, angle = 45, hjust = 1))
  p <- pa + pb + plot_layout(widths = c(1, 1.45)) + plot_annotation(tag_levels = "A")
  save_fig(p, "M1_gsa_screening", hm, w = 18, h = 9.5)
})

# ---- measured values for the waveform figures -----------------------------
pre_cmp <- rd("protocol_comparison_pre_all")
post_cmp <- rd("protocol_comparison_post_all")
if ("Dataset" %in% names(post_cmp)) post_cmp <- post_cmp[post_cmp$Dataset == "medical_record", ]
val <- function(d, pattern) {
  v <- suppressWarnings(as.numeric(d$ComparedValue[grepl(pattern, d$Variable)]))
  v <- v[is.finite(v)]
  if (length(v) == 0) stop("no measured value for ", pattern) else v[1]
}

last_cycle <- function(name) {
  d <- rd(name)
  d$t <- d$time_s - min(d$time_s)
  d
}

tmean <- function(t, y) {                       # time-weighted cycle mean (ode15s steps are uneven)
  sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2) / (max(t) - min(t))
}
wave_panel <- function(d, ycols, labels, colours, meas, title, ylab, base = NULL) {
  long <- d %>% select(t, all_of(ycols)) %>% pivot_longer(-t, names_to = "Signal", values_to = "P") %>%
    mutate(Signal = labels[Signal])
  cols <- setNames(colours, labels)
  if (!is.null(base)) {
    bl <- base %>% transmute(t, Signal = "Scaled baseline (before calibration)", P = .data[[ycols[1]]])
    long <- bind_rows(long, bl)
    cols <- c(cols, "Scaled baseline (before calibration)" = col_base)
  }
  long$Signal <- factor(long$Signal, levels = names(cols))
  p <- ggplot(long, aes(t, P, colour = Signal)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(title = title, x = "Time (s)", y = ylab) + theme_paper +
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = 7), plot.title = element_text(size = 8.5, face = "bold")) +
    guides(colour = guide_legend(nrow = 1, order = 1))
  if (nrow(meas) > 0) {
    p <- p + geom_hline(data = meas, aes(yintercept = y, linetype = What), colour = col_meas, linewidth = 0.35) +
      scale_linetype_manual(values = c("Measured systolic" = "dashed", "Measured diastolic" = "dotdash",
                                       "Measured mean" = "solid", "Measured peak" = "dashed"), name = NULL) +
      guides(linetype = guide_legend(nrow = 1, order = 2))
  }
  p
}
mk <- function(...) { x <- c(...); data.frame(What = names(x), y = unname(x)) }

# ---- M2: pre-closure waveforms against the catheter record ---------------
run_fig("M2_pre_waveforms", function() {
  d <- last_cycle("timeseries_pre_calibrated_all"); b <- last_cycle("timeseries_pre_baseline_all")
  pa <- wave_panel(d, "P_SAR_mmHg", c(P_SAR_mmHg = "Systemic arterial (model)"), col_pre,
                   mk("Measured systolic" = val(pre_cmp, "^Systolic pressure, RFA"),
                      "Measured diastolic" = val(pre_cmp, "^Diastolic pressure, RFA"),
                      "Measured mean" = val(pre_cmp, "Mean arterial pressure \\(RFA")),
                   "Systemic arterial pressure", "Pressure (mmHg)", base = b)
  pb <- wave_panel(d, c("P_LV_mmHg", "P_RV_mmHg"), c(P_LV_mmHg = "LV (model)", P_RV_mmHg = "RV (model)"),
                   c(col_pre, "#009E73"),
                   mk("Measured peak" = val(pre_cmp, "^LV systolic pressure")),
                   "Ventricular pressures", "Pressure (mmHg)")
  pb <- pb + geom_hline(yintercept = val(pre_cmp, "^RV systolic pressure"), linetype = "dashed",
                        colour = "#009E73", linewidth = 0.35) +
    labs(caption = "Dashed lines: measured LV (black) and RV (green) peak systolic pressure")
  pc <- wave_panel(d, "P_PAR_mmHg", c(P_PAR_mmHg = "Pulmonary arterial (model)"), col_pre,
                   mk("Measured systolic" = val(pre_cmp, "^PA systolic"),
                      "Measured diastolic" = val(pre_cmp, "^PA diastolic"),
                      "Measured mean" = val(pre_cmp, "^PA mean")),
                   "Pulmonary arterial pressure", "Pressure (mmHg)", base = b)
  pd <- wave_panel(d, "P_RA_mmHg", c(P_RA_mmHg = "Right atrial (model)"), col_pre,
                   mk("Measured mean" = val(pre_cmp, "^RA mean")),
                   "Right atrial pressure", "Pressure (mmHg)", base = b)
  sh <- d %>% mutate(Q = Q_VSD_mLs * 60 / 1000)
  shunt_meas <- val(pre_cmp, "^Shunt flow")
  pe <- ggplot(sh, aes(t, Q)) + geom_line(colour = col_pre, linewidth = 0.6) +
    geom_hline(yintercept = tmean(sh$t, sh$Q), colour = col_pre, linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = shunt_meas, colour = col_meas, linewidth = 0.35) +
    annotate("text", x = max(sh$t), y = shunt_meas, label = sprintf("Measured mean %.3f L/min", shunt_meas),
             hjust = 1, vjust = -0.4, size = 2.3) +
    annotate("text", x = max(sh$t), y = tmean(sh$t, sh$Q), label = sprintf("Model mean %.3f L/min", tmean(sh$t, sh$Q)),
             hjust = 1, vjust = 1.4, size = 2.3, colour = col_pre) +
    labs(title = "Shunt flow (positive = left-to-right)", x = "Time (s)", y = "Flow (L/min)") + theme_paper +
    theme(plot.title = element_text(size = 8.5, face = "bold"))
  p <- (pa + pb) / (pc + pd) / (pe + plot_spacer()) + plot_annotation(tag_levels = "A")
  save_fig(p, "M2_pre_waveforms", d, w = 18, h = 20)
})

# ---- M3: post-closure predicted waveforms --------------------------------
run_fig("M3_post_waveforms", function() {
  d <- last_cycle("timeseries_post_prediction_all")
  pa <- wave_panel(d, "P_SAR_mmHg", c(P_SAR_mmHg = "Systemic arterial (predicted)"), col_post,
                   mk("Measured systolic" = val(post_cmp, "^Systolic pressure, RFA"),
                      "Measured diastolic" = val(post_cmp, "^Diastolic pressure, RFA"),
                      "Measured mean" = val(post_cmp, "Mean arterial pressure \\(RFA")),
                   "Systemic arterial pressure", "Pressure (mmHg)")
  pb <- wave_panel(d, c("P_LV_mmHg", "P_RV_mmHg"), c(P_LV_mmHg = "LV (predicted)", P_RV_mmHg = "RV (predicted)"),
                   c(col_post, "#009E73"), data.frame(), "Ventricular pressures", "Pressure (mmHg)")
  pc <- wave_panel(d, "P_PAR_mmHg", c(P_PAR_mmHg = "Pulmonary arterial (predicted)"), col_post,
                   mk("Measured systolic" = val(post_cmp, "^PA systolic"),
                      "Measured diastolic" = val(post_cmp, "^PA diastolic"),
                      "Measured mean" = val(post_cmp, "^PA mean")),
                   "Pulmonary arterial pressure", "Pressure (mmHg)")
  pd <- wave_panel(d, "P_RA_mmHg", c(P_RA_mmHg = "Right atrial (predicted)"), col_post,
                   mk("Measured mean" = val(post_cmp, "^RA mean")),
                   "Right atrial pressure", "Pressure (mmHg)")
  p <- (pa + pb) / (pc + pd) + plot_annotation(tag_levels = "A")
  save_fig(p, "M3_post_waveforms", d, w = 18, h = 14)
})

# ---- M4: parity plots ------------------------------------------------------
run_fig("M4_parity", function() {
  keep_pre <- pre_cmp %>% filter(Role %in% c("fitted", "evaluation_only"), is.finite(ModelValue),
                                 is.finite(ComparedValue), !grepl("NIBP|form$", Variable),
                                 !grepl("^SVR|gradient", Variable)) %>%
    mutate(State = ifelse(Role == "fitted", "Pre-closure, fitted", "Pre-closure, not fitted"))
  keep_post <- post_cmp %>% filter(grepl("predicted", Role), is.finite(ModelValue), is.finite(ComparedValue),
                                   !grepl("^SVR", Variable)) %>%
    mutate(State = "Post-closure, predicted")
  cols <- c("State", "Variable", "Unit", "ComparedValue", "ModelValue", "ErrorPct")
  d <- bind_rows(keep_pre[, cols], keep_post[, cols]) %>%
    mutate(Panel = case_when(Unit == "mmHg" ~ "Pressures (mmHg)", Unit == "mL" ~ "Volumes (mL)",
                             Unit == "L/min" ~ "Flows (L/min)", TRUE ~ "Ratios (-)")) %>%
    select(State, Panel, Variable, Unit, Measured = ComparedValue, Model = ModelValue, ErrorPct)
  stats <- d %>% group_by(Panel, State) %>%
    summarise(n = n(), R2 = ifelse(n >= 5, cor(Measured, Model)^2, NA),
              RMSE = sqrt(mean((Model - Measured)^2)), .groups = "drop")
  panels <- c("Pressures (mmHg)", "Flows (L/min)", "Volumes (mL)", "Ratios (-)")
  state_levels <- c("Pre-closure, fitted", "Pre-closure, not fitted", "Post-closure, predicted")
  d <- d %>% mutate(Panel = factor(Panel, levels = panels), State = factor(State, levels = state_levels))
  lims <- d %>% group_by(Panel) %>%
    summarise(lo = min(c(Measured, Model)), hi = max(c(Measured, Model)), .groups = "drop") %>%
    mutate(pad = 0.08 * (hi - lo), lo = lo - pad, hi = hi + pad)
  corners <- bind_rows(lims %>% transmute(Panel, Measured = lo, Model = lo),
                       lims %>% transmute(Panel, Measured = hi, Model = hi))
  band <- lims %>% rowwise() %>%
    do(data.frame(Panel = .$Panel, x = seq(.$lo, .$hi, length.out = 50))) %>% ungroup() %>%
    mutate(lo = 0.9 * x, hi = 1.1 * x)
  lab_txt <- stats %>% mutate(Panel = factor(Panel, levels = panels),
                              State = factor(State, levels = state_levels)) %>% arrange(Panel, State) %>%
    group_by(Panel) %>%
    summarise(txt = paste(sprintf("%s: n = %d, RMSE = %.3g%s", State, n, RMSE,
                                  ifelse(is.na(R2), "", sprintf(", R² = %.2f", R2))), collapse = "
"),
              .groups = "drop") %>% left_join(lims, by = "Panel")
  p <- ggplot(d, aes(Measured, Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_blank(data = corners) +
    geom_point(aes(colour = State, shape = State), size = 2) +
    geom_text(data = lab_txt, aes(x = lo, y = hi, label = txt), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 2.1) +
    scale_colour_manual(values = c("Pre-closure, fitted" = col_pre, "Pre-closure, not fitted" = col_pre,
                                   "Post-closure, predicted" = col_post), name = NULL, drop = FALSE) +
    scale_shape_manual(values = c("Pre-closure, fitted" = 16, "Pre-closure, not fitted" = 1,
                                  "Post-closure, predicted" = 17), name = NULL, drop = FALSE) +
    facet_wrap(~Panel, scales = "free", ncol = 2) +
    labs(x = "Measured", y = "Model") + theme_paper + theme(aspect.ratio = 1)
  save_fig(p, "M4_parity", d, w = 18, h = 18)
  write.csv(stats, file.path(out_dir, "data", "M4_parity_stats.csv"), row.names = FALSE)
})

# ---- M5: pressure-volume loops + loop metrics ------------------------------
loop_metrics <- function(d, state) {
  one <- function(P, V, side) {
    sw <- abs(sum(P * c(diff(V), V[1] - V[length(V)])))       # mmHg*mL, closed loop
    # End-systole = maximum P/V (the upper-left corner of the loop). With an
    # open VSD the minimum volume comes later, during relaxation, because the
    # LV keeps emptying into the RV through the defect.
    edv <- max(V); esv <- min(V); ies <- which.max(P / V)
    data.frame(State = state, Ventricle = side, EDV_mL = edv, ESV_mL = esv, SV_mL = edv - esv,
               EF = (edv - esv) / edv, EDP_mmHg = P[which.max(V)], ESP_mmHg = P[ies],
               ESV_at_ESP_mL = V[ies],
               Pmax_mmHg = max(P), SW_mJ = sw * 133.322e-3)
  }
  rbind(one(d$P_LV_mmHg, d$state_V_LV, "LV"), one(d$P_RV_mmHg, d$state_V_RV, "RV"))
}
run_fig("M5_pv_loops", function() {
  b <- last_cycle("timeseries_pre_baseline_all"); c <- last_cycle("timeseries_pre_calibrated_all")
  q <- last_cycle("timeseries_post_prediction_all")
  st <- c(b = "Pre-closure, scaled baseline", c = "Pre-closure, calibrated", q = "Post-closure, predicted")
  lp <- bind_rows(
    lapply(list(b = b, c = c, q = q), function(x) x), .id = "k") %>%
    transmute(State = factor(st[k], levels = st), LV_V = state_V_LV, LV_P = P_LV_mmHg,
              RV_V = state_V_RV, RV_P = P_RV_mmHg)
  pal <- setNames(c(col_base, col_pre, col_post), st)
  meas <- data.frame(Ventricle = c("Left ventricle", "Left ventricle", "Right ventricle", "Right ventricle"),
                     V = c(val(post_cmp, "^LVEDV"), val(post_cmp, "^LVESV"),
                           val(post_cmp, "^RVEDV"), val(post_cmp, "^RVESV")),
                     What = rep(c("Measured post-closure end-diastolic volume",
                                  "Measured post-closure end-systolic volume"), 2))
  long <- bind_rows(lp %>% transmute(State, Ventricle = "Left ventricle", V = LV_V, P = LV_P),
                    lp %>% transmute(State, Ventricle = "Right ventricle", V = RV_V, P = RV_P))
  p <- ggplot(long, aes(V, P, colour = State)) + geom_path(linewidth = 0.6) +
    geom_vline(data = meas, aes(xintercept = V, linetype = What), linewidth = 0.35) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = c("dashed", "dotted"), name = NULL) +
    facet_wrap(~Ventricle, scales = "free") +
    labs(x = "Volume (mL)", y = "Pressure (mmHg)") + theme_paper +
    guides(colour = guide_legend(nrow = 1, order = 1), linetype = guide_legend(nrow = 1, order = 2)) +
    theme(legend.box = "vertical")
  m <- bind_rows(loop_metrics(b, st["b"]), loop_metrics(c, st["c"]), loop_metrics(q, st["q"]))
  save_fig(p, "M5_pv_loops", long, w = 18, h = 10)
  write.csv(m, file.path(out_dir, "data", "M5_pv_loop_metrics.csv"), row.names = FALSE)
  print(m)
})

# ---- M6: measurement-noise simulation (after the batch) --------------------
if (nzchar(noise_dir)) run_fig("M6_noise", function() {
  P <- read.csv(file.path(noise_dir, "noise_parameters_long.csv"), stringsAsFactors = FALSE)
  S <- read.csv(file.path(noise_dir, "noise_summary_parameters.csv"), stringsAsFactors = FALSE)
  O <- read.csv(file.path(noise_dir, "noise_post_long.csv"), stringsAsFactors = FALSE)
  P <- P %>% left_join(S %>% select(Parameter, Reference), by = "Parameter") %>%
    mutate(Ratio = Calibrated / Reference, Label = cand_label[Parameter])
  pa <- ggplot(P, aes(x = Label, y = Ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey92") +
    geom_jitter(width = 0.12, size = 0.8, colour = col_pre, alpha = 0.7) +
    labs(x = NULL, y = "Calibrated value / reference calibration") + theme_paper +
    coord_flip()
  post_label <- c(SAP_max = "Systolic pressure, RFA", SAP_min = "Diastolic pressure, RFA",
                  SAP_mean = "Mean arterial pressure, RFA", PAP_max = "PA systolic pressure",
                  PAP_min = "PA diastolic pressure", PAP_mean = "Mean PA pressure", RAP_mean = "Mean RA pressure",
                  CO_Lmin = "Systemic flow Qs", QpQs = "Qp/Qs", SVR = "SVR", LVEDV = "LVEDV", LVESV = "LVESV",
                  RVEDV = "RVEDV", RVESV = "RVESV", LVEF = "LVEF", RVEF = "RVEF")
  O <- O %>% mutate(Label = factor(post_label[Metric], levels = rev(post_label)))
  pb <- ggplot(O, aes(x = Error_pct, y = Label)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA, colour = col_post) +
    geom_jitter(height = 0.12, size = 0.6, colour = col_post, alpha = 0.6) +
    labs(x = "Prediction minus measurement (% of measurement)", y = NULL) + theme_paper
  p <- pa / pb + plot_layout(heights = c(0.7, 1.3)) + plot_annotation(tag_levels = "A")
  save_fig(p, "M6_noise", O, w = 16, h = 18)
})
cat("done\n")
