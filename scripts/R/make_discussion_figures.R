# make_discussion_figures.R
# -----------------------------------------------------------------------
# Discussion figures for the P1 manuscript:
#   D1_closure_response  change at closure, model against measurement
#   D2_zscores           standardized residuals before and after closure
#   D3_noise_boxplots    measurement-noise replicates: calibrated parameters
#                        and pre/post-closure errors
#
# D1 and D2 use the values of the main-text tables (reference run
# 20260920_005419, post-closure medians of repeated readings). D3 reads the
# summary of the noise batch and is regenerated when more replicates finish.
#
# Usage: Rscript make_discussion_figures.R <output_dir> <noise_summary_dir>
# -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript make_discussion_figures.R <output_dir> <noise_summary_dir>")
out_dir <- args[1]; noise_dir <- args[2]
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)

col_pre  <- "#0072B2"
col_post <- "#D55E00"
col_meas <- "grey35"
theme_paper <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5),
        legend.position = "bottom", plot.tag = element_text(face = "bold"))
save_fig <- function(p, name, data, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "cm")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  write.csv(data, file.path(out_dir, "data", paste0(name, ".csv")), row.names = FALSE)
  cat("  written:", name, "\n")
}

# ---- D1: change at closure -------------------------------------------------
closure <- tribble(
  ~Quantity,                      ~MeasPre, ~MeasPost, ~ModelPre, ~ModelPost,
  "Pulmonary flow Qp",            4.087,    NA,        3.645,     3.174,
  "Systemic flow Qs",             3.423,    2.54,      3.149,     3.170,
  "Qp/Qs",                        1.194,    0.90,      1.157,     1.000,
  "Mean PA pressure",             15,       13,        15.82,     14.86,
  "Systolic PA pressure",         20,       17,        25.10,     23.92,
  "Diastolic PA pressure",        10,       9,         9.70,      9.29,
  "Mean arterial pressure",       71.33,    75.00,     69.56,     70.19,
  "Systolic RFA pressure",        100,      89,        88.09,     90.08,
  "Diastolic RFA pressure",       57,       68,        52.75,     51.95,
  "Mean RA pressure",             5,        5,         5.23,      5.41,
  "Mean LA pressure",             NA,       NA,        7.85,      7.92,
  "LVEDV",                        NA,       36.13,     39.8,      41.10,
  "RVEDV",                        NA,       30.5,      43.3,      43.35)
d1 <- closure %>%
  mutate(Model = 100 * (ModelPost - ModelPre) / ModelPre,
         Measured = 100 * (MeasPost - MeasPre) / MeasPre) %>%
  select(Quantity, Model, Measured) %>%
  pivot_longer(-Quantity, names_to = "Source", values_to = "Change") %>%
  mutate(Quantity = factor(Quantity, levels = rev(closure$Quantity)),
         Source = factor(Source, levels = c("Model", "Measured")))
p1 <- ggplot(mapping = aes(y = Quantity)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_col(data = filter(d1, Source == "Model"), aes(x = Change, fill = "Model"), width = 0.6) +
  geom_point(data = filter(d1, Source == "Measured", !is.na(Change)), aes(x = Change, shape = "Measured"),
             size = 2.4, colour = "black") +
  scale_fill_manual(values = c(Model = col_post), name = NULL) +
  scale_shape_manual(values = c(Measured = 18), name = NULL) +
  labs(x = "Change at closure (% of pre-closure value)", y = NULL) + theme_paper
save_fig(p1, "D1_closure_response", d1, w = 14, h = 10)

# ---- D2: standardized residuals before and after closure --------------------
ztab <- tribble(
  ~Quantity,                ~Pre,   ~Post,
  "Systolic RFA pressure",  -1.19,  0.12,
  "Diastolic RFA pressure", -0.75, -2.36,
  "Mean arterial pressure", -0.50, -1.28,
  "Systolic PA pressure",    5.10,  8.14,
  "Diastolic PA pressure",  -0.60,  0.65,
  "Mean PA pressure",        1.09,  2.86,
  "Mean RA pressure",        0.93,  1.64,
  "Systemic flow Qs",       -0.55,  4.97,
  "Qp/Qs",                  -0.61,  2.23,
  "SVR",                      NA,  -1.29,
  "LVEDV",                    NA,   1.38,
  "LVESV",                    NA,  -1.58,
  "RVEDV",                    NA,   4.21,
  "RVESV",                    NA,   1.48,
  "LVEF",                     NA,   1.61,
  "RVEF",                     NA,   1.25)
d2 <- ztab %>% pivot_longer(-Quantity, names_to = "State", values_to = "z") %>%
  filter(!is.na(z)) %>%
  mutate(Quantity = factor(Quantity, levels = rev(ztab$Quantity)),
         State = factor(State, levels = c("Pre", "Post"),
                        labels = c("Pre-closure, calibrated", "Post-closure, predicted")))
p2 <- ggplot(d2, aes(x = z, y = Quantity, colour = State, shape = State)) +
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf, fill = "grey92") +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "grey82") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(col_pre, col_post), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(breaks = seq(-2, 8, 2)) +
  labs(x = expression("Standardized residual " * z == (model - measured) / sigma), y = NULL) +
  theme_paper
save_fig(p2, "D2_zscores", d2, w = 14, h = 10)

# ---- D3: measurement-noise replicates --------------------------------------
par_lab <- c("group.R_sys_scale" = "Systemic resistance scale", "R.SVEN" = "Systemic venous resistance",
             "group.R_pul_scale" = "Pulmonary resistance scale", "C.SAR" = "Systemic arterial compliance",
             "C.PAR" = "Pulmonary arterial compliance", "E.RV.EB" = "RV passive elastance",
             "vsd.Cd" = "Discharge coefficient")
met_lab <- c(SAP_max = "Systolic RFA", SAP_min = "Diastolic RFA", SAP_mean = "Mean arterial",
             PAP_max = "Systolic PA", PAP_min = "Diastolic PA", PAP_mean = "Mean PA", RAP_mean = "Mean RA",
             CO_Lmin = "Qs", QpQs = "Qp/Qs", LVEDV = "LVEDV", LVESV = "LVESV", RVEDV = "RVEDV", RVESV = "RVESV")
P <- read.csv(file.path(noise_dir, "noise_parameters_long.csv"), stringsAsFactors = FALSE)
SP <- read.csv(file.path(noise_dir, "noise_summary_parameters.csv"), stringsAsFactors = FALSE)
PRE <- read.csv(file.path(noise_dir, "noise_pre_long.csv"), stringsAsFactors = FALSE)
POST <- read.csv(file.path(noise_dir, "noise_post_long.csv"), stringsAsFactors = FALSE)
n_rep <- length(unique(P$Replicate))

pa <- P %>% left_join(SP %>% select(Parameter, Reference), by = "Parameter") %>%
  filter(Parameter %in% names(par_lab)) %>%
  mutate(Ratio = Calibrated / Reference, Label = factor(par_lab[Parameter], levels = rev(par_lab)))
pA <- ggplot(pa, aes(x = Ratio, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(fill = "grey90", outlier.shape = NA, width = 0.6) +
  geom_jitter(height = 0.12, width = 0, size = 1.1, colour = col_pre, alpha = 0.8) +
  scale_x_log10() +
  labs(x = "Calibrated value / reference calibration (log scale)", y = NULL) + theme_paper

err <- bind_rows(
  PRE %>% filter(Metric %in% names(met_lab)) %>%
    transmute(Replicate, Metric, Error = ErrorVsOriginal_pct, State = "Pre-closure, calibrated"),
  POST %>% filter(Metric %in% names(met_lab)) %>%
    transmute(Replicate, Metric, Error = Error_pct, State = "Post-closure, predicted")) %>%
  mutate(Label = factor(met_lab[Metric], levels = rev(met_lab)),
         State = factor(State, levels = c("Pre-closure, calibrated", "Post-closure, predicted")))
pB <- ggplot(err, aes(x = Error, y = Label, fill = State)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_point(aes(colour = State), position = position_jitterdodge(jitter.height = 0.1, jitter.width = 0,
             dodge.width = 0.8), size = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = c(col_pre, col_post), name = NULL) +
  scale_colour_manual(values = c(col_pre, col_post)) +
  labs(x = "Model minus measurement (% of measurement)", y = NULL) + theme_paper
p3 <- pA + pB + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A", caption = sprintf("n = %d recalibrated replicates", n_rep))
save_fig(p3, "D3_noise_boxplots", err, w = 18, h = 10)
