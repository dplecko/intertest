# MIMIC-IV interaction testing
# Produces a 3-column table: measure | difference-scale p-val | log-risk p-val
#
# Run from the ia-testing project root:
#   Rscript scripts/mimic-interactions.R

library(intertest)
library(data.table)
library(ggplot2)

# ── load data ─────────────────────────────────────────────────────────────────

root      <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
miiv_path <- file.path(root, "data", "dat-miiv-death.RData")
load(miiv_path)  # loads: dat

# apply split_elective: elective surgical cases get diag_index shifted by +20
dat <- as.data.frame(dat)
dat$death      <- as.integer(dat$death)
dat$diag_index <- ifelse(dat$diag_index >= 5,
                         dat$diag_index + 20L * dat$elective,
                         dat$diag_index)

X <- "majority"
Z <- c("age", "sex")
W <- c("charlson", "acu_24", "diag_index", "elective")
Y <- "death"

dat <- dat[complete.cases(dat[, c(X, Z, W, Y)]), ]
cat("n =", nrow(dat), "\n\n")

# ── interaction testing ───────────────────────────────────────────────────────

set.seed(2024)
fit <- inter_test(dat, X, Z, W, Y,
                  scale = c("difference", "log-risk"),
                  nthread = 4)

# ── format LaTeX table ────────────────────────────────────────────────────────
ia_measures <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
fmt_pval <- function(p) {
  stars <- cut(p,
               breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
               labels = c("***", "**", "*", ".", ""))
  paste0(formatC(p, format = "f", digits = 4), "$^{", stars, "}$")
}
pull <- function(scale) {
  sub <- fit$results[fit$results$scale == scale &
                       fit$results$measure %in% ia_measures, ]
  sub <- sub[match(ia_measures, sub$measure), ]
  fmt_pval(sub$p_value)
}
tbl <- data.frame(
  Measure      = ia_measures,
  Difference   = pull("difference"),
  `Log-risk`   = pull("log-risk"),
  check.names  = FALSE,
  stringsAsFactors = FALSE
)

to_tex <- function(tbl) {
  
  format_meas <- function(m) gsub(" x ", " $\\\\oplus$ ", m)
  
  cat("\\begin{table}[t]\n")
  cat("\\centering\n")
  cat("\\caption{Interaction effects on MIMIC-IV (majority vs. minority).}\n")
  cat("\\label{tab:mimic-interactions}\n")
  cat("\\begin{tabular}{lcc}\n")
  cat("\\toprule\n")
  cat("Measure & Difference & Log-risk \\\\\n")
  cat("\\midrule\n")
  for (i in seq_len(nrow(tbl))) {
    cat(format_meas(tbl$Measure[i]), "&", tbl$Difference[i], "&", 
        tbl$`Log-risk`[i], "\\\\\n")
  }
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\\\[2pt]\n")
  cat("{\\footnotesize Signif.: $^{***}\\,p<0.001$, $^{**}\\,p<0.01$, $^{*}\\,p<0.05$, $^{.}\\,p<0.1$.}\n")
  cat("\\end{table}\n")
}

to_tex(tbl)

# look at heterogeneity across admission groups
diag_grp <- list(
  "Medical"              = list(diag_index = 0:3),
  "Surgical (Emergency)" = list(diag_index = 5:15),
  "Surgical (Elective)"  = list(diag_index = 25:35)
)

cat("Group sizes:\n")
for (nm in names(diag_grp)) {
  n_g <- sum(dat$diag_index %in% diag_grp[[nm]]$diag_index)
  cat(" ", nm, ":", n_g, "\n")
}

# ── conditional DE ────────────────────────────────────────────────────────────

set.seed(2024)
source(file.path(rprojroot::find_root(rprojroot::has_file(".gitignore")),
                 "scripts", "helpers-cnd-de.R"))

cde <- cnd_de(dat, X, Z, W, Y, E_lst = diag_grp, nthread = 4)

cat("\nConditional direct effects:\n")
print(cde, row.names = FALSE, digits = 4)

# ── plot ──────────────────────────────────────────────────────────────────────

cde$group <- factor(cde$group, levels = names(diag_grp))

p <- ggplot(cde, aes(x = group, y = effect)) +
  geom_col(fill = "grey25", width = 0.5, colour = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, colour = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    x = NULL,
    y = "Ctf-DE (risk difference)",
    title = "Direct racial disparity in ICU mortality by admission type",
    subtitle = "MIMIC-IV  |  majority vs. minority  |  95% CI"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    axis.text.x   = element_text(size = 11)
  )

out_path <- file.path(root, "results", "mimic-de-cond.png")
if (!dir.exists(file.path(root, "results")))
  dir.create(file.path(root, "results"))

ggsave(out_path, p, width = 7, height = 4, bg = "white")
cat("\nPlot saved to", out_path, "\n")
