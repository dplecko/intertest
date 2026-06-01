
# nohup Rscript scripts/ia-synthetic.R > parallel-synth.log 2>&1 &
library(intertest)
library(data.table)
library(ggplot2)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
source(file.path(root, "scripts", "helpers-gen-scm.R"))
set.seed(2024)

scm_class <- LETTERS[1:5]
sample_grid <- c(500, 750, 1500, 3000, 5000, 8000)[1:6]
nrep <- 100
ias <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
res <- c()
for (sclass in scm_class) {

  for (sample_size in sample_grid) {

    cat("Class:", sclass, ", n = ", sample_size, "parallelized...\n")

    ia_chunk <- parallel::mclapply(
      seq_len(nrep),
      function(rep) {

        scm <- gen_from_scm(sclass, n = sample_size)
        data <- scm$data
        X <- scm$mapping$X; Z <- scm$mapping$Z
        W <- scm$mapping$W; Y <- scm$mapping$Y
        gt <- scm$gt

        ia_iter <- inter_test(data, X, Z, W, Y, scale = "difference")$results
        ia_iter <- as.data.table(ia_iter)
        ia_iter <- ia_iter[measure %in% ias]
        ia_iter <- merge(ia_iter, gt, by = "measure") # merge-in ground truth
        ia_iter[, c("scm_class","method","sample_size","rep") := 
                  .(sclass, "one-step", sample_size, rep)]
        return(ia_iter)
      }, mc.cores = parallel::detectCores() / 2
    )
    res <- rbind(res, do.call(rbind, ia_chunk))
  }
}


#' * Analysis I: distribution of p-values *
res[, measure := factor(measure, levels = c("TE x SE", "DE x IE", "DE x SE", "IE x SE",
                                            "DE x IE x SE"))]
res[, ia_ind := ifelse(gt, "Interaction", "No Interaction")]
res[, ia_ind := factor(ia_ind, levels = c("No Interaction", "Interaction"))]
ggplot(res, aes(x = p_value, color = factor(scm_class),
                linewidth = sqrt(sample_size),
                group = interaction(scm_class, sample_size))) +
  stat_ecdf() + theme_bw() +
  facet_grid(rows = vars(ia_ind), cols = vars(measure)) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed") +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_linewidth_continuous(range = c(0.3, 2),
                             name = latex2exp::TeX("\\sqrt{sample size}")) +
  scale_color_discrete(
    name = "SCM", labels = sapply(1:5, function(i) latex2exp::TeX(paste0("$M_", i, "$")))
  ) +
  ylab("Empirical Cumulative Distribution Function") +
  xlab("p-value") +
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("results/p-value-distr.png", width = 10, height = 5)

#' * Analysis II: type I and II error rates *
alpha <- 0.05
res[, h0_true := !gt]
res[, reject := p_value < alpha / 2]
res[h0_true == TRUE, err := (reject == 1)]
res[h0_true == FALSE, err := (reject == 0)]

type_err <- res[, list(err = mean(err), h0_true = h0_true),
                by = c("scm_class", "sample_size", "measure")]
type_err[, type := ifelse(h0_true == TRUE, "Type II", "Type I")]
type_err[, type := factor(type, levels = c("Type II", "Type I"))]
ggplot(type_err,
       aes(y = err, x = sample_size, color = scm_class)) +
  geom_line() + theme_bw() +
  facet_grid(rows = vars(type), cols = vars(measure),
             scales = "free") +
  xlab("Sample size") + ylab("Testing Error") +
  scale_y_continuous(labels = scales::percent) +
  geom_hline(data = subset(type_err, type == "Type II"),
             aes(yintercept = 0.05), linetype = "dashed", color = "gray") +
  scale_color_discrete(
    name = "SCM", labels = sapply(1:5, function(i) latex2exp::TeX(paste0("$M_", i, "$")))
  ) +
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave("results/test-errors.png", width = 10, height = 5)
