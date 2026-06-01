
# suppress R CMD check notes for ggplot2 column-name variables
utils::globalVariables(c("measure", "value", "type", "ci_lo", "ci_hi",
                          "scale_label"))

# IA measures shown first in every output; xspec shown below in summary()
.ia_measures    <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
.xspec_measures <- c("tv", "ctfde", "ctfie", "ctfse", "ett")

# display labels (internal key -> printed/plotted label)
.measure_labels <- c(
  "tv"           = "TV",
  "ctfde"        = "Ctf-DE",
  "ctfie"        = "Ctf-IE",
  "ctfse"        = "Ctf-SE",
  "ett"          = "ETT",
  "TE x SE"      = "TE x SE",
  "DE x IE"      = "DE x IE",
  "DE x SE"      = "DE x SE",
  "IE x SE"      = "IE x SE",
  "DE x IE x SE" = "DE x IE x SE"
)

# label for the scale column in printed output
.scale_label <- c(
  "difference" = "Difference scale",
  "log-risk"   = "Log-risk scale",
  "log-odds"   = "Log-odds scale"
)

.sig_stars <- function(p) {
  cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
      labels = c("***", "**", "*", ".", " "))
}

# order rows: ia measures first, then xspec; within each group, preserve order
.order_rows <- function(df) {

  ia_rows    <- df[df$measure %in% .ia_measures, , drop = FALSE]
  xspec_rows <- df[df$measure %in% .xspec_measures, , drop = FALSE]

  ia_rows$measure    <- factor(ia_rows$measure,    levels = .ia_measures)
  xspec_rows$measure <- factor(xspec_rows$measure, levels = .xspec_measures)

  rbind(
    ia_rows[order(ia_rows$measure), ],
    xspec_rows[order(xspec_rows$measure), ]
  )
}

# print
#' @export
print.intertest <- function(x, ...) {

  cat("Call: ")
  print(x$call)
  cat("\n")

  cat("Graph:  X =", paste(x$graph$X, collapse = ", "))
  if (length(x$graph$Z) > 0) cat(" | Z =", paste(x$graph$Z, collapse = ", "))
  if (length(x$graph$W) > 0) cat(" | W =", paste(x$graph$W, collapse = ", "))
  cat(" | Y =", paste(x$graph$Y, collapse = ", "), "\n")
  cat("n =", x$n, "\n\n")

  for (s in unique(x$results$scale)) {

    cat(.scale_label[s], "\n")
    df <- x$results[x$results$scale == s & x$results$measure %in% .ia_measures, ]
    df <- .order_rows(df)
    df$measure <- as.character(df$measure)

    # format columns
    fmt <- data.frame(
      measure  = format(df$measure, justify = "left"),
      estimate = formatC(df$value,   format = "f", digits = 4, width = 9),
      SE       = formatC(df$sd,      format = "f", digits = 4, width = 7),
      p.value  = formatC(df$p_value, format = "f", digits = 4, width = 8),
      stars    = as.character(.sig_stars(df$p_value)),
      stringsAsFactors = FALSE
    )

    header <- sprintf("  %-18s %9s %7s %8s\n",
                      "measure", "estimate", "SE", "p-value")
    cat(header)
    cat(strrep("-", nchar(trimws(header, "right")) + 4), "\n")

    for (i in seq_len(nrow(fmt))) {
      cat(sprintf("  %-18s %s %s %s %s\n",
                  fmt$measure[i], fmt$estimate[i], fmt$SE[i],
                  fmt$p.value[i], fmt$stars[i]))
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  }

  if (length(unique(x$results$scale)) > 0 &&
      any(x$results$measure %in% .xspec_measures)) {
    cat("Use summary() for counterfactual (x-specific) effects.\n")
  }

  invisible(x)
}

# summary
#' @export
summary.intertest <- function(object, ...) {

  cat("Call: ")
  print(object$call)
  cat("\n")

  cat("Graph:  X =", paste(object$graph$X, collapse = ", "))
  if (length(object$graph$Z) > 0)
    cat(" | Z =", paste(object$graph$Z, collapse = ", "))
  if (length(object$graph$W) > 0)
    cat(" | W =", paste(object$graph$W, collapse = ", "))
  cat(" | Y =", paste(object$graph$Y, collapse = ", "), "\n")
  cat("n =", object$n, "\n\n")

  for (s in unique(object$results$scale)) {

    cat(.scale_label[s], "\n\n")
    df <- object$results[object$results$scale == s, ]
    df <- .order_rows(df)
    df$measure <- as.character(df$measure)

    # section: interaction measures
    ia_df <- df[df$measure %in% .ia_measures, , drop = FALSE]
    if (nrow(ia_df) > 0) {

      cat("  Interaction effects:\n")
      .print_results_table(ia_df)
    }

    # section: counterfactual (x-specific) effects
    xs_df <- df[df$measure %in% .xspec_measures, , drop = FALSE]
    if (nrow(xs_df) > 0) {

      cat("\n  Counterfactual effects:\n")
      .print_results_table(xs_df)
    }

    cat("\n")
  }

  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  invisible(object)
}

.print_results_table <- function(df) {

  ci_lo <- df$value - 1.96 * df$sd
  ci_hi <- df$value + 1.96 * df$sd

  # apply display labels
  df$measure <- ifelse(df$measure %in% names(.measure_labels),
                       .measure_labels[df$measure], df$measure)

  fmt <- data.frame(
    measure  = format(df$measure, justify = "left"),
    estimate = formatC(df$value,   format = "f", digits = 4, width = 9),
    SE       = formatC(df$sd,      format = "f", digits = 4, width = 7),
    ci_lo    = formatC(ci_lo,      format = "f", digits = 4, width = 7),
    ci_hi    = formatC(ci_hi,      format = "f", digits = 4, width = 7),
    p.value  = formatC(df$p_value, format = "f", digits = 4, width = 8),
    stars    = as.character(.sig_stars(df$p_value)),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-18s %9s %7s  %15s  %8s\n",
              "measure", "estimate", "SE", "95% CI", "p-value"))
  cat(strrep("-", 72), "\n")

  for (i in seq_len(nrow(fmt))) {
    ci_str <- sprintf("[%s, %s]", fmt$ci_lo[i], fmt$ci_hi[i])
    cat(sprintf("  %-18s %s %s  %-17s %s %s\n",
                fmt$measure[i], fmt$estimate[i], fmt$SE[i],
                ci_str, fmt$p.value[i], fmt$stars[i]))
  }
}

# autoplot

# measure type: used for coloring
.measure_type <- function(measures) {
  dplyr_like <- function(m) {
    ifelse(m %in% .xspec_measures, "1st-order",
    ifelse(m %in% c("TE x SE", "DE x IE", "DE x SE", "IE x SE"), "2nd-order",
           "3rd-order"))
  }
  dplyr_like(measures)
}

#' @importFrom ggplot2 autoplot ggplot aes geom_col geom_errorbar coord_flip theme_bw labs scale_fill_manual scale_x_discrete facet_wrap geom_hline theme element_text
#' @export
autoplot.intertest <- function(object, scales_arg = "fixed", ...) {

  df <- object$results

  # factor level order (bottom to top after coord_flip):
  #   3rd-order | 2nd-order | 1st-order
  # so the plot reads top-to-bottom: 1st -> 2nd -> 3rd
  ia_3rd <- .ia_measures[5]
  ia_2nd <- .ia_measures[1:4]
  measure_order <- c(
    ia_3rd[ia_3rd %in% df$measure],
    rev(ia_2nd[ia_2nd %in% df$measure]),
    rev(.xspec_measures[.xspec_measures %in% df$measure])
  )
  df$measure <- factor(df$measure, levels = measure_order)
  df$type    <- .measure_type(as.character(df$measure))
  df$type    <- factor(df$type, levels = c("1st-order", "2nd-order", "3rd-order"))
  df$ci_lo   <- df$value - 1.96 * df$sd
  df$ci_hi   <- df$value + 1.96 * df$sd

  scale_labs <- stats::setNames(.scale_label[unique(df$scale)], unique(df$scale))
  df$scale_label <- scale_labs[df$scale]

  pal <- c("1st-order" = "#4E79A7",
           "2nd-order" = "#F28E2B",
           "3rd-order" = "#E15759")

  p <- ggplot(df, aes(x = measure, y = value, fill = type)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    coord_flip() +
    scale_fill_manual(
      values = pal,
      name   = "Effect order",
      drop   = FALSE
    ) +
    scale_x_discrete(labels = .measure_labels) +
    labs(x = NULL, y = "Estimate (95% CI)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11))

  if (length(unique(df$scale)) > 1) {

    p <- p + facet_wrap(~ scale_label, scales = scales_arg)
  } else {

    p <- p + labs(title = scale_labs[unique(df$scale)])
  }

  p
}
