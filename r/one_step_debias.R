#' One-step debiased estimator for potential outcome means
#'
#' Workhorse estimator for a single scale. Cross-fits nuisance functions via
#' gradient boosting and returns a data frame of effect estimates with standard
#' errors.
#'
#' @param data A numeric data frame (all columns must be numeric or integer; use
#'   `inter_test()` for automatic factor encoding).
#' @param X Name of the binary exposure variable.
#' @param Z Character vector of pre-exposure confounder names. Pass
#'   `character(0)` if absent.
#' @param W Character vector of mediator names. Pass `character(0)` if absent.
#' @param Y Name of the outcome variable.
#' @param scale One of `"difference"`, `"log-risk"`, or `"log-odds"`.
#' @param nested_mean How to estimate the nested conditional mean
#'   E\[E\[Y(x') | X=x, Z, W\] | X=x, Z\]. `"refit"` (default) regresses
#'   cross-fitted outcome predictions on (X, Z); `"wregr"` uses a weighted
#'   regression on Z only (not available for `"log-risk"`).
#' @param eps_trim Propensity score trimming threshold. Observations with
#'   P(X=x | Z) or P(X=x | Z, W) below this value are excluded. Default `0`
#'   (no trimming).
#' @param cfit Optional pre-computed cross-fit object (output of the internal
#'   `cross_fit()` function). When supplied, the cross-fitting step is skipped.
#'   Useful when `inter_test()` shares one cross-fit across multiple scales.
#' @param ... Additional arguments passed to `xgb.cv` / `xgb.train` (e.g.
#'   `nthread`).
#'
#' @return A data frame with columns `measure`, `value`, `sd`, and `scale`, one
#'   row per effect. The data frame carries a `"pw"` attribute with the fitted
#'   propensity weights.
#'
#' @export
one_step_debias <- function(data, X, Z = character(0), W = character(0), Y,
                            scale = c("difference", "log-risk", "log-odds"),
                            nested_mean = c("refit", "wregr"),
                            eps_trim = 0, cfit = NULL, ...) {

  scale <- match.arg(scale)
  nested_mean <- match.arg(nested_mean)
  sc <- switch(scale,
               "difference" = "diff", "log-risk" = "logr", "log-odds" = "logo")

  if (is.null(cfit))
    cfit <- cross_fit(data, X, Z, W, Y, nested_mean, scale = sc, ...)

  phi <- pso(cfit, data, X, Z, W, Y, scale = sc, ...)

  # trim extreme propensity observations
  extrm_pxz  <- (cfit$px_z[[1]]  < eps_trim) | (1 - cfit$px_z[[1]]  < eps_trim)
  extrm_pxzw <- (cfit$px_zw[[1]] < eps_trim) | (1 - cfit$px_zw[[1]] < eps_trim)
  extrm_idx  <- extrm_pxz | extrm_pxzw

  if (mean(extrm_idx) > 0.02) {
    message(round(100 * mean(extrm_idx), 2),
            "% extreme propensity weights at threshold = ", eps_trim,
            ". Results are for the overlap population.")
  }

  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1))
    phi[[xz+1]][[xw+1]][[xy+1]][extrm_idx] <- NA

  ias <- measure_spec()
  res <- NULL

  for (i in seq_along(ias)) {

    pseudo_out <- 0
    for (j in seq_along(ias[[i]]$sgn)) {

      xz <- ias[[i]]$spc[[j]][1]
      xw <- ias[[i]]$spc[[j]][2]
      xy <- ias[[i]]$spc[[j]][3]
      pseudo_out <- pseudo_out + ias[[i]]$sgn[j] * phi[[xz+1]][[xw+1]][[xy+1]]
    }
    psi_osd <- mean(pseudo_out, na.rm = TRUE)
    dev <- sqrt(stats::var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))

    res <- rbind(res,
      data.frame(measure = ias[[i]]$ia, value = psi_osd, sd = dev,
                 scale = scale, stringsAsFactors = FALSE))
  }

  attr(res, "pw")  <- list(px_z = cfit$px_z, px_zw = cfit$px_zw)
  attr(res, "pso") <- phi
  res
}
