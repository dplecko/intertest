#' Test for interaction effects in a causal graph
#'
#' Main entry point for the `intertest` package. Estimates and tests for
#' pairwise and higher-order interaction effects between the direct effect (DE),
#' indirect effect (IE), and spurious effect (SE) of a binary exposure X on
#' outcome Y, mediated through confounders Z and mediators W.
#'
#' @param data A data frame. Numeric and integer columns are used as-is.
#'   Logical, character, and factor columns are auto-encoded (binary factors to
#'   0/1; multi-level factors to dummy columns) and a message is emitted
#'   listing the encoded variables. Pass `preproc = FALSE` to disable this
#'   behaviour and supply a fully numeric data frame.
#' @param X Name of the binary exposure variable (character scalar).
#' @param Z Character vector of pre-exposure confounder names. Pass
#'   `character(0)` (default) if absent.
#' @param W Character vector of mediator names. Pass `character(0)` (default)
#'   if absent.
#' @param Y Name of the outcome variable (character scalar).
#' @param scale Character vector, one or more of `"difference"`, `"log-risk"`,
#'   `"log-odds"`. Nuisance estimation is shared across `"difference"` and
#'   `"log-odds"` but is run separately for `"log-risk"`. Default
#'   `"difference"`.
#' @param nested_mean How to estimate the nested conditional mean. `"refit"`
#'   (default) regresses cross-fitted outcome predictions on (X, Z); `"wregr"`
#'   uses a weighted regression on Z only (not compatible with `"log-risk"`).
#' @param eps_trim Propensity score trimming threshold. Default `0` (no
#'   trimming).
#' @param preproc Logical. If `TRUE` (default), non-numeric columns are
#'   auto-encoded and a message is emitted. Set to `FALSE` to bypass encoding
#'   entirely.
#' @param save_pso Logical. If `TRUE`, pseudo-outcomes for each potential
#'   outcome are stored in the returned object (useful for diagnostics). Default
#'   `FALSE`.
#' @param ... Additional arguments forwarded to `xgb.cv` / `xgb.train` (e.g.
#'   `nthread`).
#'
#' @return An object of class `"intertest"`, a list with components:
#'   \describe{
#'     \item{`results`}{Data frame with columns `measure`, `value`, `sd`,
#'       `scale`, and `p_value`.}
#'     \item{`pso`}{List of pseudo-outcomes, one element per scale (only
#'       populated when `save_pso = TRUE`).}
#'     \item{`call`}{The matched call.}
#'     \item{`graph`}{Named list `list(X, Z, W, Y)` giving the variable groups
#'       after any encoding.}
#'     \item{`n`}{Number of observations used (after trimming).}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500
#' Z <- rnorm(n); X <- rbinom(n, 1, plogis(Z))
#' W <- X + Z + rnorm(n); Y <- X * W + rnorm(n)
#' dat <- data.frame(Z = Z, X = X, W = W, Y = Y)
#' fit <- inter_test(dat, X = "X", Z = "Z", W = "W", Y = "Y")
#' summary(fit)
#' }
#'
#' @export
inter_test <- function(data, X, Z = character(0), W = character(0), Y,
                       scale = "difference",
                       nested_mean = c("refit", "wregr"),
                       eps_trim = 0,
                       preproc = TRUE,
                       save_pso = FALSE,
                       ...) {

  cl <- match.call()
  scale <- match.arg(scale, c("difference", "log-risk", "log-odds"),
                     several.ok = TRUE)
  nested_mean <- match.arg(nested_mean)

  if (preproc) {

    pp <- preproc_data(data, X, Z, W, Y)
    data <- pp$data
    X <- pp$sfm$X
    Z <- pp$sfm$Z
    W <- pp$sfm$W
    Y <- pp$sfm$Y

    if (length(pp$encoded) > 0) {
      message("Auto-encoded: ", paste(pp$encoded, collapse = ", "),
              " -- pass preproc = FALSE for full control.")
    }
  }

  y_vals <- data[[Y]]
  if ("log-risk" %in% scale && any(y_vals < 0, na.rm = TRUE))
    stop("scale = 'log-risk' requires Y >= 0.")
  if ("log-odds" %in% scale && any(y_vals < 0 | y_vals > 1, na.rm = TRUE))
    stop("scale = 'log-odds' requires Y in [0, 1] (binary or probability outcome).")

  results  <- NULL
  pso_list <- list()

  for (s in scale) {

    res_s <- one_step_debias(data, X, Z, W, Y, scale = s,
                             nested_mean = nested_mean,
                             eps_trim = eps_trim, ...)
    results <- rbind(results, res_s)
    if (save_pso) pso_list[[s]] <- attr(res_s, "pso")
  }

  results$p_value <- 2 * stats::pnorm(-abs(results$value / results$sd))
  results$scale   <- factor(results$scale, levels = scale)
  results         <- results[order(results$scale), ]
  results$scale   <- as.character(results$scale)

  structure(
    list(
      results = results,
      pso     = if (save_pso) pso_list else NULL,
      call    = cl,
      graph   = list(X = X, Z = Z, W = W, Y = Y),
      n       = nrow(data)
    ),
    class = "intertest"
  )
}
