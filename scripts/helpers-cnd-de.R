# Helpers for conditional direct effect analysis.
# Source this file from analysis scripts; not part of the intertest package.

# Convert a named list of allowed values to a logical row-index vector.
# E.g.: E_to_ind(list(diag_index = 0:3), dat)
E_to_ind <- function(E, data) {

  ind <- rep(TRUE, nrow(data))
  for (i in seq_along(E)) {
    var <- names(E)[i]
    ind <- ind & (data[[var]] %in% E[[i]])
  }
  ind
}

# Conditional counterfactual direct effect E[Y(1,W(0)) - Y(0,W(0)) | E]
# for a list of subgroups E_lst.
#
# Uses intertest::cross_fit + pso_diff on the full data; pseudo-outcomes are
# then restricted to each subgroup to estimate the conditional DE.
#
# Returns a data frame: group, effect, sd, ci_lo, ci_hi.
cnd_de <- function(data, X, Z = character(0), W = character(0), Y,
                   E_lst, ...) {

  cfit   <- intertest:::cross_fit(data, X, Z, W, Y,
                                  nested_mean = "refit", scale = "diff", ...)
  phi    <- intertest:::pso(cfit, data, X, Z, W, Y, scale = "diff", ...)

  phi_de <- phi[[1]][[1]][[2]] - phi[[1]][[1]][[1]]

  res <- NULL
  for (ei in seq_along(E_lst)) {

    E_ind <- E_to_ind(E_lst[[ei]], data)
    phi_e <- phi_de[E_ind]
    n_e   <- sum(!is.na(phi_e))
    est   <- mean(phi_e, na.rm = TRUE)
    se    <- sqrt(var(phi_e, na.rm = TRUE) / n_e)

    res <- rbind(res,
      data.frame(
        group  = names(E_lst)[ei],
        effect = est,
        sd     = se,
        ci_lo  = est - 1.96 * se,
        ci_hi  = est + 1.96 * se,
        stringsAsFactors = FALSE
      )
    )
  }

  res
}
